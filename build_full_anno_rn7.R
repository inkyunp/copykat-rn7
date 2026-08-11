#!/usr/bin/env Rscript

# Phase 0: Build a CopyKAT-compatible rat annotation table (full.anno.rn7)
# from an mRatBN7.2 Ensembl GTF.
#
# Output schema is IDENTICAL (name + order) to CopyKAT's internal
# full.anno.mm10 so the mm10 gene-space code path can be reused for rat:
#   1 abspos            cumulative genome coordinate (sole ordering key)
#   2 chromosome_name   INTEGER label (autosomes 1..20, X = 21)
#   3 start_position
#   4 end_position
#   5 ensembl_gene_id
#   6 mgi_symbol        <- column name kept verbatim; holds RAT gene symbol
#   7 band              NA in v1 (visualization only, unused by math)
#
# Usage:
#   Rscript build_full_anno_rn7.R [GTF] [OUT_RDS]
# Defaults: GTF = ../genes.gtf , OUT_RDS = full.anno.rn7.rds (next to this script)

args <- commandArgs(trailingOnly = TRUE)

script_dir <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())

gtf <- if (length(args) >= 1L) args[[1]] else file.path(dirname(script_dir), "genes.gtf")
out_rds <- if (length(args) >= 2L) args[[2]] else file.path(script_dir, "full.anno.rn7.rds")

if (!file.exists(gtf)) stop(sprintf("GTF not found: %s", gtf))

# Chromosomes kept, in genome order, with integer encoding (X = 21).
keep_chr <- c(as.character(1:20), "X")
chr_int  <- stats::setNames(c(1:20, 21L), keep_chr)

message(sprintf("[rn7] reading gene rows from %s", gtf))
# Prefilter to gene feature rows with awk to keep this light on RAM.
cmd <- sprintf("grep -v '^#' %s | awk -F'\\t' '$3==\"gene\"'", shQuote(gtf))
con <- pipe(cmd, open = "r")
lines <- readLines(con)
close(con)
message(sprintf("[rn7] gene feature rows: %d", length(lines)))
if (length(lines) == 0L) stop("No gene rows parsed from GTF")

f <- do.call(rbind, strsplit(lines, "\t", fixed = TRUE))
chr   <- f[, 1]
start <- suppressWarnings(as.integer(f[, 4]))
end   <- suppressWarnings(as.integer(f[, 5]))
attr9 <- f[, 9]

gid <- sub('.*gene_id "([^"]+)".*', "\\1", attr9)
gid[!grepl('gene_id "', attr9)] <- NA_character_
gname <- sub('.*gene_name "([^"]+)".*', "\\1", attr9)
gname[!grepl('gene_name "', attr9)] <- NA_character_

# Keep only standard chromosomes (1..20, X); drop Y/MT/scaffolds so that
# as.numeric(chromosome_name) and the 21-label ngene.chr filter stay valid.
keep <- chr %in% keep_chr & !is.na(start) & !is.na(end)
n_all <- length(chr)
chr <- chr[keep]; start <- start[keep]; end <- end[keep]
gid <- gid[keep]; gname <- gname[keep]
message(sprintf("[rn7] genes on 1-20,X: %d (dropped %d off-target rows)",
                length(chr), n_all - length(chr)))

# Symbol column: gene_name when present, else fall back to gene_id.
missing_name <- is.na(gname) | !nzchar(gname)
mgi_symbol <- ifelse(missing_name, gid, gname)
message(sprintf("[rn7] genes with real symbol: %d ; gene_id fallback: %d",
                sum(!missing_name), sum(missing_name)))

enc <- unname(chr_int[chr])

# abspos = per-chromosome offset (cumulative max end of preceding chromosomes)
#          + start_position. This yields a real cumulative coordinate, avoiding
#          the constant-abspos defect in the bundled full.anno.mm10.
ord0 <- order(enc, start)
chr <- chr[ord0]; enc <- enc[ord0]; start <- start[ord0]; end <- end[ord0]
gid <- gid[ord0]; mgi_symbol <- mgi_symbol[ord0]

chrmax  <- tapply(end, enc, max)
chrmax  <- chrmax[order(as.integer(names(chrmax)))]
offset  <- c(0, cumsum(as.numeric(chrmax)))
names(offset) <- c(names(chrmax), NA)
offset  <- offset[seq_along(chrmax)]
names(offset) <- names(chrmax)
abspos  <- as.numeric(offset[as.character(enc)]) + start

full.anno.rn7 <- data.frame(
  abspos          = abspos,
  chromosome_name = as.integer(enc),
  start_position  = as.integer(start),
  end_position    = as.integer(end),
  ensembl_gene_id = gid,
  mgi_symbol      = mgi_symbol,
  band            = NA_character_,
  stringsAsFactors = FALSE
)

# Order by abspos (CopyKAT's sole ordering key) and drop duplicate symbols,
# keeping the first (lowest abspos) occurrence. CopyKAT itself dedupes on
# mgi_symbol in both id.type branches, so doing it here is behavior-consistent
# and deterministic.
full.anno.rn7 <- full.anno.rn7[order(full.anno.rn7$abspos), , drop = FALSE]
dup <- duplicated(full.anno.rn7$mgi_symbol)
if (any(dup)) {
  message(sprintf("[rn7] dropping %d duplicate-symbol rows (first-keep)", sum(dup)))
  full.anno.rn7 <- full.anno.rn7[!dup, , drop = FALSE]
}
rownames(full.anno.rn7) <- NULL

# -------------------- integrity gates --------------------
n <- nrow(full.anno.rn7)
u_abspos <- length(unique(full.anno.rn7$abspos))
labels   <- sort(unique(full.anno.rn7$chromosome_name))
per_chr  <- table(full.anno.rn7$chromosome_name)
min_per_chr <- min(per_chr)
monotonic <- all(diff(full.anno.rn7$abspos) >= 0)
# Same-chromosome rows must be contiguous after abspos ordering.
contiguous <- length(rle(full.anno.rn7$chromosome_name)$lengths) == length(labels)
u_sym <- length(unique(full.anno.rn7$mgi_symbol))
u_ens <- length(unique(full.anno.rn7$ensembl_gene_id))

cat("\n==================== full.anno.rn7 integrity ====================\n")
cat(sprintf("rows (genes)              : %d\n", n))
cat(sprintf("unique abspos             : %d  (must be >> 20)\n", u_abspos))
cat(sprintf("abspos monotonic (sorted) : %s\n", monotonic))
cat(sprintf("chromosome labels         : %d  [%s]\n", length(labels),
            paste(labels, collapse = ",")))
cat(sprintf("chromosome contiguous     : %s\n", contiguous))
cat(sprintf("min genes on any chr      : %d  (must be >= ngene.chr, default 5)\n", min_per_chr))
cat(sprintf("unique mgi_symbol         : %d  (== rows: %s)\n", u_sym, u_sym == n))
cat(sprintf("unique ensembl_gene_id    : %d\n", u_ens))
cat("genes per chromosome:\n"); print(per_chr)
cat("head:\n"); print(utils::head(full.anno.rn7, 4))

stopifnot(
  u_abspos > 20,
  monotonic,
  contiguous,
  identical(colnames(full.anno.rn7),
            c("abspos","chromosome_name","start_position","end_position",
              "ensembl_gene_id","mgi_symbol","band")),
  min_per_chr >= 5L,
  u_sym == n
)

saveRDS(full.anno.rn7, out_rds)
cat(sprintf("\n[rn7] wrote %s (%s)\n", out_rds,
            format(structure(file.info(out_rds)$size, class = "object_size"),
                   units = "auto")))
cat("[rn7] ALL INTEGRITY GATES PASSED\n")
