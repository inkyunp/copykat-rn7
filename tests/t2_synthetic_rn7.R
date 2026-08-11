#!/usr/bin/env Rscript

# T2 synthetic smoke test for the rat (genome="rn7") CopyKAT fork.
#
# Builds a tiny integer gene x cell matrix whose gene symbols are drawn from
# full.anno.rn7$mgi_symbol (>= genes.per.chr genes on every one of the 21 rat
# chromosome labels, so CopyKAT's ngene.chr filter passes), with a "tumor"
# subpopulation carrying an amplified chromosome 1 so segmentation has signal.
# Then runs the MODIFIED copykat(..., genome = "rn7") end to end and checks:
#   - it returns prediction / CNAmat / hclustering
#   - CNA_results.txt (gene-space, mm10/rn7 branch) is written
#   - predicted classes are literal aneuploid / diploid / not.defined
#
# This does not assert biological correctness; it proves the rn7 code path
# executes on the real target R 4.2.2 + real deps + real rat annotation.
#
# Usage:
#   Rscript t2_synthetic_rn7.R FORK_R_FILE SYSDATA_RDA [OUTDIR]

args <- commandArgs(trailingOnly = TRUE)
fork_r  <- args[[1]]
sysdata <- args[[2]]
outdir  <- if (length(args) >= 3L) args[[3]] else tempfile("t2_rn7_")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

# ---- load the modified fork (functions + internal data) into one env ----
env <- new.env()
load(sysdata, envir = env)                 # full.anno.rn7, full.anno.mm10, ...
sys.source(fork_r, envir = env)            # copykat, annotateGenes.rn7, ...
stopifnot(is.function(env$copykat), is.function(env$annotateGenes.rn7))

fa <- env$full.anno.rn7

# ---- pick genes.per.chr real rat symbols on each of the 21 chromosomes ----
genes.per.chr <- 40L
chr_labels <- sort(unique(fa$chromosome_name))
stopifnot(length(chr_labels) == 21L)
pick <- unlist(lapply(chr_labels, function(cl) {
  rows <- fa[fa$chromosome_name == cl, ]
  rows <- rows[order(rows$abspos), ]
  head(rows$mgi_symbol, genes.per.chr)
}), use.names = FALSE)
pick <- unique(pick)
G <- length(pick)

# ---- cells: 30 normal + 30 tumor ----
n_norm <- 30L; n_tum <- 30L; N <- n_norm + n_tum
cell_ids <- sprintf("CELL%03d", seq_len(N))
norm_cells <- cell_ids[seq_len(n_norm)]
tum_cells  <- cell_ids[(n_norm + 1L):N]

lambda <- 6
mat <- matrix(rpois(G * N, lambda), nrow = G, ncol = N,
              dimnames = list(pick, cell_ids))

# amplify chromosome 1 in tumor cells (CNV-like gain signal)
chr1_syms <- fa$mgi_symbol[fa$chromosome_name == 1L]
chr1_rows <- which(rownames(mat) %in% chr1_syms)
mat[chr1_rows, tum_cells] <- mat[chr1_rows, tum_cells] +
  rpois(length(chr1_rows) * n_tum, 2 * lambda)

cat(sprintf("[T2] synthetic matrix: %d genes x %d cells (%d normal, %d tumor)\n",
            G, N, n_norm, n_tum))
cat(sprintf("[T2] amplified chr1 genes in tumor cells: %d\n", length(chr1_rows)))

# ---- run the modified copykat on the rn7 path ----
old <- setwd(outdir)
on.exit(setwd(old), add = TRUE)

res <- env$copykat(
  rawmat = mat,
  id.type = "S",
  genome = "rn7",
  norm.cell.names = norm_cells,   # known-normal path (stable on tiny data)
  sam.name = "t2rn7",
  n.cores = 1,
  plot.genes = "FALSE"
)

# ---- assertions ----
cat("\n[T2] returned names:", paste(names(res), collapse = ", "), "\n")
stopifnot(all(c("prediction", "CNAmat", "hclustering") %in% names(res)))

pred <- res$prediction
cls  <- as.character(pred[[2]])
cat("[T2] prediction class table:\n"); print(table(cls))
stopifnot(all(cls %in% c("aneuploid", "diploid", "not.defined")))

cna_file <- file.path(outdir, "t2rn7_copykat_CNA_results.txt")
cat("[T2] CNA_results.txt exists:", file.exists(cna_file), "\n")
stopifnot(file.exists(cna_file))
hdr <- strsplit(readLines(cna_file, n = 1L), "\t")[[1]]
cat("[T2] CNA_results.txt first cols:", paste(head(hdr, 7), collapse = ", "), "\n")
stopifnot(identical(head(hdr, 7),
                    c("abspos","chromosome_name","start_position","end_position",
                      "ensembl_gene_id","mgi_symbol","band")))

cat("\n[T2] rn7 code path executed; gene-space output; literal classes preserved.\n")
cat("[T2] PASS\n")
