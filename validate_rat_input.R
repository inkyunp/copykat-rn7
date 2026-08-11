#!/usr/bin/env Rscript

# Data-validation gate for the rat CopyKAT (genome="rn7") pipeline.
# Reusable on the big machine (right after sampling, to gate transfer) and
# locally (to re-check the slim RDS before running). It NEVER modifies inputs.
#
# Checks the CopyKAT input contract, verified against copykat 1.1.0 source:
#   copykat() filters by positive-count detection (sum(x>0)) and internally
#   normalizes via log(sqrt(x)+sqrt(x+1)); therefore the input MUST be raw
#   integer UMI counts, genes x cells, never normalized/scaled/SCT/integrated.
#
# Usage:
#   Rscript validate_rat_input.R \
#     --input sample.rds \
#     --anno full.anno.rn7.rds \
#     [--assay RNA] [--celltype-column COL] [--normal-labels "A,B"] \
#     [--norm-cells norm_cell_names.txt] [--ngene-chr 5] \
#     [--report input_validation_report.tsv]
#
# Exit code 0 only if every hard gate passes; non-zero otherwise.

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(args) {
  out <- list(assay = "RNA", ngene_chr = 5L, report = "input_validation_report.tsv")
  i <- 1L
  val <- c("input","anno","assay","celltype-column","normal-labels",
           "norm-cells","ngene-chr","report")
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) stop("unexpected arg: ", a)
    k <- sub("^--", "", a)
    if (!(k %in% val)) stop("unknown option: ", a)
    if (i == length(args)) stop("missing value for ", a)
    out[[gsub("-", "_", k)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  if (is.null(out$input)) stop("--input is required")
  if (is.null(out$anno))  stop("--anno is required (full.anno.rn7.rds)")
  out$ngene_chr <- as.integer(out$ngene_chr)
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

rows <- list()
add <- function(check, status, detail) {
  rows[[length(rows) + 1L]] <<- data.frame(check = check, status = status,
                                            detail = detail, stringsAsFactors = FALSE)
}
hard_fail <- FALSE
fail <- function() hard_fail <<- TRUE

message("[validate] reading ", opt$input)
obj <- readRDS(opt$input)
anno <- readRDS(opt$anno)

# ---- extract a raw gene x cell counts matrix ----
get_counts <- function(obj, assay) {
  if (inherits(obj, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE) &&
        !requireNamespace("Seurat", quietly = TRUE))
      stop("Seurat/SeuratObject needed to read a Seurat object")
    ver <- tryCatch(as.integer(substr(as.character(obj@version), 1, 1)), error = function(e) NA)
    getlayer <- function() {
      # Seurat v5 first, fall back to v4 slot
      out <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
                      error = function(e) NULL)
      if (is.null(out) || !length(out))
        out <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, slot = "counts"),
                        error = function(e) NULL)
      if (is.null(out) || !length(out))
        out <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
                        error = function(e) NULL)
      out
    }
    list(counts = getlayer(), meta = obj@meta.data, sver = ver)
  } else if (inherits(obj, c("dgCMatrix", "matrix", "Matrix"))) {
    list(counts = obj, meta = NULL, sver = NA)
  } else {
    stop("unsupported input class: ", paste(class(obj), collapse = "/"))
  }
}

ext <- get_counts(obj, opt$assay)
counts <- ext$counts
if (is.null(counts) || !length(counts)) { add("counts_extract","FAIL","could not extract counts layer"); fail() } else {
  add("counts_extract","PASS", sprintf("%d genes x %d cells (assay=%s, seurat_major=%s)",
      nrow(counts), ncol(counts), opt$assay, ext$sver %||% NA))
}

# ---- gate 1: raw integer counts, not normalized ----
if (!is.null(counts) && length(counts)) {
  xv <- if (inherits(counts, "dgCMatrix")) counts@x else as.numeric(counts[counts != 0])
  xv <- xv[is.finite(xv)]
  has_neg <- any(xv < 0)
  frac_int <- if (length(xv)) mean(abs(xv - round(xv)) < 1e-8) else NA_real_
  mx <- if (length(xv)) max(xv) else NA_real_
  detail <- sprintf("negatives=%s frac_integer=%.4f max=%.1f", has_neg, frac_int %||% NA, mx %||% NA)
  if (isTRUE(has_neg)) { add("raw_counts_not_normalized","FAIL", paste0(detail, " (negatives -> scale.data/z-scores)")); fail() }
  else if (is.finite(frac_int) && frac_int < 0.999) { add("raw_counts_not_normalized","FAIL", paste0(detail, " (non-integer -> log-normalized data)")); fail() }
  else if (is.finite(mx) && mx <= 12) { add("raw_counts_not_normalized","WARN", paste0(detail, " (low max; verify these are UMI counts, not normalized)")) }
  else { add("raw_counts_not_normalized","PASS", detail) }
}

# ---- gate 2: gene-ID intersection with full.anno.rn7 + id.type ----
id_type <- NA_character_; kept <- NULL
if (!is.null(counts) && length(counts)) {
  g <- rownames(counts)
  n_sym <- length(intersect(g, anno$mgi_symbol))
  n_ens <- length(intersect(g, anno$ensembl_gene_id))
  id_type <- if (n_ens > n_sym) "E" else "S"
  key <- if (id_type == "E") anno$ensembl_gene_id else anno$mgi_symbol
  kept <- intersect(g, key)
  detail <- sprintf("symbol_hits=%d ensembl_hits=%d -> id.type=%s kept_genes=%d/%d",
                    n_sym, n_ens, id_type, length(kept), length(g))
  if (length(kept) >= 8000) add("gene_intersection","PASS", detail)
  else if (length(kept) >= 3000) add("gene_intersection","WARN", paste0(detail, " (low; check same GTF build)"))
  else { add("gene_intersection","FAIL", paste0(detail, " (too low; likely different genome build/gene_name version)")); fail() }
}

# ---- gate 3: per-chromosome coverage of kept genes (21 labels) ----
if (!is.null(kept) && length(kept)) {
  key <- if (id_type == "E") "ensembl_gene_id" else "mgi_symbol"
  sub <- anno[anno[[key]] %in% kept, ]
  per <- table(factor(sub$chromosome_name, levels = 1:21))
  min_per <- min(per); labels_covered <- sum(per > 0)
  detail <- sprintf("labels_with_genes=%d/21 min_genes_per_label=%d", labels_covered, min_per)
  if (labels_covered == 21 && min_per >= opt$ngene_chr) add("chromosome_coverage","PASS", detail)
  else { add("chromosome_coverage","FAIL", paste0(detail, sprintf(" (need 21 labels each >= ngene.chr=%d)", opt$ngene_chr))); fail() }
}

# ---- gate 4: normal cells available ----
normal_ids <- NULL
if (!is.null(opt$norm_cells) && file.exists(opt$norm_cells)) {
  normal_ids <- readLines(opt$norm_cells)
  normal_ids <- normal_ids[nzchar(normal_ids)]
} else if (!is.null(opt$celltype_column) && !is.null(opt$normal_labels) && !is.null(ext$meta)) {
  labs <- trimws(strsplit(opt$normal_labels, ",")[[1]])
  col <- ext$meta[[opt$celltype_column]]
  if (!is.null(col)) normal_ids <- rownames(ext$meta)[as.character(col) %in% labs]
}
if (!is.null(counts) && length(counts)) {
  if (is.null(normal_ids)) {
    add("normal_cells","WARN","no --norm-cells / celltype+labels given; copykat will auto-pick a baseline")
  } else {
    present <- intersect(normal_ids, colnames(counts))
    detail <- sprintf("provided=%d present_in_matrix=%d", length(normal_ids), length(present))
    if (length(present) >= 30) add("normal_cells","PASS", detail)
    else if (length(present) >= 1) add("normal_cells","WARN", paste0(detail, " (<30; auto-baseline may be used)"))
    else { add("normal_cells","FAIL", paste0(detail, " (none present; barcodes do not match)")); fail() }
  }
}

# ---- gate 5: barcode integrity ----
if (!is.null(counts) && length(counts)) {
  cn <- colnames(counts)
  dup <- sum(duplicated(cn))
  ok_names <- !is.null(cn) && all(nzchar(cn))
  detail <- sprintf("cells=%d duplicate_barcodes=%d all_nonempty=%s", length(cn), dup, ok_names)
  if (dup == 0 && ok_names) add("barcode_integrity","PASS", detail)
  else { add("barcode_integrity","FAIL", detail); fail() }
}

report <- do.call(rbind, rows)
write.table(report, opt$report, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n==================== input validation report ====================\n")
print(report, row.names = FALSE)
cat(sprintf("\nid.type (recommended): %s\n", id_type %||% NA))
cat(sprintf("report written: %s\n", opt$report))

if (hard_fail) {
  cat("\n[validate] RESULT: FAIL - do not proceed / do not transfer.\n")
  quit(status = 1L)
}
cat("\n[validate] RESULT: PASS - input satisfies the CopyKAT rn7 contract.\n")
