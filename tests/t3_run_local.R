#!/usr/bin/env Rscript

# T3: run the rat (genome="rn7") CopyKAT fork on a REAL sampled Seurat RDS,
# loading the fork sources directly (no package install needed) so it runs on
# the extracted target R 4.2.2. Mirrors run_copykat_rat.R behavior.
#
# Usage:
#   Rscript t3_run_local.R FORK_R SYSDATA INPUT_RDS OUTDIR \
#     [CELLTYPE_COL] [NORMAL_LABELS_CSV]

args <- commandArgs(trailingOnly = TRUE)
fork_r  <- args[[1]]
sysdata <- args[[2]]
input   <- args[[3]]
outdir  <- args[[4]]
ct_col  <- if (length(args) >= 5L) args[[5]] else NA
norm_csv<- if (length(args) >= 6L) args[[6]] else NA

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

env <- new.env()
load(sysdata, envir = env)
sys.source(fork_r, envir = env)
stopifnot(is.function(env$copykat), is.function(env$annotateGenes.rn7))

suppressPackageStartupMessages(library(Seurat))
message("[T3] reading ", input)
obj <- readRDS(input)

counts <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts"),
                   error = function(e) NULL)
if (is.null(counts) || !length(counts))
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
stopifnot(!any(counts < 0), mean(abs(counts - round(counts)) < 1e-8) > 0.999)
message(sprintf("[T3] counts: %d genes x %d cells (raw integer)", nrow(counts), ncol(counts)))

# id.type auto
anno <- env$full.anno.rn7
n_sym <- length(intersect(rownames(counts), anno$mgi_symbol))
n_ens <- length(intersect(rownames(counts), anno$ensembl_gene_id))
id_type <- if (n_ens > n_sym) "E" else "S"
message(sprintf("[T3] id.type=%s (symbol=%d ensembl=%d)", id_type, n_sym, n_ens))

# known-normal barcodes
norm_cells <- ""
if (!is.na(norm_csv) && !is.na(ct_col)) {
  labs <- trimws(strsplit(norm_csv, ",")[[1]])
  col <- obj@meta.data[[ct_col]]
  nc <- rownames(obj@meta.data)[as.character(col) %in% labs]
  norm_cells <- intersect(nc, colnames(counts))
  message(sprintf("[T3] known-normal (%s in {%s}): %d cells",
                  ct_col, paste(labs, collapse=","), length(norm_cells)))
  if (!length(norm_cells)) norm_cells <- ""
}

old <- setwd(outdir); on.exit(setwd(old), add = TRUE)
t0 <- Sys.time()
ck <- env$copykat(rawmat = counts, id.type = id_type, genome = "rn7",
                  norm.cell.names = norm_cells, sam.name = "ratpilot",
                  n.cores = 1, plot.genes = "TRUE")
message(sprintf("[T3] copykat wall time: %.1f s", as.numeric(difftime(Sys.time(), t0, units="secs"))))

pred <- ck$prediction
cls <- setNames(as.character(pred[[2]]), as.character(pred[[1]]))
classes <- cls[colnames(counts)]; classes[is.na(classes)] <- "not.defined"

# join with celltype for eyeball cross-check
ctv <- if (!is.na(ct_col)) as.character(obj@meta.data[colnames(counts), ct_col]) else NA
ev <- data.frame(cell = colnames(counts), celltype = ctv,
                 copykat_class = unname(classes),
                 interpretation = "review_required", stringsAsFactors = FALSE)
setwd(old)
write.table(ev, file.path(outdir, "ratpilot.cnv_evidence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(prediction = ck$prediction, CNAmat = ck$CNAmat,
             hclustering = ck$hclustering, classes = classes,
             id_type = id_type, norm_cells_used = norm_cells,
             interpretation = "review_required"),
        file.path(outdir, "ratpilot_result.rds"))

cat("\n[T3] class table:\n"); print(table(classes))
if (!is.na(ct_col)) {
  cat("\n[T3] class x celltype (eyeball cross-check):\n")
  print(table(celltype = ctv, class = classes))
}
cat("\n[T3] artifacts in:", normalizePath(outdir), "\n")
cat("[T3] interpretation = review_required (cross-check with WGS/SNP-array/karyotype)\n")
cat("[T3] PASS\n")
