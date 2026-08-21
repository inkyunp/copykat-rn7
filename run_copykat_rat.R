#!/usr/bin/env Rscript

# Run the CopyKAT rn7 fork on a Seurat v4/v5 RDS or a plain gene x cell counts
# matrix RDS. Supports three genomes: human (hg20), mouse (mm10), and rat (rn7).
# Records candidate CNV evidence only; it never assigns a final malignant/normal
# diagnosis (results remain review_required).
#
# The input contract is enforced by validate_rat_input.R; this runner also
# re-derives id.type and (optionally) known-normal barcodes.
#
# Usage (inside the copykat-rn7 image):
#   singularity exec --cleanenv --containall copykat-rn7.sif \
#     Rscript /work/run_copykat_rat.R \
#       --input /work/rat_copykat_test_500.rds \
#       --output /work/rat_copykat_result.rds \
#       [--genome rn7|mm10|hg20] \
#       [--assay RNA] [--sample-id ratpilot] [--cores 1] \
#       [--id-type auto|S|E] \
#       [--norm-cells /work/norm_cell_names.txt] \
#       [--celltype-column broad_lineage --normal-labels "T cell,B cell"]
#
# --genome selects the species gene-space path (default rn7):
#       hg20  human (GRCh38), id.type="S" matches HGNC symbols
#       mm10  mouse (GRCm38), id.type="S" matches MGI symbols
#       rn7   rat   (mRatBN7.2), id.type="S" matches rat symbols
#
# CopyKAT tuning parameters (call-time args of copykat(); changing them does NOT
# require re-editing or rebuilding the fork/image):
#       [--ngene-chr 5]     genes required per chromosome label per cell (rat has 21 labels)
#       [--low-dr 0.05]     min detection rate to keep a gene for smoothing/segmentation
#       [--up-dr 0.1]       min detection rate to keep a gene in the final result
#       [--ks-cut 0.1]      KS test cutoff for segmentation sensitivity (lower = more segments)
#       [--win-size 25]     min genes per segmentation window

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(args) {
  out <- list(genome = "rn7", assay = "RNA", sample_id = NULL, cores = 1L,
              id_type = "auto", output = NULL,
              ngene_chr = 5L, low_dr = 0.05, up_dr = 0.1, ks_cut = 0.1,
              win_size = 25L)
  val <- c("input","output","genome","assay","sample-id","cores","id-type",
           "norm-cells","celltype-column","normal-labels",
           "ngene-chr","low-dr","up-dr","ks-cut","win-size")
  i <- 1L
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
  if (is.null(out$output)) out$output <- sub("\\.rds$", "", out$input, ignore.case = TRUE)
  out$output <- sub("\\.rds$", "", out$output, ignore.case = TRUE)
  out$genome <- tolower(out$genome)
  if (!(out$genome %in% c("hg20","mm10","rn7")))
    stop("--genome must be one of: hg20, mm10, rn7")
  out$cores <- as.integer(out$cores)
  # CopyKAT tuning parameters (call-time args of copykat(); NOT baked into the fork)
  out$ngene_chr <- as.integer(out$ngene_chr)
  out$win_size  <- as.integer(out$win_size)
  out$low_dr    <- as.numeric(out$low_dr)
  out$up_dr     <- as.numeric(out$up_dr)
  out$ks_cut    <- as.numeric(out$ks_cut)
  if (is.na(out$ngene_chr) || out$ngene_chr < 1L) stop("--ngene-chr must be a positive integer")
  if (is.na(out$win_size)  || out$win_size  < 1L) stop("--win-size must be a positive integer")
  for (nm in c("low_dr","up_dr","ks_cut"))
    if (is.na(out[[nm]]) || out[[nm]] < 0 || out[[nm]] > 1) stop(sprintf("--%s must be in [0,1]", gsub("_","-",nm)))
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages(library(copykat))
stopifnot("annotateGenes.rn7" %in% getNamespaceExports("copykat"))

# genome -> (annotation object, symbol column). copykat body references the
# symbol column name directly, so it differs per genome.
geno_cfg <- list(
  hg20 = list(anno = "full.anno",      sym = "hgnc_symbol"),
  mm10 = list(anno = "full.anno.mm10", sym = "mgi_symbol"),
  rn7  = list(anno = "full.anno.rn7",  sym = "mgi_symbol")
)[[opt$genome]]
message(sprintf("[run] genome=%s (annotation=%s, symbol column=%s)",
                opt$genome, geno_cfg$anno, geno_cfg$sym))

message("[run] reading ", opt$input)
obj <- readRDS(opt$input)

extract_counts <- function(obj, assay) {
  if (inherits(obj, "Seurat")) {
    out <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
                    error = function(e) NULL)
    if (is.null(out) || !length(out))
      out <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, slot = "counts"),
                      error = function(e) NULL)
    if (is.null(out) || !length(out))
      out <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
                      error = function(e) NULL)
    list(counts = out, meta = obj@meta.data)
  } else if (inherits(obj, c("dgCMatrix","matrix","Matrix"))) {
    list(counts = obj, meta = NULL)
  } else stop("unsupported input class: ", paste(class(obj), collapse = "/"))
}

ex <- extract_counts(obj, opt$assay)
counts <- as.matrix(ex$counts)
storage.mode(counts) <- "numeric"
if (any(counts < 0) || mean(abs(counts - round(counts)) < 1e-8) < 0.999)
  stop("input does not look like raw integer counts; run validate_rat_input.R first")
message(sprintf("[run] counts: %d genes x %d cells", nrow(counts), ncol(counts)))

# id.type
anno_ns <- get(geno_cfg$anno, asNamespace("copykat"))
id_type <- toupper(opt$id_type)
if (id_type == "AUTO") {
  n_sym <- length(intersect(rownames(counts), anno_ns[[geno_cfg$sym]]))
  n_ens <- length(intersect(rownames(counts), anno_ns$ensembl_gene_id))
  id_type <- if (n_ens > n_sym) "E" else "S"
  message(sprintf("[run] id.type auto -> %s (symbol=%d, ensembl=%d)", id_type, n_sym, n_ens))
}

# known-normal barcodes (optional)
norm_cells <- ""
if (!is.null(opt$norm_cells) && file.exists(opt$norm_cells)) {
  nc <- readLines(opt$norm_cells); nc <- nc[nzchar(nc)]
  norm_cells <- intersect(nc, colnames(counts))
  message(sprintf("[run] known-normal cells matched: %d/%d", length(norm_cells), length(nc)))
  if (!length(norm_cells)) norm_cells <- ""
} else if (!is.null(opt$celltype_column) && !is.null(opt$normal_labels) && !is.null(ex$meta)) {
  labs <- trimws(strsplit(opt$normal_labels, ",")[[1]])
  col <- ex$meta[[opt$celltype_column]]
  if (!is.null(col)) {
    nc <- rownames(ex$meta)[as.character(col) %in% labs]
    norm_cells <- intersect(nc, colnames(counts))
    message(sprintf("[run] known-normal from %s%%in%%{%s}: %d cells",
                    opt$celltype_column, paste(labs, collapse=","), length(norm_cells)))
    if (!length(norm_cells)) norm_cells <- ""
  }
}

sample_id <- opt$sample_id %||% sub("\\.rds$", "", basename(opt$input), ignore.case = TRUE)
run_dir <- paste0(opt$output, "_", opt$genome, "_run")
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
old <- setwd(run_dir); on.exit(setwd(old), add = TRUE)

message(sprintf("[run] running copykat(genome=\"%s\") ...", opt$genome))
ck <- copykat(rawmat = counts, id.type = id_type, genome = opt$genome,
              norm.cell.names = norm_cells, sam.name = sample_id,
              ngene.chr = opt$ngene_chr, LOW.DR = opt$low_dr, UP.DR = opt$up_dr,
              KS.cut = opt$ks_cut, win.size = opt$win_size,
              n.cores = opt$cores, plot.genes = "TRUE")
message(sprintf("[run] params: ngene.chr=%d LOW.DR=%.3f UP.DR=%.3f KS.cut=%.3f win.size=%d",
                opt$ngene_chr, opt$low_dr, opt$up_dr, opt$ks_cut, opt$win_size))

pred <- ck$prediction
cls <- setNames(as.character(pred[[2]]), as.character(pred[[1]]))
# align to matrix columns; cells copykat dropped are not.defined
classes <- cls[colnames(counts)]
classes[is.na(classes)] <- "not.defined"

setwd(old)
result <- list(
  sample_id = sample_id,
  genome = opt$genome,
  id_type = id_type,
  params = list(ngene.chr = opt$ngene_chr, LOW.DR = opt$low_dr, UP.DR = opt$up_dr,
                KS.cut = opt$ks_cut, win.size = opt$win_size, cores = opt$cores),
  copykat_class = classes,                 # literal: aneuploid/diploid/not.defined
  candidate_state = classes,               # candidate CNV evidence, NOT a diagnosis
  prediction = ck$prediction,
  CNAmat = ck$CNAmat,
  hclustering = ck$hclustering,
  norm_cells_used = norm_cells,
  interpretation = "review_required",      # cross-check with DNA-level evidence before use
  copykat_version = as.character(packageVersion("copykat")),
  run_dir = normalizePath(run_dir),
  note = "Candidate CNV evidence from CopyKAT; hg20 uses the upstream human path, mm10/rn7 use the gene-space (mouse) path. Do not treat aneuploid as a diagnosis."
)
out_rds <- paste0(opt$output, ".rds")
saveRDS(result, out_rds)

ev <- data.frame(cell = names(classes), copykat_class = unname(classes),
                 interpretation = "review_required", stringsAsFactors = FALSE)
write.table(ev, paste0(opt$output, ".cnv_evidence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n[run] class table:\n"); print(table(classes))
cat(sprintf("[run] wrote %s and %s.cnv_evidence.tsv\n", out_rds, opt$output))
cat(sprintf("[run] artifacts in %s\n", normalizePath(run_dir)))
cat("[run] interpretation = review_required (cross-check with WGS/SNP-array/karyotype)\n")
