# 출력 디렉토리 설정
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Output directory argument is required.")

outdir <- args[1]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", outdir, "\n")

# 1. 라이브러리 로드 (copykat-rn7 이미지에는 Seurat 없음 → base-R slot 접근)
suppressPackageStartupMessages(library(Matrix))
library(copykat)

# 2. Seurat 객체 로드 (slot 직접 접근; GetAssayData/subset 미사용)
obj <- readRDS("/data/test.obj_final_add_celltype_scp2711.rds")
meta <- obj@meta.data
counts_all <- obj@assays$RNA@counts   # RNA raw counts (dgCMatrix)

# 3. epithelial cell만 선택
epithelial_types <- c(
  "ciliated cell",
  "pulmonary alveolar type 2 cell",
  "club cell")

ct <- as.character(meta[["ref_cell_type__ontology_label"]])
epi_cells <- rownames(meta)[ct %in% epithelial_types]

# RNA raw count 추출 (genes × cells, 정규화/scale.data 금지)
counts <- as.matrix(counts_all[, epi_cells])

# 4. condition + cell type 기반 group (정상 baseline 선택용)
condition_celltype <- paste(
  meta[epi_cells, "condition"],
  meta[epi_cells, "ref_cell_type__ontology_label"],
  sep = "_"
)
names(condition_celltype) <- epi_cells
print(table(condition_celltype))

# 5. 정상 baseline = 정상 조직(IMHS78w-10) epithelial cells
#    (inferCNV의 ref_group_names ↔ copykat의 norm.cell.names)
reference_groups <- c(
  "IMHS78w-10_ciliated cell",
  "IMHS78w-10_pulmonary alveolar type 2 cell",
  "IMHS78w-10_club cell"
)
norm_cells <- names(condition_celltype)[condition_celltype %in% reference_groups]
if (length(norm_cells) == 0) norm_cells <- ""   # 없으면 copykat 자동 baseline

cat("norm.cell.names:", if (identical(norm_cells, "")) 0 else length(norm_cells), "cells\n")
print(reference_groups)
print(unique(condition_celltype))

# 6. id.type 결정 (rownames가 rat gene symbol이면 "S", Ensembl이면 "E")
anno <- get("full.anno.rn7", asNamespace("copykat"))
n_sym <- length(intersect(rownames(counts), anno$mgi_symbol))
n_ens <- length(intersect(rownames(counts), anno$ensembl_gene_id))
id_type <- if (n_ens > n_sym) "E" else "S"
cat("id.type =", id_type, " (symbol", n_sym, "/ ensembl", n_ens, ")\n")

# 7. CopyKAT 실행 (rn7 = rat mRatBN7.2 gene-space 경로; full.anno.rn7 이미지에 baked)
setwd(outdir)
options(bitmapType = "cairo", scipen = 100)

ck <- copykat(
  rawmat          = counts,
  id.type         = id_type,
  genome          = "rn7",
  norm.cell.names = norm_cells,   # 정상 baseline (""이면 자동 추정)
  sam.name        = "HN00283643_epi",
  ngene.chr       = 5,
  KS.cut          = 0.1,
  win.size        = 25,
  n.cores         = 8,
  output.seg      = "FALSE"       # "TRUE" 시 IGV용 .seg (start_position 사용)
)

# 8. 결과 저장 (리터럴 aneuploid/diploid/not.defined 보존 → review_required)
saveRDS(ck, file.path(outdir, "copykat_result.rds"))
print(table(ck$prediction$copykat.pred))
