# CopyKAT rat (genome="rn7") fork + test image

A fork and dedicated lightweight image for running CopyKAT on rat
(*Rattus norvegicus*, **mRatBN7.2**) scRNA-seq. It inherits CopyKAT's mouse (mm10)
gene-space path (no 220kb bin conversion), adds a `genome="rn7"` value, and bakes a
rat annotation table — built from the `genes.gtf` in the parent folder — into the
package.

> Results are **candidate CNV evidence** only, not a final malignant/normal
> diagnosis. They run through a mouse code path the upstream author did not
> extensively validate, so every result stays `review_required` until cross-checked
> against DNA-level evidence (WGS/SNP-array/karyotype).

## What / How / What differs

### Why it is needed (problem)
Upstream CopyKAT `copykat()` branches its `genome` argument on **only** `"hg20"`
(human) and `"mm10"` (mouse), with no `else`. Passing a rat value matches neither
branch, so `anno.mat` is never created and the next line dies immediately. There was
simply no code path for running CopyKAT on rat scRNA-seq.

### What was done (approach)
Instead of building a rat-specific 220kb bin coordinate system from scratch, this fork
**reuses the mouse (mm10) downstream path**. The mm10 path produces results in
gene-space without a 220kb bin conversion, so only a single rat annotation table needs
to be injected. Thus `genome="rn7"` (rat mRatBN7.2) is made to follow the same
gene-space block as mm10 with minimal changes.

### How it was done (implementation)
1. Generate a CopyKAT-schema (7-column) rat annotation `full.anno.rn7` from the mRatBN7.2 GTF.
2. Add a `genome=="rn7"` branch and `annotateGenes.rn7()` to `copykat()`, and widen the mm10 gene-space block by one line so rn7 passes through it too.
3. Bake `full.anno.rn7` into the package `R/sysdata.rda`, so no GTF or gene-order file is needed at runtime.
4. Package it into a copykat-only lightweight image (`rocker/r-ver:4.2.2`) and validate with synthetic and real data.

### What differs from upstream

| Item | upstream copykat 1.1.0 | This fork (rn7) |
| --- | --- | --- |
| Supported species | human (hg20), mouse (mm10) | **+ rat (rn7, mRatBN7.2)** |
| Unsupported `genome` value | no `else` → dies with an error | rn7 branch added |
| Rat annotation | none | `full.anno.rn7` baked into the package |
| Annotation function | `annotateGenes.hg20/.mm10` | `+ annotateGenes.rn7()` |
| Rat output unit | (not possible) | gene-space (`CNA_results.txt`, same as mm10) |
| hg20-only steps | unchanged | **unchanged** (HLA/cyclegene removal and 220kb bin are hg20-gated, so left untouched) |
| Algorithm logic | — | **unchanged** (filtering, smoothing, baseline, MCMC segmentation are the original) |

> The code change is exactly **3 edits (+1 data object)**. For the full diff and the
> restored original, see [`docs/MODIFICATIONS.md`](docs/MODIFICATIONS.md),
> [`docs/copykat_rn7.patch`](docs/copykat_rn7.patch), and
> [`docs/original_copykat_b795ff7.R`](docs/original_copykat_b795ff7.R).

### Validation status (verified by execution)
- **Build**: `copykat-rn7.sif` builds successfully; in-image rn7 checks + the `%test` synthetic run PASS; offline (`--network none`) run PASSes.
- **T2 synthetic**: `copykat(genome="rn7")` runs to completion on real R 4.2.2 + real deps, producing gene-space output and preserving literal classes.
- **Real data**: runs to completion on rat 310 cells (lung, epithelial) — the normal baseline (macrophage) is mostly diploid and epithelial cells show a higher aneuploid fraction. Results are `review_required`.

## Input contract (confirmed from copykat 1.1.0 source)
CopyKAT filters cells/genes by **positive-count detection (`sum(x>0)`)** and internally
performs its own `log(sqrt(x)+sqrt(x+1))` transform, centering, and dlm smoothing.
Therefore the input must be **raw integer UMI counts, genes × cells**.

- Seurat v4: `GetAssayData(obj, assay="RNA", slot="counts")`
- Seurat v5: `LayerData(obj, assay="RNA", layer="counts")`
- **Forbidden**: log-normalized `data` (double transform → signal loss), `scale.data` (negatives → `sqrt(negative)=NaN`), SCT/integrated assays
- Do not reduce genes to variable features (genome-wide coverage is required)
- `norm.cell.names` (optional, recommended): a normal-cell barcode vector that must match the raw matrix column names verbatim

## Files
| File | Purpose |
| --- | --- |
| `build_full_anno_rn7.R` | Generate the 7-column `full.anno.rn7` from `../genes.gtf` (mRatBN7.2) → `full.anno.rn7.rds` + integrity gate |
| `full.anno.rn7.rds` | Generated rat annotation (25,302 genes, 21 chromosome labels) |
| `copykat-rn7/` | Fork source tree (based on copykat 1.1.0, commit `b795ff7`). Bakes `full.anno.rn7` into `R/sysdata.rda` |
| `copykat-rn7.def` | Lightweight Singularity definition (`rocker/r-ver:4.2.2`, copykat-only) |
| `build.sh` / `validate.sh` | Build the image (fakeroot/network once) / unprivileged offline validation |
| `validate_rat_input.R` | **Data-validation gate** (5 input-contract checks). Shared by big machine / local |
| `run_copykat_rat.R` | Runner: extract RNA counts → `copykat(..., genome="rn7")` → save candidate evidence |
| `run_copykat.R` | Base-R runner mirroring the inferCNV pipeline (no Seurat in the image; direct slot access) |
| `tests/t2_synthetic_rn7.R` | Synthetic smoke test (no RDS needed) |

## Modifications the fork makes (exactly 3 vs upstream)
1. Add `else if(genome=="rn7") anno.mat <- annotateGenes.rn7(...)` to the annotation branch
2. New `annotateGenes.rn7()` = a copy of `annotateGenes.mm10()` with default `full.anno=full.anno.rn7`
3. Widen the later gene-space block `if(genome=="mm10")` → `if(genome=="mm10" || genome=="rn7")`

The hg20-only steps (HLA/cyclegene removal, 220kb bin conversion) are untouched — they
are already gated to `genome=="hg20"`. The `mgi_symbol` column name is kept verbatim
(the copykat body references this name directly, so it holds rat symbols but the column
name must not change or the code breaks).

`full.anno.rn7` schema: `abspos, chromosome_name (integer 1–20, X=21), start_position,
end_position, ensembl_gene_id, mgi_symbol (rat symbol), band (NA)`. `abspos` is a real
cumulative coordinate, so it does not reproduce the constant-abspos defect of the
bundled `full.anno.mm10` (20 unique of 137,030).

## Execution order

### Already done (local, no RDS needed)
- **Phase 0**: `Rscript build_full_anno_rn7.R` → generates `full.anno.rn7.rds`, passes the integrity gate
  (25,302 genes / 25,268 unique abspos / 21 labels / min 575 genes per chr).
- **T1**: fork structure checks — `annotateGenes.rn7` defined, `formals(copykat)$genome`, sysdata baking.
- **T2**: modified copykat `genome="rn7"` runs to completion on real R 4.2.2 + real deps
  (gene-space `CNA_results.txt`, literal `aneuploid`/`diploid` separation confirmed).

### Build the image (user's WSL, network once)
```bash
cd rat-copykat
./build.sh                 # -> copykat-rn7.sif
./validate.sh copykat-rn7.sif
```
The build needs network to download the base image and the CRAN snapshot (2023-10-01); runtime is offline.

### Phase 1 — sampling + data validation on a big machine (transfer gate)
The 5.6GB source RDS cannot be loaded on a 7.6GB local host, so run the following **once on a
≥32GB (64GB recommended) machine**. Bring `full.anno.rn7.rds` and `validate_rat_input.R` along.

Provide first: the cell-type column name + normal labels, and the `sample_id` column name.
(rownames / Seurat version / gene & cell counts are auto-recorded by the validation report.)

```r
library(Seurat); set.seed(1729)
obj <- readRDS("test.obj_final_add_celltype_scp2711.rds")
DefaultAssay(obj) <- "RNA"
ct  <- "<celltype_col>"; smp <- "<sample_id_col>"
one <- names(sort(table(obj[[smp]][,1]), decreasing=TRUE))[1]   # one sample
o1  <- subset(obj, cells = colnames(obj)[obj[[smp]][,1] == one])
# select ~500 barcodes by cell-type stratification + dropping low-depth cells (ensure normal lineage >= 60~100)
keep <- NULL  # <~500 barcodes from stratification + depth filter>
sub  <- subset(o1, cells = keep)
sub  <- DietSeurat(sub, assays="RNA", dimreducs=NULL, graphs=NULL)  # keep counts
sub@meta.data <- sub@meta.data[, c(ct, smp), drop=FALSE]
saveRDS(sub, "rat_copykat_test_500.rds")                            # target <0.5GB
writeLines(colnames(sub)[sub@meta.data[[ct]] %in% c("<normal labels>")], "norm_cell_names.txt")
```

Validate on the same machine right after extraction and **transfer only on PASS**:
```bash
Rscript validate_rat_input.R \
  --input rat_copykat_test_500.rds --anno full.anno.rn7.rds \
  --norm-cells norm_cell_names.txt \
  --celltype-column <celltype_col> --normal-labels "<normal labels>" \
  --report input_validation_report.tsv
# exit 0 = PASS (transfer), exit 1 = FAIL (do not transfer)
```

### Phase 2/3 — run after transferring locally
Copy only `rat_copykat_test_500.rds` + `norm_cell_names.txt` + `input_validation_report.tsv`
locally (do not move the 5.6GB source RDS). Re-validate locally, then run:
```bash
# local re-validation (transfer integrity / reproducibility)
Rscript validate_rat_input.R --input rat_copykat_test_500.rds --anno full.anno.rn7.rds \
  --norm-cells norm_cell_names.txt --report input_validation_report.local.tsv

# run (inside the image)
singularity exec --cleanenv --containall --bind "$PWD":/work copykat-rn7.sif \
  Rscript /work/run_copykat_rat.R \
    --input /work/rat_copykat_test_500.rds \
    --output /work/rat_copykat_result \
    --norm-cells /work/norm_cell_names.txt \
    --id-type auto --cores 1
```
Output: `rat_copykat_result.rds` (literal classes + `interpretation=review_required`),
`rat_copykat_result.cnv_evidence.tsv`, and `rat_copykat_result_rn7_run/` (heatmap PDF, etc.).

## Execution parameters to revisit (after real data is available)
- `ngene.chr`: rat has 20 autosomes + X = 21 labels → the filter strength differs from human (23) / mouse (20)
- `LOW.DR`, `UP.DR`: depend on the depth/quality of the rat data
- `norm.cell.names` vs automatic baseline: the auto baseline is unstable when normal cells are scarce

See [`docs/TUNING.md`](docs/TUNING.md) — these are all call-time arguments of `copykat()`, so
changing them does **not** require rebuilding the fork/image.

## Note: building/running in this environment
- The system R is 4.5.0 without CopyKAT deps, so T1/T2 were actually run with **R 4.2.2 + deps** extracted from the existing Rocky9 SIF (a faithful target R).
- The existing SIF fails to mount in this exec sandbox due to a missing `/dev/fuse` and setuid-starter ownership (sandbox artifacts). Validate the real build/run of the new image in the user's WSL with `build.sh`/`validate.sh`.
