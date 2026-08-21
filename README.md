# CopyKAT multi-species fork (hg20 / mm10 / rn7)

A minimal runtime fork of [navinlabcode/copykat](https://github.com/navinlabcode/copykat)
(1.1.0, commit `b795ff7`) that adds a rat (*Rattus norvegicus*, **mRatBN7.2**)
gene-space path via `genome="rn7"`, so the same package runs on **human, mouse, and rat**.

> Results are **candidate CNV evidence**, not a malignant/normal diagnosis
> (`review_required`). rat/mouse share the gene-space path; cross-check with
> DNA-level evidence before use.

## What changed vs upstream

The package keeps the **exact upstream b795ff7 file layout** (per-function R
scripts, `data/`, `man/`). Every original file is byte-for-byte identical; only
`copykat.R` gains the rn7 branch, plus one new file and one new data object:

1. `copykat()` annotation branch: added `else if(genome=="rn7") annotateGenes.rn7(...)`.
2. gene-space block: `if(genome=="mm10")` → `if(genome=="mm10" || genome=="rn7")`.
3. new `R/annotateGenes.rn7.R` (copy of `annotateGenes.mm10()`, default `full.anno=full.anno.rn7`) + `man/annotateGenes.rn7.Rd` + NAMESPACE export.
4. `full.anno.rn7` (mRatBN7.2) added to `data/sysdata.rda` — no GTF needed at runtime.

hg20-only steps (HLA/cell-cycle removal, 220kb bins) are untouched.
All other functions (hg20 human, mm10 mouse, MCMC, baselines, heatmap) are the
unmodified upstream originals.

## Install

```bash
R CMD INSTALL copykat-rn7
```

Runtime deps: `parallelDist`, `dlm`, `gplots`, `RColorBrewer`, `mixtools`, `cluster`, `MCMCpack`.

## Input contract

**Raw integer UMI counts, genes × cells.** Not log-normalized `data`, not
`scale.data`, not SCT/integrated; do not subset to variable features.
`norm.cell.names` (optional): normal-cell barcodes matching the matrix columns.

## Usage

```r
library(copykat)
res <- copykat(rawmat = counts, id.type = "S", genome = "hg20",  # or mm10 / rn7
               norm.cell.names = normal_barcodes, sam.name = "sample", n.cores = 1)
```

Or via the runner (`--genome`, default `rn7`):

```bash
Rscript run_copykat_rat.R --input in.rds --output out \
  --genome hg20 --id-type auto --cores 1     # hg20 | mm10 | rn7
```

| species | `--genome` | `id.type="S"` symbols | annotation |
| --- | --- | --- | --- |
| human | `hg20` | HGNC (`hgnc_symbol`) | `full.anno` |
| mouse | `mm10` | MGI (`mgi_symbol`) | `full.anno.mm10` |
| rat | `rn7` | rat (`mgi_symbol`) | `full.anno.rn7` |

Output: `<output>.rds`, `<output>.cnv_evidence.tsv`, `<output>_<genome>_run/`.

## Files

| File | Purpose |
| --- | --- |
| `copykat-rn7/` | Fork package (upstream b795ff7 layout): `R/*.R`, `data/sysdata.rda` (three species), `man/`. Only `copykat.R` edited + `annotateGenes.rn7.R` added |
| `run_copykat_rat.R` | Runner (RNA counts → `copykat(genome=…)` → candidate evidence) |
| `run_copykat.R` | Base-R runner (direct Seurat slot access, no Seurat dependency) |

## License

Inherits upstream [copykat](https://github.com/navinlabcode/copykat)'s **GPL-2**;
derivative work.
