# CopyKAT multi-species fork (hg20 / mm10 / rn7)

A minimal fork of [navinlabcode/copykat](https://github.com/navinlabcode/copykat)
(1.1.0, commit `b795ff7`) that adds a rat (*Rattus norvegicus*, **mRatBN7.2**)
gene-space path via `genome="rn7"`, so the same tool runs on **human, mouse, and rat**.

> Results are **candidate CNV evidence**, not a malignant/normal diagnosis
> (`review_required`). rat/mouse share the gene-space path; cross-check with
> DNA-level evidence before use.

## What changed vs upstream

3 source edits + 1 data object (`full.anno.rn7`, baked into `R/sysdata.rda`).
See [`docs/MODIFICATIONS.md`](docs/MODIFICATIONS.md) and [`docs/copykat_rn7.patch`](docs/copykat_rn7.patch).

1. `copykat()` annotation branch: added `else if(genome=="rn7") annotateGenes.rn7(...)`.
2. gene-space block: `if(genome=="mm10")` → `if(genome=="mm10" || genome=="rn7")`.
3. new `annotateGenes.rn7()` (copy of `annotateGenes.mm10()`, default `full.anno=full.anno.rn7`).

hg20-only steps (HLA/cell-cycle removal, 220kb bins) are untouched.

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
singularity exec --cleanenv --containall --bind "$PWD":/work copykat-rn7.sif \
  Rscript /work/run_copykat_rat.R --input /work/in.rds --output /work/out \
    --genome hg20 --id-type auto --cores 1     # hg20 | mm10 | rn7
```

| species | `--genome` | `id.type="S"` symbols | annotation |
| --- | --- | --- | --- |
| human | `hg20` | HGNC (`hgnc_symbol`) | `full.anno` |
| mouse | `mm10` | MGI (`mgi_symbol`) | `full.anno.mm10` |
| rat | `rn7` | rat (`mgi_symbol`) | `full.anno.rn7` |

Output: `<output>.rds`, `<output>.cnv_evidence.tsv`, `<output>_<genome>_run/`.

## Build

```bash
./build.sh                    # -> copykat-rn7.sif (needs network once)
./validate.sh copykat-rn7.sif # offline validation
```

## Files

| File | Purpose |
| --- | --- |
| `copykat-rn7/` | Fork source (baked `full.anno.rn7`) |
| `copykat-rn7.def` / `build.sh` / `validate.sh` | Image build + validation |
| `build_full_anno_rn7.R` | Build `full.anno.rn7` from `../genes.gtf` (mRatBN7.2) |
| `validate_rat_input.R` | Input-contract gate |
| `run_copykat_rat.R` | Runner (RNA counts → `copykat(genome=…)` → candidate evidence) |
| `run_copykat.R` | Base-R runner (no Seurat in image) |
| `tests/` | Synthetic smoke tests |

## Validation

All three genomes run end-to-end in-image on a synthetic smoke test
(prediction/CNAmat/hclustering returned, gene-space `CNA_results.txt` written,
literal `aneuploid`/`diploid` preserved). Results stay `review_required`.

## License

Inherits upstream [copykat](https://github.com/navinlabcode/copykat)'s **GPL-2**;
derivative work.
