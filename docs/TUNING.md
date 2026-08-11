# Parameter tuning — where and what

**Important: tuning is not a source-code change.** `ngene.chr`, `LOW.DR`, `UP.DR`,
`KS.cut`, and `win.size` are all **call-time arguments** of `copykat()` (in the fork
source these names appear only in the function signature and as variable references, not
as hardcoded constants — see `docs/MODIFICATIONS.md`). Changing their values therefore
**does not require rebuilding the fork/image.** Just pass them as flags to
`run_copykat_rat.R`.

## Where

In `run_copykat_rat.R` (run inside the image). Example:

```bash
singularity exec --cleanenv --containall --bind "$PWD":/work copykat-rn7.sif \
  Rscript /work/run_copykat_rat.R \
    --input /work/rat_copykat_test_500.rds \
    --output /work/rat_copykat_result \
    --norm-cells /work/norm_cell_names.txt \
    --ngene-chr 5 --low-dr 0.05 --up-dr 0.10 --ks-cut 0.10 --win-size 25
```

Re-run with the same image, changing only the values, as many times as needed.

## What (meaning and tuning direction per parameter)

| Flag | copykat arg | Default | What it changes | When to tune |
| --- | --- | --- | --- | --- |
| `--ngene-chr` | `ngene.chr` | 5 | minimum genes a cell must have on **every chromosome label (rat=21)** to pass | lower it (e.g. 3) if low-depth cells get filtered out en masse. Rat has a different label count than human (23)/mouse (20), so the filter strength differs |
| `--low-dr` | `LOW.DR` | 0.05 | minimum detection rate (fraction of cells) to keep a gene for smoothing/segmentation | lower it to keep more genes when the data is shallow. If `nrow(rawmat)<7000`, the code forces `UP.DR<-LOW.DR` |
| `--up-dr` | `UP.DR` | 0.10 | minimum detection rate to keep a gene in the final result | keep it larger than LOW.DR. Raise it if there is a lot of noise |
| `--ks-cut` | `KS.cut` | 0.10 | segmentation sensitivity (KS test cutoff) | **lower = finer segmentation** (more sensitive), higher = more conservative. Lower it if you see "too few breakpoints" |
| `--win-size` | `win.size` | 25 | minimum genes per segmentation window | lower it when genes are few. Too small increases noisy segments |

`--cores` (n.cores), `--id-type` (auto/S/E), and `--norm-cells` are also tunable.
`--id-type auto` picks whichever of symbol↔`mgi_symbol` or Ensembl↔`ensembl_gene_id` has
the larger intersection.

## Signal observed on this data (HN00283643, 310 cells) and tuning hints

- The log did not show `WARNING: low data quality; assigned LOW.DR to UP.DR` (8,845 genes
  passed LOW.DR), and `step 5: segmentation` proceeded normally. → the filter was not too aggressive.
- If the heatmap looks scattered rather than showing clear arm-level blocks, the signal is
  weak. Things to try, in order:
  1. lower `--ks-cut 0.05` to increase segmentation sensitivity
  2. increase the number of cells (≥500~1000) and re-run — the 310-cell pilot has a shallow baseline/clustering
  3. specify the normal reference more firmly with `--norm-cells` (currently 100 alveolar macrophage cells)
- `ngene.chr`/`LOW.DR`/`UP.DR` should be set **based on the actual data depth**. The defaults
  are not pinned as recommended values (they vary by dataset).

## When a source change *is* needed (reference — not now)

The following are code/data, not arguments, so changing them requires editing the fork +
rebuilding the image:

- `set.seed(1234)` (copykat's internal fixed seed)
- creating a rat-specific 220kb bin conversion (currently it finishes in gene-space like mm10 — a separate large effort)
- the chromosome scope of `full.anno.rn7` (currently 1–20+X; to include Y/MT, edit `build_full_anno_rn7.R` and rebuild)

Everything in scope for tuning here is an argument, so **no image rebuild is needed.**
