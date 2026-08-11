# copykat (rat rn7 fork)

A minimal fork of [navinlabcode/copykat](https://github.com/navinlabcode/copykat)
**1.1.0** (upstream commit `b795ff793522499f814f6ae282aad1aab790902f`) that adds a
rat (*Rattus norvegicus*, **mRatBN7.2**) gene-space code path via `genome="rn7"`,
reusing the mouse (mm10) downstream branch with a rat annotation table
(`full.anno.rn7`) baked into `R/sysdata.rda`.

## What changed vs upstream

Exactly three source edits plus one data object — see [`../docs/MODIFICATIONS.md`](../docs/MODIFICATIONS.md)
and [`../docs/copykat_rn7.patch`](../docs/copykat_rn7.patch):

1. `copykat()` annotation branch: added `else if(genome=="rn7") annotateGenes.rn7(...)`.
2. `copykat()` gene-space block: `if(genome=="mm10")` -> `if(genome=="mm10" || genome=="rn7")`.
3. New `annotateGenes.rn7()` (copy of `annotateGenes.mm10()`, default `full.anno=full.anno.rn7`).
4. `R/sysdata.rda`: added `full.anno.rn7` (mRatBN7.2, 25,302 genes, 21 chromosome labels).

The hg20-only steps (cell-cycle/HLA removal, 220kb bin conversion) are untouched;
they are already gated to `genome=="hg20"`, so rat inherits the gene-space mouse path.

## Usage

```r
library(copykat)
# rawmat: raw integer UMI counts, genes x cells (NOT normalized/scaled)
res <- copykat(rawmat = counts, id.type = "S", genome = "rn7",
               norm.cell.names = normal_barcodes, sam.name = "sample",
               n.cores = 1)
```

Gene identifiers must match `full.anno.rn7`: rat gene symbols (`id.type="S"`, matched
against the `mgi_symbol` column which holds rat symbols) or Ensembl IDs (`id.type="E"`).

> Results are **candidate CNV evidence only**, not a malignant/normal diagnosis, and
> run through a mouse code path the upstream author did not extensively validate.
> Cross-check with DNA-level evidence (WGS/SNP-array/karyotype) before use.

## License

Inherits upstream copykat's **GPL-2** license (see `DESCRIPTION`). This is a
derivative work of navinlabcode/copykat.
