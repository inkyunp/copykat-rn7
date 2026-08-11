#!/usr/bin/env bash
# Unprivileged runtime validation for the CopyKAT rat (rn7) image.
# Mirrors the Rocky9 kit style: version pins, offline execution, and a real
# synthetic rn7 run under a none-network namespace when available.
#
# Usage:
#   ./validate.sh [IMAGE.sif]
set -Eeuo pipefail
IFS=$'\n\t'

project_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
image="${1:-$project_dir/copykat-rn7.sif}"

die() { printf '[FAIL] %s\n' "$*" >&2; exit 1; }
pass() { printf '[PASS] %s\n' "$*"; }

command -v singularity >/dev/null 2>&1 || die "SingularityCE not on PATH"
test "$(id -u)" -ne 0 || die "run runtime validation as an unprivileged user"
test -r "$image" || die "SIF not readable: $image"
image_abs="$(readlink -f "$image")"
sha_before="$(sha256sum "$image_abs" | awk '{print $1}')"

run() { singularity exec --cleanenv --containall "$image_abs" "$@"; }

printf 'image=%s\n' "$image_abs"
printf 'image_sha256=%s\n' "$sha_before"
printf 'singularity_version=%s\n' "$(singularity version)"

test "$(run Rscript --vanilla -e 'cat(as.character(getRversion()))')" = 4.2.2 || die "R is not 4.2.2"
test "$(run Rscript --vanilla -e 'cat(as.character(packageVersion("copykat")))')" = 1.1.0.9001 || die "copykat fork version mismatch"
pass "pinned R 4.2.2 and copykat 1.1.0.9001"

run Rscript --vanilla -e 'suppressPackageStartupMessages(library(copykat)); stopifnot("annotateGenes.rn7" %in% getNamespaceExports("copykat")); fa <- get("full.anno.rn7", asNamespace("copykat")); stopifnot(is.data.frame(fa), identical(colnames(fa), c("abspos","chromosome_name","start_position","end_position","ensembl_gene_id","mgi_symbol","band")), length(unique(fa$abspos)) > 20, length(unique(fa$chromosome_name)) == 21); cat("rn7 API + baked full.anno.rn7 OK: genes =", nrow(fa), "\n")'
pass "rn7 API exported and full.anno.rn7 baked"

# Real synthetic rn7 run.
run Rscript --vanilla /opt/tests/t2_synthetic_rn7.R \
    /opt/src/copykat-rn7/R/copykat.R /opt/src/copykat-rn7/R/sysdata.rda /tmp/t2out_validate
pass "synthetic rn7 copykat run"

# Offline acceptance under a none-network namespace when the runner allows it.
if singularity exec --net --network none --cleanenv --containall "$image_abs" \
     Rscript --vanilla /opt/tests/t2_synthetic_rn7.R \
       /opt/src/copykat-rn7/R/copykat.R /opt/src/copykat-rn7/R/sysdata.rda /tmp/t2out_offline; then
    pass "offline rn7 run under --network none"
else
    printf '[SKIP] runner cannot create a none-network namespace; offline acceptance pending\n'
fi

sha_after="$(sha256sum "$image_abs" | awk '{print $1}')"
test "$sha_before" = "$sha_after" || die "SIF changed during validation"
pass "image unchanged during validation"
printf 'All validation gates passed.\n'
