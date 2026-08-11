#!/usr/bin/env bash
# Build the CopyKAT rat (rn7) image with SingularityCE (fakeroot / userns).
# Run this on a host with working Singularity + network (the build pulls the
# rocker/r-ver:4.2.2 base and CRAN packages once). Runtime is fully offline.
#
# Usage:
#   ./build.sh [OUTPUT.sif]
# Default output: ./copykat-rn7.sif
set -Eeuo pipefail
IFS=$'\n\t'

project_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
definition="$project_dir/copykat-rn7.def"
image="${1:-$project_dir/copykat-rn7.sif}"
build_log="$project_dir/build.log"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

command -v singularity >/dev/null 2>&1 || die "SingularityCE not on PATH"
test -f "$definition" || die "definition not found: $definition"
test -d "$project_dir/copykat-rn7" || die "fork source tree missing: copykat-rn7/"
test -f "$project_dir/copykat-rn7/R/sysdata.rda" || \
    die "full.anno.rn7 not baked; run build_full_anno_rn7.R then regenerate the fork"
test -f "$project_dir/tests/t2_synthetic_rn7.R" || die "tests/t2_synthetic_rn7.R missing"
test "$(uname -m)" = x86_64 || die "locked to linux x86_64"

# fakeroot preflight (matches the Rocky9 kit convention).
if singularity config fakeroot --list >/dev/null 2>&1; then
    :
elif command -v unshare >/dev/null 2>&1 && unshare -Ur true >/dev/null 2>&1; then
    :
else
    die "fakeroot preflight failed; configure /etc/subuid & /etc/subgid, or build as root"
fi

printf 'building %s from %s\n' "$image" "$definition" | tee "$build_log"
( cd "$project_dir" && singularity build --fakeroot "$image" "$definition" ) 2>&1 | tee -a "$build_log"

test -s "$image" || die "build produced no SIF"
chmod 0444 "$image" 2>/dev/null || true
sha256sum "$image" | tee "${image}.sha256"
printf 'Build OK: %s\n' "$image"
printf 'Next: ./validate.sh %s\n' "$image"
