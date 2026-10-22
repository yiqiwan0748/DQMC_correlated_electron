#!/usr/bin/env bash
set -euo pipefail
: "${SCRATCHDIR:?ISAAC scratch is unset}"
validation_root="$SCRATCHDIR/advisor-validation-20261022"
mkdir -p "$validation_root/runtime" "$validation_root/logs"
cd "$validation_root/runtime"
archive=julia-1.12.1-linux-x86_64.tar.gz
curl --fail --location --retry 2 --max-time 600 --output "$archive" "https://julialang-s3.julialang.org/bin/linux/x64/1.12/$archive"
curl --fail --location --retry 2 --max-time 60 --output julia-1.12.1.sha256 https://julialang-s3.julialang.org/bin/checksums/julia-1.12.1.sha256
awk '$2=="julia-1.12.1-linux-x86_64.tar.gz" {print}' julia-1.12.1.sha256 > selected.sha256
test -s selected.sha256
sha256sum -c selected.sha256
test ! -e julia-1.12.1
tar -xzf "$archive"
./julia-1.12.1/bin/julia --version
