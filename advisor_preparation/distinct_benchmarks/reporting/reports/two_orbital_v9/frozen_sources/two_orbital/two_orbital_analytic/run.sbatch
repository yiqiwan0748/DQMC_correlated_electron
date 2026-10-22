#!/bin/bash
#SBATCH --job-name=orbital-analytic
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:15:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${TRANCHE_ROOT:?Set TRANCHE_ROOT to this frozen release}"
module load anaconda3/2024.06
ANALYTIC_RUN_MODE=${ANALYTIC_RUN_MODE:-pilot}
case "$ANALYTIC_RUN_MODE" in pilot|full) ;; *) exit 2;; esac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
cd "$TRANCHE_ROOT"
sha256sum -c SHA256SUMS
OUT="$VALIDATION_ROOT/results/two_orbital_analytic_${ANALYTIC_RUN_MODE}_${SLURM_JOB_ID}"
mkdir "$OUT"
export ANALYTIC_TWO_ORBITAL_OUTPUT="$OUT/target.toml"
export ANALYTIC_TWO_ORBITAL_MODE="$ANALYTIC_RUN_MODE"
date -u +'%FT%TZ' > "$OUT/start_utc.txt"
if [[ "$ANALYTIC_RUN_MODE" == full ]]; then
  /usr/bin/time -v -o "$OUT/oracle.resources.txt" timeout 180 \
    "$VALIDATION_ROOT/quspin-venv/bin/python" oracle.py --output "$OUT/oracle.json"
fi
/usr/bin/time -v -o "$OUT/target.resources.txt" timeout 600 \
  "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --startup-file=no target.jl
if [[ "$ANALYTIC_RUN_MODE" == full ]]; then
  "$VALIDATION_ROOT/quspin-venv/bin/python" compare.py --target "$OUT/target.toml" --oracle "$OUT/oracle.json" --output "$OUT/comparison.json"
fi
date -u +'%FT%TZ' > "$OUT/end_utc.txt"
cd "$OUT"
sha256sum ./* > receipts.sha256
