#!/usr/bin/env bash
# Submission/scheduling belongs to root. This excluded pilot cannot complete X01.
#SBATCH --job-name=x01_timing_pilot
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:10:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${X01_PILOT_TRANCHE:?Set absolute staged pilot directory}"
module load anaconda3/2024.06
cd "$X01_PILOT_TRANCHE"
sha256sum -c SHA256SUMS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_disjoint_corrected"
X01_OUTPUT_DIR="$VALIDATION_ROOT/results/continuum_bias_pilot_${SLURM_JOB_ID}"
mkdir "$X01_OUTPUT_DIR"
export X01_PILOT_OUTPUT="$X01_OUTPUT_DIR/target.toml"
set +e
/usr/bin/time -v -o "$X01_OUTPUT_DIR/oracle.resources.txt" timeout 90 \
 "$VALIDATION_ROOT/quspin-venv/bin/python" "$X01_PILOT_TRANCHE/oracle_quspin.py" \
 --manifest "$X01_PILOT_TRANCHE/manifest.toml" --output "$X01_OUTPUT_DIR/oracle.json"
X01_ORACLE_EXIT=$?
/usr/bin/time -v -o "$X01_OUTPUT_DIR/target.resources.txt" timeout 480 \
 "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
 "$X01_PILOT_TRANCHE/target_pilot.jl"
X01_TARGET_EXIT=$?
set -e
printf 'oracle_exit=%s\ntarget_exit=%s\n' "$X01_ORACLE_EXIT" "$X01_TARGET_EXIT" > "$X01_OUTPUT_DIR/process_exit_codes.txt"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$X01_PILOT_TRANCHE/summarize_pilot.py" \
 --manifest "$X01_PILOT_TRANCHE/manifest.toml" --target "$X01_OUTPUT_DIR/target.toml" \
 --oracle "$X01_OUTPUT_DIR/oracle.json" --target-exit-code "$X01_TARGET_EXIT" --output "$X01_OUTPUT_DIR/summary.json"
if (( X01_ORACLE_EXIT != 0 || X01_TARGET_EXIT != 0 )); then exit 1; fi
