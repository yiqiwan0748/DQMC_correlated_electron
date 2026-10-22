#!/usr/bin/env bash
# Root owns submission and serial scheduling; this pilot cannot pass A18/A19/A20.
#SBATCH --job-name=bethe_l24_pilot
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=00:45:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${BETHE_TRANCHE:?Set absolute staged tranche directory}"
module load anaconda3/2024.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_disjoint_corrected"
BETHE_OUTPUT="$VALIDATION_ROOT/results/bethe_l24_pilot_${SLURM_JOB_ID}"
mkdir "$BETHE_OUTPUT"
/usr/bin/time -v -o "$BETHE_OUTPUT/oracle.resources.txt" \
 "$VALIDATION_ROOT/quspin-venv/bin/python" "$BETHE_TRANCHE/oracle_lieb_wu.py" \
 --manifest "$BETHE_TRANCHE/pilot.toml" --output "$BETHE_OUTPUT/oracle.json"
/usr/bin/time -v -o "$BETHE_OUTPUT/target.resources.txt" timeout 2400 \
 "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
 "$BETHE_TRANCHE/pilot_target.jl" --manifest "$BETHE_TRANCHE/pilot.toml" --output "$BETHE_OUTPUT/target.toml"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$BETHE_TRANCHE/summarize_pilot.py" \
 --manifest "$BETHE_TRANCHE/pilot.toml" --oracle "$BETHE_OUTPUT/oracle.json" \
 --target "$BETHE_OUTPUT/target.toml" --output "$BETHE_OUTPUT/summary.json"
