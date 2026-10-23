#!/usr/bin/env bash
# Root owns submission and serial scheduling; this pilot cannot pass A18/A19/A20.
#SBATCH --job-name=bethe_l32_pilot
#SBATCH --account=acf-utk0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3500M
#SBATCH --time=01:00:00
set -euo pipefail
: "${VALIDATION_ROOT:?Set VALIDATION_ROOT}"
: "${BETHE_TRANCHE:?Set absolute staged tranche directory}"
module load anaconda3/2024.06
cd "$BETHE_TRANCHE"
sha256sum -c SHA256SUMS
python3 - "$BETHE_TRANCHE" <<'PY'
from pathlib import Path
import hashlib,json,sys,tomllib
root=Path(sys.argv[1]);cfg=tomllib.loads((root/'pilot.toml').read_text());p=cfg['engineering_prerequisite']
audit=root/p['audit_relative'];anchor=root/p['anchor_relative']
assert hashlib.sha256(audit.read_bytes()).hexdigest()==p['audit_sha256']
assert hashlib.sha256(anchor.read_bytes()).hexdigest()==p['preexecution_engineering_anchor_sha256']
report=json.loads(audit.read_text());gate=report['next_L32_engineering_eligibility']
assert report['status']=='PASS' and report['scientific_completions']==0
assert gate['job_elapsed_seconds']==p['L24_job_elapsed_seconds']<=p['maximum_L24_wall_seconds']
assert gate['peak_rss_kib']==p['L24_peak_rss_kib']<=p['maximum_L24_peak_rss_MiB']*1024
assert gate['wall_at_most_30_minutes'] and gate['peak_at_most_2300_MiB']
assert gate['preexecution_anchor']['commit']==p['preexecution_engineering_anchor_commit']
assert cfg['L']==[32] and cfg['scientific_completion_allowed'] is False
PY
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="$VALIDATION_ROOT/depot"
export RENYI_ROOT="$VALIDATION_ROOT/source/renyi_disjoint_corrected"
BETHE_OUTPUT="$VALIDATION_ROOT/results/bethe_l32_pilot_${SLURM_JOB_ID}"
mkdir "$BETHE_OUTPUT"
/usr/bin/time -v -o "$BETHE_OUTPUT/oracle.resources.txt" \
 "$VALIDATION_ROOT/quspin-venv/bin/python" "$BETHE_TRANCHE/oracle_lieb_wu.py" \
 --manifest "$BETHE_TRANCHE/pilot.toml" --output "$BETHE_OUTPUT/oracle.json"
/usr/bin/time -v -o "$BETHE_OUTPUT/target.resources.txt" timeout 3300 \
 "$VALIDATION_ROOT/runtime/julia-1.12.1/bin/julia" --project="$RENYI_ROOT" \
 "$BETHE_TRANCHE/pilot_target.jl" --manifest "$BETHE_TRANCHE/pilot.toml" --output "$BETHE_OUTPUT/target.toml"
"$VALIDATION_ROOT/quspin-venv/bin/python" "$BETHE_TRANCHE/summarize_pilot.py" \
 --manifest "$BETHE_TRANCHE/pilot.toml" --oracle "$BETHE_OUTPUT/oracle.json" \
 --target "$BETHE_OUTPUT/target.toml" --output "$BETHE_OUTPUT/summary.json"
