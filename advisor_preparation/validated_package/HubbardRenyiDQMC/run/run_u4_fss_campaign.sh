#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
JULIA_BIN="${JULIA_BIN:-/home/yqw/.juliaup/bin/julia}"

readonly -a SIZES=(4 6 8 10 12 14 16 18 20 22)
readonly -a ED_SIZES=(4 6 8 10 12)
readonly -a DMRG_SIZES=(12 14 16 18 20 22)
readonly -a DMRG_REPEAT_SIZES=(20 22)
readonly -a ROOT_SEEDS=(1234 11234 21234 31234)
readonly DMRG_LEVEL1_SCHEDULE="20,50,100,200,400,800,1200,1600,2400,3200,4800,4800,4800,4800,5200,5200,5200,5200"
readonly DMRG_LEVEL2_SCHEDULE="20,50,100,200,400,800,1200,1600,2400,3200,4800,4800,4800,4800,5200,5200,5200,5200,5600,5600,5600,5600"
readonly DMRG_S2_SPAN_LIMIT="1e-4"
readonly DMRG_ENERGY_DENSITY_SPAN_LIMIT="1e-8"
readonly DEFAULT_TASK_NAME="HubbardRenyiDQMC-U4-FSS"
readonly DEFAULT_ALLOWED_RUN_ROOT="/mnt/c/Users/YQW/HubbardRenyiDQMC_runs"
readonly EXPECTED_NODE_PAIRS=68
readonly EXPECTED_RAW_FILES=136

run_julia() {
    # Active compute inherits the campaign lock. If the owner shell is killed
    # but Julia survives, a recovery attempt cannot overlap the same outputs.
    "$JULIA_BIN" "$@"
}

usage() {
    cat <<'EOF'
Usage:
  bash run/run_u4_fss_campaign.sh RUN_DIR [MAX_JOBS] [TASK_NAME]
  bash run/run_u4_fss_campaign.sh --preflight
  bash run/run_u4_fss_campaign.sh --status RUN_DIR

RUN_DIR must be a persistent directory below:
  /mnt/c/Users/YQW/HubbardRenyiDQMC_runs/
EOF
}

parameter_path() {
    local size="$1"
    local beta=$((4 * size))
    local subsystem=$((size / 2))
    printf '%s/parameters/L%d_U4_beta%d_nA%d_Ninc16_burnin8000.toml' \
        "$REPO_ROOT" "$size" "$beta" "$subsystem"
}

run_name() {
    local size="$1"
    local beta=$((4 * size))
    local subsystem=$((size / 2))
    printf 'L%d_U4_beta%d_nA%d_Ninc16_burnin8000' \
        "$size" "$beta" "$subsystem"
}

preflight() {
    [[ -x "$JULIA_BIN" ]] || {
        echo "Julia executable not found: $JULIA_BIN" >&2
        return 1
    }
    command -v flock >/dev/null
    command -v sha256sum >/dev/null
    command -v realpath >/dev/null

    local parameter_args=()
    local size
    for size in "${SIZES[@]}"; do
        local parameter
        parameter="$(parameter_path "$size")"
        [[ -f "$parameter" ]] || {
            echo "Missing parameter file: $parameter" >&2
            return 1
        }
        parameter_args+=("$parameter")
    done

    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        -e '
            using TOML
            expected_sizes = collect(4:2:22)
            length(ARGS) == length(expected_sizes) || error("wrong parameter count")
            for (path, L) in zip(ARGS, expected_sizes)
                cfg = TOML.parsefile(path)
                model = cfg["model"]
                sampling = cfg["sampling"]
                grid = cfg["lambda_grid"]
                subsystem = cfg["subsystem"]
                @assert model["L"] == L
                @assert model["U"] == 4.0
                @assert model["beta"] == 4L
                @assert model["dtau"] == 0.05
                @assert model["n_stab"] == 10
                @assert model["N_sweeps"] == 1
                @assert model["boundary_condition"] == "pbc"
                @assert subsystem["nA"] == L ÷ 2
                @assert grid["Ninc"] == 16
                @assert sampling["N_burnin"] == 8000
                @assert sampling["N_measurements"] == 4000
                @assert sampling["N_bins"] == 40
                @assert sampling["root_seeds"] == [1234, 11234, 21234, 31234]
            end
            println("Validated U=4 PBC parameter matrix: L=4:2:22")
        ' "${parameter_args[@]}"

    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/analysis/compute_s2_ed.jl" --help >/dev/null
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/analysis/compute_s2_itensor_dmrg.jl" --help >/dev/null
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/analysis/analyze_u4_fss.jl" --help >/dev/null
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/run/run_grid_local.jl" --help >/dev/null
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/analysis/summarize_grid.jl" --help >/dev/null
    echo "Preflight PASS: no simulation was started."
}

status_campaign() {
    local requested="$1"
    local directory
    directory="$(realpath -m -- "$requested")"
    [[ -d "$directory" ]] || {
        echo "Campaign directory does not exist: $directory" >&2
        return 1
    }

    echo "Campaign: $directory"
    if [[ -f "$directory/state/status.txt" ]]; then
        printf 'Status: '
        cat "$directory/state/status.txt"
    else
        echo "Status: not started"
    fi
    if [[ -f "$directory/state/heartbeat.txt" ]]; then
        printf 'Heartbeat: '
        cat "$directory/state/heartbeat.txt"
    fi

    printf '%-4s %10s %10s %8s\n' L summaries bins complete
    local size
    for size in "${SIZES[@]}"; do
        local raw="$directory/dqmc/L$size/raw"
        local summaries=0
        local bins=0
        [[ -d "$raw" ]] && summaries="$(find "$raw" -maxdepth 1 -type f -name '*_summary.csv' | wc -l)"
        [[ -d "$raw" ]] && bins="$(find "$raw" -maxdepth 1 -type f -name '*_bins.csv' | wc -l)"
        local complete="no"
        [[ -f "$directory/dqmc/L$size/state/dqmc.done" ]] && complete="yes"
        printf '%-4d %10d %10d %8s\n' "$size" "$summaries" "$bins" "$complete"
    done

    [[ -f "$directory/DONE" ]] && echo "Overall: DONE" || echo "Overall: running or pending"
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    usage
    exit 0
elif [[ "${1:-}" == "--preflight" ]]; then
    [[ $# -eq 1 ]] || { usage >&2; exit 2; }
    preflight
    exit 0
elif [[ "${1:-}" == "--status" ]]; then
    [[ $# -eq 2 ]] || { usage >&2; exit 2; }
    status_campaign "$2"
    exit 0
elif [[ $# -lt 1 || $# -gt 3 ]]; then
    usage >&2
    exit 2
fi

RUN_DIR="$(realpath -m -- "$1")"
MAX_JOBS="${2:-4}"
TASK_NAME="${3:-$DEFAULT_TASK_NAME}"
ALLOWED_RUN_ROOT="$(realpath -m -- "${U4_FSS_ALLOWED_RUN_ROOT:-$DEFAULT_ALLOWED_RUN_ROOT}")"

[[ "$ALLOWED_RUN_ROOT" != "/" && "$RUN_DIR" == "$ALLOWED_RUN_ROOT"/* ]] || {
    echo "RUN_DIR must be below $ALLOWED_RUN_ROOT" >&2
    exit 2
}
[[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]] || {
    echo "MAX_JOBS must be a positive integer" >&2
    exit 2
}

LOG_DIR="$RUN_DIR/logs"
STATE_DIR="$RUN_DIR/state"
REFERENCE_DIR="$RUN_DIR/references"
DQMC_DIR="$RUN_DIR/dqmc"
PROCESSED_DIR="$RUN_DIR/processed"
QUARANTINE_DIR="$RUN_DIR/quarantine"
mkdir -p "$LOG_DIR" "$STATE_DIR" "$REFERENCE_DIR/ed" \
    "$REFERENCE_DIR/dmrg" "$DQMC_DIR" "$PROCESSED_DIR" "$QUARANTINE_DIR"

exec 9>"$STATE_DIR/campaign.lock"
if ! flock -n 9; then
    echo "Campaign is already active: $RUN_DIR"
    exit 0
fi

ATTEMPT_ID="$(date -u +%Y%m%dT%H%M%SZ)_pid$$"
ATTEMPT_LOG="$LOG_DIR/campaign_attempt_$ATTEMPT_ID.log"
exec > >(tee -a "$ATTEMPT_LOG" "$LOG_DIR/campaign_master.log" 9>&-) 2>&1

write_status() {
    local message="$1"
    local temp="$STATE_DIR/status.txt.tmp.$$"
    printf '%s\n' "$message" > "$temp"
    mv -- "$temp" "$STATE_DIR/status.txt"
}

disable_task() {
    if [[ -x /mnt/c/Windows/System32/schtasks.exe ]]; then
        /mnt/c/Windows/System32/schtasks.exe /Change /TN "$TASK_NAME" /Disable \
            >/dev/null 2>&1 || true
    fi
}

heartbeat_loop() {
    local parent_pid="$1"
    while kill -0 "$parent_pid" >/dev/null 2>&1; do
        local temp="$STATE_DIR/heartbeat.txt.tmp.$$"
        printf '%s attempt=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ATTEMPT_ID" > "$temp"
        mv -- "$temp" "$STATE_DIR/heartbeat.txt"
        sleep 60
    done
}

heartbeat_loop "$$" 9>&- >/dev/null 2>&1 &
HEARTBEAT_PID=$!

cleanup() {
    kill "$HEARTBEAT_PID" >/dev/null 2>&1 || true
    wait "$HEARTBEAT_PID" >/dev/null 2>&1 || true
}
trap cleanup EXIT

record_failure() {
    local code="$1"
    trap - ERR
    local count=0
    [[ -f "$STATE_DIR/consecutive_failures.txt" ]] && \
        read -r count < "$STATE_DIR/consecutive_failures.txt"
    [[ "$count" =~ ^[0-9]+$ ]] || count=0
    count=$((count + 1))
    printf '%d\n' "$count" > "$STATE_DIR/consecutive_failures.txt.tmp.$$"
    mv -- "$STATE_DIR/consecutive_failures.txt.tmp.$$" \
        "$STATE_DIR/consecutive_failures.txt"
    write_status "FAILED attempt=$ATTEMPT_ID exit=$code consecutive_failures=$count"
    echo "Campaign attempt failed with exit code $code (consecutive failures: $count)."
    if (( count >= 3 )); then
        printf '%s repeated deterministic failure; inspect logs before resuming\n' \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$RUN_DIR/STOPPED"
        disable_task
        echo "STOPPED after three consecutive handled failures."
    fi
    exit "$code"
}
trap 'record_failure $?' ERR

fatal() {
    echo "$1" >&2
    record_failure 1
}

mark_progress() {
    printf '0\n' > "$STATE_DIR/consecutive_failures.txt.tmp.$$"
    mv -- "$STATE_DIR/consecutive_failures.txt.tmp.$$" \
        "$STATE_DIR/consecutive_failures.txt"
}

verify_overall_done() {
    local manifest="$STATE_DIR/final_artifacts.sha256"
    [[ -f "$RUN_DIR/DONE" && -f "$manifest" ]] ||
        fatal "DONE state is incomplete"
    (cd "$RUN_DIR" && sha256sum --check "$manifest" >/dev/null)
    local recorded_manifest_hash actual_manifest_hash
    recorded_manifest_hash="$(awk -F= '$1 == "final_manifest_sha256" {print $2}' "$RUN_DIR/DONE")"
    actual_manifest_hash="$(sha256sum "$manifest" | awk '{print $1}')"
    [[ -n "$recorded_manifest_hash" && "$recorded_manifest_hash" == "$actual_manifest_hash" ]] ||
        fatal "DONE manifest hash mismatch"
    grep -qx 'PASS: U=4 PBC FSS L=4:2:22 complete' "$RUN_DIR/DONE" ||
        fatal "DONE marker is missing its PASS line"
}

quarantine_path() {
    local path="$1"
    local label="$2"
    [[ -e "$path" ]] || return 0
    local bucket="$QUARANTINE_DIR/$label/$ATTEMPT_ID"
    mkdir -p "$bucket"
    mv -- "$path" "$bucket/"
    echo "Quarantined interrupted artifact: $path"
}

if [[ -f "$RUN_DIR/DONE" ]]; then
    verify_overall_done
    echo "Campaign already complete."
    disable_task
    exit 0
fi
if [[ -f "$RUN_DIR/STOPPED" ]]; then
    echo "Campaign is stopped for inspection: $RUN_DIR/STOPPED" >&2
    disable_task
    exit 2
fi

echo "============================================================"
echo "U=4 PBC FSS campaign attempt $ATTEMPT_ID"
echo "run directory = $RUN_DIR"
echo "max jobs      = $MAX_JOBS"
echo "============================================================"

preflight

SOURCE_PROVENANCE_DIR="$STATE_DIR/source_provenance"
SOURCE_COMMIT_FILE="$SOURCE_PROVENANCE_DIR/source_commit.txt"
SOURCE_MANIFEST="$SOURCE_PROVENANCE_DIR/source_manifest.sha256"
CURRENT_COMMIT="$(git -C "$REPO_ROOT" rev-parse HEAD)"
SOURCE_FILES=(
    Project.toml
    Manifest.toml
    analysis/compute_s2_ed.jl
    analysis/compute_s2_itensor_dmrg.jl
    analysis/analyze_u4_fss.jl
    analysis/summarize_grid.jl
    run/run_grid_local.jl
    run/run_lambda_node.jl
    run/run_u4_fss_campaign.sh
    run/start_u4_fss_detached.ps1
    src/HubbardRenyiDQMC.jl
    src/dqmc_kernel.jl
    src/grover_rank1_update.jl
)
for size in "${SIZES[@]}"; do
    SOURCE_FILES+=("parameters/$(basename "$(parameter_path "$size")")")
done

for orphan in "$STATE_DIR"/source_provenance.tmp.*; do
    quarantine_path "$orphan" "source_provenance_orphan"
done

if [[ -d "$SOURCE_PROVENANCE_DIR" && \
      ( ! -f "$SOURCE_COMMIT_FILE" || ! -f "$SOURCE_MANIFEST" ) ]]; then
    quarantine_path "$SOURCE_PROVENANCE_DIR" "source_provenance_incomplete"
fi

if [[ -d "$SOURCE_PROVENANCE_DIR" ]]; then
    read -r RECORDED_COMMIT < "$SOURCE_COMMIT_FILE"
    if [[ "$CURRENT_COMMIT" != "$RECORDED_COMMIT" ]]; then
        fatal "Source commit changed: expected $RECORDED_COMMIT, got $CURRENT_COMMIT"
    fi
    (cd "$REPO_ROOT" && sha256sum --check "$SOURCE_MANIFEST" >/dev/null)
else
    for source_file in "${SOURCE_FILES[@]}"; do
        if ! git -C "$REPO_ROOT" ls-files --error-unmatch -- "$source_file" >/dev/null 2>&1; then
            fatal "Production source is not tracked: $source_file"
        fi
    done
    SOURCE_STATUS="$(git -C "$REPO_ROOT" status --porcelain -- "${SOURCE_FILES[@]}")"
    if [[ -n "$SOURCE_STATUS" ]]; then
        fatal "Production sources differ from HEAD; commit them before launch"
    fi

    SOURCE_PROVENANCE_TEMP="$STATE_DIR/source_provenance.tmp.$ATTEMPT_ID"
    mkdir "$SOURCE_PROVENANCE_TEMP"
    printf '%s\n' "$CURRENT_COMMIT" > "$SOURCE_PROVENANCE_TEMP/source_commit.txt"
    (cd "$REPO_ROOT" && sha256sum "${SOURCE_FILES[@]}") \
        > "$SOURCE_PROVENANCE_TEMP/source_manifest.sha256"
    mv -- "$SOURCE_PROVENANCE_TEMP" "$SOURCE_PROVENANCE_DIR"
fi

for orphan in "$REFERENCE_DIR"/ed/*.partial.* \
              "$REFERENCE_DIR"/ed/jl_* \
              "$REFERENCE_DIR"/dmrg/*.partial.* \
              "$REFERENCE_DIR"/dmrg/jl_* \
              "$PROCESSED_DIR"/*.partial.* \
              "$PROCESSED_DIR"/.*.tmp.*; do
    quarantine_path "$orphan" "orphan_temps"
done

validate_ed() {
    local path="$1"
    local size="$2"
    local subsystem=$((size / 2))
    awk -F, -v L="$size" -v nA="$subsystem" '
        NR == 1 {
            if ($0 != "L,t,U,mu,Nup,Ndn,nA,boundary,S2_ED,energy,iterations,residual,lanczos_tol") exit 1
        }
        NR == 2 {
            good = ($1 + 0 == L && $3 + 0 == 4 && $7 + 0 == nA &&
                    $8 == "pbc" && $9 + 0 > 0 && $10 + 0 < 0 &&
                    $11 + 0 > 0 && $11 + 0 <= 300 && $12 + 0 >= 0 &&
                    $13 + 0 == 1e-11 && $12 + 0 < $13 + 0)
        }
        END { if (NR != 2 || !good) exit 1 }
    ' "$path"
}

run_ed_reference() {
    local size="$1"
    local subsystem=$((size / 2))
    local output="$REFERENCE_DIR/ed/L${size}_U4_PBC_nA${subsystem}.csv"
    if [[ -f "$output" ]]; then
        validate_ed "$output" "$size"
        echo "ED L=$size already complete."
        return
    fi

    write_status "RUNNING ED L=$size"
    local partial="$output.partial.$ATTEMPT_ID"
    local case_log="$LOG_DIR/ed_L${size}_attempt_$ATTEMPT_ID.log"
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/analysis/compute_s2_ed.jl" \
        --L "$size" --U 4 --nA "$subsystem" --boundary pbc \
        --lanczos-maxiter 300 --lanczos-tol 1e-11 \
        --output "$partial" 2>&1 | tee "$case_log" 9>&-
    validate_ed "$partial" "$size"
    mv -- "$partial" "$output"
    mark_progress
    echo "ED L=$size complete."
}

validate_dmrg() {
    local path="$1"
    local size="$2"
    local seed="$3"
    local sweeps="$4"
    local cap="$5"
    local input_schedule="$6"
    local subsystem=$((size / 2))
    local output_schedule="${input_schedule//,/;}"
    awk -F, -v L="$size" -v nA="$subsystem" -v expected_seed="$seed" \
        -v sweeps="$sweeps" -v cap="$cap" -v schedule="$output_schedule" '
        NR == 1 {
            if ($0 != "L,t,U,mu,nA,periodic,boundary,seed,energy,S2_MPS,normcheck,maxlinkdim,nsweeps,cutoff,maxdim_cap,maxdim_schedule,energy_density,final_sweep_max_truncerr") exit 1
        }
        NR == 2 {
            good = ($1 + 0 == L && $3 + 0 == 4 && $5 + 0 == nA &&
                    $6 == "true" && $7 == "pbc" && $8 + 0 == expected_seed &&
                    $9 + 0 < 0 && $10 + 0 > 0 &&
                    ($11 + 0 - 1)^2 <= 1e-20 &&
                    $12 + 0 > 0 && $12 + 0 <= cap &&
                    $13 + 0 == sweeps && $14 + 0 == 1e-12 &&
                    $15 + 0 == cap && $16 == schedule && $17 + 0 < 0 &&
                    $18 + 0 >= 0 && $18 + 0 <= 1e-8)
        }
        END { if (NR != 2 || !good) exit 1 }
    ' "$path"
}

run_dmrg_reference() {
    local size="$1"
    local seed="$2"
    local sweeps="$3"
    local cap="$4"
    local schedule="$5"
    local subsystem=$((size / 2))
    local output="$REFERENCE_DIR/dmrg/L${size}_U4_PBC_nA${subsystem}_seed${seed}_maxdim${cap}.csv"
    if [[ -f "$output" ]]; then
        validate_dmrg "$output" "$size" "$seed" "$sweeps" "$cap" "$schedule"
        echo "DMRG L=$size seed=$seed maxdim=$cap already complete."
        return
    fi

    write_status "RUNNING DMRG L=$size seed=$seed maxdim=$cap"
    local case_log="$LOG_DIR/dmrg_L${size}_seed${seed}_maxdim${cap}_attempt_$ATTEMPT_ID.log"
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/analysis/compute_s2_itensor_dmrg.jl" \
        --L "$size" --U 4 --nA "$subsystem" --boundary pbc \
        --seed "$seed" --nsweeps "$sweeps" --cutoff 1e-12 \
        --maxdim-schedule "$schedule" \
        --output "$output" 2>&1 | tee "$case_log" 9>&-
    validate_dmrg "$output" "$size" "$seed" "$sweeps" "$cap" "$schedule"
    mark_progress
    echo "DMRG L=$size seed=$seed maxdim=$cap complete."
}

check_dmrg_repeat() {
    local size="$1"
    local subsystem=$((size / 2))
    local level1="$REFERENCE_DIR/dmrg/L${size}_U4_PBC_nA${subsystem}_seed1234_maxdim5200.csv"
    local primary="$REFERENCE_DIR/dmrg/L${size}_U4_PBC_nA${subsystem}_seed1234_maxdim5600.csv"
    local repeat="$REFERENCE_DIR/dmrg/L${size}_U4_PBC_nA${subsystem}_seed11234_maxdim5600.csv"
    local record="$REFERENCE_DIR/dmrg/L${size}_convergence.csv"
    local level1_s2 primary_s2 repeat_s2 level1_energy primary_energy repeat_energy
    local s2_span energy_span
    level1_s2="$(awk -F, 'NR == 2 {print $10}' "$level1")"
    primary_s2="$(awk -F, 'NR == 2 {print $10}' "$primary")"
    repeat_s2="$(awk -F, 'NR == 2 {print $10}' "$repeat")"
    level1_energy="$(awk -F, 'NR == 2 {print $17}' "$level1")"
    primary_energy="$(awk -F, 'NR == 2 {print $17}' "$primary")"
    repeat_energy="$(awk -F, 'NR == 2 {print $17}' "$repeat")"
    read -r s2_span energy_span < <(awk \
        -v a="$level1_s2" -v b="$primary_s2" -v c="$repeat_s2" \
        -v ea="$level1_energy" -v eb="$primary_energy" -v ec="$repeat_energy" \
        -v s2_limit="$DMRG_S2_SPAN_LIMIT" \
        -v energy_limit="$DMRG_ENERGY_DENSITY_SPAN_LIMIT" 'BEGIN {
            min_s = a; if (b < min_s) min_s = b; if (c < min_s) min_s = c
            max_s = a; if (b > max_s) max_s = b; if (c > max_s) max_s = c
            min_e = ea; if (eb < min_e) min_e = eb; if (ec < min_e) min_e = ec
            max_e = ea; if (eb > max_e) max_e = eb; if (ec > max_e) max_e = ec
            ds = max_s - min_s
            de = max_e - min_e
            if (ds > s2_limit || de > energy_limit) exit 1
            printf "%.16e %.16e\n", ds, de
        }')

    if [[ -f "$record" ]]; then
        awk -F, -v L="$size" -v s2_limit="$DMRG_S2_SPAN_LIMIT" \
            -v energy_limit="$DMRG_ENERGY_DENSITY_SPAN_LIMIT" '
            NR == 1 {ok = ($0 == "L,level1_maxdim,level2_maxdim,primary_seed,repeat_seed,S2_span,S2_span_limit,energy_density_span,energy_density_span_limit,status")}
            NR == 2 {seen = ($1 + 0 == L && $2 + 0 == 5200 && $3 + 0 == 5600 &&
                              $4 + 0 == 1234 && $5 + 0 == 11234 &&
                              $6 + 0 <= s2_limit && $7 + 0 == s2_limit &&
                              $8 + 0 <= energy_limit && $9 + 0 == energy_limit &&
                              $10 == "PASS")}
            END {if (NR != 2 || !ok || !seen) exit 1}
        ' "$record"
    else
        local partial="$record.partial.$ATTEMPT_ID"
        printf '%s\n' 'L,level1_maxdim,level2_maxdim,primary_seed,repeat_seed,S2_span,S2_span_limit,energy_density_span,energy_density_span_limit,status' > "$partial"
        printf '%d,5200,5600,1234,11234,%s,%s,%s,%s,PASS\n' \
            "$size" "$s2_span" "$DMRG_S2_SPAN_LIMIT" \
            "$energy_span" "$DMRG_ENERGY_DENSITY_SPAN_LIMIT" >> "$partial"
        mv -- "$partial" "$record"
    fi
    echo "DMRG convergence L=$size PASS: S2 span=$s2_span, energy-density span=$energy_span"
}

for size in "${ED_SIZES[@]}"; do
    run_ed_reference "$size"
done

for size in "${DMRG_SIZES[@]}"; do
    run_dmrg_reference "$size" 1234 22 5600 "$DMRG_LEVEL2_SCHEDULE"
    if [[ "$size" -eq 12 ]]; then
        ED_L12="$REFERENCE_DIR/ed/L12_U4_PBC_nA6.csv"
        DMRG_L12="$REFERENCE_DIR/dmrg/L12_U4_PBC_nA6_seed1234_maxdim5600.csv"
        ED_VALUE="$(awk -F, 'NR == 2 {print $9}' "$ED_L12")"
        DMRG_VALUE="$(awk -F, 'NR == 2 {print $10}' "$DMRG_L12")"
        awk -v ed="$ED_VALUE" -v dmrg="$DMRG_VALUE" 'BEGIN {
            delta = dmrg - ed
            if (delta < 0) delta = -delta
            if (delta > 1e-6) exit 1
        }'
        echo "ED/DMRG L=12 overlap PASS: ED=$ED_VALUE DMRG=$DMRG_VALUE"
    fi
done


for size in "${DMRG_REPEAT_SIZES[@]}"; do
    run_dmrg_reference "$size" 1234 18 5200 "$DMRG_LEVEL1_SCHEDULE"
    run_dmrg_reference "$size" 11234 22 5600 "$DMRG_LEVEL2_SCHEDULE"
    check_dmrg_repeat "$size"
done

prepare_dqmc_resume() {
    local size="$1"
    local raw="$2"
    local name="$3"
    local logs="$raw/logs"
    mkdir -p "$raw" "$logs"

    local seed k padded prefix summary bins log have_summary have_bins
    for seed in "${ROOT_SEEDS[@]}"; do
        for k in $(seq 0 16); do
            printf -v padded '%03d' "$k"
            prefix="${name}_k${padded}_seed${seed}"
            summary="$raw/${prefix}_summary.csv"
            bins="$raw/${prefix}_bins.csv"
            log="$logs/${prefix}.log"
            have_summary=0
            have_bins=0
            [[ -f "$summary" ]] && have_summary=1
            [[ -f "$bins" ]] && have_bins=1

            if (( have_summary != have_bins )); then
                quarantine_path "$summary" "dqmc_L${size}_partial_${prefix}"
                quarantine_path "$bins" "dqmc_L${size}_partial_${prefix}"
                quarantine_path "$log" "dqmc_L${size}_partial_${prefix}"
            elif (( have_summary == 0 )) && [[ -f "$log" ]]; then
                quarantine_path "$log" "dqmc_L${size}_interrupted_${prefix}"
            fi
        done
    done

    local orphan
    for orphan in "$raw"/jl_* "$raw"/.jl_*; do
        quarantine_path "$orphan" "dqmc_L${size}_orphan_temps"
    done

    # A complete-looking pair can still be truncated or internally inconsistent.
    # Validate every existing pair with the launcher's own strict parser, then
    # preserve invalid pairs in quarantine so the normal resume path recomputes
    # only those nodes.
    local validation="$LOG_DIR/dqmc_L${size}_resume_validation_$ATTEMPT_ID.tsv"
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        -e '
            include(ARGS[1])
            config = load_config(ARGS[2])
            raw = ARGS[3]
            for seed in config.root_seeds, k in 0:config.Ninc
                summary = summary_path(config, raw, k, seed)
                bins = bins_path(config, raw, k, seed)
                if isfile(summary) && isfile(bins)
                    try
                        validate_node_outputs(config, summary, bins, k, seed)
                    catch err
                        message = replace(sprint(showerror, err), "\n" => " ", "\t" => " ")
                        println(k, "\t", seed, "\t", message)
                    end
                end
            end
        ' "$REPO_ROOT/run/run_grid_local.jl" "$(parameter_path "$size")" "$raw" \
        > "$validation.partial.$$"
    mv -- "$validation.partial.$$" "$validation"

    local reason
    while IFS=$'\t' read -r k seed reason; do
        [[ -n "$k" ]] || continue
        printf -v padded '%03d' "$k"
        prefix="${name}_k${padded}_seed${seed}"
        summary="$raw/${prefix}_summary.csv"
        bins="$raw/${prefix}_bins.csv"
        log="$logs/${prefix}.log"
        echo "Invalid completed DQMC node L=$size k=$k seed=$seed: $reason"
        quarantine_path "$summary" "dqmc_L${size}_invalid_${prefix}"
        quarantine_path "$bins" "dqmc_L${size}_invalid_${prefix}"
        quarantine_path "$log" "dqmc_L${size}_invalid_${prefix}"
    done < "$validation"
}

verify_completed_dqmc() {
    local size="$1"
    local directory="$DQMC_DIR/L$size"
    local raw="$directory/raw"
    local records="$directory/state/raw_manifest.sha256"
    local summary="$directory/processed/summary.csv"
    local done="$directory/state/dqmc.done"
    [[ -f "$done" && -f "$records" && -f "$summary" ]] ||
        fatal "Incomplete DQMC completion state for L=$size"
    (cd "$raw" && sha256sum --check "$records" >/dev/null)
    local record_count
    record_count="$(wc -l < "$records")"
    [[ "$record_count" -eq "$EXPECTED_RAW_FILES" ]] ||
        fatal "Wrong DQMC raw-manifest count for L=$size: $record_count"
    awk -F, 'NR == 1 {ok = ($0 == "quantity,value,error")}
        $1 == "S2_DQMC" {seen = ($2 + 0 > 0 && $3 + 0 > 0)}
        END {if (!ok || !seen) exit 1}' "$summary"

    local recorded_commit recorded_files recorded_raw_hash recorded_summary_hash
    local actual_raw_hash actual_summary_hash
    recorded_commit="$(awk -F= '$1 == "source_commit" {print $2}' "$done")"
    recorded_files="$(awk -F= '$1 == "raw_files" {print $2}' "$done")"
    recorded_raw_hash="$(awk -F= '$1 == "raw_manifest_sha256" {print $2}' "$done")"
    recorded_summary_hash="$(awk -F= '$1 == "summary_sha256" {print $2}' "$done")"
    actual_raw_hash="$(sha256sum "$records" | awk '{print $1}')"
    actual_summary_hash="$(sha256sum "$summary" | awk '{print $1}')"
    [[ "$recorded_commit" == "$CURRENT_COMMIT" ]] ||
        fatal "DQMC source commit mismatch for L=$size"
    [[ "$recorded_files" == "$EXPECTED_RAW_FILES" ]] ||
        fatal "DQMC done-file count mismatch for L=$size"
    [[ -n "$recorded_raw_hash" && "$recorded_raw_hash" == "$actual_raw_hash" ]] ||
        fatal "DQMC raw-manifest hash mismatch for L=$size"
    [[ -n "$recorded_summary_hash" && "$recorded_summary_hash" == "$actual_summary_hash" ]] ||
        fatal "DQMC summary hash mismatch for L=$size"
}

run_dqmc_size() {
    local size="$1"
    local directory="$DQMC_DIR/L$size"
    local raw="$directory/raw"
    local processed="$directory/processed"
    local state="$directory/state"
    local name
    name="$(run_name "$size")"
    local parameter
    parameter="$(parameter_path "$size")"
    mkdir -p "$raw" "$processed" "$state"

    if [[ -f "$state/dqmc.done" ]]; then
        verify_completed_dqmc "$size"
        echo "DQMC L=$size already complete and verified."
        return
    fi

    if [[ -f "$processed/summary.csv" ]]; then
        quarantine_path "$processed/summary.csv" "dqmc_L${size}_summary_without_done"
    fi
    prepare_dqmc_resume "$size" "$raw" "$name"

    write_status "RUNNING DQMC L=$size"
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/run/run_grid_local.jl" "$parameter" \
        --all-seeds --max-jobs "$MAX_JOBS" --output-dir "$raw"

    local summary_count bins_count
    summary_count="$(find "$raw" -maxdepth 1 -type f \
        -name "${name}_k???_seed*_summary.csv" | wc -l)"
    bins_count="$(find "$raw" -maxdepth 1 -type f \
        -name "${name}_k???_seed*_bins.csv" | wc -l)"
    [[ "$summary_count" -eq "$EXPECTED_NODE_PAIRS" ]]
    [[ "$bins_count" -eq "$EXPECTED_NODE_PAIRS" ]]

    local summary_partial="$processed/summary.csv.partial.$ATTEMPT_ID"
    run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
        "$REPO_ROOT/analysis/summarize_grid.jl" "$raw" \
        --parameter-file "$parameter" --output "$summary_partial"
    awk -F, 'NR == 1 {ok = ($0 == "quantity,value,error")}
        $1 == "S2_DQMC" {seen = ($2 + 0 > 0 && $3 + 0 > 0)}
        END {if (!ok || !seen) exit 1}' "$summary_partial"
    mv -- "$summary_partial" "$processed/summary.csv"

    local records="$state/raw_manifest.sha256"
    (cd "$raw" && find . -maxdepth 1 -type f \
        \( -name '*_bins.csv' -o -name '*_summary.csv' \) -print0 | \
        sort -z | xargs -0 sha256sum) > "$records.tmp.$$"
    local record_count
    record_count="$(wc -l < "$records.tmp.$$")"
    [[ "$record_count" -eq "$EXPECTED_RAW_FILES" ]]
    mv -- "$records.tmp.$$" "$records"
    (cd "$raw" && sha256sum --check "$records" >/dev/null)

    local done_temp="$state/dqmc.done.tmp.$$"
    {
        printf 'completed_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'source_commit=%s\n' "$CURRENT_COMMIT"
        printf 'raw_files=%d\n' "$EXPECTED_RAW_FILES"
        printf 'raw_manifest_sha256=%s\n' "$(sha256sum "$records" | awk '{print $1}')"
        printf 'summary_sha256=%s\n' "$(sha256sum "$processed/summary.csv" | awk '{print $1}')"
    } > "$done_temp"
    mv -- "$done_temp" "$state/dqmc.done"
    verify_completed_dqmc "$size"
    mark_progress
    echo "DQMC L=$size complete."
}

for size in "${SIZES[@]}"; do
    run_dqmc_size "$size"
done

write_status "BUILDING final U=4 comparison"
COMPARISON="$PROCESSED_DIR/u4_fss_L4_to_L22_comparison.csv"
if [[ -f "$COMPARISON" ]]; then
    quarantine_path "$COMPARISON" "comparison_without_done"
fi
COMPARISON_PARTIAL="$COMPARISON.partial.$ATTEMPT_ID"
echo 'L,nA,beta,boundary,U,DQMC_S2,DQMC_error,reference_method,reference_S2,DQMC_minus_reference,z_DQMC_only' \
    > "$COMPARISON_PARTIAL"

for size in "${SIZES[@]}"; do
    subsystem=$((size / 2))
    beta=$((4 * size))
    summary="$DQMC_DIR/L$size/processed/summary.csv"
    dqmc_value="$(awk -F, '$1 == "S2_DQMC" {print $2}' "$summary")"
    dqmc_error="$(awk -F, '$1 == "S2_DQMC" {print $3}' "$summary")"
    if (( size <= 12 )); then
        reference_method="ED"
        reference_file="$REFERENCE_DIR/ed/L${size}_U4_PBC_nA${subsystem}.csv"
        reference_value="$(awk -F, 'NR == 2 {print $9}' "$reference_file")"
    else
        reference_method="DMRG"
        reference_file="$REFERENCE_DIR/dmrg/L${size}_U4_PBC_nA${subsystem}_seed1234_maxdim5600.csv"
        reference_value="$(awk -F, 'NR == 2 {print $10}' "$reference_file")"
    fi
    read -r difference zscore < <(awk -v dqmc="$dqmc_value" -v err="$dqmc_error" \
        -v ref="$reference_value" 'BEGIN {
            difference = dqmc - ref
            printf "%.12e %.6f\n", difference, difference / err
        }')
    printf '%d,%d,%d,pbc,4,%.12e,%.12e,%s,%.12e,%s,%s\n' \
        "$size" "$subsystem" "$beta" "$dqmc_value" "$dqmc_error" \
        "$reference_method" "$reference_value" "$difference" "$zscore" \
        >> "$COMPARISON_PARTIAL"
done
mv -- "$COMPARISON_PARTIAL" "$COMPARISON"

write_status "BUILDING final U=4 FSS diagnostics"
FIT_SUMMARY="$PROCESSED_DIR/u4_fss_L4_to_L22_fits.csv"
FSS_PLOT="$PROCESSED_DIR/u4_fss_L4_to_L22_s2_vs_L.svg"
quarantine_path "$FIT_SUMMARY" "fss_diagnostics_without_done"
quarantine_path "$FSS_PLOT" "fss_diagnostics_without_done"
run_julia --project="$REPO_ROOT" --startup-file=no --history-file=no \
    "$REPO_ROOT/analysis/analyze_u4_fss.jl" \
    --input "$COMPARISON" --fit-output "$FIT_SUMMARY" --plot-output "$FSS_PLOT"
[[ "$(wc -l < "$FIT_SUMMARY")" -eq 11 ]] ||
    fatal "U=4 FSS fit summary must contain one header and ten fits"
grep -q 'c_eff is a finite-window diagnostic only' "$FIT_SUMMARY"
grep -q '<svg ' "$FSS_PLOT"

FINAL_MANIFEST="$STATE_DIR/final_artifacts.sha256"
(cd "$RUN_DIR" && find references dqmc processed -type f \
    \( -name '*.csv' -o -name '*.svg' -o -name '*.sha256' -o -name '*.done' \) -print0 | \
    sort -z | xargs -0 sha256sum) > "$FINAL_MANIFEST.tmp.$$"
mv -- "$FINAL_MANIFEST.tmp.$$" "$FINAL_MANIFEST"
(cd "$RUN_DIR" && sha256sum --check "$FINAL_MANIFEST" >/dev/null)

DONE_TEMP="$RUN_DIR/DONE.tmp.$$"
{
    printf 'PASS: U=4 PBC FSS L=4:2:22 complete\n'
    printf 'completed_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'source_commit=%s\n' "$CURRENT_COMMIT"
    printf 'final_manifest_sha256=%s\n' \
        "$(sha256sum "$FINAL_MANIFEST" | awk '{print $1}')"
    printf 'comparison=%s\n' "$COMPARISON"
    printf 'fit_summary=%s\n' "$FIT_SUMMARY"
    printf 'plot=%s\n' "$FSS_PLOT"
} > "$DONE_TEMP"
mv -- "$DONE_TEMP" "$RUN_DIR/DONE"
mark_progress
write_status "PASS complete"
echo "PASS: U=4 PBC FSS campaign complete."
echo "Comparison: $COMPARISON"
echo "Fits: $FIT_SUMMARY"
echo "Plot: $FSS_PLOT"
disable_task
