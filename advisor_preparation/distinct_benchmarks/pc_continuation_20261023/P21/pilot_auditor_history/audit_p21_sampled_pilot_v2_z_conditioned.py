#!/usr/bin/env python3
"""Independent P21 sampled-pilot receipt audit; never executes a scientific solver.

Integrity PASS certifies source/raw-grid/statistical-report fidelity only. P21
remains NOT_RUN and this excluded pilot contributes zero scientific completions.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import statistics
import subprocess
import sys
import tarfile
import tempfile
import tomllib

FREEZE_SHA256 = "be30caae8ca42dbdf7e5915e42505e757f03279a625a4d55c5ce0238774e33a6"
ARCHIVE_SHA256 = "6b314468a6cc74a7bf427400bc028b2fcaed43035904dbb16af1be2da7498214"
SCHEMA = "P21_SAMPLED_TIMING_PILOT_V1"
RELEASE = "p21_sampling_pilot_v1"
REMOTE_ROOT = "/lustre/isaac24/scratch/ywang287/advisor-validation-20261022"
CHECKS = 0
MAX_NUMERIC_DIFFERENCE = 0.0
MAX_CONDITIONED_Z_DIFFERENCE = 0.0


def need(condition, reason):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(reason)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    if isinstance(value, float):
        need(math.isfinite(value), "Finite raw/statistical receipt")
    elif isinstance(value, dict):
        for child in value.values():
            finite(child)
    elif isinstance(value, list):
        for child in value:
            finite(child)


def equal(a, b, label):
    global MAX_NUMERIC_DIFFERENCE, MAX_CONDITIONED_Z_DIFFERENCE
    if isinstance(a, dict) or isinstance(b, dict):
        need(isinstance(a, dict) and isinstance(b, dict) and set(a) == set(b), label + ": object keys")
        zkey = "difference_over_descriptive_standard_error"
        sigma_key = "combined_descriptive_standard_error"
        difference_key = "user_minus_alf"
        conditioned = {zkey, sigma_key, difference_key} <= set(a)
        for key in a:
            if not conditioned or key != zkey:
                equal(a[key], b[key], label + "." + key)
        if conditioned:
            need((a[sigma_key] == 0) == (b[sigma_key] == 0), label + ": exact zero-SE semantics")
            for name, component in (("independent", a), ("reported", b)):
                sigma, difference, z = component[sigma_key], component[difference_key], component[zkey]
                need(sigma >= 0, label + ": nonnegative descriptive SE")
                if sigma == 0:
                    need(z is None, label + ": exact None convention for zero SE")
                else:
                    need(isinstance(z, (int, float)) and math.isfinite(z), label + ": finite standardized diagnostic")
                    # Verify division using this receipt's already independently
                    # checked inputs, so tiny SE cannot conceal an arbitrary Z.
                    equal(z, difference / sigma, label + ": " + name + " exact serialized Z definition")
                    # Cross-receipt means and SEs above still face the unchanged
                    #64-epsilon rule. Compare the reconstructed diagnostic on its
                    # original mean scale rather than amplifying their roundoff.
                    equal(z * sigma, difference, label + ": " + name + " original-mean-scale Z residual")
            if a[sigma_key] > 0:
                MAX_CONDITIONED_Z_DIFFERENCE = max(MAX_CONDITIONED_Z_DIFFERENCE, abs(a[zkey] - b[zkey]))
    elif isinstance(a, list) or isinstance(b, list):
        need(isinstance(a, list) and isinstance(b, list) and len(a) == len(b), label + ": array grid")
        for i, (x, y) in enumerate(zip(a, b)):
            equal(x, y, label + "[" + str(i) + "]")
    elif isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        need(a == b, label + ": exact value")
    else:
        difference = abs(a - b)
        MAX_NUMERIC_DIFFERENCE = max(MAX_NUMERIC_DIFFERENCE, difference)
        need(difference <= 64 * sys.float_info.epsilon * max(1.0, abs(a), abs(b)), label + ": beyond64-epsilon arithmetic allowance")


def raw_native(path, width):
    """Independently parse every OBS_LEGACY vector/paired-sign record."""
    rows = []
    expression = re.compile(r"\(\s*([^(),]+)\s*,\s*([^(),]+)\s*\)")
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        matches = list(expression.finditer(line))
        need(len(matches) == width, f"{path.name}:{line_number}: exactly{width} complex entries")
        need(line[:matches[0].start()].strip() == str(width + 1), "Native vector width marker")
        for left, right in zip(matches, matches[1:]):
            need(not line[left.end():right.start()].strip(), "No unparsed native inter-element tokens")
        def number(text):
            value = float(text.strip().replace("D", "E").replace("d", "e"))
            need(math.isfinite(value), "Finite native scalar")
            return value
        vector = [complex(number(m.group(1)), number(m.group(2))) for m in matches]
        sign = number(line[matches[-1].end():])
        rows.append((vector, sign))
    need(len(rows) == 96, "All96 native bins retained, including32 declared burn bins")
    return rows


def index(tau, spin, row, column):
    return (((tau * 2 + spin) * 4 + row) * 4 + column)


def reconstruct_alf(folder, chain_id, seed):
    path = folder / ("alf_chain_" + str(chain_id))
    gt = raw_native(path / "P21_GT0_scal", 352)
    gz = raw_native(path / "P21_G0T_scal", 352)
    counts = raw_native(path / "P21_Counts_scal", 11)
    need((path / "seeds.initial").read_bytes() == (str(seed) + "\n").encode(), "Exact native startup seed bytes")
    result = []
    endpoint_error = 0.0
    for native_bin, ((a, sa), (b, sb), (count_vector, count_sign)) in enumerate(zip(gt, gz, counts), 1):
        need(sa == sb and abs(sa - 1) <= 1e-10 and count_sign == 1, "Paired native sign averages and unit count divisor")
        need(count_vector == [2 + 0j] * 11, "Exactly2 native traversals at every tau in every raw bin")
        record = {"counts": [2] * 11, "sign": [2 * sa] * 11}
        for order, values in (("GT0", a), ("G0T", b)):
            record[order + ".real"] = [2 * value.real for value in values]
            record[order + ".imag"] = [2 * value.imag for value in values]
        for spin in range(2):
            for row in range(4):
                for column in range(4):
                    initial = index(0, spin, row, column)
                    end = index(10, spin, row, column)
                    identity = float(row == column)
                    endpoint_error = max(endpoint_error, abs(b[initial] - a[initial] + identity * sa),
                        abs(a[end] + a[initial] - identity * sa), abs(b[end] + a[initial]))
        if native_bin > 32:
            result.append(record)
    need(len(result) == 64 and endpoint_error <= 5e-10, "Retained native64-bin grid and fixed-configuration endpoints")
    return result, endpoint_error


def reconstruct_user(target, release, job_id, seeds):
    need(target["schema"] == SCHEMA and target["phase"] == "PILOT_ONLY" and target["status"] == "PILOT_ONLY_COMPLETE"
         and target["scientific_completions"] == 0 and target["job_id"] == job_id, "User pilot/job identity; zero scientific completions")
    for field, name in (("manifest_sha256", "manifest.toml"), ("driver_sha256", "user_sampling.jl"), ("observer_sha256", "UserMatrixObserver.jl")):
        need(target[field] == sha(release / name), "Executed user adapter: " + name)
    need(target["julia_version"] == "1.12.1" and target["jdqmcframework_version"] == "1.2.6"
         and target["julia_threads"] == target["blas_threads"] == 1, "User package versions/single threads")
    need(target["original_kernel_sha256"] == "b72bbbf158381b0f776c1fa259d532b069094c47e361d3490c19ea0f43ada4be"
         and target["observer_kernel_sha256"] == "911a61642b23a4465639c3b3e4c63ff9d7cf05d51da288f3224f1775c0aba891", "Original/observer native kernel pins")
    chains = target["chains"]
    need(len(chains) == 2 and {(c["chain"], c["seed"]) for c in chains} == set(seeds.items()), "Both independent fixed user chains")
    result, diagnostics = {}, []
    for chain in chains:
        need(chain["status"] == "PILOT_ONLY_COMPLETE" and chain["burn_cycles"] == 64 and chain["retained_cycles"] == 128
             and chain["update_calls"] == 384, "User burn64 plus retained128 full forward/reverse cycles")
        need(0 <= chain["acceptance"] <= 1 and chain["elapsed_seconds"] > 0, "User acceptance/runtime")
        bins = chain["bins"]
        need(len(bins) == 64 and [b["bin"] for b in bins] == list(range(1, 65)), "Ordered user64-bin coverage")
        ds = [chain["cumulative_deltaG_after_burn"]]
        reconstructed = []
        endpoint = 0.0
        for record in bins:
            need(record["counts_by_tau"] == [2] * 11, "Two actual user callbacks per tau per retained bin")
            need(record["primary_stage"] == "initial(tau0) and pre_stabilization(tau1..M)"
                 and record["measurement_origin"] == "actual user maintained buffers; no dense recomputation"
                 and record["ordering"] == "tau-major,spin(up,down),row_i,column_j; all indices zero-based externally", "Actual user primary callback representation")
            sign = record["sign_denominator_sums_by_tau"]
            need(len(sign) == 11 and all(abs(value / 2 - 1) <= 1e-10 for value in sign), "Per-tau moving-cut paired sign denominators")
            row = {"counts": record["counts_by_tau"], "sign": sign}
            for order in ("GT0", "G0T"):
                for part in ("real", "imag"):
                    values = record[order + "_numerator_" + part + "_sums"]
                    need(len(values) == 352, "Complete actual-user352-entry order/part export")
                    row[order + "." + part] = values
            reconstructed.append(row)
            ds.append(record["cumulative_deltaG"])
            endpoint = max(endpoint, record["endpoint_contemporaneous_max_error"])
        need(all(0 <= value <= 1e-8 for value in ds) and all(a <= b for a, b in zip(ds, ds[1:]))
             and ds[-1] == chain["cumulative_deltaG"], "Cumulative user stability including burn")
        need(endpoint <= 5e-10, "Contemporaneous moving-cut user endpoints")
        result[chain["chain"]] = reconstructed
        diagnostics.append({"chain": chain["chain"], "count_per_tau_per_bin": 2, "retained_bins": 64,
            "update_calls": chain["update_calls"], "endpoint_max_error": endpoint, "cumulative_deltaG": ds[-1]})
    return result, diagnostics


def paired_ratio_statistics(nums, signs):
    """Independent scalar paired-ratio jackknife from retained raw sums."""
    n = len(nums)
    sn, sd = math.fsum(nums), math.fsum(signs)
    need(n in (32, 64) and len(signs) == n and sd > 0, "Complete paired ratio bin grid")
    full = sn / sd
    leave = []
    for x, y in zip(nums, signs):
        need(sd - y > 0, "All delete-one paired denominators positive")
        leave.append((sn - x) / (sd - y))
    average = statistics.mean(leave)
    variance = (n - 1) / n * math.fsum((value - average) ** 2 for value in leave)
    return {"plug_in": full, "bias_corrected": n * full - (n - 1) * average, "standard_error": math.sqrt(variance)}


def chain_statistics(bins, order, part, tau, spin, row, column):
    ix = index(tau, spin, row, column)
    nums = [b[order + "." + part][ix] for b in bins]
    signs = [b["sign"][tau] for b in bins]
    fine = paired_ratio_statistics(nums, signs)
    coarse = paired_ratio_statistics([a + b for a, b in zip(nums[::2], nums[1::2])],
                                    [a + b for a, b in zip(signs[::2], signs[1::2])])
    series = [a / b for a, b in zip(nums, signs)]
    center = statistics.mean(series)
    s2 = math.fsum((value - center) ** 2 for value in series)
    correlation = math.fsum((a - center) * (b - center) for a, b in zip(series[:-1], series[1:])) / s2 if s2 else 0.0
    ess = min(64.0, max(1.0, 64 * (1 - correlation) / max(1e-12, 1 + correlation)))
    value = {**fine, "coarsened32_standard_error": coarse["standard_error"], "observed_constant_bins": s2 == 0,
        "lag1_bin_correlation": correlation, "lag1_bin_ESS_estimate": ess,
        "first_half_mean": statistics.mean(series[:32]), "second_half_mean": statistics.mean(series[32:])}
    return value, series


def backend_statistics(chains, coordinate):
    stats, series = [], []
    for chain in chains:
        value, sequence = chain_statistics(chain, *coordinate)
        stats.append(value)
        series.append(sequence)
    means = [value["bias_corrected"] for value in stats]
    between = statistics.stdev(means) / math.sqrt(2)
    within = math.sqrt(math.fsum(max(value["standard_error"], value["coarsened32_standard_error"]) ** 2 for value in stats)) / 2
    split = [half for sequence in series for half in (sequence[:32], sequence[32:])]
    w = statistics.mean(statistics.variance(half) for half in split)
    b = 32 * statistics.variance(statistics.mean(half) for half in split)
    rhat = math.sqrt((31 * w / 32 + b / 32) / w) if w else (1.0 if b == 0 else None)
    return {"mean": statistics.mean(means), "between_chain_standard_error": between,
        "within_chain_blocked_standard_error": within, "descriptive_standard_error": max(between, within),
        "split_Rhat": rhat, "lag1_bin_ESS_estimate": math.fsum(value["lag1_bin_ESS_estimate"] for value in stats), "chains": stats}


def accounting(text, job_id):
    header = None
    for line in text.splitlines():
        fields = [x.strip() for x in line.strip().split("|")]
        if {"JobIDRaw", "State", "ExitCode"} <= set(fields):
            header = fields
        elif header and len(fields) == len(header):
            row = dict(zip(header, fields))
            if row.get("JobIDRaw") == job_id and row.get("State") == "COMPLETED" and row.get("ExitCode") == "0:0":
                need(row.get("ReqCPUS") == "1" and row.get("ReqMem") == "3500M", "Actual oneCPU3500MB allocation")
                return row
    raise ValueError("Missing exact-job COMPLETED/0:0 named-column scheduler accounting")


def inventory(folder):
    return {path.relative_to(folder).as_posix(): sha(path) for path in folder.rglob("*") if path.is_file()}


def audit(args):
    result, release = args.result_dir.resolve(), args.release_dir.resolve()
    freeze_path = (args.freeze or release.parent / "p21-sampled-timing-pilot-v1-freeze.json").resolve()
    need(sha(freeze_path) == FREEZE_SHA256, "Pinned45-file external release freeze record")
    freeze = json.loads(freeze_path.read_text())
    need(freeze["archive_sha256"] == ARCHIVE_SHA256 and freeze["members"] == len(freeze["files"]) == 45
         and freeze["scientific_completions"] == 0 and freeze["phase"] == "PILOT_ONLY", "Frozen pilot scope")
    inputs = [freeze_path]
    for rel, digest in freeze["files"].items():
        parts = PurePosixPath(rel).parts
        need(parts[0] == RELEASE and ".." not in parts, "Safe frozen release member")
        path = release.joinpath(*parts[1:])
        need(path.is_file() and not path.is_symlink() and sha(path) == digest, "Frozen release bytes: " + rel)
        inputs.append(path)
    archive = args.archive or freeze_path.parent / freeze["archive"]
    archive_verified = False
    if archive.exists():
        need(sha(archive) == ARCHIVE_SHA256, "Frozen archive bytes")
        with tarfile.open(archive) as tf:
            members = tf.getmembers()
            need(len(members) == 45 and {m.name for m in members} == set(freeze["files"]), "Exact archive source closure")
            for member in members:
                need(member.isfile() and member.mode == 0o644 and hashlib.sha256(tf.extractfile(member).read()).hexdigest() == freeze["files"][member.name], "Archive member hash/mode")
        inputs.append(archive)
        archive_verified = True
    cfg = tomllib.loads((release / "manifest.toml").read_text())
    need(cfg["schema"] == SCHEMA and cfg["phase"] == "PILOT_ONLY" and cfg["scientific_completion_allowed"] is False
         and (cfg["L"], cfg["Ltau"], cfg["t"], cfg["U"], cfg["beta"], cfg["dtau"], cfg["mu_centered"]) == (4, 10, 1, 2, 1, 0.1, 0), "Exact P21 fixture/ensemble")
    need((cfg["chains_per_backend"], cfg["retained_bins"], cfg["cycles_per_bin"], cfg["burn_cycles"],
          cfg["retained_cycles"], cfg["alf_NBin"], cfg["alf_discard_initial_bins"]) == (2, 64, 2, 64, 128, 96, 32), "Exact predeclared pilot budgets")
    closure = json.loads((release / "closure.json").read_text())
    for rel, digest in closure["assets"].items():
        need(sha(release / rel) == digest, "Native source/proof/license asset: " + rel)
    alf_proof = json.loads((release / "provenance/alf_adapter_proof.json").read_text())
    user_proof = tomllib.loads((release / "provenance/user_adapter_proof.toml").read_text())
    linkage = json.loads((release / "provenance/observer_linkage_audit.json").read_text())
    need(alf_proof["status"] == user_proof["status"] == "ADAPTER_PROOF_PASS_NOT_EXTERNAL_BENCHMARK", "Preserved native/user eligibility proof statuses")
    need(linkage["status"] == "PINNED_SINGLE_PROVIDER_PASS" and alf_proof["binary_sha256"] == linkage["binary_sha256"] == closure["binary_sha256"], "Prior binary/proof identity")
    need(linkage["provider_sha256"] == closure["provider_sha256"], "Prior provider pin")
    seeds = json.loads((release / "seeds.json").read_text())
    need(seeds["namespace"] == cfg["seed_namespace"] and seeds["production_eligible"] is False, "Fixed excluded seed namespace")
    seed_map = {(r["backend"], r["chain"]): r["seed"] for r in seeds["chains"]}
    need(len(seed_map) == len(seeds["chains"]) == 4 and set(seed_map) == {(b, c) for b in ("alf", "user") for c in (1, 2)}
         and len(set(seed_map.values())) == 4 and not set(seed_map.values()) & {20261022, 20261022021}, "All fresh fixed native/user seeds")
    before = inventory(result)
    required = [result / name for name in ("execution.json", "preflight.json", "user_sampling.toml", "summary.json")]
    need(all(path.is_file() for path in required), "Complete sampled pilot and original summary receipts")
    execution = json.loads(required[0].read_text())
    preflight = json.loads(required[1].read_text())
    target = tomllib.loads(required[2].read_text())
    summary = json.loads(required[3].read_text())
    for data in (execution, preflight, target, summary):
        finite(data)
    need(execution["schema"] == "P21_SAMPLED_PILOT_EXECUTION_V1" and execution["status"] == "PILOT_ONLY_COMPLETE"
         and execution["job_id"] == args.job_id and execution["scientific_completions"] == 0, "Actual pilot execution/job identity")
    process_map = {row["name"]: row for row in execution["processes"]}
    need(len(process_map) == len(execution["processes"]) == 3 and set(process_map) == {"alf_chain_1", "alf_chain_2", "user_sampling"}, "All sampled process roles")
    need(all(row["exit_code"] == 0 and row["proof_mode_present"] is False and row["elapsed_seconds"] > 0 for row in process_map.values()), "Sampled mode, positive timings and successful exits")
    need(preflight["schema"] == "P21_SAMPLED_PILOT_PREFLIGHT_V1" and preflight["status"] == "SOURCE_RUNTIME_AND_PRIOR_PROOFS_VERIFIED"
         and preflight["job_id"] == args.job_id and preflight["proof_mode_unset_in_sampled_children"] is True, "Prior-proof/source/runtime preflight")
    need(preflight["binary_sha256"] == closure["binary_sha256"] and preflight["provider_sha256"] == closure["provider_sha256"]
         and preflight["provider"] == closure["provider"], "Actual binary/provider receipt pins")
    need(preflight["proof_root"] == closure["remote_proof_root"] and preflight["user_source_root"] == REMOTE_ROOT + "/source/renyi_disjoint_corrected", "Actual native roots")
    pins = {str(PurePosixPath(preflight["proof_root"]) / rel): digest for rel, digest in closure["remote_proof_pins"].items()}
    pins.update({str(PurePosixPath(preflight["user_source_root"]) / rel): digest for rel, digest in closure["user_base_pins"].items()})
    need(preflight["actual_source_receipt_hashes"] == pins, "Exact fresh native/source/proof receipt closure")
    expected_files = {str(PurePosixPath(name).relative_to(RELEASE)): digest for name, digest in freeze["files"].items() if not name.endswith("/SHA256SUMS")}
    need(preflight["release_files"] == expected_files and preflight["release_sha256sums_sha256"] == sha(release / "SHA256SUMS"), "Executed entire frozen source closure")
    need("not found" not in preflight["ldd"], "No unresolved native libraries")
    provider_lines = [line for line in preflight["ldd"].splitlines() if any(x in line for x in ("libopenblas", "libblas", "liblapack"))]
    need(bool(provider_lines) and all(str(PurePosixPath(closure["provider"]).parent) + "/" in line for line in provider_lines), "Pinned provider directory in dynamic-link receipt")
    alf, endpoint_errors = [], []
    for chain in (1, 2):
        path = result / f"alf_chain_{chain}"
        need((path / "parameters").read_bytes() == (release / "parameters.in").read_bytes(), "Exact ordinary native parameters")
        process = process_map[f"alf_chain_{chain}"]
        need(process["seed"] == seed_map["alf", chain] and process["initial_seed_sha256"] == sha(path / "seeds.initial")
             and process["parameters_sha256"] == sha(path / "parameters"), "Native executed seed/parameter hashes")
        need(process["argv"][-1] == preflight["binary"], "Native sampled executable path")
        bins, endpoint = reconstruct_alf(result, chain, seed_map["alf", chain])
        alf.append(bins)
        endpoint_errors.append(endpoint)
    user, user_diagnostics = reconstruct_user(target, release, args.job_id, {c: seed_map["user", c] for c in (1, 2)})
    need(summary["schema"] == SCHEMA and summary["phase"] == "PILOT_ONLY" and summary["catalogue_status"] == "NOT_RUN"
         and summary["status"] == "SAMPLED_PILOT_NUMERICALLY_ELIGIBLE" and summary["scientific_completions"] == 0
         and summary["job_id"] == args.job_id and summary["real_components_retained"] == 1408, "Faithful excluded-pilot summary identity")
    need(summary["manifest_sha256"] == sha(release / "manifest.toml") and summary["seed_sha256"] == sha(release / "seeds.json")
         and summary["summary_tool_sha256"] == sha(release / "summarize_pilot.py"), "Summary inputs/source pins")
    need(len(summary["components"]) == 1408, "All1408 real summarized components")
    independent_components = []
    position = 0
    for order in ("GT0", "G0T"):
        for part in ("real", "imag"):
            for tau in range(11):
                for spin in range(2):
                    for row in range(4):
                        for column in range(4):
                            coordinate = order, part, tau, spin, row, column
                            left = backend_statistics(alf, coordinate)
                            right = backend_statistics([user[1], user[2]], coordinate)
                            difference = right["mean"] - left["mean"]
                            sigma = math.hypot(left["descriptive_standard_error"], right["descriptive_standard_error"])
                            record = {"order": order, "part": part, "tau_index": tau, "tau": tau * 0.1,
                                "spin": "up" if spin == 0 else "down", "row": row, "column": column,
                                "alf": left, "user": right, "user_minus_alf": difference,
                                "combined_descriptive_standard_error": sigma,
                                "difference_over_descriptive_standard_error": difference / sigma if sigma else None}
                            equal(record, summary["components"][position], "Independent raw component " + str(position))
                            independent_components.append(record)
                            position += 1
    equal(max(endpoint_errors), summary["alf_endpoint_max_error"], "Independent native endpoint maximum")
    for rel, digest in summary["input_sha256"].items():
        path = PurePosixPath(rel)
        need(not path.is_absolute() and ".." not in path.parts and sha(result.joinpath(*path.parts)) == digest, "All exact summary raw inputs: " + rel)
    with tempfile.TemporaryDirectory(prefix="p21_independent_replay_", dir=result.parent) as temp:
        output = Path(temp) / "summary_replay.json"
        subprocess.run([sys.executable, "-B", str(release / "summarize_pilot.py"), "--result-dir", str(result), "--output", str(output)], check=True, capture_output=True, text=True, timeout=60)
        replay = json.loads(output.read_text())
        equal(replay, summary, "Frozen statistical summary replay")
        replay_hash = sha(output)
        replay_byte_identical = output.read_bytes() == (result / "summary.json").read_bytes()
    need(inventory(result) == before, "Independent audit/frozen replay did not mutate any input result file")
    transfer_path = args.transfer or result.parent / "transfer.json"
    need(transfer_path.is_file(), "Remote/local transfer metadata required")
    transfer = json.loads(transfer_path.read_text())
    need(str(transfer["job_id"]) == args.job_id and transfer["remote_writes"] is False, "Read-only remote transfer job identity")
    transferred = set()
    for row in transfer["files"]:
        rel = PurePosixPath(row["local"])
        need(not rel.is_absolute() and ".." not in rel.parts, "Safe transfer member")
        if row["status"] != "RETRIEVED":
            continue
        path = transfer_path.parent.joinpath(*rel.parts)
        need(path.is_file() and path.stat().st_size == row["bytes"] and sha(path) == row["sha256"], "Exact remote/local text bytes: " + row["local"])
        transferred.add(path.resolve())
        if path.resolve().is_relative_to(result):
            expected_remote = REMOTE_ROOT + "/results/p21_sampled_pilot_v1_" + args.job_id + "/" + path.resolve().relative_to(result).as_posix()
            need(row["remote"] == expected_remote or path.name == "slurm.log", "Exact remote sampled result/job path")
        inputs.append(path)
    need(all(path.resolve() in transferred for path in required), "Every essential raw/summary receipt is remote-hash bound")
    for chain in (1, 2):
        need(all((result / f"alf_chain_{chain}" / name).resolve() in transferred for name in
            ("P21_GT0_scal", "P21_G0T_scal", "P21_Counts_scal", "parameters", "seeds.initial")), "Every native matrix/count/startup raw file transferred")
    scheduler_path = args.scheduler or result.parent / "scheduler.txt"
    need(scheduler_path.is_file(), "Scheduler accounting receipt required")
    schedule = accounting(scheduler_path.read_text(), args.job_id)
    inputs.extend([transfer_path, scheduler_path])
    inputs.extend(result / rel for rel in before)
    unique = list(dict.fromkeys(path.resolve() for path in inputs))
    history_path = Path(__file__).resolve().parent / "p21_auditor_conditioning_correction_v2.json"
    tool_history = None
    if history_path.exists():
        history = json.loads(history_path.read_text())
        need(history["new_auditor_sha256"] == sha(__file__), "Auditor correction/version receipt")
        snapshot = history_path.parent / history["preserved_original_auditor"]
        need(sha(snapshot) == history["old_auditor_sha256"], "Preserved original auditor bytes")
        old_receipt = history_path.parent / history["original_audit_receipt"]
        need(sha(old_receipt) == history["original_audit_receipt_sha256"], "Preserved original audit failure receipt")
        inputs.extend([history_path, snapshot, old_receipt])
        unique = list(dict.fromkeys(path.resolve() for path in inputs))
        tool_history = {"correction_receipt": str(history_path), "correction_receipt_sha256": sha(history_path), **history}
    return {"schema": "P21_SAMPLED_PILOT_INDEPENDENT_AUDIT_V1", "status": "PASS", "job_id": args.job_id,
        "audit_scope": "Independent source/raw-grid/statistical-report integrity; no scientific equivalence test or catalogue promotion",
        "catalogue_id": "P21", "catalogue_status": "NOT_RUN", "scientific_completions": 0,
        "real_components_independently_reconstructed": 1408, "complex_components_per_order": 352,
        "native_raw_bins_per_chain": 96, "native_burn_bins_excluded": 32, "retained_bins_per_backend_chain": 64,
        "chains_per_backend": 2, "count_per_tau_per_bin": 2, "source_members": 45,
        "archive_verified": archive_verified, "freeze_sha256": FREEZE_SHA256,
        "native_source_runtime_receipt_pins_verified": True, "binary_sha256": closure["binary_sha256"],
        "provider_sha256": closure["provider_sha256"], "scheduler_accounting": schedule,
        "native_endpoint_max_error": max(endpoint_errors), "user_callback_diagnostics": user_diagnostics,
        "independent_components": independent_components, "frozen_replay": {"byte_identical": replay_byte_identical,
            "metrics_agree_within_64_epsilon": True, "sha256": replay_hash, "no_input_mutations": True},
        "maximum_numeric_comparison_difference": MAX_NUMERIC_DIFFERENCE, "checks": CHECKS,
        "maximum_conditioning_amplified_Z_difference": MAX_CONDITIONED_Z_DIFFERENCE,
        "diagnostic_conditioning_rule": "Raw means/SEs/differences remain subject to strict64-epsilon comparison. Every Z must independently equal its own reported difference/SE within64epsilon and reconstruct that difference on the original mean scale; zero-SE/None semantics are exact. Cross-receipt Z differences are recorded separately because division amplifies independently bounded input roundoff. No scientific gate is changed.",
        "auditor_tool_history": tool_history,
        "inputs": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)} for path in unique],
        "no_scientific_solver_executed": True,
        "limitations": "Integrity PASS is not ALF-versus-user equivalence, equilibrium, precision or1408-component confidence coverage. Independent arithmetic reuses no frozen summary functions. Binary/provider and fresh native-source identities are checked against immutable runtime/transfer receipts and pinned source snapshots; no binary/library is copied or re-executed. User callback provenance is source/trajectory-proof and exported-count evidence; raw retained measurements cannot reconstruct an unrecorded HS trajectory."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("result-dir", "release-dir", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--job-id", required=True)
    for name in ("freeze", "archive", "transfer", "scheduler"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    need(args.job_id.isdecimal() and not args.output.exists(), "Numeric job ID and fresh immutable audit output required")
    try:
        report = audit(args)
    except Exception as exc:
        report = {"schema": "P21_SAMPLED_PILOT_INDEPENDENT_AUDIT_V1", "status": "UNRESOLVED",
            "job_id": args.job_id, "catalogue_id": "P21", "catalogue_status": "NOT_RUN", "scientific_completions": 0,
            "reason": type(exc).__name__ + ": " + str(exc), "checks": CHECKS, "no_scientific_solver_executed": True}
    report["auditor_sha256"] = sha(__file__)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({key: report.get(key) for key in ("status", "scientific_completions", "real_components_independently_reconstructed", "reason")}, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
