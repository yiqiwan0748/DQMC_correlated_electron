#!/usr/bin/env python3
"""Replay saved L48 receipts; no DMRG, eigensolver, oracle or sampling execution."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    base = Path(__file__).resolve().parent
    assert not a.output.exists()
    manifest = json.loads((base / "packet_manifest.json").read_text())
    for n, r in manifest["files"].items():
        p = base / n
        assert not Path(n).is_absolute() and ".." not in Path(n).parts
        assert p.is_file() and not p.is_symlink()
        assert sha(p) == r["sha256"] and p.stat().st_size == r["bytes"]
    with tempfile.TemporaryDirectory(prefix="l48_public_failure_") as temp:
        raw = Path(temp) / "raw"
        raw.mkdir()
        with tarfile.open(base / "scientific_inputs.tar.gz") as tf:
            members = tf.getmembers()
            for m in members:
                assert m.isfile() and not Path(m.name).is_absolute() and ".." not in Path(m.name).parts
                b = tf.extractfile(m).read()
                b.decode("utf-8")
                p = raw / m.name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(b)
        for name in ["source_archive.tar.gz", "public_input_inventory.json"]:
            shutil.copyfile(base / name, raw / name)
        subprocess.run([sys.executable, "-B", str(base / "public_auditor_v3.py"), "--public-inputs", "--raw", str(raw), "--output", str(a.output)], check=True)
    current = json.loads((a.output / "audit.json").read_text())
    original = json.loads((base / "original_complete_audit.json").read_text())
    keys = ["status", "audit_integrity", "pilot_quality", "quality_failures", "original_quality_reasons", "refined", "coarse_refined_errors", "resources", "source_archive_sha256", "source_freeze_sha256", "quality_failure_count", "scientific_completions", "catalogue_status"]
    assert all(current[k] == original[k] for k in keys)
    assert (a.output / "frozen_summary_replay.json").read_bytes() == (base / "original_complete_frozen_summary_replay.json").read_bytes()
    comparison = {"schema": "L48_PUBLIC_VS_ORIGINAL_REPLAY_V1", "status": "PASS", "same_scientific_fields": keys, "quality_failures": 10, "scientific_completions": 0, "private_transport_checks_replayed": False, "public_check_count": current["check_count"], "original_check_count": original["check_count"]}
    (a.output / "public_comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    fs = sorted(p for p in a.output.iterdir() if p.is_file() and p.name != "SHA256SUMS")
    (a.output / "SHA256SUMS").write_text("".join(sha(p) + "  " + p.name + "\n" for p in fs))
    print(json.dumps(comparison))


if __name__ == "__main__":
    main()
