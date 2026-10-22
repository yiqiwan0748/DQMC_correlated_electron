#!/usr/bin/env python3
"""Deterministic design metadata only: no package imports or physics execution."""
import csv
import hashlib
import io
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
NAMESPACE = "P21_SAMPLING_20261022_v1"


def tables():
    seeds = []
    for mode in ("pilot", "final"):
        for package in ("alf", "user"):
            for chain in range(1, 9):
                text = f"{NAMESPACE}|mode={mode}|package={package}|chain={chain}"
                digest = hashlib.sha256(text.encode("utf-8")).digest()
                value = 1 + int.from_bytes(digest[:8], "big") % (2**31 - 2)
                seeds.append(dict(mode=mode, package=package, chain=chain,
                                  seed=value, sha256_input=text,
                                  sha256_digest=digest.hex()))
    assert len(seeds) == len({r["seed"] for r in seeds}) == 32
    seed_record = dict(schema="p21_design_seed_list_v1", namespace=NAMESPACE,
                       algorithm="1+big_endian_uint64(first8_sha256_bytes(input_utf8))%(2^31-2)",
                       native_seed_format="one positive default Int32 integer in native serial ALF seeds file",
                       no_simulation_executed=True, seeds=seeds)
    seed_text = json.dumps(seed_record, indent=2, allow_nan=False) + "\n"
    output = io.StringIO(newline="")
    w = csv.writer(output, lineterminator="\n")
    w.writerow(["component_id", "order", "tau_index", "tau", "spin", "i", "j", "part", "structural_zero"])
    index = 0
    for order in ("GT0", "G0T"):
        for tau in range(11):
            for spin in ("up", "down"):
                for i in range(4):
                    for j in range(4):
                        for part in ("real", "imag"):
                            index += 1
                            w.writerow([f"P21-C{index:04d}", order, tau, f"{tau/10:.1f}",
                                        spin, i, j, part, str(part == "imag").lower()])
    assert index == 1408
    return {"seeds.json": seed_text, "components.csv": output.getvalue()}


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--check", action="store_true")
    args = p.parse_args()
    for name, expected in tables().items():
        path = HERE / name
        if args.check:
            assert path.read_bytes() == expected.encode("utf-8"), name
        else:
            if path.exists():
                raise FileExistsError(path)
            path.write_bytes(expected.encode("utf-8"))
    print("Design metadata verified" if args.check else "Design metadata written; no scientific execution")
