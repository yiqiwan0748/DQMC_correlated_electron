#!/usr/bin/env python3
"""Exact-rational saved-context conditioning diagnostic; never calls a renderer."""
import argparse
import hashlib
import json
import math
import sys
import tomllib
from fractions import Fraction as F
from pathlib import Path

TRACE_SHA = "fe67ce63c4a1ad0ee365e8724df4909f694b700708577f73e3e998b1b0aad033"
RENDERER_SHA = "05a77d109916913bf8139c06e37141e614fd81307daab10a552fb65582a9ae39"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rat(x):
    assert type(x) in (float, int) and math.isfinite(x)
    return F(x)


def pair(x):
    return {"numerator": str(x.numerator), "denominator": str(x.denominator), "approximate": float(x)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trace", type=Path, required=True)
    ap.add_argument("--renderer", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    assert not a.output.exists()
    assert sha(a.trace) == TRACE_SHA and sha(a.renderer) == RENDERER_SHA
    original_freeze = a.renderer.parent / "source_freeze.json"
    assert sha(original_freeze) == "371c4e6e4ac9efa0f05ec90333c5252765535cf14925b5c9cdea08cad992c3e8"
    manifest = a.renderer.parent / "manifest.toml"
    freeze = json.loads(original_freeze.read_text())
    assert sha(manifest) == freeze["files"]["manifest.toml"]
    config = tomllib.loads(manifest.read_text())
    assert config["numerics"]["stationary_bracket_width_max"] == 1e-7
    assert config["numerics"]["stationary_binary64_padding"] == 1e-13
    saved = json.loads(a.trace.read_text())
    assert saved["diagnostic_status"] == "CAPTURED_ORIGINAL_FAILURE"
    assert not saved["history_truncated"] and not saved["extra_scientific_probes"]
    fail = saved["captured_failures"][0]
    assert fail["message"] == "Endpoint padding prevents required root precision"
    root = next(s for s in fail["original_stack"] if s["function"] == "enclosed_root")
    loc = root["locals"]
    left, right, eta = map(rat, [loc["a"], loc["b"], loc["eta"]])
    poles, weights = map(lambda xs: list(map(rat, xs)), [loc["e"], loc["w"]])
    assert len(poles) == len(weights) == 11 and left < right and eta > 0
    assert all(w >= 0 for w in weights)
    events = root["existing_probe_history"]
    pads = [rat(s["existing_locals"]["eps"]) for s in events if s["call"] == "interval_bounds" and s["existing_locals"].get("order") == 1]
    assert pads and len(set(pads)) == 1
    eps = pads[0]
    M = F(0)
    B = F(0)
    terms = []
    normal = rat(sys.float_info.min)
    maximum = rat(sys.float_info.max)
    assert sys.float_info.radix == 2 and sys.float_info.mant_dig == 53
    for e, w in zip(poles, weights):
        assert e < left or e > right
        distances = [abs(left - e), abs(right - e)]
        lo, hi = min(distances), max(distances)
        assert 3 * lo * lo > eta * eta
        curvature = 2 * eta * w * (3 * hi * hi - eta * eta) / (3 * (lo * lo + eta * eta) ** 3)
        derivative = 2 * eta * w * hi / (3 * (lo * lo + eta * eta) ** 2)
        M += curvature
        B += derivative
        # pi<4 supplies a very loose positive term lower bound. Together with
        # the explicit operands, this records a large normal-range margin.
        lower = 2 * eta * w * lo / (4 * (hi * hi + eta * eta) ** 2)
        base_lo, base_hi = lo*lo + eta*eta, hi*hi + eta*eta
        intermediates = {
            "eta": (eta, eta), "distance_abs": (lo, hi),
            "twice_eta": (2*eta, 2*eta), "eta_squared": (eta*eta, eta*eta),
            "distance_squared": (lo*lo, hi*hi), "three_distance_squared": (3*lo*lo, 3*hi*hi),
            "positive_curvature_polynomial": (3*lo*lo-eta*eta, 3*hi*hi-eta*eta),
            "denominator_base": (base_lo, base_hi),
            "denominator_squared": (base_lo**2, base_hi**2),
            "denominator_cubed": (base_lo**3, base_hi**3),
            "pi_times_squared": (3*base_lo**2, 4*base_hi**2),
            "pi_times_cubed": (3*base_lo**3, 4*base_hi**3),
        }
        if w:
            intermediates.update({
                "weight": (w,w), "twice_eta_weight": (2*eta*w,2*eta*w),
                "derivative_numerator_abs": (2*eta*w*lo,2*eta*w*hi),
                "curvature_numerator_abs": (2*eta*w*(3*lo*lo-eta*eta),2*eta*w*(3*hi*hi-eta*eta)),
                "derivative_term_abs": (lower,derivative),
                "curvature_term_abs": (2*eta*w*(3*lo*lo-eta*eta)/(4*base_hi**3),curvature),
            })
        for name, (minimum, largest) in intermediates.items():
            assert normal*1024 < minimum <= largest < maximum/1024, name
        terms.append({"pole": pair(e), "weight": pair(w), "distance_min": pair(lo), "distance_max": pair(hi), "curvature_upper_pi_gt_3": pair(curvature), "absolute_derivative_upper_pi_gt_3": pair(derivative), "positive_derivative_term_lower_pi_lt_4": pair(lower), "all_center_and_curvature_intermediates_normal_with_factor1024_margin": True, "intermediate_envelopes": {name: {"lower": pair(v[0]), "upper": pair(v[1])} for name,v in intermediates.items()}})
    u = F(1, 2**53)
    gamma = 1024 * u / (1 - 1024 * u)
    rho = gamma * B
    tol = F(config["numerics"]["stationary_bracket_width_max"])
    margin = 2 * eps - (M * tol + 2 * rho)
    minimum_width = (2 * eps - 2 * rho) / M
    verdict = "FIXED_PAD_SUBBRACKET_INFEASIBLE_UNDER_DECLARED_ARITHMETIC_MODEL" if margin > 0 else "CONDITIONING_NOT_ESTABLISHED"
    report = {"schema": "X14_SAVED_CONTEXT_CONDITIONAL_CONDITIONING_V1", "diagnostic_verdict": verdict, "scientific_completions": 0, "catalogue_status": "NOT_RUN", "original_failed_job": "6466312", "original_diagnostic_job": "6466728", "trace_sha256": sha(a.trace), "original_renderer_sha256": sha(a.renderer), "analysis_source_sha256": sha(Path(__file__)), "scope": "Only subbrackets inside the retained failed bracket under the explicitly declared normal binary64 evaluation model; not hardware-verified interval arithmetic or a complete root census.", "left": pair(left), "right": pair(right), "retained_width": pair(right-left), "original_width_tolerance": pair(tol), "observed_fixed_derivative_padding": pair(eps), "curvature_upper": pair(M), "absolute_derivative_term_sum_upper": pair(B), "gamma1024": pair(gamma), "conditional_center_evaluation_error_upper": pair(rho), "strict_infeasibility_margin": pair(margin), "conditional_minimum_padded_sign_width": pair(minimum_width), "exact_rational_bounds": True, "conditional_binary64_model": True, "original_renderer_or_probe_executed": False, "source_gate_or_padding_changed": False, "new_fixture_or_eigensolver_or_sampling": False, "original_failure_preserved": True, "terms": terms}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"diagnostic_verdict": verdict, "minimum_width": float(minimum_width), "original_tolerance": float(tol), "curvature_upper": float(M), "rho": float(rho), "scientific_completions": 0}))


if __name__ == "__main__":
    main()
