#!/usr/bin/env python3
"""Controlled short C6 vs C3 optimization from mapped warm start."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from efficient_adjoint import eval_and_grad_efficient  # noqa: E402
from efficient_adjoint_generic import eval_and_grad_generic  # noqa: E402
from mode_objectives import channel_weights_A, objective_A_from_S  # noqa: E402
from multimode_opt_core import B_PLUS, F_HZ, orbits_and_masks, te0_metrics_from_S  # noqa: E402
from sixfold_orbits import expand_tied as expand_c6  # noqa: E402
from threefold_orbits import build_orbits_c3, expand_tied as expand_c3, map_c6_q_to_c3  # noqa: E402
import sixport_common as sc  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402


def load_q_c6(design, n):
    from design_io import load_design_q

    q = load_design_q(design)
    assert len(q) == n
    return q


def opt_c6(grade, q0, orbits, maxiter, f_hz):
    hist = []
    best = {"J": -1e99, "q": q0.copy()}

    def fun(q):
        bun, g, s = eval_and_grad_efficient(grade, q, orbits, f_hz=f_hz, b_tesla=B_PLUS)
        J = bun["J"]
        hist.append({"J": J, "metrics": bun["metrics"]})
        if J > best["J"]:
            best.update(J=J, q=np.asarray(q, float).copy(), metrics=bun["metrics"])
        print(f"[C6] J={J:.6f}", flush=True)
        return -J, -g

    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=True,
        bounds=[(0.05, 1.80)] * len(q0),
        options={"maxiter": maxiter, "maxfun": maxiter * 4},
    )
    return {
        "symmetry": "C6",
        "n_var": len(q0),
        "nit": int(res.nit),
        "J_best": best["J"],
        "q_best": best["q"].tolist(),
        "metrics_best": {
            k: best["metrics"][k]
            for k in ("T_des_avg", "isolation_avg_dB", "RL_avg_dB", "TE1_avg", "IL_avg_dB")
            if k in best["metrics"]
        },
        "history": [{"J": h["J"]} for h in hist],
        "message": str(res.message),
    }


def opt_c3(grade, q0, orbits, maxiter, f_hz):
    hist = []
    best = {"J": -1e99, "q": q0.copy()}

    def fun(q):
        bun, g, s = eval_and_grad_generic(
            grade,
            q,
            orbits,
            objective_A_from_S,
            channel_weights_A,
            expand_c3,
            f_hz=f_hz,
            b_tesla=B_PLUS,
        )
        J = bun["J"]
        met = te0_metrics_from_S(bun["S"])
        hist.append({"J": J})
        if J > best["J"]:
            best.update(J=J, q=np.asarray(q, float).copy(), metrics=met)
        print(f"[C3] J={J:.6f}", flush=True)
        return -J, -g

    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=True,
        bounds=[(0.05, 1.80)] * len(q0),
        options={"maxiter": maxiter, "maxfun": maxiter * 4},
    )
    return {
        "symmetry": "C3",
        "n_var": len(q0),
        "nit": int(res.nit),
        "J_best": best["J"],
        "q_best": best["q"].tolist(),
        "metrics_best": {
            k: best["metrics"][k]
            for k in ("T_des_avg", "isolation_avg_dB", "RL_avg_dB", "TE1_avg", "IL_avg_dB")
            if "metrics" in best
        },
        "history": hist,
        "message": str(res.message),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_uniform_C_warm_M")
    ap.add_argument("--maxiter", type=int, default=5)
    ap.add_argument("--which", choices=["C6", "C3", "both"], default="both")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orbits6, _, _ = orbits_and_masks()
    q6 = load_q_c6(args.design, len(orbits6))
    centers = bulb_centers_device(sc)
    orbits3, _ = build_orbits_c3(centers)
    q3 = map_c6_q_to_c3(q6, orbits6, orbits3)
    # verify map is neutral for C6 design (expand C3 → should match C6 rods if C6-symmetric)
    s6 = expand_c6(q6, orbits6)
    s3 = expand_c3(q3, orbits3)
    map_err = float(np.max(np.abs(s6 - s3)))
    print("C6→C3 map max |Δs|", map_err, "n_var C3", len(q3), flush=True)

    results = {"map_max_abs_ds": map_err, "n_var_C6": 16, "n_var_C3": len(orbits3)}
    t0 = time.perf_counter()
    if args.which in ("C6", "both"):
        results["C6"] = opt_c6(args.grade, q6, orbits6, args.maxiter, F_HZ)
        (OUT / "opt_c6_short.json").write_text(json.dumps(results["C6"], indent=2) + "\n")
    if args.which in ("C3", "both"):
        results["C3"] = opt_c3(args.grade, q3, orbits3, args.maxiter, F_HZ)
        (OUT / "opt_c3_short.json").write_text(json.dumps(results["C3"], indent=2) + "\n")
    results["elapsed_s"] = time.perf_counter() - t0
    (OUT / "c3_c6_short_opt.json").write_text(json.dumps(results, indent=2) + "\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
