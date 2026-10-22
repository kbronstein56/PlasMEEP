#!/usr/bin/env python3
"""Short controlled L-BFGS-B opts: Objective A vs B from same start/mesh/budget."""
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
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from efficient_adjoint_generic import eval_and_grad_generic  # noqa: E402
from mode_objectives import (  # noqa: E402
    OBJECTIVE_FREEZE,
    channel_weights_A,
    channel_weights_B,
    mode_resolved_metrics,
    objective_A_from_S,
    objective_B_from_S,
    write_freeze,
)
from multimode_opt_core import B_PLUS, F_HZ, expand_tied, orbits_and_masks  # noqa: E402
from sixfold_orbits import expand_tied as expand_c6  # noqa: E402


def load_q(design, n):
    p = ROOT / "outputs/circulator_inverse_design/multimode_optimization" / f"{design}.json"
    d = json.loads(p.read_text())
    q = np.asarray(d.get("q_best", d.get("q", d.get("q0"))), float)
    assert len(q) == n
    return q


def run_opt(label, obj_fn, w_fn, grade, q0, orbits, maxiter, f_hz):
    hist = []
    best = {"J": -1e99, "q": q0.copy()}

    def fun(q):
        bun, g, s = eval_and_grad_generic(
            grade, q, orbits, obj_fn, w_fn, expand_c6, f_hz=f_hz, b_tesla=B_PLUS
        )
        J = bun["J"]
        mr = mode_resolved_metrics(bun["S"])
        hist.append({"J": J, "metrics": {k: mr[k] for k in mr if k != "per_source"}})
        if J > best["J"]:
            best["J"] = J
            best["q"] = np.asarray(q, float).copy()
            best["metrics"] = hist[-1]["metrics"]
        print(f"[{label}] J={J:.6f} des0={mr['des0_avg']:.4f} des1={mr['des1_avg']:.4f}", flush=True)
        # maximize J → minimize -J
        return -J, -g

    bounds = [(0.05, 1.80)] * len(q0)
    t0 = time.perf_counter()
    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxiter": maxiter, "ftol": 1e-10, "maxfun": maxiter * 4},
    )
    elapsed = time.perf_counter() - t0
    return {
        "label": label,
        "weights": OBJECTIVE_FREEZE,
        "elapsed_s": elapsed,
        "nit": int(res.nit),
        "success": bool(res.success),
        "message": str(res.message),
        "q0": q0.tolist(),
        "q_best": best["q"].tolist(),
        "J_best": best["J"],
        "metrics_best": best.get("metrics"),
        "history": hist,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_uniform_C_warm_M")
    ap.add_argument("--maxiter", type=int, default=5)
    ap.add_argument("--which", choices=["A", "B", "both"], default="both")
    ap.add_argument("--f-ghz", type=float, default=F_HZ / 1e9)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    write_freeze(OUT)

    orbits, _, _ = orbits_and_masks()
    q0 = load_q(args.design, len(orbits))
    f_hz = args.f_ghz * 1e9
    results = {}
    if args.which in ("A", "both"):
        results["A"] = run_opt(
            "A_mode_preserving", objective_A_from_S, channel_weights_A, args.grade, q0, orbits, args.maxiter, f_hz
        )
        (OUT / "opt_mode_A.json").write_text(json.dumps(results["A"], indent=2) + "\n")
    if args.which in ("B", "both"):
        results["B"] = run_opt(
            "B_port_power", objective_B_from_S, channel_weights_B, args.grade, q0, orbits, args.maxiter, f_hz
        )
        (OUT / "opt_mode_B.json").write_text(json.dumps(results["B"], indent=2) + "\n")
    (OUT / "mode_ab_short_opt.json").write_text(json.dumps(results, indent=2) + "\n")
    print("done", {k: results[k]["J_best"] for k in results}, flush=True)


if __name__ == "__main__":
    main()
