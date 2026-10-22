#!/usr/bin/env python3
"""Production 4–6 GHz broadband optimization per BROADBAND_PRODUCTION_SPEC_FROZEN.

Frozen:
  C6, Objective A (production multimode J), soft-min τ=0.05 α=0.25,
  initial 100 MHz grid (21 freqs), FULL91-C, L-BFGS-B.
Never changes the frequency set mid line-search.
"""
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

from broadband_objective import ALPHA, TAU, eval_bb_and_grad  # noqa: E402
from design_io import load_design_q  # noqa: E402
from multimode_opt_core import S_HI, S_LO, orbits_and_masks  # noqa: E402

# Frozen initial objective grid: 100 MHz on [4.0, 6.0] → 21 points
F_OBJ_GHZ = [4.0 + 0.1 * i for i in range(21)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_uniform_C")
    ap.add_argument("--maxiter", type=int, default=30)
    ap.add_argument("--tau", type=float, default=TAU)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    # Enforce freeze
    assert abs(args.tau - 0.05) < 1e-15 and abs(args.alpha - 0.25) < 1e-15

    orbits, _, _ = orbits_and_masks()
    q0 = load_design_q(args.design)
    assert len(q0) == 16, "C6 freeze requires 16 orbit variables"
    freqs = [f * 1e9 for f in F_OBJ_GHZ]

    history = []
    best = {"J_BB": -1e300, "q": q0.copy(), "rec": None}
    cache = {}
    t0 = time.perf_counter()

    def _eval(q):
        q = np.asarray(q, float)
        key = q.tobytes()
        if key in cache:
            return cache[key]
        rec = eval_bb_and_grad(
            args.grade, q, orbits, freqs, tau=args.tau, alpha=args.alpha
        )
        cache.clear()
        cache[key] = rec
        return rec

    def fun(q):
        rec = _eval(q)
        J = rec["J_BB"]
        if J > best["J_BB"]:
            best["J_BB"] = J
            best["q"] = np.asarray(q, float).copy()
            best["rec"] = {k: rec[k] for k in rec if k != "grad_q"}
            best["rec"]["grad_norm"] = float(np.linalg.norm(rec["grad_q"]))
            print(
                f"NEW BEST J_BB={J:.6f} worstJ={rec['worst_J']:.4f} "
                f"@ {rec['worst_f_Hz']/1e9:.2f} GHz meanJ={rec['J_mean']:.4f}",
                flush=True,
            )
        history.append(
            {
                "J_BB": float(J),
                "worst_J": float(rec["worst_J"]),
                "worst_f_GHz": float(rec["worst_f_Hz"] / 1e9),
                "J_mean": float(rec["J_mean"]),
                "elapsed_s": time.perf_counter() - t0,
            }
        )
        Path(OUT / "bb_production_history.json").write_text(
            json.dumps({"history": history, "best_J_BB": best["J_BB"]}, indent=2) + "\n"
        )
        return -J

    def jac(q):
        return -_eval(q)["grad_q"]

    meta = {
        "spec": "BROADBAND_PRODUCTION_SPEC_FROZEN.md",
        "grade": args.grade,
        "design_start": args.design,
        "n_var": 16,
        "symmetry": "C6",
        "objective": "A_mode_preserving_via_production_J",
        "tau": args.tau,
        "alpha": args.alpha,
        "freqs_GHz": F_OBJ_GHZ,
        "n_freq": len(F_OBJ_GHZ),
        "maxiter": args.maxiter,
        "bounds": [S_LO, S_HI],
    }
    (OUT / "bb_production_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print("START production BB", meta, flush=True)

    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=jac,
        bounds=[(S_LO, S_HI)] * len(q0),
        options={"maxiter": args.maxiter, "maxfun": args.maxiter * 6, "ftol": 1e-12},
    )
    elapsed = time.perf_counter() - t0
    out = {
        "meta": meta,
        "elapsed_s": elapsed,
        "optimizer_success": bool(res.success),
        "optimizer_message": str(res.message),
        "optimizer_nit": int(res.nit),
        "q0": q0.tolist(),
        "q_best": best["q"].tolist(),
        "J_BB_best": best["J_BB"],
        "rec_best": best["rec"],
        "history": history,
    }
    (OUT / "opt_broadband_production_C.json").write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"DONE J_BB={best['J_BB']:.6f} nit={res.nit} elapsed={elapsed:.0f}s msg={res.message}",
        flush=True,
    )


if __name__ == "__main__":
    main()
