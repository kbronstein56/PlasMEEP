#!/usr/bin/env python3
"""Phase 8: optional free-91 refinement from best C6 multimode design."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from multimode_opt_core import (
    B_PLUS,
    F_HZ,
    OUT,
    S_HI,
    S_LO,
    WEIGHTS,
    dump_json,
    eval_design,
    orbits_and_masks,
)
from adjoint_direct_dtn import build_system, adjoint_gradient  # noqa: E402


def full_gradient_s(grade, s, b_tesla):
    st = build_system(grade, s, F_HZ, b_tesla, "projected")
    rods = list(range(len(s)))
    # adjoint_gradient returns only selected rods — do all 91 in batches of orbits
    # for memory: compute all rods
    J0, S, grads, _ = adjoint_gradient(st, s, rods)
    g = np.zeros(len(s), float)
    for k, v in grads.items():
        g[k] = v
    met = __import__("multimode_opt_core", fromlist=["te0_metrics_from_S"]).te0_metrics_from_S(S)
    return J0, g, met, S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True)
    ap.add_argument("--grade", default="C")
    ap.add_argument("--maxiter", type=int, default=8)
    ap.add_argument("--B", type=float, default=B_PLUS)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    opt = json.loads(Path(args.opt_json).read_text())
    s0 = np.asarray(opt["best"]["s"], float)
    history = []
    best = {"J": -1e300, "s": s0.copy(), "metrics": None}
    cache = {}
    t0 = time.time()

    def _eval(s):
        key = s.tobytes()
        if key in cache:
            return cache[key]
        J, g, met, S = full_gradient_s(args.grade, s, args.B)
        cache.clear()  # keep only latest
        cache[s.tobytes()] = (J, g, met, S)
        return J, g, met, S

    def fun(s):
        J, g, met, S = _eval(np.asarray(s, float))
        if J > best["J"]:
            best["J"] = J
            best["s"] = np.asarray(s, float).copy()
            best["metrics"] = met
            print(f"NEW BEST J={J:.6f} Tdes={met['T_des_avg']:.4f} TE1={met['TE1_avg']:.4f}", flush=True)
        history.append({"J": J, "T_des_avg": met["T_des_avg"], "TE1_avg": met["TE1_avg"]})
        return -J

    def jac(s):
        J, g, met, S = _eval(np.asarray(s, float))
        return -g

    J0, g0, met0, _ = _eval(s0)
    best.update({"J": J0, "s": s0.copy(), "metrics": met0})
    print(f"free91 start J={J0:.6f}", flush=True)
    res = minimize(
        fun,
        s0,
        method="L-BFGS-B",
        jac=jac,
        bounds=[(S_LO, S_HI)] * len(s0),
        options={"maxiter": args.maxiter, "ftol": 1e-9, "disp": True},
    )
    sf = np.asarray(res.x, float)
    Jf, gf, metf, Sf = _eval(sf)
    # six-source std and C6 deviation
    orbits, _, _ = orbits_and_masks()
    orbit_std = []
    for o in orbits:
        vals = [sf[i] for i in o["rods"]]
        orbit_std.append(float(np.std(vals)))
    out = {
        "warm_start": str(args.opt_json),
        "grade": args.grade,
        "weights": WEIGHTS,
        "J0": float(J0),
        "metrics0": met0,
        "best": {
            "J": float(best["J"]),
            "s": best["s"].tolist(),
            "metrics": best["metrics"],
            "s_min": float(np.min(best["s"])),
            "s_max": float(np.max(best["s"])),
        },
        "final": {"J": float(Jf), "s": sf.tolist(), "metrics": metf},
        "delta_J": float(best["J"] - J0),
        "mean_orbit_std": float(np.mean(orbit_std)),
        "max_orbit_std": float(np.max(orbit_std)),
        "J_std_sources": float(best["metrics"]["J_std"]),
        "elapsed_s": time.time() - t0,
        "optimizer_message": str(res.message),
        "history": history,
    }
    dump_json(OUT / f"opt_free91_{args.grade}.json", out)
    np.savez_compressed(OUT / f"S12_free91_{args.grade}_Bp.npz", S=Sf, s=best["s"])
    prefer_c6 = abs(out["delta_J"]) < 0.02
    lines = [
        "# FREE91_MULTIMODE refinement",
        "",
        f"Warm start: `{args.opt_json}` on grade {args.grade}.",
        f"- J: {J0:.6f} → {best['J']:.6f} (Δ={out['delta_J']:.6f})",
        f"- T_des: {met0['T_des_avg']:.4f} → {best['metrics']['T_des_avg']:.4f}",
        f"- TE1: {met0['TE1_avg']:.4f} → {best['metrics']['TE1_avg']:.4f}",
        f"- mean orbit std (symmetry break): {out['mean_orbit_std']:.4e}",
        "",
        (
            "**Retain C6 as primary** — free-91 gain is small."
            if prefer_c6
            else "**Free-91 improves J substantially** — compare carefully before adopting."
        ),
        "",
    ]
    (OUT / "FREE91_MULTIMODE_REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"delta_J": out["delta_J"], "prefer_c6": prefer_c6, "elapsed_s": out["elapsed_s"]}, indent=2))


if __name__ == "__main__":
    main()
