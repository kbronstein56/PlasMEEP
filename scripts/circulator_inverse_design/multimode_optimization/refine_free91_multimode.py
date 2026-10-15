#!/usr/bin/env python3
"""Phase 8: optional free-91 refinement from best C6 multimode design.

Uses the same frozen objective; each rod is an independent variable.
Gradient via the validated direct-DtN adjoint (singleton 'orbits').
"""
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
    eval_and_grad,
    orbits_and_masks,
)


def singleton_orbits(n=91):
    return [{"rods": [i], "radius_a": 0.0, "size": 1} for i in range(n)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True)
    ap.add_argument("--grade", default="C")
    ap.add_argument("--maxiter", type=int, default=5)
    ap.add_argument("--B", type=float, default=B_PLUS)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    opt = json.loads(Path(args.opt_json).read_text())
    s0 = np.asarray(opt["best"]["s"], float)
    orbits = singleton_orbits(len(s0))
    history = []
    best = {"J": -1e300, "s": s0.copy(), "metrics": None}
    cache = {"s": None, "J": None, "g": None, "met": None}
    t0 = time.time()

    def _eval(s):
        s = np.asarray(s, float)
        if cache["s"] is not None and np.allclose(s, cache["s"], rtol=0, atol=0):
            return cache["J"], cache["g"], cache["met"]
        # eval_and_grad expects q with len(orbits); for singletons q==s
        bun, g, s_out = eval_and_grad(
            args.grade, s, orbits, F_HZ, args.B, with_absorption=False
        )
        cache["s"] = s.copy()
        cache["J"] = float(bun["J"])
        cache["g"] = np.asarray(g, float)
        cache["met"] = bun["metrics"]
        return cache["J"], cache["g"], cache["met"]

    def fun(s):
        J, g, met = _eval(s)
        if J > best["J"]:
            best["J"] = J
            best["s"] = np.asarray(s, float).copy()
            best["metrics"] = met
            print(
                f"NEW BEST J={J:.6f} Tdes={met['T_des_avg']:.4f} TE1={met['TE1_avg']:.4f}",
                flush=True,
            )
        history.append(
            {
                "J": J,
                "T_des_avg": met["T_des_avg"],
                "TE1_avg": met["TE1_avg"],
                "isolation_avg_dB": met["isolation_avg_dB"],
            }
        )
        return -J

    def jac(s):
        J, g, met = _eval(s)
        return -g

    J0, g0, met0 = _eval(s0)
    best.update({"J": J0, "s": s0.copy(), "metrics": met0})
    print(f"free91 start J={J0:.6f}", flush=True)
    res = minimize(
        fun,
        s0,
        method="L-BFGS-B",
        jac=jac,
        bounds=[(S_LO, S_HI)] * len(s0),
        options={"maxiter": args.maxiter, "ftol": 1e-9, "gtol": 1e-6},
    )
    sf = np.asarray(res.x, float)
    Jf, gf, metf = _eval(sf)
    c6_orbits, _, _ = orbits_and_masks()
    orbit_std = []
    for o in c6_orbits:
        vals = [float(best["s"][i]) for i in o["rods"]]
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
            else "**Free-91 improves J;** compare engineering metrics before adopting."
        ),
        "",
    ]
    (OUT / "FREE91_MULTIMODE_REPORT.md").write_text("\n".join(lines) + "\n")
    print(
        json.dumps(
            {"delta_J": out["delta_J"], "prefer_c6": prefer_c6, "elapsed_s": out["elapsed_s"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
