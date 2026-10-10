#!/usr/bin/env python3
"""Optional 91-independent refinement from a sixfold warm start."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from circulator_core import B_PLUS, F_HZ, S_HI, S_LO, gradient_s, load_device, solve_six  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"


def main():
    opt = json.loads(Path(sys.argv[1]).read_text())
    grade = sys.argv[2] if len(sys.argv) > 2 else "M"
    maxiter = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    dev = load_device(grade)
    s0 = np.asarray(opt["best"]["s"], float)
    P_ref = float(opt["P_ref"])
    history = []
    best = {"J": -1e300, "s": None, "metrics": None, "power": None}
    cache = {}

    def fun(s):
        s = np.asarray(s, float)
        bun = solve_six(dev, s, F_HZ, B_PLUS, P_ref=P_ref)
        g = gradient_s(dev, bun, None)
        cache["g"] = g
        cache["s"] = s.copy()
        history.append({
            "J": bun.J,
            "desired_avg": bun.metrics["desired_avg"],
            "desired_min": bun.metrics["desired_min"],
            "isolation_avg_dB": bun.metrics["isolation_avg_dB"],
            "grad_norm": float(np.linalg.norm(g)),
            "s_min": float(s.min()),
            "s_max": float(s.max()),
        })
        if bun.J > best["J"]:
            best.update({"J": bun.J, "s": s.tolist(), "metrics": bun.metrics, "power": bun.power.tolist()})
        print(f"free91 J={bun.J:.6e} des_avg={bun.metrics['desired_avg']:.4f} "
              f"des_min={bun.metrics['desired_min']:.4f} |g|={history[-1]['grad_norm']:.3e}", flush=True)
        return -bun.J

    def jac(s):
        s = np.asarray(s, float)
        if cache.get("s") is not None and np.allclose(s, cache["s"]):
            return -cache["g"]
        bun = solve_six(dev, s, F_HZ, B_PLUS, P_ref=P_ref)
        g = gradient_s(dev, bun, None)
        cache["g"] = g
        cache["s"] = s.copy()
        return -g

    t0 = time.perf_counter()
    res = minimize(
        fun, s0, method="L-BFGS-B", jac=jac,
        bounds=[(S_LO, S_HI)] * len(s0),
        options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-8},
    )
    out = {
        "warm_J": opt["best"]["J"],
        "grade": grade,
        "P_ref": P_ref,
        "success": bool(res.success),
        "message": str(res.message),
        "nit": int(res.nit),
        "wall_s": time.perf_counter() - t0,
        "history": history,
        "best": best,
        "final_s": np.asarray(res.x, float).tolist(),
        # rotational inequivalence: std of desired across sources
        "desired_std_best": float(np.std([r["desired"] for r in best["metrics"]["per_source"]])) if best["metrics"] else None,
    }
    path = OUT / f"opt_free91_{grade}_Bp0.0500.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE", path.name, "best", best["J"], flush=True)


if __name__ == "__main__":
    main()
