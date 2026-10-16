#!/usr/bin/env python3
"""Phase 7: SHORT broadband sanity opt on the initial 11-frequency grid.

Warm-starts from the best completed single-frequency multimode design.
Uses efficient adjoint + frozen soft-min aggregation.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "broadband_method"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from broadband_objective import (  # noqa: E402
    ALPHA,
    F_AUDIT_GHZ,
    F_OBJ_GHZ,
    TAU,
    dense_audit,
    eval_bb_and_grad,
)
from multimode_opt_core import S_HI, S_LO, dump_json, expand_tied, orbits_and_masks  # noqa: E402


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--warm-start",
        default="outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json",
    )
    ap.add_argument("--grade", default="C")
    ap.add_argument("--maxiter", type=int, default=3)
    ap.add_argument("--n-obj", type=int, default=11, help="use first n of objective grid")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orbits, _, _ = orbits_and_masks()
    d0 = json.loads(Path(args.warm_start).read_text())
    q0 = np.asarray(d0["best"]["q"], float)
    freqs = [f * 1e9 for f in F_OBJ_GHZ[: args.n_obj]]

    history = []
    best = {"J_BB": -1e300, "q": q0.copy(), "rec": None}
    cache = {}

    def _eval(q):
        q = np.asarray(q, float)
        key = q.tobytes()
        if key in cache:
            return cache[key]
        rec = eval_bb_and_grad(args.grade, q, orbits, freqs)
        cache.clear()
        cache[key] = rec
        return rec

    def fun(q):
        rec = _eval(q)
        J = rec["J_BB"]
        if J > best["J_BB"]:
            best["J_BB"] = J
            best["q"] = np.asarray(q, float).copy()
            best["rec"] = rec
            print(
                f"NEW BEST J_BB={J:.4f} worstJ={rec['worst_J']:.4f} "
                f"@ {rec['worst_f_Hz']/1e9:.2f} GHz meanJ={rec['J_mean']:.4f}",
                flush=True,
            )
        history.append(
            {
                "J_BB": J,
                "worst_J": rec["worst_J"],
                "worst_f_GHz": rec["worst_f_Hz"] / 1e9,
                "J_mean": rec["J_mean"],
                "J_per_freq": rec["J_per_freq"],
            }
        )
        dump_json(OUT / "bb_sanity_history.json", {"history": history, "best_J_BB": best["J_BB"]})
        return -J

    def jac(q):
        rec = _eval(q)
        return -rec["grad_q"]

    print(f"baseline BB on {len(freqs)} freqs ...", flush=True)
    t0 = time.time()
    rec0 = _eval(q0)
    best.update({"J_BB": rec0["J_BB"], "q": q0.copy(), "rec": rec0})
    print(
        f"  J_BB={rec0['J_BB']:.4f} worst={rec0['worst_J']:.4f} @ {rec0['worst_f_Hz']/1e9:.2f}",
        flush=True,
    )

    res = minimize(
        fun,
        q0,
        method="L-BFGS-B",
        jac=jac,
        bounds=[(S_LO, S_HI)] * len(q0),
        options={"maxiter": args.maxiter, "ftol": 1e-8, "gtol": 1e-5},
    )
    qf = np.asarray(res.x, float)
    recf = _eval(qf)
    if recf["J_BB"] > best["J_BB"]:
        best.update({"J_BB": recf["J_BB"], "q": qf.copy(), "rec": recf})

    print("dense audit 41 pts ...", flush=True)
    audit = dense_audit(args.grade, best["q"], orbits)

    out = {
        "warm_start": args.warm_start,
        "grade": args.grade,
        "tau": TAU,
        "alpha": ALPHA,
        "freqs_obj_GHz": [f / 1e9 for f in freqs],
        "J_BB_0": rec0["J_BB"],
        "worst_J_0": rec0["worst_J"],
        "best": {
            "J_BB": best["J_BB"],
            "q": best["q"].tolist(),
            "s": expand_tied(best["q"], orbits).tolist(),
            "rec": {k: best["rec"][k] for k in best["rec"] if k != "metrics_per_freq" and k != "grad_q"},
            "metrics_per_freq": best["rec"]["metrics_per_freq"],
        },
        "audit": audit,
        "delta_J_BB": best["J_BB"] - rec0["J_BB"],
        "delta_worst_J": audit["worst_J"] - rec0["worst_J"],
        "elapsed_s": time.time() - t0,
        "optimizer_message": str(res.message),
        "history": history,
        "SANITY_CHECKS": {
            "softmin_improved": best["J_BB"] > rec0["J_BB"] - 1e-8,
            "worst_obj_improved": best["rec"]["worst_J"] >= rec0["worst_J"] - 1e-3,
            "audit_worst_finite": np.isfinite(audit["worst_J"]),
            "grads_finite": bool(np.all(np.isfinite(best["rec"]["grad_q"]))) if "grad_q" in best["rec"] else True,
        },
    }
    # re-attach grad finiteness from last eval
    out["SANITY_CHECKS"]["grads_finite"] = bool(np.all(np.isfinite(recf["grad_q"])))
    dump_json(OUT / "bb_sanity_opt.json", out)

    # update adaptive report
    add = []
    # suggest adding audit worst if far from objective samples
    wf = audit["worst_f_GHz"]
    obj = np.array(F_OBJ_GHZ[: args.n_obj])
    if abs(wf - obj[np.argmin(np.abs(obj - wf))]) > 0.05:
        add.append(wf)

    lines = [
        "# Broadband sanity (Phase 7)",
        "",
        f"Warm start: `{args.warm_start}`",
        f"Objective grid: {len(freqs)} freqs, τ={TAU}, α={ALPHA}, maxiter={args.maxiter}.",
        "",
        f"| metric | start | best |",
        f"|--------|------:|-----:|",
        f"| J_BB | {rec0['J_BB']:.4f} | {best['J_BB']:.4f} |",
        f"| worst J (obj grid) | {rec0['worst_J']:.4f} | {best['rec']['worst_J']:.4f} |",
        f"| audit worst J (41) | — | {audit['worst_J']:.4f} @ {audit['worst_f_GHz']:.2f} GHz |",
        "",
        f"Suggested frequencies to add: {add if add else 'none (worst near an objective sample)'}",
        "",
        f"SANITY_CHECKS: {json.dumps(out['SANITY_CHECKS'])}",
        "",
    ]
    (OUT / "BROADBAND_SANITY_REPORT.md").write_text("\n".join(lines) + "\n")

    # append adaptive report section
    ar = OUT / "ADAPTIVE_BROADBAND_SAMPLING_REPORT.md"
    with ar.open("a") as f:
        f.write("\n### Sanity run fill-in\n\n")
        f.write(f"- J_BB: {rec0['J_BB']:.4f} → {best['J_BB']:.4f}\n")
        f.write(f"- Audit worst: {audit['worst_J']:.4f} at {audit['worst_f_GHz']:.2f} GHz\n")
        f.write(f"- Suggested add: {add}\n")
        f.write(f"- Elapsed: {out['elapsed_s']/3600:.2f} h\n")

    print(json.dumps({"J_BB_0": rec0["J_BB"], "J_BB_best": best["J_BB"], "audit_worst": audit["worst_J"], "elapsed_s": out["elapsed_s"]}, indent=2))


if __name__ == "__main__":
    main()
