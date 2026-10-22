#!/usr/bin/env python3
"""End-to-end multi-frequency soft-min gradient vs centered FD.

Uses a small 3-frequency set. Requires existing gradient tolerance 1e-3 rel.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

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

from broadband_objective import ALPHA, TAU, eval_bb_and_grad, softmin_bb  # noqa: E402
from multimode_opt_core import B_PLUS, orbits_and_masks  # noqa: E402

FD_REL_BAR = 1e-3


def load_q(design):
    from design_io import load_design_q

    return load_design_q(design)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_uniform_C_warm_M")
    ap.add_argument("--freqs-ghz", default="4.0,5.0,6.0")
    ap.add_argument("--orbits", default="1,5")
    ap.add_argument("--h", type=float, default=3e-4)
    ap.add_argument("--tau", type=float, default=TAU)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orbits, _, _ = orbits_and_masks()
    q0 = load_q(args.design)
    freqs = [float(x) * 1e9 for x in args.freqs_ghz.split(",")]
    orb_ids = [int(x) for x in args.orbits.split(",")]

    # Analytic multi-freq gradient
    print("analytic BB grad...", flush=True)
    bun = eval_bb_and_grad(
        args.grade, q0, orbits, freqs, b_tesla=B_PLUS, tau=args.tau, alpha=args.alpha
    )
    g_an = bun["grad_q"]
    print("J_BB", bun["J_BB"], "weights", bun["weights"], flush=True)

    # Soft-min weight FD (cheap, no Maxwell)
    Jf = np.asarray(bun["J_per_freq"], float)
    _, w, *_ = softmin_bb(Jf, tau=args.tau, alpha=args.alpha)
    w_fd = np.zeros_like(w)
    hJ = 1e-6
    for i in range(len(Jf)):
        Jp = Jf.copy()
        Jm = Jf.copy()
        Jp[i] += hJ
        Jm[i] -= hJ
        Jbp, *_ = softmin_bb(Jp, tau=args.tau, alpha=args.alpha)
        Jbm, *_ = softmin_bb(Jm, tau=args.tau, alpha=args.alpha)
        # softmin_bb returns (J_bb, w, ...)
        Jbp = softmin_bb(Jp, tau=args.tau, alpha=args.alpha)[0]
        Jbm = softmin_bb(Jm, tau=args.tau, alpha=args.alpha)[0]
        w_fd[i] = (Jbp - Jbm) / (2 * hJ)
    w_rel = float(np.linalg.norm(w_fd - w) / max(np.linalg.norm(w), 1e-30))

    fd_rows = []
    for m in orb_ids:
        qp = q0.copy()
        qm = q0.copy()
        qp[m] += args.h
        qm[m] -= args.h
        print(f"FD orbit {m} +h...", flush=True)
        bp = eval_bb_and_grad(
            args.grade, qp, orbits, freqs, b_tesla=B_PLUS, tau=args.tau, alpha=args.alpha
        )
        print(f"FD orbit {m} -h...", flush=True)
        bm = eval_bb_and_grad(
            args.grade, qm, orbits, freqs, b_tesla=B_PLUS, tau=args.tau, alpha=args.alpha
        )
        g_fd = (bp["J_BB"] - bm["J_BB"]) / (2 * args.h)
        g_a = float(g_an[m])
        rel = abs(g_fd - g_a) / max(abs(g_a), 1e-12)
        fd_rows.append(
            {
                "orbit": m,
                "g_analytic": g_a,
                "g_fd": float(g_fd),
                "rel_err": float(rel),
                "pass": rel < FD_REL_BAR,
                "J_BB_plus": bp["J_BB"],
                "J_BB_minus": bm["J_BB"],
            }
        )
        print(f"orbit {m} rel={rel:.3e} pass={rel < FD_REL_BAR}", flush=True)

    report = {
        "grade": args.grade,
        "design": args.design,
        "freqs_Hz": freqs,
        "tau": args.tau,
        "alpha": args.alpha,
        "J_BB": bun["J_BB"],
        "J_per_freq": bun["J_per_freq"],
        "weights": bun["weights"],
        "softmin_weight_fd_rel": w_rel,
        "softmin_weight_fd_pass": w_rel < 1e-8,
        "fd_h": args.h,
        "FD_REL_BAR": FD_REL_BAR,
        "orbit_fd": fd_rows,
        "pass": all(r["pass"] for r in fd_rows) and w_rel < 1e-8,
    }
    (OUT / "softmin_e2e_fd.json").write_text(json.dumps(report, indent=2) + "\n")
    print("PASS" if report["pass"] else "FAIL", flush=True)


if __name__ == "__main__":
    main()
