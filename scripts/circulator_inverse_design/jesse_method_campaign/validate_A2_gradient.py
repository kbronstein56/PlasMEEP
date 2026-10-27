#!/usr/bin/env python3
"""FD validation of Objective A2 gradient via existing frequency-domain adjoint.

1) Single-frequency centered FD on one C6 orbit
2) 3-frequency soft-min E2E FD on one orbit

Bar: relative error <= 1e-3 (existing project bar).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

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

from adjoint_direct_dtn import build_system  # noqa: E402
from broadband_objective import softmin_bb  # noqa: E402
from design_io import load_design_q  # noqa: E402
from diff_fem_core import assemble_dA_rod  # noqa: E402
from direct_dtn_s import build_Q_columns  # noqa: E402
from efficient_adjoint import S_from_U, _pack_C, _solve_multi  # noqa: E402
from multimode_opt_core import B_PLUS, F_HZ, expand_tied, orbits_and_masks  # noqa: E402
from objective_A2 import A2_FREEZE, EPS, objective_A2_from_S, scalar_adjoint_rhs_A2  # noqa: E402
from sixfold_orbits import expand_tied as expand_c6  # noqa: E402

FD_REL_BAR = 1e-3


def eval_and_grad_A2(grade, q, orbits, f_hz, b_tesla=B_PLUS, eps=EPS):
    s = expand_c6(np.asarray(q, float), orbits)
    st = build_system(grade, s, f_hz, b_tesla, "projected")
    o_masks = []
    for o in orbits:
        mask = np.zeros(len(st["tris"]), dtype=bool)
        for i in o["rods"]:
            mask |= st["rod_masks"][i]
        o_masks.append(mask)
    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    lu = splu(A.tocsc())
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    S = S_from_U(C, U, betas, te0_cols)
    J = objective_A2_from_S(S, eps=eps)
    R = np.column_stack([scalar_adjoint_rhs_A2(C, betas, S, j, eps=eps) for j in range(6)])
    Lam = _solve_multi(lu, R, trans="T")
    g = np.zeros(len(orbits), float)
    for m, o in enumerate(orbits):
        s_m = float(s[o["rods"][0]])
        dA = assemble_dA_rod(
            points,
            st["tris"],
            o_masks[m],
            st["f_ord"],
            st["fp_ref_ord"],
            st["gamma"],
            st["fc"],
            st["k0"],
            s=s_m,
        )
        acc = 0.0
        for j_loc in range(6):
            dAu = dA @ U[:, j_loc]
            acc += 2.0 * np.real(np.dot(Lam[:, j_loc], dAu))
        g[m] = float(acc / 6.0)
    return {"J": J, "S": S, "timings": None}, g, s


def eval_A2_only(grade, q, orbits, f_hz, b_tesla=B_PLUS, eps=EPS):
    bun, g, s = eval_and_grad_A2(grade, q, orbits, f_hz, b_tesla=b_tesla, eps=eps)
    return bun["J"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_broadband_production_C")
    ap.add_argument("--orbit", type=int, default=1)
    ap.add_argument("--h", type=float, default=3e-4)
    ap.add_argument("--f-ghz", type=float, default=F_HZ / 1e9)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orbits, _, _ = orbits_and_masks()
    # load from production BB json
    p = OUT / f"{args.design}.json"
    if p.exists():
        d = json.loads(p.read_text())
        q0 = np.asarray(d["q_best"], float)
    else:
        q0 = load_design_q(args.design)

    f_hz = args.f_ghz * 1e9
    m = args.orbit
    print(f"SF A2 analytic grad @ {args.f_ghz} GHz orbit {m}...", flush=True)
    t0 = time.perf_counter()
    bun, g, s = eval_and_grad_A2(args.grade, q0, orbits, f_hz)
    print(f"J_A2={bun['J']:.6f} g[{m}]={g[m]:.6e} dt={time.perf_counter()-t0:.1f}s", flush=True)

    qp, qm = q0.copy(), q0.copy()
    qp[m] += args.h
    qm[m] -= args.h
    print("FD +h...", flush=True)
    Jp = eval_A2_only(args.grade, qp, orbits, f_hz)
    print("FD -h...", flush=True)
    Jm = eval_A2_only(args.grade, qm, orbits, f_hz)
    g_fd = (Jp - Jm) / (2 * args.h)
    rel = abs(g_fd - g[m]) / max(abs(g[m]), 1e-12)
    sf = {
        "f_GHz": args.f_ghz,
        "orbit": m,
        "h": args.h,
        "J_A2": bun["J"],
        "g_analytic": float(g[m]),
        "g_fd": float(g_fd),
        "rel_err": float(rel),
        "pass": bool(rel < FD_REL_BAR),
        "A2_FREEZE": A2_FREEZE,
        "FD_REL_BAR": FD_REL_BAR,
    }
    print(f"SF FD rel={rel:.3e} pass={sf['pass']}", flush=True)

    # 3-frequency soft-min E2E
    freqs = [4.0e9, 5.0e9, 5.625e9]  # include the pathological enrichment freq
    print("BB soft-min A2 analytic...", flush=True)
    Js, gs = [], []
    for f in freqs:
        b, gi, _ = eval_and_grad_A2(args.grade, q0, orbits, f)
        Js.append(b["J"])
        gs.append(gi)
    Js = np.asarray(Js, float)
    J_bb, w, *_ = softmin_bb(Js, tau=0.05, alpha=0.25)
    g_bb = sum(wi * gi for wi, gi in zip(w, gs))
    print(f"J_BB={J_bb:.6f} weights={w.tolist()} Js={Js.tolist()}", flush=True)

    print("BB FD +h...", flush=True)
    Jsp = []
    for f in freqs:
        Jsp.append(eval_A2_only(args.grade, qp, orbits, f))
    print("BB FD -h...", flush=True)
    Jsm = []
    for f in freqs:
        Jsm.append(eval_A2_only(args.grade, qm, orbits, f))
    Jbp = softmin_bb(np.asarray(Jsp, float), tau=0.05, alpha=0.25)[0]
    Jbm = softmin_bb(np.asarray(Jsm, float), tau=0.05, alpha=0.25)[0]
    g_bb_fd = (Jbp - Jbm) / (2 * args.h)
    rel_bb = abs(g_bb_fd - g_bb[m]) / max(abs(g_bb[m]), 1e-12)
    bb = {
        "freqs_Hz": freqs,
        "tau": 0.05,
        "alpha": 0.25,
        "J_per_freq": Js.tolist(),
        "J_BB": float(J_bb),
        "weights": w.tolist(),
        "orbit": m,
        "g_analytic": float(g_bb[m]),
        "g_fd": float(g_bb_fd),
        "rel_err": float(rel_bb),
        "pass": bool(rel_bb < FD_REL_BAR),
    }
    print(f"BB FD rel={rel_bb:.3e} pass={bb['pass']}", flush=True)

    report = {
        "single_frequency": sf,
        "softmin_3freq": bb,
        "pass": sf["pass"] and bb["pass"],
    }
    (OUT / "A2_GRADIENT_FD_VALIDATION.json").write_text(json.dumps(report, indent=2) + "\n")
    print("OVERALL", "PASS" if report["pass"] else "FAIL", flush=True)


if __name__ == "__main__":
    main()
