#!/usr/bin/env python3
"""Benchmark forward-mode (direct sensitivity) vs efficient adjoint for C6.

After one LU:
  adjoint: 6 TE0 fwd + 6 scalar adj RHS + dA contractions over 16 orbits
  forward: 6 TE0 fwd + (N_param × 6) sensitivity RHS + contractions

Reports wall time, LU/RHS counts, memory proxy, gradient agreement.
Does NOT change production.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
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
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from adjoint_direct_dtn import channel_weights, objective_from_S  # noqa: E402
from diff_fem_core import assemble_dA_rod  # noqa: E402
from efficient_adjoint import (  # noqa: E402
    S_from_U,
    _pack_C,
    _solve_multi,
    eval_and_grad_efficient,
    scalar_adjoint_rhs,
)
from multimode_opt_core import (  # noqa: E402
    B_PLUS,
    F_HZ,
    expand_tied,
    load_start,
    orbits_and_masks,
    te0_metrics_from_S,
)


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def forward_mode_grad(st, lu, C, U_te0, betas, S, orbits, o_masks, s):
    """Direct sensitivity gradient for C6 orbits (scalar J)."""
    te0_cols = [2 * j for j in range(6)]
    n_orb = len(orbits)
    g = np.zeros(n_orb, float)
    n_sens_rhs = 0
    t_rhs = 0.0
    t_dA = 0.0
    for m, o in enumerate(orbits):
        s_m = float(s[o["rods"][0]])
        t0 = time.perf_counter()
        dA = assemble_dA_rod(
            st["points"],
            st["tris"],
            o_masks[m],
            st["f_ord"],
            st["fp_ref_ord"],
            st["gamma"],
            st["fc"],
            st["k0"],
            s=s_m,
        )
        t_dA += time.perf_counter() - t0
        acc = 0.0
        for j_loc, j_src in enumerate(range(6)):
            # A dx = -dA U
            rhs = -(dA @ U_te0[:, j_loc])
            t0 = time.perf_counter()
            dU = lu.solve(rhs)
            t_rhs += time.perf_counter() - t0
            n_sens_rhs += 1
            # dS_mn = 2i √(βm βn) C_m · dU  (bilinear)
            n = te0_cols[j_loc]
            w = channel_weights(j_src)
            dJj = 0.0
            for ch in range(12):
                if abs(w[ch]) < 1e-30:
                    continue
                sbeta = np.sqrt(betas[ch] * betas[n])
                dS = 2j * sbeta * complex(np.dot(C[:, ch], dU))
                dJj += w[ch] * 2.0 * np.real(np.conj(S[ch, n]) * dS)
            acc += dJj
        g[m] = float(acc / 6.0)
    return g, {"n_sens_rhs": n_sens_rhs, "t_sens_rhs": t_rhs, "t_dA": t_dA}


def estimate_scaling(t_lu, t_fwd6, t_adj6, t_dA_adj, t_sens_per_param, t_dA_fwd):
    """Rough wall-time model: shared LU + RHS + dA."""
    rows = []
    for n_param, label in [(16, "C6"), (31, "C3"), (91, "free91")]:
        t_adj = t_lu + t_fwd6 + t_adj6 + t_dA_adj * (n_param / 16.0)
        t_fwd = t_lu + t_fwd6 + t_sens_per_param * n_param + t_dA_fwd * (n_param / 16.0)
        rows.append(
            {
                "label": label,
                "n_param": n_param,
                "est_adjoint_s": t_adj,
                "est_forward_s": t_fwd,
                "preferred": "adjoint" if t_adj <= t_fwd else "forward",
                "ratio_fwd_over_adj": t_fwd / max(t_adj, 1e-9),
            }
        )
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_uniform_C_warm_M")
    ap.add_argument("--skip-forward", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orbits, _, _ = orbits_and_masks()
    # load_start expects tags; fall back to JSON
    design_json = (
        ROOT
        / "outputs"
        / "circulator_inverse_design"
        / "multimode_optimization"
        / f"{args.design}.json"
    )
    if design_json.exists():
        d = json.loads(design_json.read_text())
        q = np.asarray(d.get("q_best", d.get("q", d.get("q0"))), float)
    else:
        q, _, _ = load_start(args.design, orbits)

    print("RSS start MB", rss_mb(), flush=True)
    t_wall0 = time.perf_counter()
    bun, g_adj, s = eval_and_grad_efficient(
        args.grade, q, orbits, f_hz=F_HZ, b_tesla=B_PLUS, profile=True
    )
    t_adj_total = time.perf_counter() - t_wall0
    tim = bun["timings"]
    print("adjoint done", t_adj_total, "J", bun["J"], flush=True)

    report = {
        "grade": args.grade,
        "design": args.design,
        "f_Hz": F_HZ,
        "B": B_PLUS,
        "n_param": len(orbits),
        "J": bun["J"],
        "adjoint": {
            "wall_s": t_adj_total,
            "timings": tim,
            "n_LU": 1,
            "n_forward_rhs": 6,
            "n_adjoint_rhs": 6,
            "rss_mb_after": rss_mb(),
            "grad_norm": float(np.linalg.norm(g_adj)),
        },
    }

    if not args.skip_forward:
        st = bun["st"]
        _, _, o_masks = orbits_and_masks(rod_masks=st["rod_masks"])
        t0 = time.perf_counter()
        g_fwd, fwd_prof = forward_mode_grad(
            st,
            bun["lu"],
            bun["C"],
            bun["U_te0"],
            bun["betas"],
            bun["S"],
            orbits,
            o_masks,
            s,
        )
        t_fwd_only = time.perf_counter() - t0
        # Total forward-mode path ≈ assemble+LU+fwd already paid + sensitivity
        t_fwd_full = (
            tim["assemble_s"]
            + tim["factor_lu"]
            + tim["forward_te0_rhs"]
            + tim["S_and_objective"]
            + t_fwd_only
        )
        rel = float(np.linalg.norm(g_fwd - g_adj) / max(np.linalg.norm(g_adj), 1e-30))
        report["forward"] = {
            "wall_sensitivity_only_s": t_fwd_only,
            "wall_full_est_s": t_fwd_full,
            "n_LU": 1,
            "n_forward_field_rhs": 6,
            "n_sensitivity_rhs": fwd_prof["n_sens_rhs"],
            "timings_extra": fwd_prof,
            "rss_mb_after": rss_mb(),
            "grad_norm": float(np.linalg.norm(g_fwd)),
            "rel_err_vs_adjoint": rel,
            "agree_bar": 1e-3,
            "pass": rel < 1e-3,
        }
        t_sens_per = fwd_prof["t_sens_rhs"] / max(len(orbits), 1)
        report["scaling_estimate"] = estimate_scaling(
            tim["factor_lu"],
            tim["forward_te0_rhs"],
            tim["adjoint_te0_rhs"],
            tim["dA_contraction"],
            t_sens_per,
            fwd_prof["t_dA"],
        )
        print("forward rel_err", rel, "pass", rel < 1e-3, flush=True)

    (OUT / "forward_mode_bench.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in report if k != "adjoint"}, indent=2)[:2000])
    print("wrote", OUT / "forward_mode_bench.json", flush=True)


if __name__ == "__main__":
    main()
