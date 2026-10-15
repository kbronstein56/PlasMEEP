#!/usr/bin/env python3
"""Phase 1: cheap C6 tied-variable FD spot check of the direct-DtN adjoint."""
from __future__ import annotations

import argparse
import time

import numpy as np

from multimode_opt_core import (
    B_PLUS,
    F_HZ,
    OUT,
    S_HI,
    S_LO,
    dump_json,
    eval_and_grad,
    expand_tied,
    load_start,
    orbits_and_masks,
    objective_from_S,
)
from scipy.sparse.linalg import splu

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
]
from adjoint_direct_dtn import build_system  # noqa: E402
from direct_dtn_s import build_Q_columns, direct_S_from_factor  # noqa: E402

FD_REL_BAR = 1e-3


def fd_orbit(grade, q0, orbits, orbit_idx, h, b_tesla):
    qp = q0.copy()
    qm = q0.copy()
    qp[orbit_idx] = min(S_HI, qp[orbit_idx] + h)
    qm[orbit_idx] = max(S_LO, qm[orbit_idx] - h)
    hp = float(qp[orbit_idx] - q0[orbit_idx])
    hm = float(q0[orbit_idx] - qm[orbit_idx])
    denom = hp + hm
    if denom < 1e-18:
        raise RuntimeError("FD step collapsed")

    def J_of(q):
        s = expand_tied(q, orbits)
        st = build_system(grade, s, F_HZ, b_tesla, "projected")
        S = direct_S_from_factor(
            splu(st["A"].tocsc()),
            len(st["points"]),
            build_Q_columns(st["points"], st["ports"], st["modes"], 2),
        )
        return objective_from_S(S)

    Jp, Jm = J_of(qp), J_of(qm)
    return float((Jp - Jm) / denom), Jp, Jm, denom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_single_uniform_M_Bp0.0500")
    ap.add_argument("--orbits", default="0,1,5")  # center + two six-rod
    ap.add_argument("--hs", default="1e-3,3e-4")
    ap.add_argument("--B", type=float, default=B_PLUS)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    orbits, _, _ = orbits_and_masks()
    q0, s0, tag = load_start(args.design, orbits)
    oids = [int(x) for x in args.orbits.split(",")]
    hs = [float(x) for x in args.hs.split(",")]
    t0 = time.time()
    print("computing adjoint g_q ...", flush=True)
    bun, gq, _ = eval_and_grad(
        args.grade, q0, orbits, F_HZ, args.B, with_absorption=False
    )
    rows = []
    ok = True
    for oid in oids:
        best = None
        sweeps = []
        for h in hs:
            fd, Jp, Jm, denom = fd_orbit(args.grade, q0, orbits, oid, h, args.B)
            g = float(gq[oid])
            rel = abs(fd - g) / max(abs(fd), abs(g), 1e-30)
            rec = {
                "h": h,
                "h_eff": denom,
                "adjoint": g,
                "fd": fd,
                "rel_err": float(rel),
                "Jp": Jp,
                "Jm": Jm,
            }
            sweeps.append(rec)
            if best is None or rel < best["rel_err"]:
                best = rec
            print(
                f"orbit {oid} size={orbits[oid]['size']} h={h:g} "
                f"adj={g:.6e} fd={fd:.6e} rel={rel:.3e}",
                flush=True,
            )
        pass_o = best["rel_err"] <= FD_REL_BAR
        ok = ok and pass_o
        rows.append(
            {
                "orbit": oid,
                "rods": orbits[oid]["rods"],
                "size": orbits[oid]["size"],
                "best": best,
                "h_sweep": sweeps,
                "PASS": pass_o,
            }
        )
    report = {
        "grade": args.grade,
        "design": tag,
        "B": args.B,
        "J0": bun["J"],
        "fd_rel_bar": FD_REL_BAR,
        "orbits": rows,
        "C6_TIED_FD": "PASS" if ok else "FAIL",
        "elapsed_s": time.time() - t0,
    }
    dump_json(OUT / f"c6_tied_fd_{args.grade}.json", report)
    lines = [
        "# C6 tied-variable FD spot check (direct DtN)",
        "",
        f"Grade {args.grade}, design `{tag}`, B={args.B}.",
        f"J0={bun['J']:.6e}. FD rel bar={FD_REL_BAR}.",
        "",
        "| orbit | size | g_adj | best FD | rel | PASS |",
        "|------:|-----:|------:|--------:|----:|:----:|",
    ]
    for r in rows:
        b = r["best"]
        lines.append(
            f"| {r['orbit']} | {r['size']} | {b['adjoint']:.4e} | {b['fd']:.4e} | "
            f"{b['rel_err']:.3e} | {r['PASS']} |"
        )
    lines += ["", f"**C6_TIED_FD: {report['C6_TIED_FD']}**", ""]
    (OUT / "C6_TIED_FD_CHECK.md").write_text("\n".join(lines) + "\n")
    print(json_dumps_safe(report), flush=True)
    return 0 if ok else 1


def json_dumps_safe(report):
    import json

    return json.dumps(
        {
            "C6_TIED_FD": report["C6_TIED_FD"],
            "J0": report["J0"],
            "elapsed_s": report["elapsed_s"],
            "orbits": [
                {
                    "orbit": r["orbit"],
                    "rel": r["best"]["rel_err"],
                    "PASS": r["PASS"],
                }
                for r in report["orbits"]
            ],
        },
        indent=2,
    )


if __name__ == "__main__":
    raise SystemExit(main())
