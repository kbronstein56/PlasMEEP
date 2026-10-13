#!/usr/bin/env python3
"""Phase 5: TE1 power-balance residual on C vs M (projected measure)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_ports"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from boundary_ports import add_incident_rhs_from_c, modal_c_vector  # noqa: E402
from direct_dtn_s import build_Q_columns, direct_S_from_factor  # noqa: E402
from full91_direct_dtn import assemble_device, load_s  # noqa: E402
from modes import power_scale  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402

BAR = 0.05


def power_all(grade: str, design: str, measure: str = "projected"):
    t0 = time.time()
    s = load_s(design)
    points, tris, rho, A, k0, ports, modes, n_ext = assemble_device(
        grade, 3.85e9, 0.05, s, measure=measure
    )
    lu = splu(A.tocsc())
    cols = build_Q_columns(points, ports, modes, n_prop=2)
    S = direct_S_from_factor(lu, len(points), cols)
    rows = []
    for n in range(12):
        p, m = divmod(n, 2)
        idx, c = modal_c_vector(points, ports[p], m)
        scale = power_scale(modes[m]["beta"], k0)
        a_field = 1.0 / max(abs(scale), 1e-30)
        b = np.zeros(len(points), dtype=np.complex128)
        b = add_incident_rhs_from_c(b, idx, c, modes[m]["beta"], a_field)
        x = lu.solve(b)
        P_mat = float(material_absorption(points, tris, x, rho, k0))
        P_out = float(np.sum(np.abs(S[:, n]) ** 2))
        res = abs(1.0 - P_out - P_mat)
        rows.append({
            "channel": n,
            "port": p,
            "mode": m,
            "P_out": P_out,
            "P_material": P_mat,
            "residual": res,
            "PASS": bool(res <= BAR),
        })
    np.savez_compressed(OUT / f"S12_direct_{design}_Bp_{grade}.npz", S=S)
    return {
        "grade": grade,
        "n_nodes": int(len(points)),
        "n_exterior": n_ext,
        "elapsed_s": time.time() - t0,
        "channels": rows,
        "TE1_residuals": [r["residual"] for r in rows if r["mode"] == 1],
        "max_residual": max(r["residual"] for r in rows),
        "ch7_residual": rows[7]["residual"],
        "MULTIMODE_POWER_BALANCE": "PASS" if all(r["PASS"] for r in rows) else "FAIL",
    }


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="opt_single_uniform_M_Bp0.0500")
    ap.add_argument("--grades", default="C,M")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"design": args.design, "measure": "projected", "bar": BAR, "grades": {}}
    for g in args.grades.split(","):
        g = g.strip()
        print(f"=== grade {g} ===", flush=True)
        report["grades"][g] = power_all(g, args.design)
        print(json.dumps({
            "grade": g,
            "ch7": report["grades"][g]["ch7_residual"],
            "max": report["grades"][g]["max_residual"],
            "status": report["grades"][g]["MULTIMODE_POWER_BALANCE"],
            "elapsed_s": report["grades"][g]["elapsed_s"],
        }, indent=2), flush=True)
    (OUT / "te1_mesh_convergence.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# TE1 mesh convergence (projected measure)",
        "",
        f"Design `{args.design}`, bar {BAR}.",
        "",
        "| grade | nodes | ch7 residual | max residual | status |",
        "|-------|------:|-------------:|-------------:|--------|",
    ]
    for g, r in report["grades"].items():
        lines.append(
            f"| {g} | {r['n_nodes']} | {r['ch7_residual']:.4f} | "
            f"{r['max_residual']:.4f} | {r['MULTIMODE_POWER_BALANCE']} |"
        )
    lines += [
        "",
        "If residuals decrease C→M, staircase/measure error is converging.",
        "F skipped unless needed (memory).",
    ]
    (OUT / "TE1_MESH_CONVERGENCE.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
