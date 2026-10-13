#!/usr/bin/env python3
"""Phase 1–2: localize TE1 power-balance failure + staircase measure audit."""
from __future__ import annotations

import json
import sys
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
from feed_constants import clear_width  # noqa: E402
from full91_direct_dtn import assemble_device, load_s  # noqa: E402
from modes import power_scale, transverse_mode  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402

BAR = 0.05
W = float(clear_width)


def port_measure_stats(points, port) -> dict:
    u = port["outward"] / np.linalg.norm(port["outward"])
    t = port["tangent"] / np.linalg.norm(port["tangent"])
    c = port["center"]
    ds_sum = 0.0
    dproj_sum = 0.0
    alongs = []
    for i, j in port["edges"]:
        evec = points[j] - points[i]
        ds = float(np.linalg.norm(evec))
        dproj = float(abs(evec @ t))
        ds_sum += ds
        dproj_sum += dproj
        mid = 0.5 * (points[i] + points[j])
        alongs.append(float((mid - c) @ u))
    alongs = np.asarray(alongs, float) if alongs else np.array([0.0])
    # Modal norm integrals with ds and projected measure
    nrm2_ds = 0.0
    nrm2_proj = 0.0
    for m in (0, 1):
        raw = {}
        for i, j in port["edges"]:
            for n in (i, j):
                if n not in raw:
                    y = float((points[n] - c) @ t)
                    raw[n] = float(transverse_mode(np.array([y]), m, W)[0])
        n2ds = n2pr = 0.0
        for i, j in port["edges"]:
            evec = points[j] - points[i]
            L = float(np.linalg.norm(evec))
            Lp = float(abs(evec @ t))
            a, b = raw[i], raw[j]
            n2ds += (L / 3.0) * (a * a + a * b + b * b)
            n2pr += (Lp / 3.0) * (a * a + a * b + b * b)
        if m == 1:
            nrm2_ds, nrm2_proj = n2ds, n2pr
    idx0, c0 = modal_c_vector(points, port, 0)
    idx1, c1 = modal_c_vector(points, port, 1)
    return {
        "n_edges": int(port["n_edges"]),
        "n_nodes": int(port["n_nodes"]),
        "sum_ds": float(ds_sum),
        "sum_dproj": float(dproj_sum),
        "W": W,
        "y_span": float(port.get("y_span", 0.0)),
        "along_mean": float(np.mean(alongs)),
        "along_std": float(np.std(alongs)),
        "along_min": float(np.min(alongs)),
        "along_max": float(np.max(alongs)),
        "ds_over_W": float(ds_sum / W),
        "dproj_over_W": float(dproj_sum / W),
        "TE1_nrm2_ds": float(nrm2_ds),
        "TE1_nrm2_proj": float(nrm2_proj),
        "c0_dot_c0": float(np.real(np.dot(c0, c0))) if len(c0) else 0.0,
        "c1_dot_c1": float(np.real(np.dot(c1, c1))) if len(c1) else 0.0,
        "|c1|": float(np.linalg.norm(c1)) if len(c1) else 0.0,
    }


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_single_uniform_M_Bp0.0500")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    s = load_s(args.design)
    points, tris, rho, A, k0, ports, modes, _ = assemble_device(args.grade, 3.85e9, 0.05, s)
    lu = splu(A.tocsc())
    Sp = np.load(OUT / f"S12_direct_{args.design}_Bp.npz")["S"]

    geo = [port_measure_stats(points, p) for p in ports]
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
        P_out = float(np.sum(np.abs(Sp[:, n]) ** 2))
        P_te0 = float(sum(abs(Sp[2 * i, n]) ** 2 for i in range(6)))
        P_te1 = float(sum(abs(Sp[2 * i + 1, n]) ** 2 for i in range(6)))
        residual = abs(1.0 - P_out - P_mat)
        # Discrete anti-Hermitian power identity: for A x = f,
        #  Im(x^H A x) related to sources/loss. Use 0.5 Im(x^H (A-A^H) x) / norm.
        Ax = A @ x
        xHAx = complex(np.vdot(x, Ax))
        # Input power proxy from RHS: for our drive f = 2i β a_f c,
        #  Re or Im of x^H f related to injected power.
        f = b
        xHf = complex(np.vdot(x, f))
        rows.append({
            "channel": n,
            "port": p,
            "mode": m,
            "P_out": P_out,
            "P_TE0": P_te0,
            "P_TE1": P_te1,
            "P_material": P_mat,
            "residual": residual,
            "PASS": bool(residual <= BAR),
            "Re_xHAx": float(xHAx.real),
            "Im_xHAx": float(xHAx.imag),
            "Re_xHf": float(xHf.real),
            "Im_xHf": float(xHf.imag),
            "geometry": geo[p],
        })

    fail = [r for r in rows if not r["PASS"]]
    report = {
        "design": args.design,
        "grade": args.grade,
        "bar": BAR,
        "port_geometry": geo,
        "channels": rows,
        "failing": fail,
        "TE1_residuals": [r["residual"] for r in rows if r["mode"] == 1],
        "TE0_residuals": [r["residual"] for r in rows if r["mode"] == 0],
    }
    (OUT / "te1_power_balance_diagnostic.json").write_text(json.dumps(report, indent=2) + "\n")

    # Markdown summary
    lines = [
        "# TE1 power balance diagnostic",
        "",
        f"Design `{args.design}`, grade {args.grade}, bar {BAR}.",
        "",
        "## Per-channel residuals",
        "",
        "| ch | port | mode | Σ|S|² | TE0 | TE1 | P_mat | residual | PASS |",
        "|----|------|------|-------|-----|-----|-------|----------|------|",
    ]
    for r in rows:
        lines.append(
            f"| {r['channel']} | {r['port']} | {r['mode']} | {r['P_out']:.4f} | "
            f"{r['P_TE0']:.4f} | {r['P_TE1']:.4f} | {r['P_material']:.4f} | "
            f"{r['residual']:.4f} | {r['PASS']} |"
        )
    lines += ["", "## Port geometry (staircase audit)", "",
              "| port | n_edges | Σds | Σ|e·t| | Σds/W | Σproj/W | y_span | along_std | TE1 nrm²_ds | TE1 nrm²_proj | |c1| |",
              "|------|---------|-----|--------|-------|---------|--------|-----------|-------------|---------------|-----|"]
    for i, g in enumerate(geo):
        lines.append(
            f"| {i} | {g['n_edges']} | {g['sum_ds']:.4f} | {g['sum_dproj']:.4f} | "
            f"{g['ds_over_W']:.4f} | {g['dproj_over_W']:.4f} | {g['y_span']:.4f} | "
            f"{g['along_std']:.5f} | {g['TE1_nrm2_ds']:.4f} | {g['TE1_nrm2_proj']:.4f} | {g['|c1|']:.5f} |"
        )
    if fail:
        f0 = fail[0]
        lines += [
            "",
            f"## Failing channel: ch={f0['channel']} port={f0['port']} TE{f0['mode']}",
            "",
            f"- residual = **{f0['residual']:.4f}** (bar {BAR})",
            f"- Σ|S|² = {f0['P_out']:.4f} (exceeds 1 by {f0['P_out']-1:.4f})",
            f"- P_material = {f0['P_material']:.4f}",
            f"- Im(xᴴ A x) = {f0['Im_xHAx']:.6e}, Im(xᴴ f) = {f0['Im_xHf']:.6e}",
            "",
            "Compare this port's geometry row above to the other five ports.",
        ]
    (OUT / "TE1_POWER_BALANCE_DIAGNOSTIC.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({
        "failing": [{"ch": f["channel"], "port": f["port"], "residual": f["residual"]} for f in fail],
        "TE1_residuals": report["TE1_residuals"],
        "max_ds_over_W": max(g["ds_over_W"] for g in geo),
        "max_dproj_over_W": max(g["dproj_over_W"] for g in geo),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
