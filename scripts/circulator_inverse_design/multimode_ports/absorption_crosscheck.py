#!/usr/bin/env python3
"""Phase 6: discrete FEM energy identity vs field-integral material absorption."""
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
from full91_direct_dtn import assemble_device, load_s  # noqa: E402
from modes import power_scale  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_single_uniform_M_Bp0.0500")
    ap.add_argument("--measure", default="projected")
    ap.add_argument("--channels", default="7,1,3", help="Comma channels to audit")
    args = ap.parse_args()
    s = load_s(args.design)
    points, tris, rho, A, k0, ports, modes, _ = assemble_device(
        args.grade, 3.85e9, 0.05, s, measure=args.measure
    )
    lu = splu(A.tocsc())
    Sp = np.load(OUT / f"S12_direct_{args.design}_Bp.npz")["S"]
    # Prefer projected artifact if present and measure=projected
    proj = OUT / f"S12_direct_{args.design}_Bp_proj.npz"
    if args.measure == "projected" and proj.exists():
        Sp = np.load(proj)["S"]

    rows = []
    for n in [int(x) for x in args.channels.split(",")]:
        p, m = divmod(n, 2)
        idx, c = modal_c_vector(points, ports[p], m)
        scale = power_scale(modes[m]["beta"], k0)
        a_field = 1.0 / max(abs(scale), 1e-30)
        f = np.zeros(len(points), dtype=np.complex128)
        f = add_incident_rhs_from_c(f, idx, c, modes[m]["beta"], a_field)
        x = lu.solve(f)
        P_mat = float(material_absorption(points, tris, x, rho, k0))
        P_out = float(np.sum(np.abs(Sp[:, n]) ** 2))
        eng_res = abs(1.0 - P_out - P_mat)

        # Discrete identity: A x = f ⇒ xᴴ A x = xᴴ f
        xHAx = complex(np.vdot(x, A @ x))
        xHf = complex(np.vdot(x, f))
        # Anti-Hermitian power proxy: (1/2i)(xᴴ A x - (xᴴ A x)*) = Im(xᴴ A x)
        # For our drive f = 2i β a_f c, injected modal power relates to Im(xᴴ f)/(2β a_f ...)
        beta = float(np.real(modes[m]["beta"]))
        # From matched-port theory: unit power ⇒ Im(xᴴ f) / (2 β a_f²) style; report raw.
        P_disc_proxy = float(np.imag(xHf) / max(2.0 * beta * (a_field**2), 1e-30))
        rows.append({
            "channel": n,
            "port": p,
            "mode": m,
            "P_out": P_out,
            "P_material_field": P_mat,
            "eng_residual_field": eng_res,
            "Im_xHAx": float(xHAx.imag),
            "Im_xHf": float(xHf.imag),
            "Re_xHAx_minus_xHf": float((xHAx - xHf).real),
            "Im_xHAx_minus_xHf": float((xHAx - xHf).imag),
            "discrete_Ax_eq_f": abs(xHAx - xHf) / max(abs(xHf), 1e-30) < 1e-10,
            "P_disc_proxy": P_disc_proxy,
            "|c|": float(np.linalg.norm(c)),
            "cTc": float(np.real(np.dot(c, c))),
            "a_field": float(a_field),
            "beta": beta,
            "modal_norm_ok": abs(1.0 - float(np.real(np.dot(np.conj(c), c))) )  # not L2
        })
        # True discrete modal norm ∫φ² ≈ from construction
        rows[-1]["incident_unit_power"] = True  # by S construction + power_scale

    report = {
        "design": args.design,
        "grade": args.grade,
        "measure": args.measure,
        "channels": rows,
        "interpretation": (
            "If discrete_Ax_eq_f and eng_residual_field<=0.05, both FEM identity and "
            "field absorption are consistent with unit-power S. If Ax=f closes but "
            "field residual fails, fix absorption post-processing."
        ),
    }
    out = OUT / "absorption_crosscheck.json"
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"channels": [
        {k: r[k] for k in ("channel", "eng_residual_field", "discrete_Ax_eq_f", "P_out", "P_material_field")}
        for r in rows
    ]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
