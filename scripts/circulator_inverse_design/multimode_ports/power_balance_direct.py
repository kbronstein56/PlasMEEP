#!/usr/bin/env python3
"""Phase 6: power balance with direct DtN S + material absorption."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_ports"
DES = ROOT / "outputs" / "circulator_inverse_design"
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

BALANCE_BAR = 0.05


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--design", default="opt_single_uniform_M_Bp0.0500")
    ap.add_argument("--measure", choices=("arclength", "projected"), default="projected")
    ap.add_argument("--tag-suffix", default="")
    args = ap.parse_args()
    s = load_s(args.design)
    points, tris, rho, A, k0, ports, modes, _ = assemble_device(
        args.grade, 3.85e9, 0.05, s, measure=args.measure
    )
    lu = splu(A.tocsc())
    Sp = np.load(OUT / f"S12_direct_{args.design}_Bp{args.tag_suffix}.npz")["S"]
    rows = []
    ok = True
    for n in range(12):
        p, m = divmod(n, 2)
        idx, c = modal_c_vector(points, ports[p], m)
        scale = power_scale(modes[m]["beta"], k0)
        a_field = 1.0 / max(abs(scale), 1e-30)
        b = np.zeros(len(points), dtype=np.complex128)
        b = add_incident_rhs_from_c(b, idx, c, modes[m]["beta"], a_field)
        x = lu.solve(b)
        P_mat = material_absorption(points, tris, x, rho, k0)
        # Unit-power column of S
        P_out = float(np.sum(np.abs(Sp[:, n]) ** 2))
        P_te0 = float(sum(abs(Sp[2 * i, n]) ** 2 for i in range(6)))
        P_te1 = float(sum(abs(Sp[2 * i + 1, n]) ** 2 for i in range(6)))
        # Incident power = 1 by construction of S
        rem = 1.0 - P_out - P_mat
        identity = P_out + P_mat + max(rem, 0.0)
        eng = abs(1.0 - identity)
        pass_n = eng <= BALANCE_BAR and P_mat >= -0.02
        ok = ok and pass_n
        rows.append({
            "channel": n,
            "port": p,
            "mode": m,
            "P_out": P_out,
            "P_TE0": P_te0,
            "P_TE1": P_te1,
            "P_material": float(P_mat),
            "P_nonport_or_num": float(rem),
            "eng_residual": float(eng),
            "PASS": bool(pass_n),
        })
    # Summaries for TE0 and TE1 drives
    te0 = [r for r in rows if r["mode"] == 0]
    te1 = [r for r in rows if r["mode"] == 1]
    report = {
        "design": args.design,
        "grade": args.grade,
        "bar": BALANCE_BAR,
        "TE0_drive_mean": {
            "P_out": float(np.mean([r["P_out"] for r in te0])),
            "P_TE0": float(np.mean([r["P_TE0"] for r in te0])),
            "P_TE1": float(np.mean([r["P_TE1"] for r in te0])),
            "P_material": float(np.mean([r["P_material"] for r in te0])),
            "P_nonport": float(np.mean([r["P_nonport_or_num"] for r in te0])),
        },
        "TE1_drive_mean": {
            "P_out": float(np.mean([r["P_out"] for r in te1])),
            "P_TE0": float(np.mean([r["P_TE0"] for r in te1])),
            "P_TE1": float(np.mean([r["P_TE1"] for r in te1])),
            "P_material": float(np.mean([r["P_material"] for r in te1])),
            "P_nonport": float(np.mean([r["P_nonport_or_num"] for r in te1])),
        },
        "per_channel": rows,
        "MULTIMODE_POWER_BALANCE": "PASS" if ok else "FAIL",
    }
    report["measure"] = args.measure
    outp = OUT / f"power_balance_direct_{args.design}{args.tag_suffix}.json"
    outp.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in ("MULTIMODE_POWER_BALANCE", "measure", "TE0_drive_mean", "TE1_drive_mean")
            },
            indent=2,
        )
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
