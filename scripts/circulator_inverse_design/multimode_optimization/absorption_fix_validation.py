#!/usr/bin/env python3
"""Phase 2: validate fixed material_absorption vs discrete identity + power gate.

Cases (as available on disk / cheap to assemble):
  1. homogeneous lossy Voigt (analytic dissipation operator identity)
  2. one-cylinder (if mesh/helpers available; else skip with note)
  3. uniform FULL91 on C
  4. optimized FULL91 (opt_uniform_C) on C

Requires agreement of volume absorption with |1 - P_out - P_mat| ≤ 0.05
and Ax=f discrete closure. Does not loosen the gate.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_optimization"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from adjoint_direct_dtn import build_system  # noqa: E402
from boundary_ports import add_incident_rhs_from_c, modal_c_vector  # noqa: E402
from direct_dtn_s import build_Q_columns, direct_S_from_factor  # noqa: E402
from gyrotropic_tensor import dissipation_matrix, matrix_of, tensor_ordinary  # noqa: E402
from modes import power_scale  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402

BALANCE_BAR = 0.05


def case_voigt_operator_identity():
    """Homogeneous lossy Voigt: D=(ε−ε†)/(2i) must be PSD; elementwise Im fails off-diag."""
    f, fp, gamma, fc = 0.128, 0.4, 0.01, 0.05
    exx, exy, eyx, eyy, _ezz, _eta = tensor_ordinary(f, fp, gamma, fc)
    eps = matrix_of(exx, exy, eyx, eyy)
    D = dissipation_matrix(eps)
    herm = float(np.linalg.norm(D - D.conj().T))
    evals = np.linalg.eigvalsh(0.5 * (D + D.conj().T))
    Im = np.array(
        [[eps[0, 0].imag, eps[0, 1].imag], [eps[1, 0].imag, eps[1, 1].imag]],
        dtype=np.complex128,
    )
    off_diff = float(abs(Im[0, 1] - D[0, 1]) + abs(Im[1, 0] - D[1, 0]))
    return {
        "case": "homogeneous_lossy_Voigt_operator",
        "dissipation_hermitian_resid": herm,
        "dissipation_min_eig": float(np.min(evals)),
        "Im_eps_vs_D_offdiag": off_diff,
        "PASS": bool(herm < 1e-12 and np.min(evals) >= -1e-10 and off_diff > 1e-6),
        "note": "PASS means D is PSD/Hermitian AND differs from elementwise Im off-diag (bug was real)",
    }


def te0_balance(grade: str, s: np.ndarray, label: str):
    t0 = time.time()
    st = build_system(grade, s, 3.85e9, 0.05, "projected")
    points, tris, rho, A, k0 = st["points"], st["tris"], st["rho"], st["A"], st["k0"]
    ports, modes = st["ports"], st["modes"]
    lu = splu(A.tocsc())
    cols = build_Q_columns(points, ports, modes, n_prop=2)
    S = direct_S_from_factor(lu, len(points), cols)
    rows = []
    for j in range(6):
        n = 2 * j
        idx, c = modal_c_vector(points, ports[j], 0)
        scale = power_scale(modes[0]["beta"], k0)
        a_field = 1.0 / max(abs(scale), 1e-30)
        f = np.zeros(len(points), dtype=np.complex128)
        f = add_incident_rhs_from_c(f, idx, c, modes[0]["beta"], a_field)
        x = lu.solve(f)
        P_out = float(np.sum(np.abs(S[:, n]) ** 2))
        P_mat = float(material_absorption(points, tris, x, rho, k0))
        xHAx = complex(np.vdot(x, A @ x))
        xHf = complex(np.vdot(x, f))
        ax_ok = abs(xHAx - xHf) / max(abs(xHf), 1e-30) < 1e-9
        res = abs(1.0 - P_out - P_mat)
        rows.append(
            {
                "port": j,
                "P_out": P_out,
                "P_material": P_mat,
                "residual": res,
                "Ax_eq_f": ax_ok,
                "PASS": bool(res <= BALANCE_BAR and ax_ok),
            }
        )
    max_res = float(max(r["residual"] for r in rows))
    return {
        "case": label,
        "grade": grade,
        "max_residual": max_res,
        "P_out_mean": float(np.mean([r["P_out"] for r in rows])),
        "P_material_mean": float(np.mean([r["P_material"] for r in rows])),
        "Ax_eq_f_all": all(r["Ax_eq_f"] for r in rows),
        "GATE": "PASS" if max_res <= BALANCE_BAR and all(r["Ax_eq_f"] for r in rows) else "FAIL",
        "elapsed_s": time.time() - t0,
        "per_TE0": rows,
        "BALANCE_BAR": BALANCE_BAR,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"cases": []}
    report["cases"].append(case_voigt_operator_identity())
    print(json.dumps(report["cases"][-1], indent=2), flush=True)

    # uniform FULL91
    s_u = np.ones(91)
    u = te0_balance("C", s_u, "uniform_FULL91_C")
    report["cases"].append(u)
    print(
        json.dumps({k: u[k] for k in ("case", "GATE", "max_residual", "P_material_mean", "elapsed_s")}, indent=2),
        flush=True,
    )

    # optimized
    d = json.loads((OUT / "opt_uniform_C.json").read_text())
    s_o = np.asarray(d["best"]["s"], float)
    o = te0_balance("C", s_o, "optimized_FULL91_C")
    report["cases"].append(o)
    print(
        json.dumps({k: o[k] for k in ("case", "GATE", "max_residual", "P_material_mean", "elapsed_s")}, indent=2),
        flush=True,
    )

    report["cases"].append(
        {
            "case": "one_cylinder",
            "PASS": None,
            "note": "Skipped in this pass — FULL91 device is multi-rod; Voigt + uniform + optimized cover the operator bug.",
        }
    )

    report["ABSORPTION_FIX"] = (
        "PASS"
        if all(
            c.get("GATE") == "PASS" or c.get("PASS") is True
            for c in report["cases"]
            if c["case"] != "one_cylinder"
        )
        else "FAIL"
    )
    (OUT / "absorption_fix_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    md = [
        "# Absorption fix validation",
        "",
        "Root cause: `material_absorption` used elementwise `Im(ε)`; corrected to",
        "`D=(ε−ε†)/(2i)`. Gate remains **≤5%**. Frozen DtN S package unchanged.",
        "",
        f"**ABSORPTION_FIX: {report['ABSORPTION_FIX']}**",
        "",
        "| case | GATE / PASS | max residual | P_mat mean |",
        "|------|-------------|-------------:|-----------:|",
    ]
    for c in report["cases"]:
        g = c.get("GATE", c.get("PASS"))
        md.append(
            f"| {c['case']} | {g} | {c.get('max_residual', '—')} | {c.get('P_material_mean', '—')} |"
        )
    (OUT / "ABSORPTION_FIX_VALIDATION.md").write_text("\n".join(md) + "\n")
    print("ABSORPTION_FIX", report["ABSORPTION_FIX"])


if __name__ == "__main__":
    main()
