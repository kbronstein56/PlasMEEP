#!/usr/bin/env python3
"""Phase 1: discrete FEM energy identity vs volume P_material for optimized design.

Does NOT modify the frozen DtN S package. Compares:
  A) Σ|S|^2 (propagating outgoing modal power)
  B) material_absorption (current multimode post-processing)
  C) volume_absorption with (ε−ε†)/(2i) (validated FULL91 control-volume form)
  D) discrete Im(xᴴ A_vol x) / Im(xᴴ f) identities from A x = f
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
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from modes import power_scale  # noqa: E402
from validate_modal_s import material_absorption  # noqa: E402


def eps_from_rho(rho):
    rxx, rxy, ryx, ryy = rho
    det = rxx * ryy - rxy * ryx
    det_safe = np.where(np.abs(det) > 1e-30, det, 1.0)
    exx = ryy / det_safe
    exy = -rxy / det_safe
    eyx = -ryx / det_safe
    eyy = rxx / det_safe
    return exx, exy, eyx, eyy


def volume_absorption_correct(points, tris, x, rho, k0) -> float:
    """∫ (ω/2) E† ((ε−ε†)/(2i)) E with E=(i/ω) ρ curl H — gyrotropic-safe."""
    omega = float(k0)
    px = points[tris, 0]
    py = points[tris, 1]
    twice = (px[:, 1] - px[:, 0]) * (py[:, 2] - py[:, 0]) - (px[:, 2] - px[:, 0]) * (
        py[:, 1] - py[:, 0]
    )
    area = 0.5 * np.abs(twice)
    good = area > 1e-18
    bx = np.stack(
        [py[:, 1] - py[:, 2], py[:, 2] - py[:, 0], py[:, 0] - py[:, 1]], axis=1
    ) / twice[:, None]
    by = np.stack(
        [px[:, 2] - px[:, 1], px[:, 0] - px[:, 2], px[:, 1] - px[:, 0]], axis=1
    ) / twice[:, None]
    un = x[tris]
    gx = np.sum(bx * un, axis=1)
    gy = np.sum(by * un, axis=1)
    # curl H = (∂y H, −∂x H)
    cx = gy
    cy = -gx
    rxx, rxy, ryx, ryy = rho
    Ex = (1j / omega) * (rxx * cx + rxy * cy)
    Ey = (1j / omega) * (ryx * cx + ryy * cy)
    exx, exy, eyx, eyy = eps_from_rho(rho)
    # D = (ε − ε†)/(2i)
    # ε†_ij = conj(ε_ji)
    dxx = (exx - np.conj(exx)) / (2j)  # = Im(exx) for scalar diagonal
    # For off-diagonal: D_xy = (ε_xy - conj(ε_yx))/(2i)
    dxy = (exy - np.conj(eyx)) / (2j)
    dyx = (eyx - np.conj(exy)) / (2j)
    dyy = (eyy - np.conj(eyy)) / (2j)
    # Only elements with nontrivial dissipation
    lossy = good & (
        (np.abs(dxx) > 1e-18)
        | (np.abs(dyy) > 1e-18)
        | (np.abs(dxy) > 1e-18)
        | (np.abs(dyx) > 1e-18)
    )
    if not np.any(lossy):
        return 0.0
    Ex_l, Ey_l = Ex[lossy], Ey[lossy]
    Ix = dxx[lossy] * Ex_l + dxy[lossy] * Ey_l
    Iy = dyx[lossy] * Ex_l + dyy[lossy] * Ey_l
    dens = np.real(np.conj(Ex_l) * Ix + np.conj(Ey_l) * Iy)
    return float(np.sum(area[lossy] * (omega / 2.0) * dens))


def volume_absorption_im_eps(points, tris, x, rho, k0) -> float:
    """Current multimode formula: E† Im(ε) E (elementwise Im)."""
    return float(material_absorption(points, tris, x, rho, k0))


def diagnose_grade(grade: str, s_vec: np.ndarray, f_hz: float = 3.85e9, b_tesla: float = 0.05):
    t0 = time.time()
    st = build_system(grade, s_vec, f_hz, b_tesla, measure="projected")
    points, tris, rho, A, k0 = st["points"], st["tris"], st["rho"], st["A"], st["k0"]
    ports, modes = st["ports"], st["modes"]
    # Volume-only operator (no DtN Robin) — includes PML stretch + materials
    A_vol = assemble_anisotropic(
        points, tris, *rho, k0, np.zeros(len(points), dtype=bool)
    )
    lu = splu(A.tocsc())
    cols = build_Q_columns(points, ports, modes, n_prop=2)
    S = direct_S_from_factor(lu, len(points), cols)

    rows = []
    for j in range(6):  # TE0 sources only for primary diagnostic
        n = 2 * j
        p, m = j, 0
        idx, c = modal_c_vector(points, ports[p], m)
        scale = power_scale(modes[m]["beta"], k0)
        a_field = 1.0 / max(abs(scale), 1e-30)
        f = np.zeros(len(points), dtype=np.complex128)
        f = add_incident_rhs_from_c(f, idx, c, modes[m]["beta"], a_field)
        x = lu.solve(f)

        P_out = float(np.sum(np.abs(S[:, n]) ** 2))
        P_im = volume_absorption_im_eps(points, tris, x, rho, k0)
        P_corr = volume_absorption_correct(points, tris, x, rho, k0)

        xHAx = complex(np.vdot(x, A @ x))
        xHf = complex(np.vdot(x, f))
        xHAvolx = complex(np.vdot(x, A_vol @ x))
        # Port-only contribution
        A_port_x = A @ x - A_vol @ x
        xHAportx = complex(np.vdot(x, A_port_x))

        # Discrete Ax=f check
        ax_eq = abs(xHAx - xHf) / max(abs(xHf), 1e-30)

        # Heuristic normalizations of Im(xᴴ A_vol x) to unit-power scale.
        # Empirically for Pref: Im(xHf)~-10, unit power. Report raw + ratios.
        Im_vol = float(xHAvolx.imag)
        Im_f = float(xHf.imag)
        Im_port = float(xHAportx.imag)
        Im_A = float(xHAx.imag)

        # Candidate discrete absorptions (sign chosen so absorption ≥ 0 when lossy)
        # For e^{-iωt} with our Robin, Im(xHf) is typically negative for unit drive.
        # Physical absorbed power often ~ -Im(xᴴ A_vol x) / (2ω) or similar.
        omega = float(k0)
        cand = {
            "neg_Im_xHAvolx": -Im_vol,
            "neg_Im_xHAvolx_over_2omega": -Im_vol / (2.0 * omega),
            "neg_Im_xHAvolx_over_omega": -Im_vol / omega,
            "Im_xHf_minus_Im_xHAportx": Im_f - Im_port,  # = Im_vol if Ax=f
        }

        rows.append(
            {
                "channel": n,
                "port": p,
                "mode": m,
                "P_out": P_out,
                "P_material_Im_eps": P_im,
                "P_material_dissipation": P_corr,
                "residual_Im_eps": abs(1.0 - P_out - P_im),
                "residual_dissipation": abs(1.0 - P_out - P_corr),
                "Im_xHAx": Im_A,
                "Im_xHf": Im_f,
                "Im_xHAvolx": Im_vol,
                "Im_xHAportx": Im_port,
                "Ax_eq_f_rel": float(ax_eq),
                "discrete_Ax_eq_f": bool(ax_eq < 1e-9),
                "discrete_candidates": cand,
                "sum_out_Imeps": P_out + P_im,
                "sum_out_diss": P_out + P_corr,
            }
        )

    # Pick best discrete candidate by matching mean residual with Pref-like scale
    # Also compare P_corr vs candidates via correlation with P_corr across ports
    def mean_key(key):
        return float(np.mean([r[key] for r in rows]))

    summary = {
        "grade": grade,
        "P_out_mean": mean_key("P_out"),
        "P_Im_eps_mean": mean_key("P_material_Im_eps"),
        "P_diss_mean": mean_key("P_material_dissipation"),
        "residual_Im_eps_max": float(max(r["residual_Im_eps"] for r in rows)),
        "residual_diss_max": float(max(r["residual_dissipation"] for r in rows)),
        "Ax_eq_f_all": all(r["discrete_Ax_eq_f"] for r in rows),
        "elapsed_s": time.time() - t0,
        "per_TE0": rows,
    }
    # Which absorption closes the 5% gate?
    summary["GATE_Im_eps"] = "PASS" if summary["residual_Im_eps_max"] <= 0.05 else "FAIL"
    summary["GATE_dissipation"] = (
        "PASS" if summary["residual_diss_max"] <= 0.05 else "FAIL"
    )
    return summary


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--design-json",
        default="outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json",
    )
    ap.add_argument("--grades", default="C,M")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    d = json.loads(Path(args.design_json).read_text())
    s = np.asarray(d["best"]["s"], float)

    report = {"design": args.design_json, "grades": {}}
    for g in args.grades.split(","):
        g = g.strip()
        print(f"=== grade {g} ===", flush=True)
        try:
            summary = diagnose_grade(g, s)
            report["grades"][g] = summary
            print(
                json.dumps(
                    {
                        "grade": g,
                        "P_out": summary["P_out_mean"],
                        "P_Im_eps": summary["P_Im_eps_mean"],
                        "P_diss": summary["P_diss_mean"],
                        "res_Im_eps_max": summary["residual_Im_eps_max"],
                        "res_diss_max": summary["residual_diss_max"],
                        "GATE_Im_eps": summary["GATE_Im_eps"],
                        "GATE_dissipation": summary["GATE_dissipation"],
                        "Ax_eq_f": summary["Ax_eq_f_all"],
                        "elapsed_s": summary["elapsed_s"],
                    },
                    indent=2,
                ),
                flush=True,
            )
        except Exception as e:
            report["grades"][g] = {"error": repr(e)}
            print(f"FAILED {g}: {e}", flush=True)

    (OUT / "power_balance_rootcause.json").write_text(json.dumps(report, indent=2) + "\n")

    # Markdown
    lines = [
        "# OPTIMIZED_POWER_BALANCE_DIAGNOSTIC",
        "",
        f"Design: `{args.design_json}` (transferred C optimum).",
        "Frozen DtN S package **not** modified.",
        "",
        "## Method",
        "",
        "For each TE0 unit-power drive:",
        "",
        "- **A.** `P_out = Σ|S_{:,n}|²` from direct DtN S",
        "- **B.** `P_Im_eps` = current `material_absorption` (uses elementwise `Im(ε)`)",
        "- **C.** `P_diss` = `E† ((ε−ε†)/(2i)) E` (gyrotropic-safe dissipation)",
        "- **D.** Discrete: `A x = f` ⇒ `xᴴ A x = xᴴ f`; split `A = A_vol + A_port`",
        "",
        "## Results",
        "",
    ]
    for g, sry in report["grades"].items():
        if "error" in sry:
            lines += [f"### Grade {g}", "", f"ERROR: `{sry['error']}`", ""]
            continue
        lines += [
            f"### Grade {g}",
            "",
            f"| quantity | value |",
            f"|----------|------:|",
            f"| P_out mean | {sry['P_out_mean']:.4f} |",
            f"| P_Im_eps mean (current) | {sry['P_Im_eps_mean']:.4f} |",
            f"| P_diss mean (ε−ε†)/(2i) | {sry['P_diss_mean']:.4f} |",
            f"| residual max (Im_eps) | {sry['residual_Im_eps_max']:.4f} |",
            f"| residual max (diss) | {sry['residual_diss_max']:.4f} |",
            f"| GATE Im_eps (≤0.05) | **{sry['GATE_Im_eps']}** |",
            f"| GATE dissipation (≤0.05) | **{sry['GATE_dissipation']}** |",
            f"| discrete Ax=f all TE0 | {sry['Ax_eq_f_all']} |",
            f"| elapsed | {sry['elapsed_s']:.1f} s |",
            "",
            "| ch | P_out | P_Im_eps | P_diss | res_Im | res_diss | Im(xᴴA_vol x) |",
            "|---:|------:|---------:|-------:|-------:|---------:|-------------:|",
        ]
        for r in sry["per_TE0"]:
            lines.append(
                f"| {r['channel']} | {r['P_out']:.4f} | {r['P_material_Im_eps']:.4f} | "
                f"{r['P_material_dissipation']:.4f} | {r['residual_Im_eps']:.4f} | "
                f"{r['residual_dissipation']:.4f} | {r['Im_xHAvolx']:.4f} |"
            )
        lines.append("")

    # Interpretation
    lines += [
        "## Interpretation",
        "",
    ]
    # auto interpret from C if present
    if "C" in report["grades"] and "GATE_Im_eps" in report["grades"]["C"]:
        c = report["grades"]["C"]
        if c["GATE_Im_eps"] == "FAIL" and c["GATE_dissipation"] == "PASS":
            lines += [
                "**Root cause: volume post-processing.** Discrete `Ax=f` closes and",
                "`P_diss = E†((ε−ε†)/(2i))E` restores the ≤5% gate, while the current",
                "`material_absorption` uses elementwise `Im(ε)`, which is incorrect for",
                "a non-Hermitian gyrotropic tensor.",
                "",
                "**Action:** fix absorption post-processing only; do **not** change the",
                "frozen direct-DtN S package.",
                "",
            ]
        elif c["GATE_Im_eps"] == "FAIL" and c["GATE_dissipation"] == "FAIL":
            lines += [
                "Both volume formulas fail the 5% gate. Inspect discrete port power",
                "normalization / PML bookkeeping further.",
                "",
            ]
        else:
            lines += [
                f"C gate Im_eps={c['GATE_Im_eps']}, diss={c['GATE_dissipation']}.",
                "",
            ]

    (OUT / "OPTIMIZED_POWER_BALANCE_DIAGNOSTIC.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUT / "OPTIMIZED_POWER_BALANCE_DIAGNOSTIC.md")


if __name__ == "__main__":
    main()
