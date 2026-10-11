#!/usr/bin/env python3
"""Validate modal a,b / S_ij formulation before design re-evaluation.

Cases
-----
0. Synthetic bidirectional field → recover a,b and power identity
1. Straight matched guide: S21→1, S11→0
2. Dielectric slab: |S11|²,|S21|² vs analytic Fresnel
3. Full91 B=0: S≈Sᵀ
4. Full91 magnetized: S(+B)≈Sᵀ(−B)
5. Passive power balance: Σ(|a|²−|b|²) ≈ material absorption ≥ 0
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from circulator_core import load_device, solve_six  # noqa: E402
from modal_s_core import (  # noqa: E402
    complex_matrix_to_json,
    dump_json,
    extract_ab,
    modal_port_powers,
    onsager_error,
    reciprocity_error,
    scattering_matrix,
)
from straight_guide_fem import make_port_line, solve_straight_guide  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design" / "modal_s"


def case0_synthetic(k0: float) -> dict:
    """Plant known a,b on a fake sampler line and recover them."""
    W = float(sc.clear_width)
    span = 0.92 * W
    n_points = 41
    line = make_port_line(0.0, np.array([-1.0, 0.0]), W, n_points=n_points)
    # Build a tiny 1-triangle dummy — we bypass sampler by monkeypatching extract
    # Instead: direct formula test with synthetic Hz, dn arrays.
    a_true, b_true = 0.7 + 0.1j, 0.2 - 0.05j
    beta = k0
    omega = k0
    phi = 1.0 / np.sqrt(span)
    # n̂ = −x̂; plant fields consistent with decomposition
    Hz = (a_true + b_true) * phi * np.ones(n_points, dtype=np.complex128)
    dn = 1j * beta * (b_true - a_true) * phi * np.ones(n_points, dtype=np.complex128)
    weights = line["weights"]
    c_H = complex(np.sum(weights * phi * Hz))  # wait phi already in Hz
    # Hz was built as (a+b)*phi, so ∫ φ* Hz = (a+b) ∫ φ² = a+b
    c_H = complex(np.sum(weights * (1.0 / np.sqrt(span)) * Hz))
    c_d = complex(np.sum(weights * (1.0 / np.sqrt(span)) * dn))
    sum_ab = c_H
    diff_ba = c_d / (1j * beta)
    a_f = 0.5 * (sum_ab - diff_ba)
    b_f = 0.5 * (sum_ab + diff_ba)
    scale = np.sqrt(beta / (2.0 * omega))
    a_hat, b_hat = a_f * scale, b_f * scale
    a_true_p = a_true * scale
    b_true_p = b_true * scale
    # Power identity from Poynting
    alpha = -1j / omega
    P_poynt = float(np.sum(weights * 0.5 * np.real(alpha * dn * np.conj(Hz))))
    P_modal = float(abs(b_hat) ** 2 - abs(a_hat) ** 2)
    err_a = abs(a_hat - a_true_p) / max(abs(a_true_p), 1e-30)
    err_b = abs(b_hat - b_true_p) / max(abs(b_true_p), 1e-30)
    err_p = abs(P_poynt - P_modal) / max(abs(P_modal), 1e-30)
    ok = err_a < 1e-10 and err_b < 1e-10 and err_p < 1e-8
    return {
        "case": 0,
        "name": "synthetic_bidirectional",
        "err_a": err_a,
        "err_b": err_b,
        "err_power_identity": err_p,
        "P_modal": P_modal,
        "P_poynting": P_poynt,
        "PASS": bool(ok),
        "threshold": {"err_a": 1e-10, "err_b": 1e-10, "err_p": 1e-8},
    }


def case1_matched_guide() -> dict:
    r = solve_straight_guide(eps_slab=None, drive_port=0, h=0.10, L=8.0)
    S11 = abs(complex(*r["S_col"][0])) ** 2
    S21 = abs(complex(*r["S_col"][1])) ** 2
    # Matched guide: |S11|≪1, |S21|≈1, power |S11|²+|S21|²≈1
    ok = S11 < 0.02 and abs(S21 - 1.0) < 0.05 and abs(S11 + S21 - 1.0) < 0.05
    # Also check modal vs poynting at driven port
    p_err = abs(r["P_net_modal"][0] - r["P_net_poynting"][0]) / max(abs(r["P_net_poynting"][0]), 1e-30)
    return {
        "case": 1,
        "name": "straight_matched_guide",
        "result": r,
        "|S11|^2": S11,
        "|S21|^2": S21,
        "power_sum": S11 + S21,
        "modal_vs_poynting_rel": p_err,
        "PASS": bool(ok and p_err < 0.05),
        "threshold": {"|S11|^2": 0.02, "|S21|^2-1": 0.05, "power_sum-1": 0.05},
    }


def case2_slab() -> dict:
    """Quarter-wave dielectric slab (substantial R, not a Fabry–Perot null)."""
    eps = 3.8 + 0j
    k0 = 2 * np.pi * sc.fs_a
    n = np.sqrt(complex(eps))
    d = float(np.pi / (2.0 * k0 * n.real))  # λ/4 in slab
    r = solve_straight_guide(
        eps_slab=eps, drive_port=0, h=0.08, L=8.0, slab_x0=3.5, slab_x1=3.5 + d
    )
    R = r["|S11|^2"]
    T = r["|S21|^2"]
    Ra = r["analytic"]["R"]
    Ta = r["analytic"]["T"]
    ok = abs(R - Ra) < 0.03 and abs(T - Ta) < 0.03 and abs(R + T - 1.0) < 0.02
    return {
        "case": 2,
        "name": "quarter_wave_slab_fresnel",
        "slab_d": d,
        "result": {k: r[k] for k in r if k != "result"},
        "FEM_R": R,
        "FEM_T": T,
        "analytic_R": Ra,
        "analytic_T": Ta,
        "dR": abs(R - Ra),
        "dT": abs(T - Ta),
        "PASS": bool(ok),
        "threshold": {"dR": 0.03, "dT": 0.03, "R+T-1": 0.02},
    }


def material_absorption(points, tris, x, rho, k0) -> float:
    """Vectorized ∫ (ω/2) E† Im(ε) E with E=(i/ω) ρ curl H, ε=ρ⁻¹.

    Counts only material loss (plasma/quartz). PML loss appears in the modal
    balance P_inc−ΣP_out but not here.
    """
    omega = float(k0)
    px = points[tris, 0]
    py = points[tris, 1]
    twice = (px[:, 1] - px[:, 0]) * (py[:, 2] - py[:, 0]) - (px[:, 2] - px[:, 0]) * (py[:, 1] - py[:, 0])
    area = 0.5 * np.abs(twice)
    good = area > 1e-18
    bx = np.stack([py[:, 1] - py[:, 2], py[:, 2] - py[:, 0], py[:, 0] - py[:, 1]], axis=1) / twice[:, None]
    by = np.stack([px[:, 2] - px[:, 1], px[:, 0] - px[:, 2], px[:, 1] - px[:, 0]], axis=1) / twice[:, None]
    un = x[tris]
    gx = np.sum(bx * un, axis=1)
    gy = np.sum(by * un, axis=1)
    rxx, rxy, ryx, ryy = rho[0], rho[1], rho[2], rho[3]
    det = rxx * ryy - rxy * ryx
    good &= np.abs(det) > 1e-18
    det_safe = np.where(good, det, 1.0)
    exx = ryy / det_safe
    exy = -rxy / det_safe
    eyx = -ryx / det_safe
    eyy = rxx / det_safe
    # Skip lossless elements
    lossy = good & (
        (np.abs(exx.imag) > 1e-18)
        | (np.abs(eyy.imag) > 1e-18)
        | (np.abs(exy.imag) > 1e-18)
        | (np.abs(eyx.imag) > 1e-18)
    )
    if not np.any(lossy):
        return 0.0
    cx = gy[lossy]
    cy = -gx[lossy]
    Ex = (1j / omega) * (rxx[lossy] * cx + rxy[lossy] * cy)
    Ey = (1j / omega) * (ryx[lossy] * cx + ryy[lossy] * cy)
    # E† Im(ε) E
    Ix = exx[lossy].imag * Ex + exy[lossy].imag * Ey
    Iy = eyx[lossy].imag * Ex + eyy[lossy].imag * Ey
    dens = np.real(np.conj(Ex) * Ix + np.conj(Ey) * Iy)
    return float(np.sum(area[lossy] * (omega / 2.0) * dens))


def sixport_modal(dev, s_vec, f_hz, b_tesla, P_ref=None):
    bun = solve_six(dev, s_vec, f_hz, b_tesla, P_ref=P_ref)
    A, B, details = modal_port_powers(dev.sampler, bun.fields, bun.rho, bun.k0, dev.lines)
    S = scattering_matrix(A, B)
    return bun, A, B, S, details


def case3_reciprocity(grade: str, f_hz: float) -> dict:
    dev = load_device(grade)
    s = np.ones(91)
    # Use production B=0 with plasma present (lossy → not unitary, but reciprocal)
    bun, A, B, S, _ = sixport_modal(dev, s, f_hz, 0.0)
    err = reciprocity_error(S)
    ok = err < 0.05
    return {
        "case": 3,
        "name": "B0_reciprocity",
        "grade": grade,
        "f_Hz": f_hz,
        "reciprocity_rel_F": err,
        "S": complex_matrix_to_json(S),
        "PASS": bool(ok),
        "threshold": {"reciprocity_rel_F": 0.05},
        "P_ref_frozen": bun.P_ref,
    }


def case4_onsager(grade: str, f_hz: float, B: float) -> dict:
    dev = load_device(grade)
    s = np.ones(91)
    _, _, _, Sp, _ = sixport_modal(dev, s, f_hz, +B)
    _, _, _, Sm, _ = sixport_modal(dev, s, f_hz, -B)
    err = onsager_error(Sp, Sm)
    ok = err < 0.05
    return {
        "case": 4,
        "name": "Onsager_S(+B)=S^T(-B)",
        "grade": grade,
        "f_Hz": f_hz,
        "B_T": B,
        "onsager_rel_F": err,
        "S_plus": complex_matrix_to_json(Sp),
        "S_minus": complex_matrix_to_json(Sm),
        "PASS": bool(ok),
        "threshold": {"onsager_rel_F": 0.05},
    }


def case5_power_balance(grade: str, f_hz: float, B: float) -> dict:
    dev = load_device(grade)
    s = np.ones(91)
    bun, A, Bmat, S, _ = sixport_modal(dev, s, f_hz, B)
    rows = []
    ok_all = True
    for j in range(6):
        P_inc = float(abs(A[j, j]) ** 2)
        P_out = float(np.sum(np.abs(Bmat[:, j]) ** 2))
        P_abs_modal = P_inc - P_out
        P_abs_mat = material_absorption(dev.points, dev.tris, bun.fields[j], bun.rho, bun.k0)
        # Passive: absorbed ≥ −tol; modal and material should agree loosely
        rel = abs(P_abs_modal - P_abs_mat) / max(P_inc, 1e-30)
        ok = P_abs_modal > -0.05 * P_inc and rel < 0.25
        ok_all = ok_all and ok
        rows.append({
            "source": j,
            "P_inc": P_inc,
            "P_out_sum": P_out,
            "P_abs_modal": P_abs_modal,
            "P_abs_material": P_abs_mat,
            "rel_diff_over_Pinc": rel,
            "PASS": bool(ok),
        })
    return {
        "case": 5,
        "name": "passive_power_balance",
        "grade": grade,
        "f_Hz": f_hz,
        "B_T": B,
        "per_source": rows,
        "PASS": bool(ok_all),
        "threshold": {
            "P_abs_modal_over_Pinc_min": -0.05,
            "rel_diff_max": 0.25,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--f-hz", type=float, default=3.85e9)
    ap.add_argument("--B", type=float, default=0.05)
    ap.add_argument("--skip-sixport", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    report = {"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "cases": []}

    print("=== Case 0: synthetic ===", flush=True)
    k0 = 2 * np.pi * sc.fs_a
    c0 = case0_synthetic(k0)
    report["cases"].append(c0)
    dump_json(OUT / "case0_synthetic.json", c0)
    print(" PASS" if c0["PASS"] else " FAIL", c0, flush=True)

    print("=== Case 1: matched guide ===", flush=True)
    c1 = case1_matched_guide()
    report["cases"].append({k: c1[k] for k in c1 if k != "result"})
    dump_json(OUT / "case1_matched_guide.json", c1)
    print(
        f" PASS={c1['PASS']} |S11|^2={c1['|S11|^2']:.4f} |S21|^2={c1['|S21|^2']:.4f}",
        flush=True,
    )

    print("=== Case 2: dielectric slab ===", flush=True)
    c2 = case2_slab()
    # strip bulky result nesting
    c2_save = {k: c2[k] for k in c2}
    report["cases"].append({k: c2[k] for k in ("case", "name", "FEM_R", "FEM_T", "analytic_R", "analytic_T", "dR", "dT", "PASS", "threshold")})
    dump_json(OUT / "case2_slab.json", c2_save)
    print(
        f" PASS={c2['PASS']} R={c2['FEM_R']:.4f} (an {c2['analytic_R']:.4f}) "
        f"T={c2['FEM_T']:.4f} (an {c2['analytic_T']:.4f})",
        flush=True,
    )

    if not args.skip_sixport:
        print(f"=== Case 3: B=0 reciprocity grade={args.grade} ===", flush=True)
        t0 = time.perf_counter()
        c3 = case3_reciprocity(args.grade, args.f_hz)
        c3["wall_s"] = time.perf_counter() - t0
        report["cases"].append({k: c3[k] for k in c3 if k != "S"})
        dump_json(OUT / "case3_reciprocity.json", c3)
        print(f" PASS={c3['PASS']} reciprocity_rel_F={c3['reciprocity_rel_F']:.4e}", flush=True)

        print(f"=== Case 4: Onsager ±B={args.B} ===", flush=True)
        t0 = time.perf_counter()
        c4 = case4_onsager(args.grade, args.f_hz, args.B)
        c4["wall_s"] = time.perf_counter() - t0
        report["cases"].append({k: c4[k] for k in c4 if k not in ("S_plus", "S_minus")})
        dump_json(OUT / "case4_onsager.json", c4)
        print(f" PASS={c4['PASS']} onsager_rel_F={c4['onsager_rel_F']:.4e}", flush=True)

        print(f"=== Case 5: power balance ===", flush=True)
        t0 = time.perf_counter()
        c5 = case5_power_balance(args.grade, args.f_hz, args.B)
        c5["wall_s"] = time.perf_counter() - t0
        report["cases"].append(c5)
        dump_json(OUT / "case5_power_balance.json", c5)
        print(f" PASS={c5['PASS']}", flush=True)

    report["all_PASS"] = all(c.get("PASS", False) for c in report["cases"])
    report["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    dump_json(OUT / "VALIDATION_SUMMARY.json", report)
    print("=== SUMMARY all_PASS =", report["all_PASS"], flush=True)
    return 0 if report["all_PASS"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
