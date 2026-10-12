#!/usr/bin/env python3
"""Phase 1: adversarially audit modal-S Onsager on the Pref sixfold peak design.

Frozen acceptance criteria: ../modal_s_optimization/ONSAGER_CRITERION_FROZEN.md
(read before judging PASS/FAIL; do not loosen after seeing numbers).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "modal_s_opt"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from circulator_core import F_HZ, B_PLUS, load_device, solve_six  # noqa: E402
from modal_s_core import (  # noqa: E402
    complex_matrix_to_json,
    dump_json,
    modal_port_powers,
    onsager_error,
    scattering_matrix,
)
from port_line_variants import make_port_lines  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design" / "modal_s_optimization"
DESIGNS = ROOT / "outputs" / "circulator_inverse_design"

# Frozen thresholds (must match ONSAGER_CRITERION_FROZEN.md)
R_F_MAX = 5e-3
ELEM_ABS_MAX = 0.02
ELEM_MAG_MAX = 0.02
SIG_THRESH = 0.10
MON_MAG_MAX = 0.015
NPTS_MAG_MAX = 0.005


def load_s_vec(tag: str) -> np.ndarray:
    if tag == "uniform":
        return np.ones(91)
    d = json.loads((DESIGNS / f"{tag}.json").read_text())
    return np.asarray(d["best"]["s"], float)


def S_from_bundle(bun, lines):
    A, B, details = modal_port_powers(bun.sampler if hasattr(bun, "sampler") else None, bun.fields, bun.rho, bun.k0, lines)
    # sampler lives on device — pass via lines context; fix below
    return A, B, details


def compute_S(dev, bun, lines):
    A, B, details = modal_port_powers(dev.sampler, bun.fields, bun.rho, bun.k0, lines)
    S = scattering_matrix(A, B)
    return S, A, B, details


def elementwise_onsager(Sp: np.ndarray, Sm: np.ndarray, sig_thresh: float = SIG_THRESH) -> dict:
    """Compare S(+B)_ij vs S(-B)_ji."""
    SmT = Sm.T
    dS = Sp - SmT
    rows = []
    worst_abs = 0.0
    worst_mag = 0.0
    n_sig = 0
    n_fail = 0
    for i in range(6):
        for j in range(6):
            mag = abs(Sp[i, j])
            dab = abs(dS[i, j])
            dmag = abs(abs(Sp[i, j]) - abs(SmT[i, j]))
            sig = mag >= sig_thresh
            if sig:
                n_sig += 1
                worst_abs = max(worst_abs, dab)
                worst_mag = max(worst_mag, dmag)
                fail = dab > ELEM_ABS_MAX or dmag > ELEM_MAG_MAX
                if fail:
                    n_fail += 1
            else:
                fail = False
            if sig or dab > 0.01:
                rows.append({
                    "i": i, "j": j,
                    "S_plus": [Sp[i, j].real, Sp[i, j].imag],
                    "S_minus_T": [SmT[i, j].real, SmT[i, j].imag],
                    "abs_diff": float(dab),
                    "mag_diff": float(dmag),
                    "|S_plus|": float(mag),
                    "significant": bool(sig),
                    "fail_elem": bool(fail and sig),
                })
    rows.sort(key=lambda r: -r["abs_diff"])
    return {
        "worst_abs_diff_sig": float(worst_abs),
        "worst_mag_diff_sig": float(worst_mag),
        "n_significant": n_sig,
        "n_fail_significant": n_fail,
        "PASS_elementwise": bool(n_fail == 0 and n_sig > 0),
        "top_residuals": rows[:20],
    }


def te1_summary(details) -> dict:
    fracs = []
    pouts = []
    for j, row in enumerate(details):
        for i, d in enumerate(row):
            fracs.append(d.get("te1_power_frac", 0.0))
            pouts.append(d.get("te1_P_out", 0.0))
    return {
        "te1_frac_max": float(max(fracs) if fracs else 0.0),
        "te1_frac_mean": float(np.mean(fracs) if fracs else 0.0),
        "te1_P_out_max": float(max(pouts) if pouts else 0.0),
        "te1_P_out_mean": float(np.mean(pouts) if pouts else 0.0),
    }


def run_pair(dev, s, f_hz, B, lines, tag: str):
    t0 = time.perf_counter()
    bun_p = solve_six(dev, s, f_hz, +B, P_ref=1.0)
    Sp, Ap, Bp, dp = compute_S(dev, bun_p, lines)
    bun_m = solve_six(dev, s, f_hz, -B, P_ref=1.0)
    Sm, Am, Bm, dm = compute_S(dev, bun_m, lines)
    rF = onsager_error(Sp, Sm)
    elem = elementwise_onsager(Sp, Sm)
    te1p = te1_summary(dp)
    te1m = te1_summary(dm)
    return {
        "tag": tag,
        "wall_s": time.perf_counter() - t0,
        "factor_s_plus": bun_p.factor_s,
        "factor_s_minus": bun_m.factor_s,
        "r_F": rF,
        "elementwise": elem,
        "te1_plus": te1p,
        "te1_minus": te1m,
        "S_plus": complex_matrix_to_json(Sp),
        "S_minus": complex_matrix_to_json(Sm),
        "T_des_plus": float(np.mean([abs(Sp[(j + 1) % 6, j]) ** 2 for j in range(6)])),
        "overlap_min_plus": float(min(d["overlap_quality"] for row in dp for d in row)),
    }


def monitor_invariance(dev, bun_p, bun_m, shifts, n_points=41):
    """Reuse ±B fields; only re-extract at shifted monitors."""
    base_lines = make_port_lines(n_points=n_points, monitor_shift_outward=0.0)
    Sp0, _, _, _ = compute_S(dev, bun_p, base_lines)
    results = {"shifts": {}, "max_sig_mag_delta": 0.0}
    for dz in shifts:
        lines = make_port_lines(n_points=n_points, monitor_shift_outward=dz)
        Sp, _, _, _ = compute_S(dev, bun_p, lines)
        # magnitude deltas vs baseline for significant channels
        worst = 0.0
        for i in range(6):
            for j in range(6):
                if abs(Sp0[i, j]) >= SIG_THRESH:
                    worst = max(worst, abs(abs(Sp[i, j]) - abs(Sp0[i, j])))
        # phase check: Δφ ≈ −β Δζ for transmission (outgoing at receive)
        # report raw for desired channels
        des_phase = []
        for j in range(6):
            i = (j + 1) % 6
            if abs(Sp0[i, j]) > 1e-6 and abs(Sp[i, j]) > 1e-6:
                dphi = np.angle(Sp[i, j] / Sp0[i, j])
                des_phase.append(float(dphi))
        results["shifts"][str(dz)] = {
            "worst_sig_mag_delta": float(worst),
            "desired_dphi_rad": des_phase,
            "desired_dphi_mean": float(np.mean(des_phase)) if des_phase else 0.0,
            "r_F_vs_minus_baseline_extract": None,  # filled if Sm provided
        }
        results["max_sig_mag_delta"] = max(results["max_sig_mag_delta"], worst)
        # Onsager at this station
        Sm, _, _, _ = compute_S(dev, bun_m, lines)
        results["shifts"][str(dz)]["r_F"] = onsager_error(Sp, Sm)
    results["PASS_monitor_invariance"] = bool(results["max_sig_mag_delta"] <= MON_MAG_MAX)
    # expected phase for receive-port outward travel: monitor shift +outward moves
    # receive monitor toward exterior; outgoing b phase ~ e^{+i β ζ_local} with ζ along outward
    # Δζ_outward > 0 ⇒ path from device to monitor shorter ⇒ phase of b decreases by β|Δζ|
    k0 = bun_p.k0
    results["expected_dphi_per_plus_shift"] = float(-k0)  # rough for +1 outward on receive
    return results


def npoints_invariance(dev, bun_p, n_list):
    base = make_port_lines(n_points=41)
    Sp0, _, _, _ = compute_S(dev, bun_p, base)
    out = {"by_n": {}, "max_sig_mag_delta": 0.0}
    for n in n_list:
        lines = make_port_lines(n_points=n)
        Sp, _, _, _ = compute_S(dev, bun_p, lines)
        worst = 0.0
        for i in range(6):
            for j in range(6):
                if abs(Sp0[i, j]) >= SIG_THRESH:
                    worst = max(worst, abs(abs(Sp[i, j]) - abs(Sp0[i, j])))
        out["by_n"][str(n)] = {"worst_sig_mag_delta_vs_41": float(worst)}
        if n != 41:
            out["max_sig_mag_delta"] = max(out["max_sig_mag_delta"], worst)
    out["PASS_npoints"] = bool(out["max_sig_mag_delta"] <= NPTS_MAG_MAX)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grades", default="C,M,F")
    ap.add_argument("--design", default="opt_single_uniform_M_Bp0.0500")
    ap.add_argument("--f-hz", type=float, default=F_HZ)
    ap.add_argument("--B", type=float, default=B_PLUS)
    ap.add_argument("--skip-F", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    s = load_s_vec(args.design)
    grades = [g.strip() for g in args.grades.split(",") if g.strip()]
    if args.skip_F:
        grades = [g for g in grades if g != "F"]

    report = {
        "design": args.design,
        "f_Hz": args.f_hz,
        "B_T": args.B,
        "frozen_criteria": {
            "r_F_max": R_F_MAX,
            "elem_abs_max": ELEM_ABS_MAX,
            "elem_mag_max": ELEM_MAG_MAX,
            "sig_thresh": SIG_THRESH,
            "monitor_mag_max": MON_MAG_MAX,
            "npoints_mag_max": NPTS_MAG_MAX,
            "doc": "ONSAGER_CRITERION_FROZEN.md",
        },
        "feed_modes": {
            "note": "TE1 propagates at 3.85 GHz in clear_width=2 feed",
            "clear_width": 2.0,
            "kc1": float(np.pi / 2.0),
            "kc1_vs_k0": "kc1=1.571 < k0=1.614 at 3.85 GHz → TE1 propagates (β1≈0.37)",
        },
        "grades": {},
        "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    mesh_rF = {}
    for grade in grades:
        print(f"=== grade {grade} design={args.design} ===", flush=True)
        t0 = time.perf_counter()
        dev = load_device(grade)
        lines = make_port_lines(n_points=41, monitor_shift_outward=0.0)
        # Override device lines for extraction only; RHS still from production inject
        pair = run_pair(dev, s, args.f_hz, args.B, lines, tag=f"{grade}_n41_dz0")
        print(f"  r_F={pair['r_F']:.6e} elem_pass={pair['elementwise']['PASS_elementwise']} "
              f"worst_abs={pair['elementwise']['worst_abs_diff_sig']:.4f} "
              f"te1_frac_max={pair['te1_plus']['te1_frac_max']:.4f} "
              f"T_des={pair['T_des_plus']:.4f}", flush=True)

        # Reuse fields for monitor / npoints studies: re-solve once pair already has cost;
        # solve again to keep fields (run_pair doesn't return bundles). Re-solve +B/-B once more for invariance.
        print(f"  monitor/npoints invariance...", flush=True)
        bun_p = solve_six(dev, s, args.f_hz, +args.B, P_ref=1.0)
        bun_m = solve_six(dev, s, args.f_hz, -args.B, P_ref=1.0)
        mon = monitor_invariance(dev, bun_p, bun_m, shifts=[-0.3, 0.0, 0.3], n_points=41)
        npt = npoints_invariance(dev, bun_p, n_list=[41, 81])
        print(f"  mon_mag_delta={mon['max_sig_mag_delta']:.4f} pass={mon['PASS_monitor_invariance']} "
              f"npts_delta={npt['max_sig_mag_delta']:.4f} pass={npt['PASS_npoints']}", flush=True)

        g_rec = {
            "baseline": {k: pair[k] for k in pair if k not in ("S_plus", "S_minus")},
            "S_plus": pair["S_plus"],
            "S_minus": pair["S_minus"],
            "monitor_invariance": mon,
            "npoints_invariance": npt,
            "wall_s_total": time.perf_counter() - t0,
            "n_nodes": int(len(dev.points)),
        }
        # local PASS pieces
        g_rec["PASS_rF"] = bool(pair["r_F"] <= R_F_MAX)
        g_rec["PASS_elem"] = bool(pair["elementwise"]["PASS_elementwise"])
        g_rec["PASS_monitor"] = bool(mon["PASS_monitor_invariance"])
        g_rec["PASS_npoints"] = bool(npt["PASS_npoints"])
        report["grades"][grade] = g_rec
        mesh_rF[grade] = pair["r_F"]
        dump_json(OUT / f"onsager_audit_{args.design}_{grade}.json", g_rec)

    # Mesh trend
    mesh_pass = True
    if "C" in mesh_rF and "M" in mesh_rF:
        mesh_pass = (mesh_rF["M"] < mesh_rF["C"]) or (mesh_rF["M"] <= R_F_MAX and mesh_rF["C"] <= R_F_MAX)
    if "F" in mesh_rF and "M" in mesh_rF:
        mesh_pass = mesh_pass and (mesh_rF["F"] <= mesh_rF["M"] + 1e-3)

    # Overall: use finest available among M,F else C
    finest = "F" if "F" in report["grades"] else ("M" if "M" in report["grades"] else grades[-1])
    g = report["grades"][finest]
    overall = (
        g["PASS_rF"]
        and g["PASS_elem"]
        and g["PASS_monitor"]
        and g["PASS_npoints"]
        and mesh_pass
    )
    report["mesh_rF"] = mesh_rF
    report["PASS_mesh_trend"] = bool(mesh_pass)
    report["finest_grade"] = finest
    report["MODAL_S_OPTIMIZED_ONSAGER"] = "PASS" if overall else "FAIL"
    report["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    dump_json(OUT / "onsager_audit_summary.json", report)
    print("=== OVERALL", report["MODAL_S_OPTIMIZED_ONSAGER"], "finest", finest, "mesh_rF", mesh_rF, flush=True)
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
