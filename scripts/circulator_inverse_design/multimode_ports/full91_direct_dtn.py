#!/usr/bin/env python3
"""FULL91 direct DtN S (true boundary edges, S = I + 2i Q)."""
from __future__ import annotations

import json
import sys
import time
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

import sixport_common as sc  # noqa: E402
from boundary_ports import build_port_from_boundary  # noqa: E402
from diff_fem_core import per_rod_plasma_masks, rho_from_scales  # noqa: E402
from direct_dtn_s import (  # noqa: E402
    assemble_port_operator,
    build_Q_columns,
    direct_S_from_factor,
    onsager_metrics,
)
from dtn_port_geometry import port_outer_center  # noqa: E402
from fem_validated_solver import assemble_anisotropic  # noqa: E402
from full91_dtn_solve import exterior_beyond_ports_mask  # noqa: E402
from full91_solve import load_mesh, material_masks  # noqa: E402
from modes import mode_constants  # noqa: E402

ONSAGER_BAR = 5e-3


def load_s(tag: str) -> np.ndarray:
    if tag == "uniform":
        return np.ones(91)
    d = json.loads((DES / f"{tag}.json").read_text())
    return np.asarray(d["best"]["s"], float)


def assemble_device(grade: str, f_hz: float, b_tesla: float, s_vec: np.ndarray):
    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    plasma, quartz, walls, centers = masks
    from full91_geometry import bulb_centers_device

    centers = bulb_centers_device(sc)
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    rho, f_ord, gamma, fc = rho_from_scales(
        len(tris), rod_masks, s_vec, quartz, walls, f_hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz
    )
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    # Temporary port centers for excision (geometry only)
    rough = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        rough.append({"center": c, "outward": n_hat})
    ext = exterior_beyond_ports_mask(points, tris, rough)
    for arr in rho:
        arr[ext] = 0.0
    k0 = 2.0 * np.pi * (f_hz * sc.a / 299792458.0)
    A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))

    ports = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        ports.append(
                build_port_from_boundary(
                points,
                tris,
                rho,
                c,
                n_hat,
                tang,
                W=float(sc.clear_width),
                plane_tol=0.12,
                port=p,
            )
        )
    modes = mode_constants(k0, float(sc.clear_width), n_modes=4)
    A = assemble_port_operator(A, points, ports, modes, n_prop=2, n_evan=2)
    return points, tris, rho, A, k0, ports, modes, int(ext.sum())


def full12_direct(grade: str, s_vec: np.ndarray, b_tesla: float):
    points, tris, rho, A, k0, ports, modes, n_ext = assemble_device(grade, 3.85e9, b_tesla, s_vec)
    lu = splu(A.tocsc())
    cols = build_Q_columns(points, ports, modes, n_prop=2)
    S = direct_S_from_factor(lu, len(points), cols)
    geo = [{"port": p["port"], "n_edges": p["n_edges"], "n_nodes": p["n_nodes"], "length": p["total_length"]} for p in ports]
    return S, {"k0": k0, "n_exterior": n_ext, "ports": geo, "grade": grade, "B": b_tesla}


def te0_summary(S: np.ndarray) -> dict:
    rows = []
    for j in range(6):
        col = 2 * j
        des = (j + 1) % 6
        rev = (j - 1) % 6
        thr = (j + 3) % 6
        Tdes = abs(S[2 * des, col]) ** 2
        Trev = abs(S[2 * rev, col]) ** 2
        Tthr = abs(S[2 * thr, col]) ** 2
        R = abs(S[col, col]) ** 2
        te1 = float(sum(abs(S[2 * i + 1, col]) ** 2 for i in range(6)))
        rows.append({
            "T_des": float(Tdes),
            "T_rev": float(Trev),
            "T_thr": float(Tthr),
            "R": float(R),
            "TE1": te1,
            "IL_dB": float(-10 * np.log10(max(Tdes, 1e-30))),
            "iso_dB": float(10 * np.log10(max(Tdes, 1e-30) / max(Trev, 1e-30))),
        })
    mean = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    return {"per_source": rows, "mean": mean}


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument(
        "--designs",
        default="uniform,opt_single_uniform_M_Bp0.0500,opt_broadband_C_Bp0.0500",
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"grade": args.grade, "onsager_bar": ONSAGER_BAR, "designs": {}}
    t0 = time.time()
    for tag in args.designs.split(","):
        tag = tag.strip()
        print(f"=== {tag} ===", flush=True)
        s = load_s(tag)
        Sp, meta = full12_direct(args.grade, s, 0.05)
        print("  ports", meta["ports"], flush=True)
        Sm, _ = full12_direct(args.grade, s, -0.05)
        S0, _ = full12_direct(args.grade, s, 0.0)
        np.savez_compressed(OUT / f"S12_direct_{tag}_Bp.npz", S=Sp)
        np.savez_compressed(OUT / f"S12_direct_{tag}_Bm.npz", S=Sm)
        m_pm = onsager_metrics(Sp, Sm)
        m_0 = onsager_metrics(S0, S0)
        mets = te0_summary(Sp)
        rem = float(np.mean([1.0 - np.sum(np.abs(Sp[:, 2 * j]) ** 2) for j in range(6)]))
        rec = {
            "meta": meta,
            "onsager_pm": m_pm,
            "onsager_B0": m_0,
            "FULL12_ONSAGER": "PASS" if m_pm["PASS"] else "FAIL",
            "te0": mets,
            "mean_power_remaining": rem,
        }
        report["designs"][tag] = rec
        print(
            f"  r_F={m_pm['r_F']:.4e} worst={m_pm['worst_abs_entry']:.4e} "
            f"B0={m_0['r_F']:.4e} Tdes={mets['mean']['T_des']:.4f} "
            f"TE1={mets['mean']['TE1']:.4f} rem={rem:.4f} => {rec['FULL12_ONSAGER']}",
            flush=True,
        )
    report["elapsed_s"] = time.time() - t0
    all_pass = all(d["FULL12_ONSAGER"] == "PASS" for d in report["designs"].values())
    report["DIRECT_DTN_FULL12_ONSAGER"] = "PASS" if all_pass else "FAIL"
    (OUT / "direct_dtn_full12_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"DIRECT_DTN_FULL12_ONSAGER": report["DIRECT_DTN_FULL12_ONSAGER"], "elapsed_s": report["elapsed_s"]}, indent=2))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
