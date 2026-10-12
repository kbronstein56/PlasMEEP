#!/usr/bin/env python3
"""FULL91 multimode DtN port solves (12 unit incident modes)."""
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

import sixport_common as sc  # noqa: E402
from diff_fem_core import per_rod_plasma_masks, rho_from_scales  # noqa: E402
from dtn_assembly import add_modal_incident_rhs, add_modal_robin  # noqa: E402
from dtn_port_geometry import mesh_port_nodes, port_outer_center  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from full91_solve import load_mesh, material_masks  # noqa: E402
from modal_extract import extract_mode_m  # noqa: E402
from modes import mode_constants, normalize_on_line, power_scale  # noqa: E402


def port_lines_from_horns(n_points: int = 41):
    """Monitor lines at legacy monitor centers (for cross-check)."""
    lines = []
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    span = 0.92 * sc.clear_width
    offsets = np.linspace(-span / 2.0, span / 2.0, n_points)
    ds = span / (n_points - 1)
    weights = np.full(n_points, ds)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    for port in range(6):
        horn = sc.full_horns[port]
        n_hat = np.asarray(sc.port_dirs[port], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        mon = np.asarray(horn["monitor_center"], float) + shift
        lines.append({
            "port": port,
            "outward": n_hat,
            "tangent": tang,
            "monitor_xy": np.vstack([mon + float(s) * tang for s in offsets]),
            "weights": weights,
            "offsets": offsets,
        })
    return lines


def build_port_structs(points, shift):
    W = float(sc.clear_width)
    k0 = 2 * np.pi * sc.fs_a
    modes = mode_constants(k0, W, n_modes=4)
    structs = []
    for p in range(6):
        horn = sc.full_horns[p]
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        center = port_outer_center(horn, n_hat) + shift
        idx, y, w = mesh_port_nodes(points, center, n_hat, tang, W=W)
        structs.append({
            "port": p,
            "outward": n_hat,
            "tangent": tang,
            "center": center,
            "nodes": idx,
            "y": y,
            "w": w,
            "n_nodes": len(idx),
        })
    return structs, modes


def assemble_base(grade: str, f_hz: float, b_tesla: float, s_vec: np.ndarray):
    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    plasma, quartz, walls, centers = masks
    from full91_geometry import bulb_centers_device

    centers = bulb_centers_device(sc)
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    rho, _, fc, f_ord = rho_from_scales(
        len(tris),
        rod_masks,
        s_vec,
        quartz,
        walls,
        f_hz,
        b_tesla,
        sc.a,
        sc.fp_Hz,
        sc.gamma_Hz,
    )
    pec = np.zeros(len(points), dtype=bool)
    k0 = 2 * np.pi * f_ord
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    return points, tris, rho, A, k0, f_ord


def add_all_robin(A, port_structs, modes, n_prop=2, n_evan=2):
    W = float(sc.clear_width)
    for ps in port_structs:
        if ps["n_nodes"] < 3:
            continue
        for md in modes[: n_prop + n_evan]:
            m = md["m"]
            phi, _ = normalize_on_line(ps["y"], ps["w"], m, W)
            A = add_modal_robin(A, ps["nodes"], ps["y"], ps["w"], phi, md["beta"])
    return A


def incident_rhs(b, port_structs, modes, drive_port: int, drive_mode: int):
    W = float(sc.clear_width)
    ps = port_structs[drive_port]
    md = modes[drive_mode]
    phi, _ = normalize_on_line(ps["y"], ps["w"], drive_mode, W)
    scale = power_scale(md["beta"], 2 * np.pi * sc.fs_a)
    a_field = 1.0 / max(abs(scale), 1e-30)
    # side label unused when coeff is symmetric; left/right both use +2j in current frozen form
    return add_modal_incident_rhs(b, ps["nodes"], ps["y"], ps["w"], phi, md["beta"], a_field, side="left")


def extract_S(fields, port_structs, lines, sampler, tris, rho, k0, drive_col: int):
    W = float(sc.clear_width)
    n_ch = 12
    S_col = np.zeros(n_ch, dtype=np.complex128)
    a_in = None
    for p, line in enumerate(lines):
        for m in (0, 1):
            a, bb = extract_mode_m(sampler, fields, tris, line, k0, m, W)
            ch = 2 * p + m
            if ch == drive_col:
                a_in = a
            S_col[ch] = bb
    if a_in is None or abs(a_in) < 1e-30:
        return S_col, complex(np.nan)
    return S_col / a_in, a_in


def full12_s(grade: str, s_vec: np.ndarray, f_hz: float = 3.85e9, b_tesla: float = 0.05):
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    points, tris, rho, A0, k0, f_ord = assemble_base(grade, f_hz, b_tesla, s_vec)
    port_structs, modes = build_port_structs(points, shift)
    lines = port_lines_from_horns()
    sampler = ElementSampler(points, tris)
    S = np.zeros((12, 12), dtype=np.complex128)
    port_counts = [ps["n_nodes"] for ps in port_structs]
    for n in range(12):
        p, m = divmod(n, 2)
        A = add_all_robin(A0.copy(), port_structs, modes)
        b = np.zeros(len(points), dtype=np.complex128)
        b = incident_rhs(b, port_structs, modes, p, m)
        lu = splu(A.tocsc())
        x = lu.solve(b)
        S[:, n], ain = extract_S(x, port_structs, lines, sampler, tris, rho, k0, n)
        if n == 0:
            print(f"port nodes {port_counts} drive0 ain={abs(ain):.4f}", flush=True)
    return S, {"port_nodes": port_counts, "k0": k0, "grade": grade}


def onsager_residual(Sp, Sm):
    num = np.linalg.norm(Sp - Sm.T)
    den = max(np.linalg.norm(Sp), 1e-30)
    return float(num / den)


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument("--uniform", action="store_true")
    ap.add_argument("--json", type=Path, default=OUT / "full12_dtn_probe.json")
    args = ap.parse_args()
    s_vec = np.ones(91) if args.uniform else np.ones(91)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    Sp, meta = full12_s(args.grade, s_vec, b_tesla=0.05)
    Sm, _ = full12_s(args.grade, s_vec, b_tesla=-0.05)
    S0, _ = full12_s(args.grade, s_vec, b_tesla=0.0)
    r = {
        "grade": args.grade,
        "r_F_plus": onsager_residual(Sp, Sm),
        "r_F_zero": onsager_residual(S0, S0.T),
        "meta": meta,
        "elapsed_s": time.time() - t0,
    }
    args.json.write_text(json.dumps(r, indent=2) + "\n")
    print(json.dumps(r, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
