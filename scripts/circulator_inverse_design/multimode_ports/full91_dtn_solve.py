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


def exterior_beyond_ports_mask(points, tris, port_structs, margin: float = 1e-6) -> np.ndarray:
    """Triangles outside the device∪feeds polyhedron (beyond all outer port planes).

    Interior = ∩_p {(x−c_p)·n̂_p ≤ margin}. Excising exterior makes each port
    plane a true ∂Ω segment so modal Robin is a boundary condition, not an
    interior penalty on DOFs still coupled through PML.
    """
    cents = points[tris].mean(axis=1)
    exterior = np.zeros(len(tris), dtype=bool)
    for ps in port_structs:
        c = np.asarray(ps["center"], float)
        n = np.asarray(ps["outward"], float)
        n = n / np.linalg.norm(n)
        exterior |= (cents - c) @ n > margin
    return exterior


def assemble_base(grade: str, f_hz: float, b_tesla: float, s_vec: np.ndarray):
    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    plasma, quartz, walls, centers = masks
    from full91_geometry import bulb_centers_device

    centers = bulb_centers_device(sc)
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    # rho_from_scales returns (rho, f_ord, gamma, fc) — not (rho, _, fc, f_ord).
    rho, f_ord, gamma, fc = rho_from_scales(
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
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    port_structs, _modes = build_port_structs(points, shift)
    ext = exterior_beyond_ports_mask(points, tris, port_structs)
    # Excise like PEC walls: ρ→0 drops stiffness+mass (natural Neumann on cut).
    for arr in rho:
        arr[ext] = 0.0
    pec = np.zeros(len(points), dtype=bool)
    # Assembly k0 from physical f_Hz (matches diff_fem_core.solve_forward), not fc.
    k0 = 2.0 * np.pi * (f_hz * sc.a / 299792458.0)
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    return points, tris, rho, A, k0, f_ord, {"n_exterior_tris": int(ext.sum()), "port_structs": port_structs}


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


def port_line_from_struct(ps: dict) -> dict:
    """Monitor line colocated with the DtN outer-feed plane."""
    center = np.asarray(ps["center"], float)
    tang = np.asarray(ps["tangent"], float)
    tang = tang / np.linalg.norm(tang)
    y = np.asarray(ps["y"], float)
    return {
        "port": ps["port"],
        "outward": np.asarray(ps["outward"], float),
        "tangent": tang,
        "monitor_xy": np.vstack([center + float(yi) * tang for yi in y]),
        "weights": np.asarray(ps["w"], float),
        "offsets": y,
    }


def extract_ab_column(hz, port_structs, sampler, tris, k0) -> tuple[np.ndarray, np.ndarray]:
    """Power-normalized a,b for all 12 channels at DtN planes."""
    W = float(sc.clear_width)
    a = np.zeros(12, dtype=np.complex128)
    b = np.zeros(12, dtype=np.complex128)
    for p, ps in enumerate(port_structs):
        line = port_line_from_struct(ps)
        for m in (0, 1):
            aa, bb = extract_mode_m(sampler, hz, tris, line, k0, m, W)
            a[2 * p + m] = aa
            b[2 * p + m] = bb
    return a, b


def full12_s(grade: str, s_vec: np.ndarray, f_hz: float = 3.85e9, b_tesla: float = 0.05):
    points, tris, rho, A0, k0, f_ord, base_meta = assemble_base(grade, f_hz, b_tesla, s_vec)
    port_structs = base_meta["port_structs"]
    modes = mode_constants(k0, float(sc.clear_width), n_modes=4)
    sampler = ElementSampler(points, tris)
    Ain = np.zeros((12, 12), dtype=np.complex128)
    Bout = np.zeros((12, 12), dtype=np.complex128)
    port_counts = [ps["n_nodes"] for ps in port_structs]
    # One factorization for all 12 RHS (Robin independent of drive)
    A = add_all_robin(A0.copy(), port_structs, modes)
    lu = splu(A.tocsc())
    for n in range(12):
        p, m = divmod(n, 2)
        b = np.zeros(len(points), dtype=np.complex128)
        b = incident_rhs(b, port_structs, modes, p, m)
        x = lu.solve(b)
        a_col, b_col = extract_ab_column(x, port_structs, sampler, tris, k0)
        Ain[:, n] = a_col
        Bout[:, n] = b_col
        if n == 0:
            print(
                f"port nodes {port_counts} drive0 |a|={abs(a_col[0]):.4f} "
                f"diag_dom={abs(a_col[0])/max(np.linalg.norm(a_col),1e-30):.3f}",
                flush=True,
            )
    # Proper multimode S: outgoing = S @ incoming
    S = Bout @ np.linalg.inv(Ain)
    diag_dom = [float(abs(Ain[n, n]) / max(np.linalg.norm(Ain[:, n]), 1e-30)) for n in range(12)]
    return S, {
        "port_nodes": port_counts,
        "k0": k0,
        "grade": grade,
        "n_exterior_tris": base_meta["n_exterior_tris"],
        "diag_dominance_min": float(min(diag_dom)),
        "diag_dominance_mean": float(np.mean(diag_dom)),
        "cond_Ain": float(np.linalg.cond(Ain)),
    }


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
    np.savez_compressed(OUT / "S_full12_Bp.npz", S=Sp)
    Sm, _ = full12_s(args.grade, s_vec, b_tesla=-0.05)
    np.savez_compressed(OUT / "S_full12_Bm.npz", S=Sm)
    r_plus = onsager_residual(Sp, Sm)
    S0, _ = full12_s(args.grade, s_vec, b_tesla=0.0)
    np.savez_compressed(OUT / "S_full12_B0.npz", S=S0)
    r = {
        "grade": args.grade,
        "r_F_plus": r_plus,
        "r_F_zero": onsager_residual(S0, S0),
        "onsager_bar": 5e-3,
        "FULL12_ONSAGER": "PASS" if r_plus <= 5e-3 else "FAIL",
        "meta": meta,
        "elapsed_s": time.time() - t0,
    }
    # Block residuals: TE0↔TE0, TE0↔TE1, TE1↔TE1
    def block(S, mi, mj):
        return S[mi::2, mj::2] if False else S[np.ix_([2 * i + mi for i in range(6)], [2 * j + mj for j in range(6)])]

    Sp00, Sm00 = block(Sp, 0, 0), block(Sm, 0, 0)
    Sp01, Sm01 = block(Sp, 0, 1), block(Sm, 0, 1)
    Sp10, Sm10 = block(Sp, 1, 0), block(Sm, 1, 0)
    Sp11, Sm11 = block(Sp, 1, 1), block(Sm, 1, 1)
    r["blocks"] = {
        "r_TE0_TE0": onsager_residual(Sp00, Sm00),
        "r_TE0_TE1_vs_TE1_TE0T": float(
            np.linalg.norm(Sp01 - Sm10.T) / max(np.linalg.norm(Sp01), 1e-30)
        ),
        "r_TE1_TE0_vs_TE0_TE1T": float(
            np.linalg.norm(Sp10 - Sm01.T) / max(np.linalg.norm(Sp10), 1e-30)
        ),
        "r_TE1_TE1": onsager_residual(Sp11, Sm11),
        "mean_TE0_to_TE1_power": float(np.mean(np.abs(Sp01) ** 2)),
        "mean_TE1_to_TE0_power": float(np.mean(np.abs(Sp10) ** 2)),
    }
    args.json.write_text(json.dumps(r, indent=2) + "\n")
    print(json.dumps(r, indent=2))
    return 0 if r["FULL12_ONSAGER"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
