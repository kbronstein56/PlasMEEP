#!/usr/bin/env python3
"""
Six PEC horns, no bulbs. FEM consumes meep_horn_polygons.json (the Add_Prism vertices).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from fem_validated_solver import (  # noqa: E402
    assemble_anisotropic,
    fem_incident_reference,
    fields_from_hz,
    guide_normal_flux,
    inject_numerical_mode_rhs,
    load_mesh,
    load_numerical_mode,
    solve_system,
)
import sixport_common as sc  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_localization"
POLY = OUT / "meep_horn_polygons.json"


def quads_mesh(doc):
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    out = []
    for horn in doc["horns"]:
        for wall in horn["walls"]:
            out.append(np.asarray(wall["vertices_a"], dtype=float) + shift)
    return out


def metal_mask(points, tris, quads):
    cents = points[tris].mean(axis=1)
    mask = np.zeros(len(tris), dtype=bool)
    for q in quads:
        mask |= MPath(q).contains_points(cents, radius=1e-12)
    return mask


def boundary_element_size(points, tris, mask):
    """Minimum edge length among triangles marked metal."""
    edges = []
    mt = tris[mask]
    for a, b in ((0, 1), (1, 2), (2, 0)):
        edges.append(np.linalg.norm(points[mt[:, a]] - points[mt[:, b]], axis=1))
    if not edges:
        return float("nan")
    return float(np.min(np.concatenate(edges)))


def sample_line(points, Hz, Ex, Ey, center, tangent, span, n=41):
    s = np.linspace(-span / 2, span / 2, n)
    hz = np.zeros(n, dtype=np.complex128)
    ex = np.zeros(n, dtype=np.complex128)
    ey = np.zeros(n, dtype=np.complex128)
    for i, si in enumerate(s):
        xy = center + si * tangent
        k = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
        hz[i], ex[i], ey[i] = Hz[k], Ex[k], Ey[k]
    sn = 0.5 * np.real(ey * np.conj(hz))  # placeholder, overwritten by caller normal
    return s, hz, ex, ey


def main() -> int:
    doc = json.loads(POLY.read_text())
    quads = quads_mesh(doc)
    src = load_numerical_mode()
    grades = ["FEM-L", "FEM-M", "FEM-H"]
    results = []
    for grade in grades:
        print(f"=== {grade} horns only, Meep polygons ===", flush=True)
        points, tris = load_mesh(grade)
        mask = metal_mask(points, tris, quads)
        T = len(tris)
        rho = [
            np.ones(T, dtype=np.complex128),
            np.zeros(T, dtype=np.complex128),
            np.zeros(T, dtype=np.complex128),
            np.ones(T, dtype=np.complex128),
        ]
        for arr in rho:
            arr[mask] = 0.0
        pec = np.zeros(len(points), dtype=bool)
        k0 = 2 * np.pi * float(sc.fs_a)
        t0 = time.perf_counter()
        A = assemble_anisotropic(points, tris, *rho, k0, pec)
        t_asm = time.perf_counter() - t0
        b = inject_numerical_mode_rhs(points, src, pec)
        x, _lu, t_fac, t_sol, resid = solve_system(A, b, pec)
        Ex, Ey = fields_from_hz(points, tris, x, *rho, k0)
        raw = {p: float(guide_normal_flux(points, x, Ex, Ey, p)) for p in range(6)}
        inc = fem_incident_reference(grade, src)
        P_inc = float(inc["P_inc"])
        norm = {p: raw[p] / P_inc for p in range(6)}
        db = {p: (float(10 * np.log10(v)) if v > 0 else None) for p, v in norm.items()}
        rec = {
            "grade": grade,
            "n_nodes": int(len(points)),
            "n_tris": int(T),
            "metal_tri_frac": float(mask.mean()),
            "min_metal_edge_a": boundary_element_size(points, tris, mask),
            "assemble_s": t_asm,
            "factor_s": t_fac,
            "solve_s": t_sol,
            "residual": resid,
            "P_inc": P_inc,
            "raw_flux": raw,
            "normalized_power": norm,
            "power_dB": db,
        }
        print(json.dumps(rec, indent=2), flush=True)
        results.append(rec)
        if grade == "FEM-M":
            # Field lines in mesh coordinates for the divergence diagnostic.
            span = float(0.96 * sc.clear_width)
            shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
            h0 = doc["horns"][0]
            n_hat = np.asarray(h0["outward"], dtype=float)
            tang = np.asarray(h0["tangent"], dtype=float)
            lines = {
                "source": np.asarray(h0["source_center_a"], dtype=float) + shift,
                "throat": np.asarray(h0["throat_a"], dtype=float) + shift,
                "monitor": np.asarray(h0["monitor_center_a"], dtype=float) + shift,
                "center": shift.copy(),
            }
            saved = {}
            for name, c in lines.items():
                if name == "center":
                    s = np.linspace(-4, 4, 81)
                    tang_l = np.array([1.0, 0.0])
                    c0 = c
                else:
                    s = np.linspace(-span / 2, span / 2, 41)
                    tang_l = tang
                    c0 = c
                hz = np.zeros(len(s), dtype=np.complex128)
                ex = np.zeros(len(s), dtype=np.complex128)
                ey = np.zeros(len(s), dtype=np.complex128)
                for i, si in enumerate(s):
                    xy = c0 + float(si) * tang_l
                    k = int(np.argmin(np.sum((points - xy) ** 2, axis=1)))
                    hz[i], ex[i], ey[i] = x[k], Ex[k], Ey[k]
                saved[f"{name}_s"] = s
                saved[f"{name}_Hz"] = hz
                saved[f"{name}_Ex"] = ex
                saved[f"{name}_Ey"] = ey
                saved[f"{name}_outward"] = n_hat
            np.savez_compressed(OUT / "fem_horns_only_fields_FEM-M.npz", **saved)
    (OUT / "fem_horns_only.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
