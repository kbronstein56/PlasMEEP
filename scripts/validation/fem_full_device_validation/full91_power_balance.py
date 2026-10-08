#!/usr/bin/env python3
"""Closed control-volume balance for the 91-bulb device.

The PML-inner rectangle contains the horn sources, so its outward flux is not
a source-free balance. Each source segment is enclosed in a thin slot and
removed from the control volume. The control surface is the rectangle plus
the slot perimeters. Outward normals on a slot point into that slot.

Predeclared acceptance, frozen before these device runs:
    |F_rectangle + F_slots + P_abs| / (|F_rectangle| + |F_slots| + |P_abs|) <= 0.02
on mesh F at B = 0 and at both production signs. Mesh C is a debug solve and
does not loosen that threshold.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from full91_solve import (  # noqa: E402
    OUT,
    free_gib,
    inject,
    load_mesh,
    material_masks,
    port_lines,
    rectangle_flux,
    rho_of,
    volume_absorption,
)

# Frozen in the module docstring before the device numbers were computed.
REL_TOL = 0.02
SLOT_HALF_WIDTH = 0.04


def _polyline_flux(sampler, uh, rho, k0, points, normals) -> float:
    hz, ex, ey = sampler.fields(uh, *rho, k0, points)
    sx = 0.5 * np.real(ey * np.conj(hz))
    sy = -0.5 * np.real(ex * np.conj(hz))
    seg = points[1:] - points[:-1]
    mid_s = 0.5 * (sx[:-1] + sx[1:])
    mid_sy = 0.5 * (sy[:-1] + sy[1:])
    mid_n = 0.5 * (normals[:-1] + normals[1:])
    return float(np.sum((mid_s * mid_n[:, 0] + mid_sy * mid_n[:, 1]) * np.linalg.norm(seg, axis=1)))


def slot_flux(sampler, uh, rho, k0, line) -> float:
    """Outward-from-CV flux through one source slot. The slot contains the source."""
    xy = np.asarray(line["source_xy"], float)
    n_hat = np.asarray(line["outward"], float)
    p0, p1 = xy[0], xy[-1]
    hw = SLOT_HALF_WIDTH
    corners = np.vstack([
        p0 - hw * n_hat,
        p1 - hw * n_hat,
        p1 + hw * n_hat,
        p0 + hw * n_hat,
        p0 - hw * n_hat,
    ])
    center = 0.5 * (p0 + p1)
    samples = []
    normals = []
    for a, b in zip(corners[:-1], corners[1:]):
        length = float(np.linalg.norm(b - a))
        n = max(8, int(np.ceil(length / 0.01)))
        t = np.linspace(0.0, 1.0, n, endpoint=False)
        pts = a + t[:, None] * (b - a)
        edge = b - a
        normal = np.array([edge[1], -edge[0]], float)
        normal /= np.linalg.norm(normal)
        mid = 0.5 * (a + b)
        if np.dot(normal, center - mid) < 0.0:
            normal = -normal
        samples.append(pts)
        normals.append(np.repeat(normal[None, :], len(pts), axis=0))
    # Close the last edge samples up to the first corner by appending the start point
    # to the sample list of the last edge.
    pts = np.vstack(samples + [corners[:1]])
    nrm = np.vstack(normals + [normals[-1][:1]])
    return _polyline_flux(sampler, uh, rho, k0, pts, nrm)


def main():
    grade = sys.argv[1]
    b_values = [float(x) for x in sys.argv[2:]]
    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    lines = port_lines()
    rows = []
    for b_tesla in b_values:
        sc.fs_a = sc.fs_Hz * sc.a / 299792458.0
        k0 = 2 * np.pi * sc.fs_a
        rho, eps, _, _ = rho_of(points, tris, sc.fs_Hz, b_tesla, masks)
        plasma = masks[0]
        print(f"balance grade={grade} B={b_tesla} nodes={len(points)} free_GiB={free_gib():.1f}", flush=True)
        t0 = time.perf_counter()
        A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
        lu = splu(A.tocsc())
        print(f"factor_s={time.perf_counter() - t0:.1f} free_GiB={free_gib():.1f}", flush=True)
        sampler = ElementSampler(points, tris)
        for j, line in enumerate(lines):
            bvec = inject(points, tris, sampler, line["source_xy"], line["weights"])
            uh = lu.solve(bvec)
            f_rect = rectangle_flux(sampler, uh, rho, k0)
            f_slot = slot_flux(sampler, uh, rho, k0, line)
            p_abs = volume_absorption(points, tris, uh, rho, k0, eps, plasma)
            residual = f_rect + f_slot + p_abs
            denom = abs(f_rect) + abs(f_slot) + abs(p_abs)
            rel = abs(residual) / denom
            row = {
                "grade": grade,
                "B_T": b_tesla,
                "port": j,
                "F_rectangle": f_rect,
                "F_slots": f_slot,
                "P_abs": p_abs,
                "residual": residual,
                "relative": rel,
                "tol": REL_TOL,
                "pass": bool(rel <= REL_TOL),
            }
            rows.append(row)
            print(
                f"port {j} Frect={f_rect:.6e} Fslot={f_slot:.6e} Pabs={p_abs:.6e} "
                f"resid={residual:.6e} rel={rel:.4e} {'PASS' if row['pass'] else 'FAIL'}",
                flush=True,
            )
    tag = f"full91_balance_{grade}.json"
    (OUT / tag).write_text(json.dumps(rows, indent=2) + "\n")
    print("WROTE", tag, flush=True)


if __name__ == "__main__":
    main()
