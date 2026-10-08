#!/usr/bin/env python3
"""Port-observable spread from monitor placement and quadrature.

One factorization. The source line stays fixed. Variants change only the
receiving line: station along the feed, sample count, and span.
"""
from __future__ import annotations

import json
import sys
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
from full91_solve import inject, load_mesh, material_masks, port_lines, rho_of  # noqa: E402
from full91_solve import OUT  # noqa: E402


def receiving_line(port: int, shift: float, n_points: int, span_frac: float):
    base = port_lines(n_points=41)[port]
    n_hat = base["outward"]
    tang = base["tangent"]
    horn = sc.full_horns[port]
    shift_xy = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    origin = np.asarray(horn["monitor_center"], float) + shift_xy + shift * n_hat
    span = span_frac * sc.clear_width
    offsets = np.linspace(-span / 2.0, span / 2.0, n_points)
    ds = span / (n_points - 1)
    weights = np.full(n_points, ds)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    xy = np.vstack([origin + float(s) * tang for s in offsets])
    return xy, weights, n_hat


def observe(sampler, uh, rho, k0, port, shift, n_points, span_frac):
    xy, weights, n_hat = receiving_line(port, shift, n_points, span_frac)
    hz, ex, ey = sampler.fields(uh, *rho, k0, xy)
    amp = complex(np.dot(weights, hz))
    sx = 0.5 * np.real(ey * np.conj(hz))
    sy = -0.5 * np.real(ex * np.conj(hz))
    flux = float(np.dot(weights, sx * n_hat[0] + sy * n_hat[1]))
    return amp, flux


def main():
    grade = sys.argv[1]
    b_tesla = float(sys.argv[2])
    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    sc.fs_a = sc.fs_Hz * sc.a / 299792458.0
    k0 = 2 * np.pi * sc.fs_a
    rho, _, _, _ = rho_of(points, tris, sc.fs_Hz, b_tesla, masks)
    A = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    lu = splu(A.tocsc())
    sampler = ElementSampler(points, tris)
    lines = port_lines()
    rows = []
    # Source P1 only. The sixfold layout repeats the feed.
    bvec = inject(points, tris, sampler, lines[0]["source_xy"], lines[0]["weights"])
    uh = lu.solve(bvec)
    variants = []
    for port in (1, 3, 4):
        for shift in (-0.30, -0.15, 0.0, 0.15, 0.30):
            variants.append((port, shift, 41, 0.92))
        for n_points in (21, 41, 81):
            variants.append((port, 0.0, n_points, 0.92))
        for span in (0.80, 0.92, 0.98):
            variants.append((port, 0.0, 41, span))
    for port, shift, n_points, span in variants:
        amp, flux = observe(sampler, uh, rho, k0, port, shift, n_points, span)
        rows.append({
            "port": port,
            "shift_a": shift,
            "n_points": n_points,
            "span_frac": span,
            "amp_re": amp.real,
            "amp_im": amp.imag,
            "amp_abs": abs(amp),
            "amp_phase_deg": float(np.angle(amp) * 180 / np.pi),
            "flux": flux,
        })
        print(rows[-1], flush=True)
    payload = {"grade": grade, "B_T": b_tesla, "rows": rows}
    name = f"full91_port_uncertainty_{grade}_B{b_tesla:+.4f}.json".replace("+", "p").replace("-", "m")
    (OUT / name).write_text(json.dumps(payload, indent=2) + "\n")
    print("WROTE", name, flush=True)


if __name__ == "__main__":
    main()
