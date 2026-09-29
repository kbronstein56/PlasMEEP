#!/usr/bin/env python3
"""Normal and 60-degree PML reflection: Meep default PML vs the FEM stretch."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import meep as mp  # noqa: E402
import sixport_common as sc  # noqa: E402
from boundary_diagnosis import K0, MEEP_SIGMA_MAX, assemble, solve_system  # noqa: E402
from fem_validated_solver import fields_from_hz  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "boundary_pml"
FS = float(sc.fs_a)


def reflection_from_line(s, hz, beta):
    """Fit A e^{i β s} + B e^{-i β s}. Return power reflection |B/A|^2 in dB."""
    cols = np.column_stack([np.exp(1j * beta * s), np.exp(-1j * beta * s)])
    coef, *_ = np.linalg.lstsq(cols, hz, rcond=None)
    A, B = coef
    ratio = abs(B) / (abs(A) + 1e-30)
    return {
        "Gamma_field": float(ratio),
        "R_power_dB": float(20 * np.log10(ratio + 1e-30)),
        "|A|": float(abs(A)),
        "|B|": float(abs(B)),
    }


def fem_guide(length, height, h, dp, sigma_max, angle_deg, width):
    xs = np.arange(0, length + 0.5 * h, h)
    ys = np.arange(0, height + 0.5 * h, h)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    _, keep = np.unique(np.round(pts, 8), axis=0, return_index=True)
    pts = pts[np.sort(keep)]
    tris = Delaunay(pts).simplices
    cents = pts[tris].mean(1)
    th = 0.25
    ang = np.deg2rad(angle_deg)
    direction = np.array([np.cos(ang), np.sin(ang)])
    normal = np.array([-direction[1], direction[0]])
    origin = np.array([2.5, height / 2.0])
    rel = cents - origin
    along = rel @ direction
    across = rel @ normal
    wall = (np.abs(np.abs(across) - (width / 2 + th / 2)) <= th / 2 + 0.05) & (along > -1.0) & (along < length + 2.0)
    keep_el = ~wall
    A, empty = assemble(pts, tris, K0, keep_el, np.ones(len(tris)), length, height, dp, sigma_max)
    # Source across the guide.
    b = np.zeros(len(pts), dtype=np.complex128)
    s = np.linspace(-0.45 * width, 0.45 * width, 25)
    amp = np.cos(np.pi * s / width)
    amp = np.clip(amp, 0, None)
    amp /= np.sqrt(np.sum(amp**2))
    src = origin + 1.5 * direction
    for si, ai in zip(s, amp):
        xy = src + si * normal
        d2 = np.sum((pts - xy) ** 2, axis=1)
        d2 = d2.copy()
        d2[empty] = np.inf
        b[int(np.argmin(d2))] += ai
    x = solve_system(A, b, np.zeros(len(pts), dtype=bool))[0]
    ones = np.ones(len(tris), dtype=np.complex128)
    z = np.zeros(len(tris), dtype=np.complex128)
    # Sample before the PML along the guide.
    pml_hit = (length - dp - origin[0]) / max(direction[0], 0.15)
    s0 = 3.0
    s1 = min(pml_hit - 1.2, s0 + 8.0)
    ss = np.linspace(s0, s1, 80)
    hz = np.zeros(len(ss), dtype=np.complex128)
    for i, sv in enumerate(ss):
        xy = origin + sv * direction
        hz[i] = x[int(np.argmin(np.sum((pts - xy) ** 2, axis=1)))]
    beta = float(np.sqrt(max(K0**2 - (np.pi / width) ** 2, 0.05**2)))
    fit = reflection_from_line(ss, hz, beta)
    # Decay inside the PML: samples past the inner face, along x at mid-height if angle is 0.
    return fit | {"s0": float(s0), "s1": float(s1), "beta": beta, "n": int(len(pts))}


def meep_guide(length, height, res, dp, angle_deg, width):
    ang = np.deg2rad(angle_deg)
    direction = np.array([np.cos(ang), np.sin(ang)])
    normal = np.array([-direction[1], direction[0]])
    # Meep cell is centered. Shift geometry by -L/2.
    origin = np.array([2.5, height / 2.0])
    th = 0.25

    def block_at(along0, along1, side):
        c0 = origin + 0.5 * (along0 + along1) * direction + side * (width / 2 + th / 2) * normal
        # Axis-aligned block is wrong for a rotated guide. Use a prism.
        p0 = origin + along0 * direction + side * (width / 2) * normal
        p1 = origin + along0 * direction + side * (width / 2 + th) * normal
        p2 = origin + along1 * direction + side * (width / 2 + th) * normal
        p3 = origin + along1 * direction + side * (width / 2) * normal
        verts = []
        for p in (p0, p1, p2, p3):
            verts.append(mp.Vector3(float(p[0] - length / 2), float(p[1] - height / 2), 0))
        return mp.Prism(vertices=verts, height=mp.inf, axis=mp.Vector3(0, 0, 1), material=mp.perfect_electric_conductor)

    along_end = length + 2.0
    geom = [block_at(0.0, along_end, +1), block_at(0.0, along_end, -1)]
    sources = []
    s = np.linspace(-0.45 * width, 0.45 * width, 25)
    amp = np.clip(np.cos(np.pi * s / width), 0, None)
    amp /= np.sqrt(np.sum(amp**2))
    src = origin + 1.5 * direction
    for si, ai in zip(s, amp):
        xy = src + si * normal
        sources.append(
            mp.Source(
                mp.GaussianSource(frequency=FS, fwidth=0.05 * FS),
                component=mp.Hz,
                center=mp.Vector3(float(xy[0] - length / 2), float(xy[1] - height / 2), 0),
                amplitude=float(ai),
            )
        )
    sim = mp.Simulation(
        cell_size=mp.Vector3(length, height, 0),
        geometry=geom,
        sources=sources,
        resolution=res,
        boundary_layers=[mp.PML(thickness=dp)],
        default_material=mp.Medium(epsilon=1),
    )
    pml_hit = (length - dp - origin[0]) / max(direction[0], 0.15)
    s0, s1 = 3.0, min(pml_hit - 1.2, 11.0)
    ss = np.linspace(s0, s1, 80)
    # One DFT point region as a line of fluxes is heavier; sample with a thin volume along the guide
    # by adding a set of point volumes. A single rotated volume is awkward, so use many point DFTs.
    dfts = []
    for sv in ss:
        xy = origin + sv * direction
        dfts.append(
            sim.add_dft_fields(
                [mp.Hz],
                FS,
                0,
                1,
                where=mp.Volume(
                    center=mp.Vector3(float(xy[0] - length / 2), float(xy[1] - height / 2), 0),
                    size=mp.Vector3(0, 0, 0),
                ),
            )
        )
    sim.run(until_after_sources=40)
    hz = np.array([complex(np.array(sim.get_dft_array(d, mp.Hz, 0)).reshape(-1)[0]) for d in dfts])
    beta = float(np.sqrt(max(K0**2 - (np.pi / width) ** 2, 0.05**2)))
    return reflection_from_line(ss, hz, beta) | {"s0": s0, "s1": float(s1), "beta": beta}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cases = []
    specs = [
        ("normal", 0.0, 8.0, 18.0, 12.0),
        ("oblique60", 60.0, 6.0, 20.0, 24.0),
    ]
    for name, ang, width, length, height in specs:
        print("FEM", name, flush=True)
        for sig, label in ((2.0, "fem_sigma2"), (MEEP_SIGMA_MAX, "fem_meep_sigma")):
            fit = fem_guide(length, height, 0.1, 2.0, sig, ang, width)
            cases.append({"case": name, "solver": label, **fit})
            print(label, fit["R_power_dB"], flush=True)
        print("Meep", name, flush=True)
        fit = meep_guide(length, height, 16, 2.0, ang, width)
        cases.append({"case": name, "solver": "meep", **fit})
        print("meep", fit["R_power_dB"], flush=True)
    (OUT / "pml_reflection.json").write_text(json.dumps(cases, indent=2) + "\n")
    print(json.dumps(cases, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
