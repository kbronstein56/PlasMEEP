#!/usr/bin/env python3
"""
Straight parallel-plate guide: Meep FDTD vs custom FEM Poynting flux.

Same clear width, same cosine Hz line source, monitors at two stations.
Compare normalized transmission (P_down / |P_ref|) — absolute Gaussian
amplitude is removed by the ratio. Target: both ~ lossless and within 0.1 dB.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import meep as mp  # noqa: E402
import sixport_common as sc  # noqa: E402
from fem_validated_solver import assemble_anisotropic, fields_from_hz  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "stage3"

# Geometry in Meep a-units. Clear width matches the circulator feed.
CLEAR = float(sc.clear_width)  # 2.0
WALL = float(sc.wall_thickness)  # 0.2
FS = float(sc.fs_a)
# Short guide: source at x=3, ref monitor x=5, down monitor x=9, cell length 14
X0, X_SRC, X_REF, X_DOWN, X1 = 0.0, 3.0, 5.0, 9.0, 14.0
Y_MID = 0.0
DPML = 1.0
RES = 40  # cheaper than res50; still ~0.5 mm


def cosine_offsets(n=41):
    span = 0.96 * CLEAR
    s = np.linspace(-span / 2, span / 2, n)
    amp = np.cos(np.pi * s / span)
    amp = np.clip(amp, 0, None)
    amp = amp / np.sqrt(np.sum(np.abs(amp) ** 2))
    return s, amp


def run_meep() -> dict:
    s, amp = cosine_offsets()
    cell = mp.Vector3(X1 - X0, CLEAR + 4 * WALL + 2 * DPML, 0)
    # center the guide at y=0; cell y from -(half) 
    # Meep cell is centered at origin by default? cell_size is centered.
    # Place walls as PEC blocks.
    y_half = CLEAR / 2 + WALL / 2
    sources = []
    for si, ai in zip(s, amp):
        sources.append(
            mp.Source(
                mp.GaussianSource(frequency=FS, fwidth=0.10 * FS),
                component=mp.Hz,
                center=mp.Vector3(X_SRC - X1 / 2, si, 0),
                amplitude=float(ai),
            )
        )
    geom = [
        mp.Block(
            size=mp.Vector3(X1, WALL, mp.inf),
            center=mp.Vector3(0, CLEAR / 2 + WALL / 2, 0),
            material=mp.perfect_electric_conductor,
        ),
        mp.Block(
            size=mp.Vector3(X1, WALL, mp.inf),
            center=mp.Vector3(0, -(CLEAR / 2 + WALL / 2), 0),
            material=mp.perfect_electric_conductor,
        ),
    ]
    sim = mp.Simulation(
        cell_size=mp.Vector3(X1, CLEAR + 2 * WALL + 2 * DPML, 0),
        geometry=geom,
        sources=sources,
        resolution=RES,
        boundary_layers=[mp.PML(DPML)],
        symmetries=[],
    )
    # flux lines across the guide, size in y
    span = 0.96 * CLEAR

    def add(x):
        # cell is centered: physical x_meep = x_phys - X1/2
        return sim.add_flux(
            FS,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(x - X1 / 2, 0, 0),
                size=mp.Vector3(0, span, 0),
                direction=mp.X,
            ),
        )

    f_ref = add(X_REF)
    f_down = add(X_DOWN)
    # also a monitor upstream of source (should be opposite / smaller if directional)
    f_up = add(X_SRC - 1.5)
    t0 = time.perf_counter()
    sim.run(until_after_sources=40)
    wall = time.perf_counter() - t0
    pr, pd, pu = [float(mp.get_fluxes(f)[0]) for f in (f_ref, f_down, f_up)]
    return {
        "P_ref": pr,
        "P_down": pd,
        "P_upstream": pu,
        "T_down_over_ref": pd / pr if pr != 0 else None,
        "wall_s": wall,
        "res": RES,
        "fwidth": 0.10 * FS,
    }


def fem_mesh():
    h = 1.0 / RES
    xs = np.arange(0, X1 + h, h)
    ys = np.arange(-(CLEAR / 2 + WALL + DPML), CLEAR / 2 + WALL + DPML + h, h)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    tri = Delaunay(pts)
    return pts, tri.simplices


def run_fem() -> dict:
    points, tris = fem_mesh()
    T = len(tris)
    rho = [
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    ]
    cents = points[tris].mean(1)
    # PEC walls: |y| between CLEAR/2 and CLEAR/2+WALL
    y = np.abs(cents[:, 1])
    wall = (y >= CLEAR / 2) & (y <= CLEAR / 2 + WALL)
    for a in rho:
        a[wall] = 0.0
    # PML via stretch inside assemble: reuse solver pml which keys off nx_ports.
    # For this small domain, implement a local stretch by temporarily... 
    # assemble_anisotropic uses sc.nx_ports. We'll pass a custom assembly copy
    # by monkeypatching is messy. Instead pad and use sigma manually here.
    from fem_validated_solver import assemble_anisotropic as asm

    # Patch: call internal by building A with modified pml — duplicate stretch
    # by setting points that are in PML to use the existing function's domain
    # which won't match. Write a local assemble using the imported function
    # after shifting coordinates into a fake domain? Simplest: no PML, use
    # a long enough guide and measure between source and before the end
    # reflection. Run time-harmonic; open ends radiate into vacuum and reflect.
    # Add a crude PML by zeroing? Better call assemble then... 
    # Use the existing assemble which stretches near x<dpml of the BIG domain.
    # Our mesh x is 0..14, sc.dpml_ports=2, sc.nx=30 so only x<2 and x>28 get PML.
    # x in [0,2] is PML. Put source at 4 which is outside. End x=14 is NOT in PML
    # (nx-dp=28). So the right end reflects. Put down monitor at 8 and source at 4,
    # and extend mesh to x=30 so the right PML of the production stretcher engages.
    raise RuntimeError("replaced by run_fem_aligned")


def run_fem_aligned() -> dict:
    """Mesh x in [0, nx_ports] so production PML stretch applies at both ends."""
    h = 0.05  # 1 mm, cheap
    nx, ny = float(sc.nx_ports), 8.0
    y0 = sc.ny_ports / 2
    xs = np.arange(0, nx + h, h)
    ys = np.arange(y0 - 3, y0 + 3 + h, h)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    tris = Delaunay(pts).simplices
    T = len(tris)
    rho = [
        np.ones(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.zeros(T, dtype=np.complex128),
        np.ones(T, dtype=np.complex128),
    ]
    cents = pts[tris].mean(1)
    # walls parallel to x at y = y0 ± CLEAR/2
    yabs = np.abs(cents[:, 1] - y0)
    wall = (yabs >= CLEAR / 2) & (yabs <= CLEAR / 2 + WALL)
    for a in rho:
        a[wall] = 0.0
    pec = np.zeros(len(pts), dtype=bool)
    k0 = 2 * np.pi * FS
    A = assemble_anisotropic(pts, tris, *rho, k0, pec)
    s, amp = cosine_offsets(61)
    b = np.zeros(len(pts), dtype=np.complex128)
    x_src = 8.0  # inside, away from left PML (dpml=2)
    for si, ai in zip(s, amp):
        xy = np.array([x_src, y0 + si])
        i = int(np.argmin(np.sum((pts - xy) ** 2, axis=1)))
        b[i] += ai
    lu = splu(A.tocsc())
    x = lu.solve(b)
    Ex, Ey = fields_from_hz(pts, tris, x, *rho, k0)

    def flux_at(xq):
        span = 0.96 * CLEAR
        offs = np.linspace(-span / 2, span / 2, 31)
        ds = span / 30
        w = np.full(31, ds)
        w[0] *= 0.5
        w[-1] *= 0.5
        # outward for +x power: n = +x
        total = 0.0
        for si, wi in zip(offs, w):
            xy = np.array([xq, y0 + si])
            i = int(np.argmin(np.sum((pts - xy) ** 2, axis=1)))
            Hz = x[i]
            Sx = 0.5 * np.real(Ey[i] * np.conj(Hz))
            total += Sx * float(wi)
        return total

    p_ref = flux_at(x_src + 2.0)
    p_down = flux_at(x_src + 6.0)
    p_up = flux_at(x_src - 2.0)
    return {
        "P_ref": float(p_ref),
        "P_down": float(p_down),
        "P_upstream": float(p_up),
        "T_down_over_ref": float(p_down / p_ref) if p_ref != 0 else None,
        "n_nodes": int(len(pts)),
        "h_a": h,
    }


def db(r):
    if r is None or r <= 0:
        return float("nan")
    return 10 * np.log10(r)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print("FEM straight guide...", flush=True)
    fem = run_fem_aligned()
    print(fem, flush=True)
    print("Meep straight guide...", flush=True)
    meep = run_meep()
    print(meep, flush=True)
    # Compare transmission ratios in dB (each solver normalized to its own ref)
    d_fem = db(fem["T_down_over_ref"])
    d_meep = db(meep["T_down_over_ref"])
    out = {
        "fem": fem,
        "meep": meep,
        "fem_T_dB": d_fem,
        "meep_T_dB": d_meep,
        "fem_minus_meep_T_dB": d_fem - d_meep,
        "note": (
            "Each T is P_down/P_ref within that solver. "
            "Ideal lossless guide: T_dB ~ 0. Difference isolates port/solver physics "
            "from absolute source amplitude."
        ),
    }
    (OUT / "straight_guide_compare.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ("fem_T_dB", "meep_T_dB", "fem_minus_meep_T_dB")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
