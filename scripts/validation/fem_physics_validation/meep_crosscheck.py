#!/usr/bin/env python3
"""Short Meep cross-checks against analytic references.

Nondispersive cases are run to a fixed post-source time. Dispersive plasma
is not re-run here; the earlier long DFT campaign is the evidence for that.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import meep as mp
from mpi4py import MPI

import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, solve_clusters, stack_response  # noqa: E402
from scipy.special import hankel1

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
RANK = MPI.COMM_WORLD.Get_rank()


def point_dft(sim, center, fs):
    dft = sim.add_dft_fields([mp.Hz], fs, 0, 1, center=mp.Vector3(*center), size=mp.Vector3())
    return dft


def hz_of(sim, dft):
    arr = sim.get_dft_array(dft, mp.Hz, 0)
    # Point monitor is a length-1 array on every rank.
    return complex(np.asarray(arr).ravel()[0])


def quartz_slab(resolution: int) -> dict:
    fs = float(sc.fs_a)
    cell = mp.Vector3(12.0, 0.6)
    pml = [mp.PML(1.0, direction=mp.X)]
    src = mp.Source(
        mp.GaussianSource(fs, fwidth=0.08 * fs),
        component=mp.Hz,
        center=mp.Vector3(-3.2, 0),
        size=mp.Vector3(0, 0.6),
    )
    quartz = mp.Medium(epsilon=3.8)
    geom = [mp.Block(size=mp.Vector3(0.40, mp.inf, mp.inf), center=mp.Vector3(0.20, 0), material=quartz)]
    probe = (1.6, 0.0)
    vals = {}
    for tag, geometry in (("vac", []), ("slab", geom)):
        sim = mp.Simulation(cell_size=cell, resolution=resolution, boundary_layers=pml, sources=[src], geometry=geometry, force_complex_fields=False)
        dft = point_dft(sim, probe, fs)
        sim.run(until_after_sources=60)
        vals[tag] = hz_of(sim, dft)
        sim.reset_meep()
    ratio = vals["slab"] / vals["vac"]
    analytic = stack_response(2 * np.pi * fs, [(0.40, 3.8 + 0j)])["t_over_vacuum"]
    return {
        "resolution": resolution,
        "ratio": [ratio.real, ratio.imag],
        "analytic_t": [analytic.real, analytic.imag],
        "db": float(20 * np.log10(abs(ratio) / abs(analytic))),
        "phase_deg": float(np.angle(ratio / analytic) * 180 / np.pi),
        "abs_err": float(abs(ratio - analytic)),
    }


def dielectric_cylinder(resolution: int) -> dict:
    fs = float(sc.fs_a)
    k0 = 2 * np.pi * fs
    cell = mp.Vector3(12.0, 10.0)
    pml = [mp.PML(1.0)]
    src_xy = (-3.6, 0.0)
    src = mp.Source(mp.GaussianSource(fs, fwidth=0.08 * fs), component=mp.Hz, center=mp.Vector3(*src_xy))
    cyl = [mp.Cylinder(radius=0.23, center=mp.Vector3(), material=mp.Medium(epsilon=3.8))]
    probes = {"forward": (2.4, 0.0), "backward": (-2.2, 0.0), "side_90": (0.0, 2.4)}
    vals = {k: {} for k in probes}
    for tag, geometry in (("vac", []), ("obj", cyl)):
        sim = mp.Simulation(cell_size=cell, resolution=resolution, boundary_layers=pml, sources=[src], geometry=geometry)
        dfts = {name: point_dft(sim, xy, fs) for name, xy in probes.items()}
        sim.run(until_after_sources=80)
        for name, dft in dfts.items():
            vals[name][tag] = hz_of(sim, dft)
        sim.reset_meep()
    b = solve_clusters(np.zeros((1, 2)), k0, 0.23, 3.8 + 0j, np.array(src_xy), 10)
    out = {"resolution": resolution, "until_after_sources": 80, "probes": {}}
    for name, xy in probes.items():
        ratio = vals[name]["obj"] / vals[name]["vac"]
        pt = np.array(xy, float)[None, :]
        analytic = cluster_field(pt, np.zeros((1, 2)), k0, 0.23, 3.8 + 0j, np.array(src_xy), b, 10)[0]
        vac = hankel1(0, k0 * np.linalg.norm(pt - np.array(src_xy)))
        ar = analytic / vac
        out["probes"][name] = {
            "db": float(20 * np.log10(abs(ratio) / abs(ar))),
            "phase_deg": float(np.angle(ratio / ar) * 180 / np.pi),
            "abs_err": float(abs(ratio - ar)),
        }
    return out


def main():
    rows = {"slab": [], "dielectric_cylinder": []}
    for res in (24, 40):
        rec = quartz_slab(res)
        rows["slab"].append(rec)
        if RANK == 0:
            print("slab", rec, flush=True)
    for res in (24, 40):
        rec = dielectric_cylinder(res)
        rows["dielectric_cylinder"].append(rec)
        if RANK == 0:
            print("cyl", rec, flush=True)
    if RANK == 0:
        (OUT / "meep_crosscheck.json").write_text(json.dumps(rows, indent=2) + "\n")
        print("MEEP_DONE", flush=True)


if __name__ == "__main__":
    main()
