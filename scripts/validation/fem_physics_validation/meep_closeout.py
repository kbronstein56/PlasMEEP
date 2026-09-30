#!/usr/bin/env python3
"""Meep cross-checks against the analytic field. Registrations are not averaged.

Select a task with MEEP_TASK: guide, coated, pair, three, seven.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import meep as mp
from mpi4py import MPI

import sixport_common as sc  # noqa: E402
from analytic_maxwell import cluster_field, kz_of, solve_clusters  # noqa: E402
from analytic_sweeps import R_CORE, R_GAP, R_SHELL, plasma  # noqa: E402
from scipy.special import hankel1

OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"
RANK = MPI.COMM_WORLD.Get_rank()


def hz_of(sim, dft):
    return complex(np.asarray(sim.get_dft_array(dft, mp.Hz, 0)).ravel()[0])


def plasma_medium():
    return mp.Medium(
        epsilon=1.0,
        E_susceptibilities=[mp.DrudeSusceptibility(frequency=float(sc.fp_a), gamma=float(sc.gamma_a), sigma=1.0)],
    )


def run_case(tag, resolution, shift, until, geometry_fn, centers, radius, eps, coated, probes, src):
    fs = float(sc.fs_a)
    k0 = 2 * np.pi * fs
    sh = np.asarray(shift, float)
    cell = mp.Vector3(12.0, 10.0)
    pml = [mp.PML(1.0)]
    src_xy = np.asarray(src, float) + sh
    source = mp.Source(mp.GaussianSource(fs, fwidth=0.05 * fs), component=mp.Hz, center=mp.Vector3(*src_xy))
    geom = geometry_fn(sh)
    vals = {name: {} for name in probes}
    for label, geometry in (("vac", []), ("obj", geom)):
        sim = mp.Simulation(cell_size=cell, resolution=resolution, boundary_layers=pml, sources=[source], geometry=geometry)
        dfts = {}
        for name, xy in probes.items():
            p = np.asarray(xy, float) + sh
            dfts[name] = sim.add_dft_fields([mp.Hz], fs, 0, 1, center=mp.Vector3(*p), size=mp.Vector3())
        sim.run(until_after_sources=until)
        for name, dft in dfts.items():
            vals[name][label] = hz_of(sim, dft)
        sim.reset_meep()
    coat = ([R_CORE, R_GAP, R_SHELL], [complex(eps), 1.0 + 0j, 3.8 + 0j]) if coated else None
    outer = R_SHELL if coated else radius
    b = solve_clusters(centers + sh, k0, radius, complex(eps), src_xy, 10, coated=coat)
    out = {"tag": tag, "resolution": resolution, "shift": sh.tolist(), "until_after_sources": until, "probes": {}}
    for name, xy in probes.items():
        p = (np.asarray(xy, float) + sh)[None, :]
        ratio = vals[name]["obj"] / vals[name]["vac"]
        analytic = cluster_field(p, centers + sh, k0, outer, complex(eps), src_xy, b, 10)[0]
        vac = hankel1(0, k0 * np.linalg.norm(p - src_xy))
        ar = analytic / vac
        out["probes"][name] = {
            "db": float(20 * np.log10(abs(ratio) / abs(ar))) if abs(ar) > 1e-8 else None,
            "phase_deg": float(np.angle(ratio / ar) * 180 / np.pi),
            "abs_err": float(abs(ratio - ar)),
        }
    return out


def guide():
    """Propagating parallel-plate phase advance. Metallic y walls, PML in x."""
    width = 1.0
    k0 = 5.0
    fs = k0 / (2 * np.pi)
    # A y-uniform source on Neumann plates excites the m=0 mode, beta = k0.
    beta = k0 + 0j
    resolution = int(os.environ.get("GUIDE_RES", "40"))
    until = float(os.environ.get("GUIDE_UNTIL", "80"))
    cell = mp.Vector3(10.0, width)
    pml = [mp.PML(1.0, direction=mp.X)]
    src = mp.Source(mp.GaussianSource(fs, fwidth=0.04 * fs), component=mp.Hz, center=mp.Vector3(-3.0, 0), size=mp.Vector3(0, width))
    x1, x2 = 0.5, 2.5
    sim = mp.Simulation(cell_size=cell, resolution=resolution, boundary_layers=pml, sources=[src], geometry=[], dimensions=2)
    d1 = sim.add_dft_fields([mp.Hz], fs, 0, 1, center=mp.Vector3(x1, 0), size=mp.Vector3())
    d2 = sim.add_dft_fields([mp.Hz], fs, 0, 1, center=mp.Vector3(x2, 0), size=mp.Vector3())
    sim.run(until_after_sources=until)
    z1, z2 = hz_of(sim, d1), hz_of(sim, d2)
    ratio = z2 / z1
    analytic = np.exp(1j * beta * (x2 - x1))
    rec = {
        "resolution": resolution,
        "until_after_sources": until,
        "beta": [float(beta.real), float(beta.imag)],
        "db": float(20 * np.log10(abs(ratio) / abs(analytic))),
        "phase_deg": float(np.angle(ratio / analytic) * 180 / np.pi),
        "abs_err": float(abs(ratio - analytic)),
    }
    if RANK == 0:
        print("guide", rec, flush=True)
        prev = []
        path = OUT / "meep_closeout.json"
        if path.exists():
            prev = json.loads(path.read_text())
        prev.append({"task": "guide", **rec})
        path.write_text(json.dumps(prev, indent=2) + "\n")
    return rec


def save(rec):
    if RANK != 0:
        return
    path = OUT / "meep_closeout.json"
    prev = json.loads(path.read_text()) if path.exists() else []
    prev.append(rec)
    path.write_text(json.dumps(prev, indent=2) + "\n")
    print(rec["tag"], rec["probes"], flush=True)


def scatter_tasks(task):
    fs = float(sc.fs_a)
    eps = plasma(fs)
    quartz = mp.Medium(epsilon=3.8)
    plasma_m = plasma_medium()
    probes = {"forward": (2.6, 0.0), "backward": (-2.2, 0.0), "side_90": (0.0, 2.4)}
    src = (-3.6, 0.0)

    def one_disk(sh):
        c = mp.Vector3(*sh)
        return [mp.Cylinder(radius=R_CORE, center=c, material=plasma_m)]

    def one_coated(sh):
        c = mp.Vector3(*sh)
        return [
            mp.Cylinder(radius=float(sc.r_bulb_outer), center=c, material=quartz),
            mp.Cylinder(radius=float(sc.r_bulb_inner), center=c, material=mp.Medium(epsilon=1.0)),
            mp.Cylinder(radius=R_CORE, center=c, material=plasma_m),
        ]

    def centers_of(n, angle=0.0):
        if n == 1:
            return np.zeros((1, 2))
        if n == 2:
            ang = np.deg2rad(angle)
            rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
            arm = rot @ np.array([0.5, 0.0])
            return np.vstack([-arm, arm])
        if n == 3:
            return np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])
        return np.array([
            [0.0, 0.0], [-0.86602540378, 0.5], [-0.86602540378, -0.5],
            [0.86602540378, 0.5], [0.86602540378, -0.5], [0.0, 1.0], [0.0, -1.0],
        ])

    def cluster_geom(centers):
        def build(sh):
            geom = []
            for c in centers:
                cc = mp.Vector3(float(c[0] + sh[0]), float(c[1] + sh[1]))
                geom.extend([
                    mp.Cylinder(radius=R_CORE, center=cc, material=plasma_m),
                ])
            return geom
        return build

    jobs = []
    if task == "coated":
        for res, ox in ((24, 0.0), (24, 0.25), (40, 0.0), (40, 0.25)):
            jobs.append((f"coated_res{res}_ox{ox}", res, (ox / res, 0.0), 140, one_coated, np.zeros((1, 2)), R_CORE, eps, True))
    elif task == "pair":
        for angle in (0.0, 90.0, 30.0):
            centers = centers_of(2, angle)
            for res, ox in ((30, 0.0), (30, 0.25)):
                jobs.append((f"pair{angle:.0f}_res{res}_ox{ox}", res, (ox / res, 0.0), 160, cluster_geom(centers), centers, R_CORE, eps, False))
    elif task == "three":
        centers = centers_of(3)
        for res, ox in ((30, 0.0), (30, 0.25), (40, 0.0)):
            jobs.append((f"three_res{res}_ox{ox}", res, (ox / res, 0.0), 180, cluster_geom(centers), centers, R_CORE, eps, False))
    elif task == "seven":
        centers = centers_of(7)
        for res, ox in ((24, 0.0), (24, 0.25), (36, 0.0), (36, 0.25)):
            jobs.append((f"seven_res{res}_ox{ox}", res, (ox / res, 0.0), 200, cluster_geom(centers), centers, R_CORE, eps, False))
    else:
        raise SystemExit(task)
    for tag, res, shift, until, geom, centers, radius, eps_c, coated in jobs:
        rec = run_case(tag, res, shift, until, geom, centers, radius, eps_c, coated, probes, src)
        rec["task"] = task
        save(rec)


def main():
    task = os.environ.get("MEEP_TASK", "guide")
    OUT.mkdir(parents=True, exist_ok=True)
    if task == "guide":
        guide()
    else:
        scatter_tasks(task)
    if RANK == 0:
        print("MEEP_CLOSEOUT_DONE", task, flush=True)


if __name__ == "__main__":
    main()
