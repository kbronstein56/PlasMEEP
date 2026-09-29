#!/usr/bin/env python3
"""Bare-disk Drude boundary-ring campaign.

Resolution is points per centimeter, the same quantity previously called
25 / 35 / 50 in the interface diagnosis. It is not cells per free-space
wavelength. See resolution_record().

Tasks (argv):
  resolution     write the resolution table and exit
  replay-fem     confirm saved FEM JSON pairs convert to complex numbers
  mie-scan       point-source Mie spectrum of the bare disk
  meep           RING_TASK=harminv|dft  (MPI)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402
from physical_units import physical_resolution_report  # noqa: E402
from plasma_interface import (  # noqa: E402
    DPML,
    LX,
    LY,
    R_DISK,
    SRC_XY,
    mie_coefficients,
)

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "plasma_ring_campaign"
C_MPS = 2.99792458e8
# Hotspot azimuth measured on the 50 points/cm disk DFT.
HOT_XY = np.array([-0.205, -0.095])


def resolution_record(points_per_cm: float) -> dict:
    rec = physical_resolution_report(points_per_cm=points_per_cm, a_m=float(sc.a), fs_Hz=float(sc.fs_Hz))
    rec["meaning"] = (
        "points_per_cm is the physical grid density. "
        "dx_mm = 10/points_per_cm. "
        "meep_resolution = a_m/dx_m = 2*points_per_cm when a_m=0.020 m. "
        "pixels_per_free_space_lambda0 is cells per vacuum wavelength at 3.85 GHz, "
        "and is not the 25/35/50 label."
    )
    return rec


def freq_GHz(f_a: float) -> float:
    return float(f_a) * C_MPS / float(sc.a) / 1e9


def drude_eps(f_a: float) -> complex:
    eps, _eta = gyrotropic_drude_eps_eta(float(f_a), float(sc.fp_a), float(sc.gamma_a), 0.0)
    return complex(eps)


def write_resolution_table() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [resolution_record(p) for p in (25, 35, 40, 45, 50, 60)]
    (OUT / "resolution.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(rows, indent=2), flush=True)


def replay_fem() -> None:
    """Exercise the list-to-complex conversion that used to crash the FEM print."""
    src = ROOT / "outputs" / "validation" / "fem_meep_validation" / "plasma_interface"
    for kind in ("disk", "square", "bulb"):
        rows = json.loads((src / f"fem_{kind}.json").read_text())
        assert len(rows) >= 2, kind
        for key in ("T_forward", "R_back", "inside_over_vac", "just_out_over_vac"):
            a = complex(*rows[0][key])
            b = complex(*rows[-1][key])
            if abs(a) == 0:
                raise RuntimeError(f"{kind} {key} is zero")
            db = float(20 * np.log10(abs(b) / abs(a)))
            print(f"{kind} {rows[0]['level']}->{rows[-1]['level']} {key} {db:+.4f} dB", flush=True)
        for row in rows:
            npz = np.load(src / f"fem_{kind}_{row['level']}.npz")
            if not np.isfinite(npz["Hz_disk"]).all():
                raise RuntimeError(f"{kind} {row['level']} field is not finite")
    print("REPLAY_FEM_OK", flush=True)


def mie_forward_ratio(f_a: float, nmax: int = 30) -> complex:
    """Point-source Mie ratio Hz(disk)/Hz(vacuum) at x=+3, source at x=-6."""
    from scipy.special import hankel1, jv

    k0 = 2 * np.pi * float(f_a)
    eps = drude_eps(f_a)
    ns, coeff = mie_coefficients(k0, eps, R_DISK, nmax=nmax)
    kp = k0 * np.sqrt(complex(eps))
    rs = 6.0
    x = 3.0
    cn = np.array([hankel1(int(n), k0 * rs) * ((-1) ** int(n)) for n in ns])
    vac = hankel1(0, k0 * abs(x + rs))
    tot = 0j
    for n, a, c in zip(ns, coeff, cn):
        tot += c * (jv(int(n), k0 * x) + a * hankel1(int(n), k0 * x))
    return complex(tot / vac)


def mie_scan() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # Ordinary-frequency grid covering fs, the source band, and the plasma frequency.
    freqs = np.unique(
        np.concatenate(
            [
                np.linspace(0.15, 0.90, 76),
                np.array([float(sc.fs_a), float(sc.fp_a)]),
            ]
        )
    )
    rows = []
    for f in freqs:
        eps = drude_eps(float(f))
        k0 = 2 * np.pi * float(f)
        _ns, coeff = mie_coefficients(k0, eps, R_DISK, nmax=20)
        ratio = mie_forward_ratio(float(f), nmax=24)
        rows.append(
            {
                "f_a": float(f),
                "f_GHz": freq_GHz(float(f)),
                "eps_re": eps.real,
                "eps_im": eps.imag,
                "max_abs_a": float(np.max(np.abs(coeff))),
                "T_re": ratio.real,
                "T_im": ratio.imag,
                "T_abs": float(abs(ratio)),
                "T_deg": float(np.angle(ratio, deg=True)),
            }
        )
    (OUT / "mie_scan.json").write_text(json.dumps(rows, indent=2) + "\n")
    peak = max(rows, key=lambda r: r["max_abs_a"])
    at_fs = min(rows, key=lambda r: abs(r["f_a"] - float(sc.fs_a)))
    print(
        "MIE_PEAK",
        json.dumps({"f_GHz": peak["f_GHz"], "max_abs_a": peak["max_abs_a"], "eps_re": peak["eps_re"]}),
        flush=True,
    )
    print(
        "MIE_FS",
        json.dumps({"f_GHz": at_fs["f_GHz"], "T_abs": at_fs["T_abs"], "max_abs_a": at_fs["max_abs_a"]}),
        flush=True,
    )


def _probes(dx: float, shift: np.ndarray) -> dict:
    direction = HOT_XY / np.linalg.norm(HOT_XY)
    specs = {
        "center": np.array([0.0, 0.0]),
        "inside": np.array([R_DISK - 2 * dx, 0.0]),
        "outside": np.array([R_DISK + 2 * dx, 0.0]),
        "outside_far": np.array([R_DISK + 6 * dx, 0.0]),
        "forward": np.array([3.0, 0.0]),
        "vacuum": np.array([0.0, 2.0]),
        "inside_hot": direction * (R_DISK - 2 * dx),
        "outside_hot": direction * (R_DISK + 2 * dx),
    }
    return {name: (xy + shift).tolist() for name, xy in specs.items()}


def _modes_to_json(modes) -> list:
    out = []
    for m in modes:
        out.append(
            {
                "f_a": float(m.freq),
                "f_GHz": freq_GHz(m.freq),
                "decay_a": float(m.decay),
                "Q": float(m.Q),
                "amp_re": float(np.real(m.amp)),
                "amp_im": float(np.imag(m.amp)),
                "amp_abs": float(abs(m.amp)),
                "err": float(np.abs(m.err)),
            }
        )
    out.sort(key=lambda r: r["amp_abs"], reverse=True)
    return out


def meep_main() -> None:
    import meep as mp
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.Get_rank()
    task = os.environ.get("RING_TASK", "harminv")
    ppc = float(os.environ.get("RING_PPC", "50"))
    until = float(os.environ.get("RING_UNTIL", "120"))
    res_info = resolution_record(ppc)
    res = int(res_info["meep_resolution"])
    dx = float(res_info["grid_spacing_a"])
    ox = float(os.environ.get("RING_OX", "0"))
    oy = float(os.environ.get("RING_OY", "0"))
    shift = np.array([ox * dx, oy * dx])
    kind = os.environ.get("RING_KIND", "disk")
    fs = float(sc.fs_a)
    fp = float(sc.fp_a)
    gamma = float(sc.gamma_a)
    plasma = mp.Medium(
        epsilon=1.0,
        E_susceptibilities=[mp.DrudeSusceptibility(frequency=fp, gamma=gamma, sigma=1.0)],
    )
    quartz = mp.Medium(epsilon=3.8)
    core = core_name()
    center = mp.Vector3(float(shift[0]), float(shift[1]))

    def plasma_shape(cc):
        if core == "circle":
            return mp.Cylinder(radius=R_DISK, center=cc, material=plasma)
        half = inscribed_half_side()
        if core == "square":
            return mp.Block(size=mp.Vector3(2.0 * half, 2.0 * half, mp.inf), center=cc, material=plasma)
        c45 = float(np.cos(np.pi / 4.0))
        s45 = float(np.sin(np.pi / 4.0))
        return mp.Block(
            size=mp.Vector3(2.0 * half, 2.0 * half, mp.inf),
            center=cc,
            e1=mp.Vector3(c45, s45, 0.0),
            e2=mp.Vector3(-s45, c45, 0.0),
            material=plasma,
        )

    if kind == "disk":
        geom = [plasma_shape(center)]
    elif kind == "bulb":
        geom = [
            mp.Cylinder(radius=float(sc.r_bulb_outer), center=center, material=quartz),
            mp.Cylinder(radius=float(sc.r_bulb_inner), center=center, material=mp.Medium(epsilon=1.0)),
            plasma_shape(center),
        ]
    elif kind.startswith("cluster"):
        nbulbs = int(kind.replace("cluster", ""))
        vac_med = mp.Medium(epsilon=1.0)
        geom = []
        for c in cluster_centers(nbulbs):
            cc = mp.Vector3(float(c[0] + shift[0]), float(c[1] + shift[1]))
            geom.extend(
                [
                    mp.Cylinder(radius=float(sc.r_bulb_outer), center=cc, material=quartz),
                    mp.Cylinder(radius=float(sc.r_bulb_inner), center=cc, material=vac_med),
                    plasma_shape(cc),
                ]
            )
    else:
        raise ValueError(kind)
    probes = _probes(dx, shift)
    if kind.startswith("cluster"):
        # 2 dx off the plasma face, so the point is not on the material boundary.
        inside_off, outside_off = _boundary_probe_offsets(dx, core)
        for i, c in enumerate(cluster_centers(int(kind.replace("cluster", "")))):
            base = np.asarray(c, float) + shift
            probes[f"b{i}_in"] = (base + inside_off).tolist()
            probes[f"b{i}_out"] = (base + outside_off).tolist()
    src = SRC_XY + shift
    core_tag = "" if core == "circle" else f"_{core}"
    tag = f"{task}_{kind}_ppc{ppc:g}_ox{ox:g}_oy{oy:g}_u{until:g}{core_tag}"

    def one(geometry, with_harminv: bool, checkpoints: bool, run_until: float):
        sim = mp.Simulation(
            cell_size=mp.Vector3(LX, LY),
            resolution=res,
            boundary_layers=[mp.PML(DPML)],
            geometry=geometry,
            sources=[
                mp.Source(
                    mp.GaussianSource(fs, fwidth=0.20 * fs),
                    component=mp.Hz,
                    center=mp.Vector3(float(src[0]), float(src[1])),
                )
            ],
            default_material=mp.Medium(epsilon=1, mu=1),
            eps_averaging=True,
        )
        vol = mp.Volume(center=mp.Vector3(float(shift[0]), float(shift[1])), size=mp.Vector3(8.0, 0.0))
        dft = sim.add_dft_fields([mp.Hz], fs, 0, 1, where=vol) if checkpoints or task == "dft" else None
        harms = {}
        if with_harminv:
            harm_names = [n for n in probes if n.startswith("b") and n.endswith("_out")]
            if not harm_names:
                harm_names = ["outside_hot", "outside", "center", "inside"]
            for name in harm_names:
                h = mp.Harminv(
                    mp.Hz,
                    mp.Vector3(*probes[name]),
                    fcen=0.55,
                    df=0.90,
                    mxbands=40,
                )
                h.Q_thresh = 0.0
                h.err_thresh = 10.0
                h.rel_err_thresh = float("inf")
                h.amp_thresh = -1.0
                h.rel_amp_thresh = -1.0
                harms[name] = h
        names = list(probes)
        times: list[float] = []
        samples: list[list[float]] = []
        ckpt_t: list[float] = []
        ckpt_x: list = []
        ckpt_hz: list = []
        trace_on = os.environ.get("RING_TRACE", "1") != "0" and bool(geometry)
        every = float(os.environ.get("RING_EVERY", "40"))
        side = "vac" if not geometry else "disk"
        ckpt_path = OUT / f"{tag}_{side}_ckpt.jsonl"
        trace_path = OUT / f"{tag}_{side}_trace.jsonl"
        trace_buf: list[dict] = []
        if rank == 0 and do_dft:
            OUT.mkdir(parents=True, exist_ok=True)
            if ckpt_path.exists():
                ckpt_path.unlink()
        if rank == 0 and trace_on:
            OUT.mkdir(parents=True, exist_ok=True)
            if trace_path.exists():
                trace_path.unlink()

        def flush_trace() -> None:
            if rank != 0 or not trace_buf:
                return
            with trace_path.open("a") as fh:
                for item in trace_buf:
                    fh.write(json.dumps(item) + "\n")
            trace_buf.clear()

        def collect(sim_):
            vals = [sim_.get_field_point(mp.Hz, mp.Vector3(*probes[n])) for n in names]
            if rank == 0:
                t_now = float(sim_.meep_time())
                times.append(t_now)
                row = [float(np.real(v)) for v in vals]
                samples.append(row)
                if trace_on:
                    trace_buf.append({"t": t_now, "hz": row})
                    if t_now - trace_buf[0]["t"] >= 2.0:
                        flush_trace()

        def grab_dft(sim_):
            if dft is None:
                return
            raw = np.asarray(sim_.get_dft_array(dft, mp.Hz, 0))
            xs, _ys, *_rest = sim_.get_array_metadata(dft_cell=dft)
            xs = np.asarray(xs, float).ravel()
            field = np.asarray(raw).ravel()
            if field.size != xs.size and raw.ndim == 2:
                # A zero-thickness volume can come back as (nx, 1) or (1, nx).
                field = raw[:, 0] if raw.shape[0] == xs.size else raw[0, :]
            if rank != 0:
                return
            ckpt_t.append(float(sim_.meep_time()))
            ckpt_x.append(xs)
            ckpt_hz.append(np.asarray(field))
            stations = {}
            for x0 in (-3.0, 0.0, 0.25, 0.30, 3.0):
                z = complex(np.interp(x0, xs, np.real(field)) + 1j * np.interp(x0, xs, np.imag(field)))
                stations[str(x0)] = [z.real, z.imag]
            line = json.dumps({"t": float(sim_.meep_time()), "x0": float(xs[0]), "x1": float(xs[-1]), "stations": stations})
            with ckpt_path.open("a") as fh:
                fh.write(line + "\n")
            print("CKPT", line, flush=True)

        steps = [mp.at_every(0.25, collect)] if trace_on else []
        if harms:
            steps.extend(mp.after_sources(h) for h in harms.values())
        if dft is not None:
            steps.append(mp.at_every(every, grab_dft))
        t0 = time.perf_counter()
        sim.run(*steps, until_after_sources=run_until)
        wall = time.perf_counter() - t0
        flush_trace()
        if dft is not None:
            grab_dft(sim)
        mode_json = {}
        if rank == 0 and harms:
            try:
                mode_json = {name: _modes_to_json(h.modes) for name, h in harms.items()}
            except Exception as exc:
                print("MODE_SERIALIZE_FAIL", type(exc).__name__, exc, flush=True)
        return {
            "wall_s": wall,
            "times": np.asarray(times, float),
            "samples": np.asarray(samples, float) if samples else np.zeros((0, len(names))),
            "names": names,
            "modes": mode_json,
            "ckpt_t": ckpt_t,
            "ckpt_x": ckpt_x,
            "ckpt_hz": ckpt_hz,
            "source_end_guess": None,
        }

    if rank == 0:
        OUT.mkdir(parents=True, exist_ok=True)
        print(
            "RING",
            json.dumps(
                {
                    "task": task,
                    "kind": kind,
                    "core": core,
                    "resolution": res_info,
                    "shift_a": shift.tolist(),
                    "until_after_sources": until,
                    "probes": probes,
                }
            ),
            flush=True,
        )
    do_harm = task in ("harminv", "both")
    do_dft = task in ("dft", "both")
    # Vacuum reference is only needed for the DFT ratio.
    vac = one([], False, do_dft, min(until, 50.0)) if do_dft else None
    disk = one(geom, do_harm, do_dft, until)
    if rank != 0:
        return
    rec = {
        "task": task,
        "kind": kind,
        "resolution": res_info,
        "offset_cells": [ox, oy],
        "shift_a": shift.tolist(),
        "until_after_sources": until,
        "probes": probes,
        "fwidth_a": 0.20 * fs,
        "fs_GHz": freq_GHz(fs),
        "fp_GHz": freq_GHz(fp),
        "wall_disk_s": disk["wall_s"],
        "modes": disk["modes"],
    }
    if vac is not None:
        rec["wall_vac_s"] = vac["wall_s"]
    (OUT / f"{tag}.json").write_text(json.dumps(rec, indent=2) + "\n")
    np.savez_compressed(
        OUT / f"{tag}_traces.npz",
        times=disk["times"],
        samples=disk["samples"],
        names=np.array(disk["names"]),
    )
    if do_dft and disk["ckpt_t"]:
        # Store the last vacuum line and every disk checkpoint.
        np.savez_compressed(
            OUT / f"{tag}_dft.npz",
            vac_x=vac["ckpt_x"][-1],
            vac_hz=vac["ckpt_hz"][-1],
            vac_t=np.array(vac["ckpt_t"]),
            disk_t=np.array(disk["ckpt_t"]),
            **{f"x_{i}": disk["ckpt_x"][i] for i in range(len(disk["ckpt_t"]))},
            **{f"hz_{i}": disk["ckpt_hz"][i] for i in range(len(disk["ckpt_t"]))},
        )
    print(f"{tag} EXIT modes={ {k: len(v) for k, v in disk['modes'].items()} }", flush=True)


def core_name() -> str:
    """circle, axis-aligned inscribed square, or the same square rotated 45 degrees."""
    name = os.environ.get("RING_CORE", "circle")
    if name not in ("circle", "square", "diamond"):
        raise ValueError(name)
    return name


def inscribed_half_side() -> float:
    """Half-side of the square whose corners lie on the production plasma circle.

    Side length L = sqrt(2) * R, so the corner radius is R and the square
    does not leave the original plasma disk.
    """
    return float(R_DISK) / np.sqrt(2.0)


def core_geometry_record() -> dict:
    r = float(R_DISK)
    half = inscribed_half_side()
    side = 2.0 * half
    corner = half * np.sqrt(2.0)
    r_in = float(sc.r_bulb_inner)
    r_out = float(sc.r_bulb_outer)
    clearance = r_in - corner
    rec = {
        "plasma_radius_a": r,
        "plasma_radius_mm": r * float(sc.a) * 1e3,
        "square_side_a": side,
        "square_side_mm": side * float(sc.a) * 1e3,
        "square_half_side_a": half,
        "circle_area_a2": float(np.pi * r * r),
        "square_area_a2": float(side * side),
        "square_over_circle_area": float((side * side) / (np.pi * r * r)),
        "corner_radius_a": float(corner),
        "quartz_inner_a": r_in,
        "quartz_outer_a": r_out,
        "min_clearance_to_quartz_inner_a": float(clearance),
        "min_clearance_to_quartz_inner_mm": float(clearance * float(sc.a) * 1e3),
        "side_reduced": False,
        "corners_on_plasma_circle": bool(abs(corner - r) < 1e-12),
        "overlaps_quartz": bool(clearance <= 0.0),
        "topology": "square plasma inside the original disk, then the production vacuum gap, then the 1 mm quartz shell",
    }
    if rec["overlaps_quartz"]:
        raise RuntimeError(f"inscribed square reaches the quartz, clearance {clearance}")
    return rec


def _boundary_probe_offsets(dx: float, core: str) -> tuple[np.ndarray, np.ndarray]:
    """Points 2 dx inside and outside the plasma face, not on the interface."""
    if core == "circle":
        return np.array([R_DISK - 2.0 * dx, 0.0]), np.array([R_DISK + 2.0 * dx, 0.0])
    half = inscribed_half_side()
    if core == "square":
        return np.array([half - 2.0 * dx, 0.0]), np.array([half + 2.0 * dx, 0.0])
    # Midpoint of the diamond edge that faces +x/+y, offset along its outward normal.
    normal = np.array([1.0, 1.0]) / np.sqrt(2.0)
    midpoint = np.array([float(R_DISK) / 2.0, float(R_DISK) / 2.0])
    return midpoint - 2.0 * dx * normal, midpoint + 2.0 * dx * normal


def _plasma_mask(origin: np.ndarray, center: np.ndarray, core: str) -> np.ndarray:
    delta = origin - center
    if core == "circle":
        return np.linalg.norm(delta, axis=1) <= float(R_DISK)
    half = inscribed_half_side()
    if core == "square":
        return (np.abs(delta[:, 0]) <= half) & (np.abs(delta[:, 1]) <= half)
    return np.abs(delta[:, 0]) + np.abs(delta[:, 1]) <= float(R_DISK)


def cluster_centers(n: int) -> np.ndarray:
    """Production hex lattice, origin at the central bulb. Distances are in units of a."""
    locs = np.load(OUT / "bulb_locs.npz")["locs"]
    dist = np.linalg.norm(locs, axis=1)
    if n == 1:
        return np.zeros((1, 2))
    if n == 7:
        return np.asarray(locs[dist <= 1.0 + 1e-8], float)
    if n == 19:
        return np.asarray(locs[dist <= 2.0 + 1e-8], float)
    if n != 3:
        raise ValueError(n)
    ring = locs[np.abs(dist - 1.0) < 1e-6]
    ang = np.arctan2(ring[:, 1], ring[:, 0])
    i = int(np.argmin(np.abs(np.angle(np.exp(1j * (ang - np.pi / 2))))))
    target = ang[i] - np.pi / 3.0
    j = int(np.argmin(np.abs(np.angle(np.exp(1j * (ang - target))))))
    centers = np.vstack([np.zeros(2), ring[i], ring[j]])
    pairs = [np.linalg.norm(centers[a] - centers[b]) for a, b in ((0, 1), (0, 2), (1, 2))]
    if any(abs(p - 1.0) > 1e-6 for p in pairs):
        raise RuntimeError(f"3-bulb spacing is not the lattice pitch: {pairs}")
    return centers


def fem_cluster(n: int) -> None:
    """Body-fitted FEM for N production bulbs. Same point source as the bare disk."""
    from fem_validated_solver import (
        LAST_OPERATOR,
        ElementSampler,
        assemble_anisotropic,
        eps_tensor_at_bias,
        meep_sigma_max,
        rho_from_eps,
        solve_system,
    )
    from plasma_interface import (
        _tri,
        field_ratios,
        mesh_to_origin,
        origin_to_mesh,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    core = core_name()
    geom_rec = core_geometry_record()
    centers = cluster_centers(n)
    stem = f"cluster{n}" if core == "circle" else f"cluster{n}_{core}"
    (OUT / f"{stem}_centers.json").write_text(
        json.dumps({"n": n, "core": core, "centers_a": centers.tolist(), "geometry": geom_rec}, indent=2) + "\n"
    )
    (OUT / "square_core_geometry.json").write_text(json.dumps(geom_rec, indent=2) + "\n")
    tr = _tri()
    r_out = float(sc.r_bulb_outer)
    r_in = float(sc.r_bulb_inner)
    eps_p = complex(eps_tensor_at_bias(0.0)[0])
    rho_p = rho_from_eps(eps_p, 0j, 0j, eps_p)
    if os.environ.get("RING_FEM_LEVEL"):
        levels = [(
            os.environ["RING_FEM_LEVEL"],
            float(os.environ["RING_FEM_H"]),
            float(os.environ["RING_FEM_HE"]),
        )]
        prior = OUT / f"fem_{stem}.json"
        rows = json.loads(prior.read_text()) if prior.exists() else []
    else:
        levels = [("FEM-C", 0.030, 0.012), ("FEM-F", 0.018, 0.007)]
        rows = []
    for name, h, h_edge in levels:
        print(f"=== FEM cluster{n} {name} ===", flush=True)
        verts: list[list[float]] = []
        segs: list[list[int]] = []
        index: dict[tuple[float, float], int] = {}

        def add_vertex(p) -> int:
            key = (round(float(p[0]), 8), round(float(p[1]), 8))
            if key in index:
                return index[key]
            index[key] = len(verts)
            verts.append([key[0], key[1]])
            return index[key]

        def add_chain(pts, closed: bool) -> None:
            ids = [add_vertex(p) for p in pts]
            pairs = list(zip(ids, ids[1:]))
            if closed:
                pairs.append((ids[-1], ids[0]))
            for a, b in pairs:
                if a != b:
                    segs.append([a, b])

        def linspace_edge(p0, p1, spacing):
            dist = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
            m = max(2, int(np.ceil(dist / spacing)))
            t = np.linspace(0.0, 1.0, m)
            return (1 - t)[:, None] * np.asarray(p0, float) + t[:, None] * np.asarray(p1, float)

        corners = [(0.0, 0.0), (LX, 0.0), (LX, LY), (0.0, LY)]
        for i in range(4):
            add_chain(linspace_edge(corners[i], corners[(i + 1) % 4], h), closed=False)
        half = inscribed_half_side()
        for c in centers:
            mesh_c = origin_to_mesh(c)
            if core == "circle":
                radii = (R_DISK, r_in, r_out)
            else:
                radii = (r_in, r_out)
            for radius in radii:
                m = max(32, int(np.ceil(2 * np.pi * radius / h_edge)))
                ang = np.linspace(0.0, 2 * np.pi, m, endpoint=False)
                ring = mesh_c + radius * np.column_stack([np.cos(ang), np.sin(ang)])
                add_chain(ring, closed=True)
            if core != "circle":
                local = np.array([[-half, -half], [half, -half], [half, half], [-half, half]], float)
                if core == "diamond":
                    rot = np.array([[np.cos(np.pi / 4), -np.sin(np.pi / 4)], [np.sin(np.pi / 4), np.cos(np.pi / 4)]])
                    local = local @ rot.T
                corners_sq = mesh_c + local
                loop = []
                for i in range(4):
                    edge = linspace_edge(corners_sq[i], corners_sq[(i + 1) % 4], h_edge)
                    loop.append(edge[:-1])
                add_chain(np.vstack(loop), closed=True)
        amax = 0.45 * h * h
        mesh = tr.triangulate(
            {"vertices": np.asarray(verts, float), "segments": np.asarray(segs, np.int32)},
            f"pq20a{amax:.8f}",
        )
        pts = np.asarray(mesh["vertices"], float)
        tris = np.asarray(mesh["triangles"], int)
        print(f"  nodes {len(pts)} tris {len(tris)}", flush=True)
        sc.nx_ports, sc.ny_ports, sc.dpml_ports = LX, LY, DPML
        k0 = 2 * np.pi * float(sc.fs_a)
        pec = np.zeros(len(pts), dtype=bool)
        sampler = ElementSampler(pts, tris)
        src = origin_to_mesh(SRC_XY)
        bvec = np.zeros(len(pts), dtype=np.complex128)
        t_src = int(sampler.locate(src[None, :])[0])
        w = sampler._bary(t_src, src)
        for k in range(3):
            bvec[int(tris[t_src, k])] += complex(w[k])
        out = {}
        for filled, tag in ((False, "vac"), (True, "disk")):
            T = len(tris)
            rho = [
                np.ones(T, dtype=np.complex128),
                np.zeros(T, dtype=np.complex128),
                np.zeros(T, dtype=np.complex128),
                np.ones(T, dtype=np.complex128),
            ]
            if filled:
                origin = mesh_to_origin(pts[tris].mean(1))
                for c in centers:
                    rad = np.linalg.norm(origin - c, axis=1)
                    plasma = _plasma_mask(origin, c, core)
                    quartz = (rad > r_in) & (rad <= r_out)
                    if np.any(plasma & quartz):
                        raise RuntimeError("plasma square overlaps the quartz shell")
                    rho[0][plasma] = rho_p[0]
                    rho[3][plasma] = rho_p[3]
                    rho[0][quartz] = 1.0 / 3.8
                    rho[3][quartz] = 1.0 / 3.8
            A = assemble_anisotropic(pts, tris, *rho, k0, pec)
            x, _lu, t_fac, t_sol, resid = solve_system(A, bvec, pec)
            if tag == "disk":
                op = dict(LAST_OPERATOR)
                if abs(op["sigma_max"] - meep_sigma_max(DPML)) > 1e-9:
                    raise RuntimeError(f"sigma_max mismatch {op}")
            xy = np.column_stack([np.linspace(-4.0, 4.0, 801), np.zeros(801)])
            hz, ex, ey = sampler.fields(x, *rho, k0, origin_to_mesh(xy))
            out[f"Hz_{tag}"] = hz
            out[f"Ex_{tag}"] = ex
            out[f"Ey_{tag}"] = ey
            out["factor_s"] = t_fac
            out["solve_s"] = t_sol
            out["residual"] = resid
        ratios = field_ratios(xy[:, 0], out["Hz_vac"], out["Hz_disk"])
        rec = {
            "n": n,
            "core": core,
            "level": name,
            "h": h,
            "h_edge": h_edge,
            "n_nodes": int(len(pts)),
            "n_tris": int(len(tris)),
            "centers_a": centers.tolist(),
            "geometry": geom_rec,
            "eps": [eps_p.real, eps_p.imag],
            **ratios,
        }
        print(name, "DOFs", rec["n_nodes"], "T", rec["T_forward"], flush=True)
        rows.append(rec)
        np.savez_compressed(
            OUT / f"fem_{stem}_{name}.npz",
            x=xy[:, 0],
            Hz_vac=out["Hz_vac"],
            Hz_disk=out["Hz_disk"],
        )
        (OUT / f"fem_{stem}.json").write_text(json.dumps(rows, indent=2) + "\n")
    if len(rows) == 2:
        a = complex(*rows[0]["T_forward"])
        b = complex(*rows[1]["T_forward"])
        print(f"cluster{n} T fine/coarse dB", float(20 * np.log10(abs(b) / abs(a))), flush=True)


def staircase() -> None:
    """Yee-point inclusion of each plasma core, using Meep's own cylinder test."""
    import meep as mp

    ppc = float(os.environ.get("RING_PPC", "35"))
    info = resolution_record(ppc)
    res = int(info["meep_resolution"])
    dx = float(info["grid_spacing_a"])
    centers0 = cluster_centers(7)
    offsets = [(0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.25, 0.25)]
    sim = mp.Simulation(
        cell_size=mp.Vector3(LX, LY),
        resolution=res,
        boundary_layers=[mp.PML(DPML)],
        geometry=[],
        default_material=mp.Medium(epsilon=1.0),
    )
    sim.init_sim()
    grids = {}
    region = mp.Vector3(LX - 1e-6, LY - 1e-6)
    for comp, name in ((mp.Ex, "Ex"), (mp.Ey, "Ey"), (mp.Hz, "Hz")):
        arr = np.asarray(sim.get_array(center=mp.Vector3(), size=region, component=comp))
        _dims, lo, hi = sim.get_array_slice_dimensions(comp, center=mp.Vector3(), size=region)
        # Meep returns the 2D slice as (nx, ny) in this build.
        nx, ny = int(arr.shape[0]), int(arr.shape[1])
        xs = np.linspace(float(lo.x), float(hi.x), nx)
        ys = np.linspace(float(lo.y), float(hi.y), ny)
        grids[name] = (xs, ys)
        print(
            f"GRID {name} nx={nx} ny={ny} x0={xs[0]:.6f} y0={ys[0]:.6f} "
            f"dx={float(np.median(np.diff(xs))):.6f} dy={float(np.median(np.diff(ys))):.6f}",
            flush=True,
        )
    true_area = float(np.pi * R_DISK * R_DISK)
    rows = []
    for ox, oy in offsets:
        shift = np.array([ox * dx, oy * dx])
        for comp, (xs, ys) in grids.items():
            xx, yy = np.meshgrid(xs, ys, indexing="xy")
            for i, c0 in enumerate(centers0):
                c = np.asarray(c0, float) + shift
                window = (np.abs(xx - c[0]) < R_DISK + 4 * dx) & (np.abs(yy - c[1]) < R_DISK + 4 * dx)
                pts = np.column_stack([xx[window], yy[window]])
                cyl = mp.Cylinder(radius=R_DISK, center=mp.Vector3(float(c[0]), float(c[1])), material=mp.Medium(epsilon=1))
                inside_mask = np.array(
                    [bool(mp.is_point_in_object(mp.Vector3(float(p[0]), float(p[1])), cyl)) for p in pts],
                    dtype=bool,
                )
                inside = pts[inside_mask]
                area = float(inside.shape[0] * dx * dx)
                if inside.shape[0] == 0:
                    raise RuntimeError(f"no interior points bulb {i} {comp} offset {(ox, oy)}")
                centroid = inside.mean(axis=0) - c
                rel = inside - c
                radius = np.linalg.norm(rel, axis=1)
                # A boundary cell is an interior sample with an exterior grid neighbor.
                occupied = {(round(p[0] / dx), round(p[1] / dx)) for p in inside}
                boundary = []
                for p, r in zip(inside, radius):
                    ix = round(p[0] / dx)
                    iy = round(p[1] / dx)
                    if any((ix + a, iy + b) not in occupied for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1))):
                        boundary.append((p[0] - c[0], p[1] - c[1], r - R_DISK))
                boundary = np.asarray(boundary, float)
                theta = np.arctan2(boundary[:, 1], boundary[:, 0])
                # First angular moment of the radial error. Zero for a centered staircase.
                moment = np.mean((boundary[:, 2]) * np.exp(1j * theta))
                rows.append(
                    {
                        "ppc": ppc,
                        "dx_a": dx,
                        "offset_cells": [ox, oy],
                        "component": comp,
                        "bulb": i,
                        "center_a": c.tolist(),
                        "n_interior": int(inside.shape[0]),
                        "n_boundary": int(boundary.shape[0]),
                        "area_a2": area,
                        "area_error": area - true_area,
                        "area_error_frac": (area - true_area) / true_area,
                        "centroid_shift_a": centroid.tolist(),
                        "centroid_shift_mm": (centroid * float(sc.a) * 1e3).tolist(),
                        "radial_err_mean_a": float(np.mean(boundary[:, 2])),
                        "radial_err_rms_a": float(np.sqrt(np.mean(boundary[:, 2] ** 2))),
                        "asym_moment_re": float(np.real(moment)),
                        "asym_moment_im": float(np.imag(moment)),
                        "asym_moment_abs": float(abs(moment)),
                    }
                )
                print(
                    f"STAIR ox={ox:g} oy={oy:g} {comp} b{i} "
                    f"area_frac={rows[-1]['area_error_frac']:+.4f} "
                    f"cen_mm={rows[-1]['centroid_shift_mm']} "
                    f"nB={rows[-1]['n_boundary']} asym={rows[-1]['asym_moment_abs']:.4e}",
                    flush=True,
                )
    path = OUT / f"staircase_cluster7_ppc{ppc:g}.json"
    path.write_text(json.dumps({"true_area_a2": true_area, "rows": rows}, indent=2) + "\n")
    print("STAIR_DONE", path, flush=True)


def analyze_traces(path: str) -> None:
    """FFT the saved probe traces and rank late-time peaks."""
    data = np.load(path)
    times = np.asarray(data["times"], float)
    samples = np.asarray(data["samples"], float)
    names = [str(n) for n in data["names"]]
    if len(times) < 8:
        raise RuntimeError(f"no samples in {path}")
    # Late window: last 80 time units, after the source has died in previous runs.
    t1 = float(times[-1])
    late = times >= max(times[0], t1 - 80.0)
    tl = times[late]
    dt = float(np.median(np.diff(tl)))
    print(f"TRACE {path} n={len(times)} t={times[0]:.1f}..{t1:.1f} late_dt={dt:.3f}", flush=True)
    rows = []
    for i, name in enumerate(names):
        y = samples[late, i]
        y = y - np.mean(y)
        spec = np.fft.rfft(y)
        freq = np.fft.rfftfreq(len(y), d=dt)
        amp = np.abs(spec) / max(len(y), 1)
        # Skip the DC bin.
        k = int(np.argmax(amp[1:]) + 1) if len(amp) > 1 else 0
        # Half-window decay of the RMS.
        mid = len(y) // 2
        rms1 = float(np.sqrt(np.mean(y[:mid] ** 2))) if mid else 0.0
        rms2 = float(np.sqrt(np.mean(y[mid:] ** 2))) if mid else 0.0
        rows.append(
            {
                "probe": name,
                "f_a": float(freq[k]),
                "f_GHz": freq_GHz(float(freq[k])),
                "fft_amp": float(amp[k]),
                "rms_early": rms1,
                "rms_late": rms2,
                "rms_full": float(np.sqrt(np.mean(y**2))),
            }
        )
        print(
            f"  {name:14s} peak {freq_GHz(float(freq[k])):7.3f} GHz  "
            f"fft {amp[k]:.4e}  rms {rms1:.4e}->{rms2:.4e}",
            flush=True,
        )
    out = Path(path).with_suffix(".fft.json")
    out.write_text(json.dumps(rows, indent=2) + "\n")


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "resolution"
    if cmd == "resolution":
        write_resolution_table()
        return 0
    if cmd == "replay-fem":
        replay_fem()
        return 0
    if cmd == "mie-scan":
        mie_scan()
        return 0
    if cmd == "meep":
        meep_main()
        return 0
    if cmd == "analyze":
        analyze_traces(sys.argv[2])
        return 0
    if cmd == "staircase":
        staircase()
        return 0
    if cmd == "fem-cluster":
        fem_cluster(int(sys.argv[2]))
        return 0
    if cmd == "core-geometry":
        rec = core_geometry_record()
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "square_core_geometry.json").write_text(json.dumps(rec, indent=2) + "\n")
        print(json.dumps(rec, indent=2), flush=True)
        return 0
    raise SystemExit(f"unknown command {cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
