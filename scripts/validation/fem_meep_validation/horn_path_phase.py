#!/usr/bin/env python3
"""Complex Hz along shared cuts: converged FEM-VH vs the 25 ppc Meep DFT grid.

One complex scale is fixed on the P1 feed, near the source, and applied to every cut.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from b0_campaign import build_symmetric_mesh, horn_prism_quads, horn_rho  # noqa: E402
from fem_validated_solver import (  # noqa: E402
    ElementSampler,
    assemble_anisotropic,
    inject_mode_consistent,
    load_numerical_mode,
    solve_system,
)
from sixport_common import set_geometry_context  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_continuum"
GRID = ROOT / "outputs" / "validation" / "fem_meep_validation" / "boundary_pml" / "meep_hz_grid.npz"


def sample_meep(x, y, Hz, pts):
    ix = np.interp(pts[:, 0], x, np.arange(len(x)))
    iy = np.interp(pts[:, 1], y, np.arange(len(y)))
    i0 = np.clip(np.floor(ix).astype(int), 0, len(x) - 2)
    j0 = np.clip(np.floor(iy).astype(int), 0, len(y) - 2)
    tx = np.clip(ix - i0, 0.0, 1.0)
    ty = np.clip(iy - j0, 0.0, 1.0)
    f00 = Hz[i0, j0]
    f10 = Hz[i0 + 1, j0]
    f01 = Hz[i0, j0 + 1]
    f11 = Hz[i0 + 1, j0 + 1]
    return (1 - tx) * (1 - ty) * f00 + tx * (1 - ty) * f10 + (1 - tx) * ty * f01 + tx * ty * f11


def inside_metal(pts, quads):
    mask = np.zeros(len(pts), dtype=bool)
    for q in quads:
        mask |= MPath(q).contains_points(pts, radius=1e-9)
    return mask


def line(p0, p1, n):
    t = np.linspace(0.0, 1.0, n)[:, None]
    xy = np.asarray(p0, float) * (1 - t) + np.asarray(p1, float) * t
    s = np.linspace(0.0, float(np.linalg.norm(p1 - p0)), n)
    return xy, s


def metrics(s, fem, meep, air):
    use = air & np.isfinite(fem.real) & np.isfinite(meep.real)
    if use.sum() < 8:
        return {"n": int(use.sum())}
    f = fem[use]
    m = meep[use]
    ss = s[use]
    rel = float(np.linalg.norm(f - m) / (np.linalg.norm(m) + 1e-30))
    mag = float(np.linalg.norm(np.abs(f) - np.abs(m)) / (np.linalg.norm(np.abs(m)) + 1e-30))
    amp = np.maximum(np.abs(m), np.abs(f))
    strong = amp > 0.08 * np.max(amp)
    dphi = np.angle(f * np.conj(m))
    mean_phase = float(np.mean(dphi[strong])) if np.any(strong) else float(np.mean(dphi))
    # Local wavenumber from unwrapped phase, only where the field is not in a null.
    def slope(z):
        ph = np.unwrap(np.angle(z))
        ph = np.where(strong, ph, np.nan)
        # Fit on the longest finite run.
        ok = np.isfinite(ph)
        if ok.sum() < 8:
            return None
        coef = np.polyfit(ss[ok], ph[ok], 1)
        return float(coef[0])

    k_f = slope(f)
    k_m = slope(m)
    # Where the complex error leaves the feed floor: rolling window from s=0.
    win = max(12, int(0.4 / max(ss[1] - ss[0], 1e-6)))
    roll = np.full(len(ss), np.nan)
    for i in range(len(ss) - win):
        sl = slice(i, i + win)
        roll[i] = np.linalg.norm(f[sl] - m[sl]) / (np.linalg.norm(m[sl]) + 1e-30)
    feed = roll[np.isfinite(roll)]
    floor = float(np.median(feed[: max(3, len(feed) // 10)])) if len(feed) else None
    diverge_s = None
    if floor is not None:
        for i, val in enumerate(roll):
            if np.isfinite(val) and val > max(0.10, 3.0 * floor) and ss[i] > ss[0] + 0.4:
                diverge_s = float(ss[i])
                break
    grad = None
    if k_f is not None and k_m is not None:
        ph_f = np.gradient(np.unwrap(np.angle(f)), ss)
        ph_m = np.gradient(np.unwrap(np.angle(m)), ss)
        grad = {
            "median_dphi_ds_fem": float(np.nanmedian(ph_f[strong])) if np.any(strong) else None,
            "median_dphi_ds_meep": float(np.nanmedian(ph_m[strong])) if np.any(strong) else None,
        }
    return {
        "n_air": int(use.sum()),
        "complex_rel_l2": rel,
        "magnitude_rel_l2": mag,
        "mean_phase_offset_rad": mean_phase,
        "mean_phase_offset_deg": float(np.degrees(mean_phase)),
        "phase_slope_fem": k_f,
        "phase_slope_meep": k_m,
        "phase_slope_difference": None if k_f is None or k_m is None else float(k_f - k_m),
        "feed_window_rel_l2": floor,
        "diverge_s_a": diverge_s,
        "local_gradient": grad,
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    set_geometry_context(res=50, horn_walls="prism", grid_offset_cells=(0, 0), coord_rotation_deg=0.0)
    print("building FEM-VH", flush=True)
    points, tris = build_symmetric_mesh(0.04, 0.012, 0.022, 0.028)
    k0 = 2 * np.pi * float(sc.fs_a)
    pec = np.zeros(len(points), dtype=bool)
    rho, _mask = horn_rho(points, tris)
    sampler = ElementSampler(points, tris)
    A = assemble_anisotropic(points, tris, *rho, k0, pec)
    src = load_numerical_mode()
    b = inject_mode_consistent(points, tris, src, sampler)
    hz, *_ = solve_system(A, b, pec)
    print("solved", len(points), flush=True)

    grid = np.load(GRID)
    mx, my, MH = grid["x"], grid["y"], grid["Hz"]
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    quads = [q - shift for q in horn_prism_quads(sc.nx_ports, sc.ny_ports)]

    h1 = sc.horn_for_port(0, 50)
    src_xy = np.asarray(h1["source_center"], float)[:2]
    throat = np.asarray(h1["throat_center"], float)[:2]
    mon = np.asarray(sc.monitor_center_for_port(0, 50), float)[:2]
    mon4 = np.asarray(sc.monitor_center_for_port(3, 50), float)[:2]
    u2 = np.asarray(sc.effective_port_dir(1), float)
    u2 /= np.linalg.norm(u2)
    u6 = np.asarray(sc.effective_port_dir(5), float)
    u6 /= np.linalg.norm(u6)

    # s = 0 at the P1 source for the axis. Other cuts use their own arc length.
    axis_a = src_xy + np.array([1.2, 0.0])
    axis_b = mon4 - np.array([0.8, 0.0])
    paths = {
        "P1_axis": line(axis_a, axis_b, 900),
        "diag_plus60_P2": line(-7.5 * u2, 8.5 * u2, 700),
        "diag_minus60_P6": line(-7.5 * u6, 8.5 * u6, 700),
        "transverse_x0": line(np.array([0.0, -6.0]), np.array([0.0, 6.0]), 500),
    }
    print("axis ends", axis_a, axis_b, "P2 dir", u2, "throat", throat, "monitor", mon, flush=True)

    # One scale, near the source, on the axis, in the feed air.
    xy_sc, _ = paths["P1_axis"]
    feed = (xy_sc[:, 0] > src_xy[0] - 0.35) & (xy_sc[:, 0] < src_xy[0] + 0.15) & (np.abs(xy_sc[:, 1]) < 0.7)
    feed &= ~inside_metal(xy_sc, quads)
    fem_sc, _, _ = sampler.fields(hz, *rho, k0, xy_sc + shift)
    meep_sc = sample_meep(mx, my, MH, xy_sc)
    den = np.vdot(fem_sc[feed], fem_sc[feed])
    scale = np.vdot(fem_sc[feed], meep_sc[feed]) / den
    print("scale", scale, "feed points", int(feed.sum()), flush=True)

    records = {"scale": {"real": float(scale.real), "imag": float(scale.imag), "abs": float(abs(scale)), "n_feed": int(feed.sum())}}
    fig, axes = plt.subplots(4, 3, figsize=(12.2, 11.0), sharex=False)
    for row, (name, (xy, s)) in enumerate(paths.items()):
        fem, _, _ = sampler.fields(hz, *rho, k0, xy + shift)
        fem = fem * scale
        meep = sample_meep(mx, my, MH, xy)
        air = ~inside_metal(xy, quads)
        rec = metrics(s, fem, meep, air)
        rec["s"] = s.tolist()
        rec["x"] = xy[:, 0].tolist()
        rec["y"] = xy[:, 1].tolist()
        records[name] = {k: v for k, v in rec.items() if k not in ("s", "x", "y")}
        # Store curves separately so the summary json stays small.
        np.savez_compressed(
            OUT / f"path_{name}.npz",
            s=s,
            x=xy[:, 0],
            y=xy[:, 1],
            Hz_fem=fem,
            Hz_meep=meep,
            air=air,
        )
        show = air
        ax = axes[row]
        ax[0].plot(s[show], np.abs(meep[show]), lw=1.0, label="Meep")
        ax[0].plot(s[show], np.abs(fem[show]), lw=1.0, label="FEM")
        ax[0].set_ylabel(name, fontsize=8)
        ax[1].plot(s[show], np.unwrap(np.angle(meep[show])), lw=1.0)
        ax[1].plot(s[show], np.unwrap(np.angle(fem[show])), lw=1.0)
        dphi = np.angle(fem * np.conj(meep))
        ax[2].plot(s[show], np.degrees(dphi[show]), lw=1.0, color="C2")
        if row == 0:
            ax[0].set_title("|Hz|")
            ax[1].set_title("unwrapped phase (rad)")
            ax[2].set_title("phase FEM−Meep (deg)")
            ax[0].legend(fontsize=7)
        ax[2].set_ylim(-40, 40)
        print(name, {k: records[name][k] for k in ("complex_rel_l2", "magnitude_rel_l2", "mean_phase_offset_deg", "phase_slope_difference", "diverge_s_a")}, flush=True)

    for a in axes[-1]:
        a.set_xlabel("arc length s (a)")
    fig.tight_layout()
    fig.savefig(OUT / "path_phase.png", dpi=140)
    (OUT / "path_phase.json").write_text(json.dumps(records, indent=2) + "\n")
    print("wrote", OUT / "path_phase.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
