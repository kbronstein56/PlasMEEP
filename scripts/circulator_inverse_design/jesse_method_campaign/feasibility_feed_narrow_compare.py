#!/usr/bin/env python3
"""Forward-only original vs narrowed-feed comparison (fixed q).

No optimization. Stage-6 NOT STARTED.

Compares production clear_width feed vs TE1-suppressing narrowed feed
(f_c,TE1 = 6.5 GHz) for:
  - Band-B C6 Stage-5 optimum, B=0.05 T
  - Band-C B=0.10 Band-B-start optimum
"""
from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
FEED = FEAS / "feed_narrow_compare"

sys_paths = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sys

sys.path[:0] = sys_paths

from feed_narrow_geometry import (  # noqa: E402
    W_NARROW_LU,
    W_ORIG_LU,
    feed_geometry_document,
    narrowed_feed_wall_mask,
    write_geometry_doc,
)

DESIGNS = {
    "B_stage5": {
        "opt": FEAS / "B_4p75_5p25_C6_B0.0500_stage5_opt.json",
        "b_tesla": 0.05,
        "f0": 4.75,
        "f1": 5.25,
        "watch": [4.909375, 5.164062],
        "label": "Band-B C6 Stage-5 B=0.05",
    },
    "C_B010_Bs5": {
        "opt": FEAS / "C_4p5_5p5_C6_B0.1000_Bs5_opt.json",
        "b_tesla": 0.10,
        "f0": 4.5,
        "f1": 5.5,
        "watch": [4.609375, 4.70625, 4.996875],
        "label": "Band-C B=0.10 from Band-B",
    },
}


def _init_worker():
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    os.environ.setdefault("PLASMEEP_A2_FAST_KERNEL", "1")


def freq_grid(f0, f1, df):
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]


def build_system_feed(grade, s_vec, f_hz, b_tesla, clear_width_lu, measure="projected"):
    """Production build_system with optional narrowed straight-feed aperture.

    When TE1 is cut off (narrowed feed), use n_prop=1 and keep TE1 only as an
    evanescent Robin load (n_evan=1). Extracting TE1 as a 'propagating' S
    channel with imaginary β produces unphysical |S|^2 ≫ 1 on this mesh.
    """
    import sixport_common as sc
    from boundary_ports import build_port_from_boundary
    from diff_fem_core import per_rod_plasma_masks, rho_from_scales
    from direct_dtn_s import assemble_port_operator
    from dtn_port_geometry import port_outer_center
    from fem_validated_solver import assemble_anisotropic
    from full91_dtn_solve import exterior_beyond_ports_mask
    from full91_geometry import bulb_centers_device
    from full91_solve import load_mesh, material_masks
    from modes import mode_constants

    points, tris = load_mesh(grade)
    masks = material_masks(points, tris, sc)
    plasma, quartz, walls, _ = masks
    centers = bulb_centers_device(sc)
    r_p = 4.6 * sc.r_bulb_inner / 6.5
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    rho, f_ord, gamma, fc = rho_from_scales(
        len(tris), rod_masks, s_vec, quartz, walls, f_hz, b_tesla, sc.a, sc.fp_Hz, sc.gamma_Hz
    )
    W = float(clear_width_lu)
    if abs(W - float(sc.clear_width)) > 1e-12:
        wmask = narrowed_feed_wall_mask(points, tris, W_filter=W)
        for arr in rho:
            arr[wmask] = 0.0
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    rough = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        rough.append({"center": c, "outward": n_hat})
    ext = exterior_beyond_ports_mask(points, tris, rough)
    for arr in rho:
        arr[ext] = 0.0
    k0 = 2.0 * np.pi * (f_hz * sc.a / 299792458.0)
    A0 = assemble_anisotropic(points, tris, *rho, k0, np.zeros(len(points), dtype=bool))
    ports = []
    for p in range(6):
        n_hat = np.asarray(sc.port_dirs[p], float)
        n_hat /= np.linalg.norm(n_hat)
        tang = np.array([-n_hat[1], n_hat[0]])
        c = port_outer_center(sc.full_horns[p], n_hat) + shift
        ps = build_port_from_boundary(
            points, tris, rho, c, n_hat, tang, W=W, plane_tol=0.12, port=p
        )
        ps["measure"] = measure
        ports.append(ps)
    modes = mode_constants(k0, W, n_modes=4)
    te1_prop = bool(modes[1]["propagating"])
    if te1_prop:
        n_prop, n_evan = 2, 2
    else:
        # TE0 propagating + TE1 as single evanescent Robin terminator
        n_prop, n_evan = 1, 1
    A = assemble_port_operator(A0, points, ports, modes, n_prop=n_prop, n_evan=n_evan)
    return {
        "points": points,
        "tris": tris,
        "rho": rho,
        "A": A,
        "k0": k0,
        "ports": ports,
        "modes": modes,
        "n_prop": n_prop,
        "n_evan": n_evan,
        "clear_width_lu": W,
        "TE1_propagating": te1_prop,
        "TE1_f_cut_GHz": None
        if modes[1].get("f_cut_Hz") is None
        else float(modes[1]["f_cut_Hz"]) / 1e9,
    }


def s_matrix_te0_drive(st):
    """Return 12×12 S for TE0 drives (TE1 rows filled by projection if needed)."""
    from boundary_ports import modal_c_vector
    from direct_dtn_s import build_Q_columns
    from efficient_adjoint import S_from_U, _pack_C, _solve_multi
    from scipy.sparse.linalg import splu

    points, ports, modes = st["points"], st["ports"], st["modes"]
    n_prop = int(st["n_prop"])
    cols = build_Q_columns(points, ports, modes, n_prop=n_prop)
    C, betas = _pack_C(points, cols)
    lu = splu(st["A"].tocsc())
    if n_prop == 2:
        te0_cols = [2 * j for j in range(6)]
        U = _solve_multi(lu, C[:, te0_cols], trans="N")
        S = S_from_U(C, U, betas, te0_cols)
        return S, U, te0_cols
    # n_prop == 1: 6×6 TE0 S, then pad TE1 rows via DtN projection
    te0_cols = list(range(6))
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    S6 = S_from_U(C, U, betas, te0_cols)
    S = np.zeros((12, 12), dtype=np.complex128)
    for j in range(6):
        for i in range(6):
            S[2 * i, 2 * j] = S6[i, j]
    # TE1←TE0: S_{1p,0j} = 2i √(β1 β0) c1pᵀ Uj  (same frozen DtN convention)
    beta0 = complex(modes[0]["beta"])
    beta1 = complex(modes[1]["beta"])
    scale = 2j * np.sqrt(beta1 * beta0)
    for p in range(6):
        idx, c = modal_c_vector(points, ports[p], 1)
        if len(idx) == 0:
            continue
        c = np.asarray(c, np.complex128)
        for j in range(6):
            S[2 * p + 1, 2 * j] = scale * complex(np.dot(c, U[idx, j]))
    return S, U, [2 * j for j in range(6)]  # logical TE0 columns in 12×12


def diagnose_S(S: np.ndarray) -> dict:
    from mode_objectives import _port_powers
    from multimode_opt_core import te0_metrics_from_S
    from objective_A2 import objective_A2_from_S

    keys = ["des0", "des1", "rev0", "rev1", "thr0", "thr1", "other0", "other1", "R0", "R1"]
    acc = {k: 0.0 for k in keys}
    for j in range(6):
        p = _port_powers(S, j)
        for k in keys:
            acc[k] += float(p[k])
    for k in keys:
        acc[k] /= 6.0
    te1 = acc["des1"] + acc["rev1"] + acc["thr1"] + acc["other1"] + acc["R1"]
    wrong0 = acc["rev0"] + acc["thr0"] + acc["other0"]
    rtot = acc["R0"] + acc["R1"]
    pout = sum(acc.values())
    met = te0_metrics_from_S(S)
    return {
        "J_A2": float(objective_A2_from_S(S)),
        "TE0_des": acc["des0"],
        "TE0_rev": acc["rev0"],
        "TE0_thr": acc["thr0"],
        "TE0_other": acc["other0"],
        "TE0_wrong": wrong0,
        "R0": acc["R0"],
        "R1": acc["R1"],
        "R_tot": rtot,
        "TE1_total": te1,
        "TE1_at_des": acc["des1"],
        "P_out": pout,
        "missing": 1.0 - pout,
        "IL_avg_dB": float(met["IL_avg_dB"]),
        "isolation_avg_dB": float(met["isolation_avg_dB"]),
        "T_des_avg": float(met["T_des_avg"]),
        "T_rev_avg": float(met["T_rev_avg"]),
        "R_avg_te0": float(met["R_avg"]),
        "TE1_avg": float(met["TE1_avg"]),
    }


def _worker(payload: dict) -> dict:
    _init_worker()
    from multimode_opt_core import expand_tied, orbits_and_masks

    q = np.asarray(payload["q"], float)
    orbits, _, _ = orbits_and_masks()
    s = expand_tied(q, orbits)
    st = build_system_feed(
        payload["grade"],
        s,
        float(payload["f_GHz"]) * 1e9,
        float(payload["b_tesla"]),
        float(payload["clear_width_lu"]),
    )
    S, _U, _te0_cols = s_matrix_te0_drive(st)
    row = diagnose_S(S)
    row.update(
        {
            "design": payload["design"],
            "feed": payload["feed"],
            "f_GHz": float(payload["f_GHz"]),
            "b_tesla": float(payload["b_tesla"]),
            "clear_width_lu": float(payload["clear_width_lu"]),
            "TE1_propagating": bool(st["TE1_propagating"]),
            "TE1_f_cut_GHz": st["TE1_f_cut_GHz"],
            "n_prop": int(st["n_prop"]),
            "S_note": (
                "12x12 production S"
                if st["TE1_propagating"]
                else "TE0 6x6 S + TE1 rows by DtN projection (TE1 cut off)"
            ),
        }
    )
    return row


def build_jobs(df_MHz: float) -> list[dict]:
    jobs = []
    for dkey, meta in DESIGNS.items():
        opt = json.loads(meta["opt"].read_text())
        q = np.asarray(opt["q_best"], float).tolist()
        freqs = set(freq_grid(meta["f0"], meta["f1"], df_MHz * 1e-3))
        # densify around watch frequencies
        for fw in meta["watch"]:
            for d in np.arange(-25e-3, 25e-3 + 1e-9, 3.125e-3):
                freqs.add(round(fw + float(d), 6))
        for feed, W in [("original", W_ORIG_LU), ("narrowed", W_NARROW_LU)]:
            for f in sorted(freqs):
                if f < meta["f0"] - 1e-9 or f > meta["f1"] + 1e-9:
                    continue
                jobs.append(
                    {
                        "design": dkey,
                        "feed": feed,
                        "q": q,
                        "f_GHz": f,
                        "b_tesla": meta["b_tesla"],
                        "clear_width_lu": W,
                        "grade": "C",
                    }
                )
    return jobs


def run_jobs(jobs, n_workers=12):
    FEED.mkdir(parents=True, exist_ok=True)
    results = []
    t0 = time.time()
    print(f"FEED NARROW COMPARE: {len(jobs)} evals, n_workers={n_workers}", flush=True)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(_worker, j): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            done += 1
            j = futs[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = {
                    "design": j["design"],
                    "feed": j["feed"],
                    "f_GHz": j["f_GHz"],
                    "error": repr(e),
                }
                print(f"FAIL {j['design']} {j['feed']} f={j['f_GHz']}: {e}", flush=True)
            results.append(r)
            if done % 10 == 0 or done == len(jobs):
                print(
                    f"  [{done}/{len(jobs)}] {r.get('design')} {r.get('feed')} "
                    f"f={r.get('f_GHz')} A2={r.get('J_A2', 'err')}",
                    flush=True,
                )
                (FEED / "feed_narrow_partial.json").write_text(
                    json.dumps({"rows": results, "elapsed_s": time.time() - t0}, indent=2) + "\n"
                )
    results = sorted(
        results,
        key=lambda r: (r.get("design") or "", r.get("feed") or "", r.get("f_GHz") or 0),
    )
    return results, time.time() - t0


def classify_resonance(orig: dict, narr: dict) -> dict:
    """Classify outcome of feed narrowing at a watch frequency."""
    if orig is None or narr is None or "error" in (orig or {}) or "error" in (narr or {}):
        return {"class": "missing", "reason": "missing row"}
    d_te1 = narr["TE1_total"] - orig["TE1_total"]
    d_des = narr["TE0_des"] - orig["TE0_des"]
    d_wrong = narr["TE0_wrong"] - orig["TE0_wrong"]
    d_r = narr["R_tot"] - orig["R_tot"]
    d_a2 = narr["J_A2"] - orig["J_A2"]
    # frequency shift estimated from local spectra by caller; here pointwise
    weakened = d_a2 > 2.0 and narr["J_A2"] > -3.0
    little = (
        abs(d_te1) < 0.05
        and abs(d_des) < 0.05
        and abs(d_wrong) < 0.05
        and abs(d_r) < 0.05
        and abs(d_a2) < 1.0
    )
    if little:
        label = "6_little_change"
        reason = "pointwise deltas small"
    elif weakened and d_te1 < -0.05:
        label = "5_weakened_with_TE1_drop"
        reason = f"A2 {orig['J_A2']:.2f}→{narr['J_A2']:.2f}, TE1 {orig['TE1_total']:.3f}→{narr['TE1_total']:.3f}"
    elif d_te1 < -0.08 and d_des > 0.05 and d_des >= max(d_wrong, d_r, 0):
        label = "1_TE1_to_desired_TE0"
        reason = f"ΔTE1={d_te1:.3f}, Δdes0={d_des:.3f}"
    elif d_te1 < -0.08 and d_wrong > 0.05 and d_wrong >= max(d_des, d_r, 0):
        label = "2_TE1_to_wrong_TE0"
        reason = f"ΔTE1={d_te1:.3f}, Δwrong0={d_wrong:.3f}"
    elif d_te1 < -0.08 and d_r > 0.05 and d_r >= max(d_des, d_wrong, 0):
        label = "3_TE1_to_reflection"
        reason = f"ΔTE1={d_te1:.3f}, ΔR={d_r:.3f}"
    elif weakened:
        label = "5_resonance_weakened"
        reason = f"A2 improved {orig['J_A2']:.2f}→{narr['J_A2']:.2f}"
    else:
        label = "mixed_or_persistent"
        reason = (
            f"ΔA2={d_a2:.2f}, ΔTE1={d_te1:.3f}, Δdes={d_des:.3f}, "
            f"Δwrong={d_wrong:.3f}, ΔR={d_r:.3f}"
        )
    return {
        "class": label,
        "reason": reason,
        "delta": {
            "A2": d_a2,
            "TE1": d_te1,
            "des0": d_des,
            "wrong0": d_wrong,
            "R_tot": d_r,
        },
        "original": {k: orig[k] for k in ("J_A2", "TE0_des", "TE0_wrong", "R_tot", "TE1_total", "P_out")},
        "narrowed": {k: narr[k] for k in ("J_A2", "TE0_des", "TE0_wrong", "R_tot", "TE1_total", "P_out")},
    }


def nearest_row(rows, f, feed):
    cands = [r for r in rows if r.get("feed") == feed and "error" not in r]
    if not cands:
        return None
    return min(cands, key=lambda r: abs(r["f_GHz"] - f))


def find_local_min(rows, feed, f_center, window=0.04):
    cands = [
        r
        for r in rows
        if r.get("feed") == feed
        and "error" not in r
        and abs(r["f_GHz"] - f_center) <= window
    ]
    if not cands:
        return None
    return min(cands, key=lambda r: r["J_A2"])


def plot_spectra(design, rows):
    feeds = ["original", "narrowed"]
    colors = {"original": "#1f4e79", "narrowed": "#b85c38"}
    keys = [
        ("TE0_des", "Desired TE0"),
        ("TE0_wrong", "Wrong TE0"),
        ("R_tot", "Reflection R0+R1"),
        ("TE1_total", "TE1 total"),
        ("J_A2", "A2"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(10, 9), constrained_layout=True)
    axes = axes.ravel()
    for ax, (key, title) in zip(axes, keys):
        for feed in feeds:
            seq = sorted(
                [r for r in rows if r.get("feed") == feed and "error" not in r],
                key=lambda r: r["f_GHz"],
            )
            if not seq:
                continue
            ax.plot(
                [r["f_GHz"] for r in seq],
                [r[key] for r in seq],
                "-",
                lw=1.2,
                color=colors[feed],
                label=feed,
            )
        ax.set_title(title)
        ax.set_xlabel("f [GHz]")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    axes[-1].axis("off")
    # watch lines on A2 panel
    ax = axes[4]
    for fw in DESIGNS[design]["watch"]:
        ax.axvline(fw, color="k", ls=":", lw=0.8, alpha=0.5)
    fig.suptitle(f"{DESIGNS[design]['label']} — original vs narrowed feed (fixed q)")
    fig.savefig(FEED / f"{design}_spectra_orig_vs_narrow.png", dpi=150)
    plt.close(fig)


def export_field_pair(design="C_B010_Bs5", f_ghz=4.609375):
    """|Hz| maps at Band-C notch for original and narrowed feeds."""
    from multimode_opt_core import expand_tied, orbits_and_masks, te0_metrics_from_S
    from matplotlib.tri import LinearTriInterpolator, Triangulation
    from objective_A2 import objective_A2_from_S

    _init_worker()
    meta = DESIGNS[design]
    q = np.asarray(json.loads(meta["opt"].read_text())["q_best"], float)
    orbits, _, _ = orbits_and_masks()
    s = expand_tied(q, orbits)
    grids = []
    for feed, W in [("original", W_ORIG_LU), ("narrowed", W_NARROW_LU)]:
        st = build_system_feed("C", s, f_ghz * 1e9, meta["b_tesla"], W)
        S, U, _te0_cols = s_matrix_te0_drive(st)
        met = te0_metrics_from_S(S)
        j_src = int(max(met["per_source"], key=lambda r: r["IL_dB"])["source"])
        # U columns are always the 6 TE0 drives
        Hz = np.asarray(U[:, j_src], complex)
        points, tris = st["points"], st["tris"]
        tri = Triangulation(points[:, 0], points[:, 1], tris)
        xmin, xmax = float(points[:, 0].min()), float(points[:, 0].max())
        ymin, ymax = float(points[:, 1].min()), float(points[:, 1].max())
        x = np.linspace(xmin, xmax, 420)
        y = np.linspace(ymin, ymax, 400)
        X, Y = np.meshgrid(x, y)
        Mag = np.ma.filled(LinearTriInterpolator(tri, np.abs(Hz))(X, Y), np.nan).astype(np.float32)
        grids.append(
            {
                "feed": feed,
                "Mag": Mag,
                "x": x,
                "y": y,
                "source": j_src + 1,
                "J_A2": float(objective_A2_from_S(S)),
                "met": {
                    "T_des": float(met["T_des_avg"]),
                    "T_rev": float(met["T_rev_avg"]),
                    "R": float(met["R_avg"]),
                    "TE1": float(met["TE1_avg"]),
                },
                "TE1_propagating": st["TE1_propagating"],
            }
        )
        print(f"FIELD {feed} A2={grids[-1]['J_A2']:.3f} TE1={grids[-1]['met']['TE1']:.3f}", flush=True)

    vmax = float(
        np.nanpercentile(
            np.concatenate([g["Mag"][~np.isnan(g["Mag"])].ravel() for g in grids]), 99.0
        )
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), constrained_layout=True)
    for ax, g in zip(axes, grids):
        im = ax.imshow(
            g["Mag"],
            origin="lower",
            extent=[g["x"][0], g["x"][-1], g["y"][0], g["y"][-1]],
            cmap="magma",
            vmin=0,
            vmax=vmax,
            aspect="equal",
        )
        ax.set_title(
            f"{g['feed']} feed @ {f_ghz:.6f} GHz P{g['source']}\n"
            f"A2={g['J_A2']:.2f} Tdes={g['met']['T_des']:.3f} Trev={g['met']['T_rev']:.3f} "
            f"R={g['met']['R']:.3f} TE1={g['met']['TE1']:.3f}\n"
            f"TE1_prop={g['TE1_propagating']}",
            fontsize=9,
        )
        ax.set_xlabel("x [cm]")
        ax.set_ylabel("y [cm]")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="|Hz|")
    fig.suptitle(f"Band-C notch field — original vs narrowed feed (identical scale, vmax={vmax:.3g})")
    png = FEED / "C_B010_field_4p609_orig_vs_narrow.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)
    (FEED / "C_B010_field_4p609_meta.json").write_text(
        json.dumps({"f_GHz": f_ghz, "vmax": vmax, "frames": grids, "png": str(png)}, indent=2, default=lambda o: None)
        + "\n"
    )
    # strip Mag from json meta — rewrite compact
    meta_out = {
        "f_GHz": f_ghz,
        "vmax": vmax,
        "png": str(png),
        "frames": [
            {
                "feed": g["feed"],
                "source": g["source"],
                "J_A2": g["J_A2"],
                "met": g["met"],
                "TE1_propagating": g["TE1_propagating"],
            }
            for g in grids
        ],
    }
    (FEED / "C_B010_field_4p609_meta.json").write_text(json.dumps(meta_out, indent=2) + "\n")
    return meta_out


def write_report(results, geom_doc, elapsed_s):
    by = {}
    for r in results:
        if "error" in r:
            continue
        by.setdefault(r["design"], {}).setdefault(r["feed"], []).append(r)

    classifications = {}
    lines = [
        "# Forward-only narrowed-feed comparison",
        "",
        "**Stage-6: NOT STARTED / AWAITING AUTHORIZATION**",
        "**No optimization.** Fixed plasma-density vectors q.",
        "",
        "## Feed geometry and analytic cutoffs",
        "",
        f"- Original clear width: **{geom_doc['original_feed']['W_mm']:.3f} mm** "
        f"({geom_doc['original_feed']['W_lu']:.4f} lu) → TE1 cutoff "
        f"**{geom_doc['original_feed']['TE1_cutoff_GHz']:.4f} GHz**",
        f"- Narrowed clear width: **{geom_doc['narrowed_feed']['W_mm']:.3f} mm** "
        f"({geom_doc['narrowed_feed']['W_lu']:.4f} lu) → TE1 cutoff "
        f"**{geom_doc['narrowed_feed']['TE1_cutoff_GHz']:.4f} GHz**",
        f"- Formula: `{geom_doc['original_feed']['formula']}`",
        f"- Target: {geom_doc['target']['rationale']}",
        "",
        f"Runtime: {elapsed_s/60:.1f} min for {len(results)} forward evals.",
        "",
        "## Watch-frequency table",
        "",
        "| design | f_watch | feed | A2 | des0 | wrong0 | Rtot | TE1 | P_out | TE1_prop |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|---|",
    ]

    for dkey, meta in DESIGNS.items():
        rows_all = [r for r in results if r.get("design") == dkey]
        classifications[dkey] = {}
        for fw in meta["watch"]:
            o = nearest_row(rows_all, fw, "original")
            n = nearest_row(rows_all, fw, "narrowed")
            # also locate local A2 minima near watch for shift estimate
            o_min = find_local_min(rows_all, "original", fw)
            n_min = find_local_min(rows_all, "narrowed", fw)
            for feed, row in [("original", o), ("narrowed", n)]:
                if row is None:
                    continue
                lines.append(
                    f"| {dkey} | {fw:.6f} | {feed} | {row['J_A2']:.3f} | {row['TE0_des']:.3f} | "
                    f"{row['TE0_wrong']:.3f} | {row['R_tot']:.3f} | {row['TE1_total']:.3f} | "
                    f"{row['P_out']:.3f} | {row.get('TE1_propagating')} |"
                )
            cls = classify_resonance(o, n)
            # frequency shift
            if o_min and n_min:
                df = n_min["f_GHz"] - o_min["f_GHz"]
                cls["local_min_original_f"] = o_min["f_GHz"]
                cls["local_min_narrowed_f"] = n_min["f_GHz"]
                cls["local_min_shift_MHz"] = df * 1e3
                cls["local_min_A2_orig"] = o_min["J_A2"]
                cls["local_min_A2_narr"] = n_min["J_A2"]
                if abs(df) > 0.01 and abs(n_min["J_A2"] - o_min["J_A2"]) < 2.0:
                    # substantial shift with similar depth
                    if "little" not in cls["class"] and "weakened" not in cls["class"]:
                        cls["class"] = "4_resonance_shifted"
                        cls["reason"] += f"; local-min shift {df*1e3:.1f} MHz"
            classifications[dkey][f"{fw:.6f}"] = cls

    lines += ["", "## Resonance classifications under feed narrowing", ""]
    for dkey, cmap in classifications.items():
        lines.append(f"### {dkey}")
        for fw, cls in cmap.items():
            lines.append(
                f"- **{fw} GHz → {cls['class']}**: {cls['reason']} "
                f"(ΔA2={cls['delta']['A2']:.2f}, ΔTE1={cls['delta']['TE1']:.3f}, "
                f"Δdes={cls['delta']['des0']:.3f}, Δwrong={cls['delta']['wrong0']:.3f}, "
                f"ΔR={cls['delta']['R_tot']:.3f}"
                + (
                    f", shift={cls.get('local_min_shift_MHz'):.1f} MHz"
                    if cls.get("local_min_shift_MHz") is not None
                    else ""
                )
                + ")"
            )
        lines.append("")

    # Scientific interpretation scaffold — filled with actual classes
    all_classes = [c["class"] for m in classifications.values() for c in m.values()]
    te1_bottleneck = any("TE1_to_desired" in c or "weakened_with_TE1" in c for c in all_classes)
    te1_outlet = any("TE1_to_wrong" in c or "TE1_to_reflection" in c for c in all_classes)
    little = sum(1 for c in all_classes if "little" in c)
    weakened = sum(1 for c in all_classes if "weakened" in c)

    lines += [
        "## Scientific interpretation",
        "",
        "### Distinguishing hypotheses",
        "",
    ]
    if weakened and not te1_outlet:
        lines.append(
            "- Evidence leans toward **TE1 as a true bandwidth bottleneck** and/or "
            "**feed geometry as a contributor**: narrowing substantially weakens notches "
            "without simply dumping TE1 into wrong-port/reflection."
        )
    elif te1_outlet and not weakened:
        lines.append(
            "- Evidence leans toward **TE1 as an outlet for a deeper cavity resonance**: "
            "TE1 drops but power reappears in wrong TE0 or reflection; A2 stays poor."
        )
    elif little >= len(all_classes) // 2:
        lines.append(
            "- Evidence leans toward **feed-TE1 not being the controlling knob** for these "
            "particular resonances (little change at watch frequencies)."
        )
    else:
        lines.append(
            "- Mixed evidence across watch frequencies — see per-resonance classes above. "
            "TE1 suppression alone does not uniformly heal all notch families."
        )

    lines += [
        "",
        "### Is a narrowed-feed optimization justified?",
        "",
    ]
    # Decision rule
    if weakened >= 1 and te1_outlet == 0:
        verdict = (
            "**Conditionally yes** — as a *small* next design study: re-optimize q under the "
            "narrowed feed on a manageable band (e.g. Band B), still not Stage-6 full-band."
        )
    elif te1_outlet and weakened:
        verdict = (
            "**Partially** — narrowing helps some channels but can redistribute failures; "
            "a short Band-B reopt under narrowed feed is justified only if TE1 drop coincides "
            "with A2 improvement on the soft floor. Do not jump to Stage-6."
        )
    elif te1_outlet and not weakened:
        verdict = (
            "**Not yet as the primary path** — TE1 appears to be an outlet for a deeper "
            "resonance; changing feed width alone redistributes modal power. Prefer geometry "
            "changes that target wrong-port/cavity content, or a combined feed+array study."
        )
    else:
        verdict = (
            "**Not justified as the next major effort** based on this forward test — "
            "narrowing did not materially heal the known notches. Seek a different physical lever."
        )
    lines.append(verdict)
    lines += [
        "",
        "No optimization was launched after this comparison.",
        "",
        "Stage-6 remains: **NOT STARTED / AWAITING AUTHORIZATION.**",
        "",
        "## Artifacts",
        "",
        f"- `{FEED/'feed_geometry_cutoffs.json'}`",
        f"- `{FEED/'feed_narrow_results.json'}`",
        f"- `{FEED/'B_stage5_spectra_orig_vs_narrow.png'}`",
        f"- `{FEED/'C_B010_Bs5_spectra_orig_vs_narrow.png'}`",
        f"- `{FEED/'C_B010_field_4p609_orig_vs_narrow.png'}`",
        "",
    ]
    path = FEED / "FEED_NARROW_COMPARE_REPORT.md"
    path.write_text("\n".join(lines) + "\n")
    (FEED / "feed_narrow_classifications.json").write_text(json.dumps(classifications, indent=2) + "\n")
    return path, classifications


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--df-MHz", type=float, default=12.5, help="base grid; watch freqs densified")
    ap.add_argument("--skip-fields", action="store_true")
    ap.add_argument("--fields-only", action="store_true")
    args = ap.parse_args()
    FEED.mkdir(parents=True, exist_ok=True)

    geom = write_geometry_doc(FEED / "feed_geometry_cutoffs.json")
    print("GEOMETRY DOC:", json.dumps({k: geom[k] for k in ("original_feed", "narrowed_feed", "target")}, indent=2), flush=True)

    results = []
    elapsed = 0.0
    if not args.fields_only:
        jobs = build_jobs(args.df_MHz)
        (FEED / "feed_narrow_jobs.json").write_text(json.dumps({"n": len(jobs), "df_MHz": args.df_MHz}, indent=2) + "\n")
        results, elapsed = run_jobs(jobs, n_workers=args.n_workers)
        payload = {
            "stage6_status": "NOT STARTED / AWAITING AUTHORIZATION",
            "no_optimization": True,
            "df_MHz": args.df_MHz,
            "n_jobs": len(jobs),
            "elapsed_s": elapsed,
            "geometry": geom,
            "rows": results,
        }
        (FEED / "feed_narrow_results.json").write_text(json.dumps(payload, indent=2) + "\n")
        for dkey in DESIGNS:
            plot_spectra(dkey, [r for r in results if r.get("design") == dkey])
        report, classifications = write_report(results, geom, elapsed)
        print(f"wrote {report}", flush=True)
        print("\n=== WATCH FREQUENCY SUMMARY ===", flush=True)
        for dkey, cmap in classifications.items():
            for fw, cls in cmap.items():
                print(
                    f"{dkey} @{fw}: {cls['class']} | "
                    f"A2 {cls['original']['J_A2']:.3f}→{cls['narrowed']['J_A2']:.3f} | "
                    f"TE1 {cls['original']['TE1_total']:.3f}→{cls['narrowed']['TE1_total']:.3f} | "
                    f"wrong {cls['original']['TE0_wrong']:.3f}→{cls['narrowed']['TE0_wrong']:.3f} | "
                    f"R {cls['original']['R_tot']:.3f}→{cls['narrowed']['R_tot']:.3f}",
                    flush=True,
                )

    if not args.skip_fields:
        export_field_pair("C_B010_Bs5", 4.609375)

    print("STAGE-6: NOT STARTED / AWAITING AUTHORIZATION", flush=True)
    print("NO OPTIMIZATION LAUNCHED", flush=True)


if __name__ == "__main__":
    main()
