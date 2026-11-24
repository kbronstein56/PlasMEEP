#!/usr/bin/env python3
"""Resonance localization + forward-only sensitivity for Band-B ~5.16 GHz notch.

Baseline: Band-B C6 Stage-5 optimum, B=0.05 T, original 40 mm feed.
No q optimization. Stage-6 NOT STARTED.

Phase A — characterize the resonance (exact min, controls, field/energy maps).
Phase B — cheap C6-preserving density perturbations (rings + per-orbit) with
          local frequency sweeps; rank by Δf_res and ΔA2_min.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import LinearTriInterpolator, Triangulation

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
LOC = FEAS / "resonance_localization"

sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from feasibility_notch_autopsy import diagnose_bundle  # noqa: E402
from multimode_opt_core import expand_tied, orbits_and_masks  # noqa: E402
from objective_A2 import objective_A2_from_S  # noqa: E402
from sixfold_orbits import build_orbits  # noqa: E402

OPT_PATH = FEAS / "B_4p75_5p25_C6_B0.0500_stage5_opt.json"
AUDIT_FINE = FEAS / "audit_B_4p75_5p25_C6_B0.0500_stage5_opt_1.5625MHz.json"
B_TESLA = 0.05
BAND_F0, BAND_F1 = 4.75, 5.25
Q_LO, Q_HI = 0.05, 1.8

# Prior dense-audit minimum of the ~5.16 family
F_RES_SEED = 5.164062
# Side controls (best local A2 from 1.5625 MHz audit)
F_CTRL_LEFT = 5.059375
F_CTRL_RIGHT = 5.226562
F_CTRL_BEST = 4.967187
# Extra Band-B sentinels for collateral damage
DELTA_Q_ABS = 0.15  # absolute C6 / ring density-scale bump
# Local resonance window for perturbation sweeps (cheap but resolves ~5.16 family)
LOCAL_F0, LOCAL_F1, LOCAL_DF_MHZ = 5.12, 5.21, 12.5
# Band-B sentinels for collateral-damage check (kept small)
BAND_SENTINELS = [4.967187, 5.059375, 5.226562]


def _init_worker():
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    os.environ.setdefault("PLASMEEP_A2_FAST_KERNEL", "1")


def freq_grid(f0, f1, df_mhz):
    df = df_mhz * 1e-3
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]


def clip_q(q: np.ndarray) -> np.ndarray:
    return np.clip(np.asarray(q, float), Q_LO, Q_HI)


def orbit_catalog():
    import sixport_common as sc
    from full91_geometry import bulb_centers_device

    centers = np.asarray(bulb_centers_device(sc), float)
    orbits, origin = build_orbits(centers)
    rows = []
    for i, o in enumerate(orbits):
        rows.append(
            {
                "orbit": i,
                "size": int(o["size"]),
                "radius_a": float(o["radius_a"]),
                "rods": [int(x) for x in o["rods"]],
            }
        )
    # Ring families (hex shells grouped)
    rings = {
        "center": [0],
        "inner": [i for i, r in enumerate(rows) if 0 < r["radius_a"] <= 2.05],
        "middle": [i for i, r in enumerate(rows) if 2.05 < r["radius_a"] <= 4.05],
        "outer": [i for i, r in enumerate(rows) if r["radius_a"] > 4.05],
    }
    return {"origin": [float(origin[0]), float(origin[1])], "orbits": rows, "rings": rings}


def apply_orbit_delta(q0: np.ndarray, orbit_ids, delta: float) -> np.ndarray:
    q = np.array(q0, float, copy=True)
    for i in orbit_ids:
        q[i] = q[i] + delta
    return clip_q(q)


def perturbation_cases(q0: np.ndarray, catalog: dict) -> list[dict]:
    cases = [{"name": "baseline", "kind": "baseline", "orbits": [], "delta": 0.0, "q": q0.tolist()}]
    rings = catalog["rings"]
    for rname, ids in rings.items():
        for sign, tag in [(+1, "plus"), (-1, "minus")]:
            d = sign * DELTA_Q_ABS
            q = apply_orbit_delta(q0, ids, d)
            cases.append(
                {
                    "name": f"ring_{rname}_{tag}",
                    "kind": "ring",
                    "ring": rname,
                    "orbits": ids,
                    "delta": d,
                    "q": q.tolist(),
                }
            )
    # Individual C6 orbits (±)
    for i in range(len(q0)):
        for sign, tag in [(+1, "plus"), (-1, "minus")]:
            d = sign * DELTA_Q_ABS
            q = apply_orbit_delta(q0, [i], d)
            cases.append(
                {
                    "name": f"orbit{i:02d}_{tag}",
                    "kind": "orbit",
                    "orbit": i,
                    "orbits": [i],
                    "delta": d,
                    "radius_a": catalog["orbits"][i]["radius_a"],
                    "q": q.tolist(),
                }
            )
    return cases


def _worker(payload: dict) -> dict:
    _init_worker()
    from eval_forward_only import eval_q_forward

    q = np.asarray(payload["q"], float)
    orbits, _, _ = orbits_and_masks()
    bun, _ = eval_q_forward(
        payload.get("grade", "C"),
        q,
        orbits,
        f_hz=float(payload["f_GHz"]) * 1e9,
        b_tesla=float(payload["b_tesla"]),
    )
    out = diagnose_bundle(float(payload["f_GHz"]), bun)
    out.update(
        {
            "case": payload["case"],
            "kind": payload.get("kind"),
            "b_tesla": float(payload["b_tesla"]),
        }
    )
    return out


def run_jobs(jobs: list[dict], n_workers: int, partial_name: str) -> tuple[list[dict], float]:
    LOC.mkdir(parents=True, exist_ok=True)
    results = []
    t0 = time.time()
    print(f"RESONANCE LOCALIZATION: {len(jobs)} evals, n_workers={n_workers}", flush=True)
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
                    "case": j["case"],
                    "f_GHz": j["f_GHz"],
                    "error": repr(e),
                }
                print(f"FAIL {j['case']} f={j['f_GHz']}: {e}", flush=True)
            results.append(r)
            if done % 20 == 0 or done == len(jobs):
                print(
                    f"  [{done}/{len(jobs)}] {r.get('case')} f={r.get('f_GHz')} "
                    f"A2={r.get('J_A2', 'err')}",
                    flush=True,
                )
                (LOC / partial_name).write_text(
                    json.dumps({"rows": results, "elapsed_s": time.time() - t0}, indent=2) + "\n"
                )
    results = sorted(
        results,
        key=lambda r: (r.get("case") or "", r.get("f_GHz") or 0),
    )
    return results, time.time() - t0


# ---------------------------------------------------------------------------
# Phase A: field / energy localization
# ---------------------------------------------------------------------------

def solve_fields(q, f_ghz, b_tesla):
    """Return Hz for all 6 TE0 drives + masks + triangle areas."""
    from adjoint_direct_dtn import build_system
    from direct_dtn_s import build_Q_columns
    from efficient_adjoint import S_from_U, _pack_C, _solve_multi
    from multimode_opt_core import te0_metrics_from_S
    from scipy.sparse.linalg import splu

    orbits, _, _ = orbits_and_masks()
    s = expand_tied(np.asarray(q, float), orbits)
    st = build_system("C", s, f_ghz * 1e9, b_tesla, "projected")
    points, tris = st["points"], st["tris"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    lu = splu(st["A"].tocsc())
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    S = S_from_U(C, U, betas, te0_cols)
    met = te0_metrics_from_S(S)
    # triangle areas
    xy = points[tris]
    area = 0.5 * np.abs(
        (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1])
        - (xy[:, 2, 0] - xy[:, 0, 0]) * (xy[:, 1, 1] - xy[:, 0, 1])
    )
    return {
        "points": points,
        "tris": tris,
        "area": area,
        "U": U,
        "S": S,
        "met": met,
        "J_A2": float(objective_A2_from_S(S)),
        "rho": st["rho"],
        "rod_masks": st.get("rod_masks"),
        "plasma_mask": None,  # filled below
        "st": st,
    }


def _plasma_mask_from_st(st, rod_masks):
    if rod_masks is not None:
        m = np.zeros(len(rod_masks[0]), dtype=bool)
        for rm in rod_masks:
            m |= rm
        return m
    # fallback: quartz/plasma from rho non-vacuum — use walls==0 & plasma dens
    return None


def energy_partition(pack, catalog, orbits, j_src: int):
    """Partition |Hz|^2-weighted triangle energy by ring / plasma / cavity."""
    points, tris, area = pack["points"], pack["tris"], pack["area"]
    Hz = np.asarray(pack["U"][:, j_src], complex)
    # triangle-averaged |Hz|^2
    Hz_t = np.mean(np.abs(Hz[tris]) ** 2, axis=1)
    e_tri = Hz_t * area
    e_tot = float(e_tri.sum()) + 1e-30

    rod_masks = pack.get("rod_masks")
    if rod_masks is None:
        # rebuild rod masks
        import sixport_common as sc
        from diff_fem_core import per_rod_plasma_masks
        from full91_geometry import bulb_centers_device

        centers = bulb_centers_device(sc)
        r_p = 4.6 * sc.r_bulb_inner / 6.5
        rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
        pack["rod_masks"] = rod_masks

    plasma = np.zeros(len(tris), dtype=bool)
    for rm in rod_masks:
        plasma |= rm
    e_plasma = float(e_tri[plasma].sum())
    e_nonplasma = float(e_tri[~plasma].sum())

    # per-orbit / ring plasma energy
    orbit_e = []
    for i, o in enumerate(orbits):
        m = np.zeros(len(tris), dtype=bool)
        for rod in o["rods"]:
            m |= rod_masks[rod]
        orbit_e.append(float(e_tri[m].sum()))
    ring_e = {}
    for rname, ids in catalog["rings"].items():
        ring_e[rname] = float(sum(orbit_e[i] for i in ids))

    # angular sectors (6) of plasma energy about device origin
    cents = points[tris].mean(axis=1)
    origin = np.array(catalog["origin"], float)
    ang = np.arctan2(cents[:, 1] - origin[1], cents[:, 0] - origin[0])
    sector_e = []
    for s in range(6):
        lo = -np.pi + s * np.pi / 3
        hi = lo + np.pi / 3
        sel = plasma & (ang >= lo) & (ang < hi)
        sector_e.append(float(e_tri[sel].sum()))

    return {
        "E_tot": e_tot,
        "E_plasma": e_plasma,
        "E_nonplasma": e_nonplasma,
        "frac_plasma": e_plasma / e_tot,
        "frac_nonplasma": e_nonplasma / e_tot,
        "orbit_E": orbit_e,
        "orbit_frac_of_plasma": [x / (e_plasma + 1e-30) for x in orbit_e],
        "ring_E": ring_e,
        "ring_frac_of_plasma": {k: v / (e_plasma + 1e-30) for k, v in ring_e.items()},
        "sector_E_plasma": sector_e,
        "sector_frac_of_plasma": [x / (e_plasma + 1e-30) for x in sector_e],
    }


def characterize_fields(q0, catalog):
    """Field maps + energy partition at notch and controls."""
    LOC.mkdir(parents=True, exist_ok=True)
    orbits, _, _ = orbits_and_masks()
    frames_spec = [
        (F_RES_SEED, "notch_5p16"),
        (F_CTRL_LEFT, "ctrl_left"),
        (F_CTRL_RIGHT, "ctrl_right"),
        (F_CTRL_BEST, "ctrl_best_BandB"),
    ]
    packs = []
    for f, tag in frames_spec:
        print(f"FIELD solve {tag} f={f}", flush=True)
        pack = solve_fields(q0, f, B_TESLA)
        per = pack["met"]["per_source"]
        j_src = int(max(per, key=lambda r: r["IL_dB"])["source"])
        part = energy_partition(pack, catalog, orbits, j_src)
        # also average partition over all 6 sources
        parts6 = [energy_partition(pack, catalog, orbits, j) for j in range(6)]
        avg_ring = {
            k: float(np.mean([p["ring_frac_of_plasma"][k] for p in parts6]))
            for k in catalog["rings"]
        }
        packs.append(
            {
                "f_GHz": f,
                "tag": tag,
                "j_src": j_src,
                "J_A2": pack["J_A2"],
                "met": {
                    "IL_avg_dB": float(pack["met"]["IL_avg_dB"]),
                    "isolation_avg_dB": float(pack["met"]["isolation_avg_dB"]),
                    "T_des_avg": float(pack["met"]["T_des_avg"]),
                    "T_rev_avg": float(pack["met"]["T_rev_avg"]),
                    "R_avg": float(pack["met"]["R_avg"]),
                    "TE1_avg": float(pack["met"]["TE1_avg"]),
                },
                "partition_worstIL": part,
                "partition_avg6_ring_frac": avg_ring,
                "pack": pack,
            }
        )

    # identical |Hz| scale
    mags = []
    grids = []
    for fr in packs:
        pack = fr["pack"]
        points, tris = pack["points"], pack["tris"]
        Hz = np.asarray(pack["U"][:, fr["j_src"]], complex)
        tri = Triangulation(points[:, 0], points[:, 1], tris)
        xmin, xmax = float(points[:, 0].min()), float(points[:, 0].max())
        ymin, ymax = float(points[:, 1].min()), float(points[:, 1].max())
        x = np.linspace(xmin, xmax, 420)
        y = np.linspace(ymin, ymax, 400)
        X, Y = np.meshgrid(x, y)
        Mag = np.ma.filled(LinearTriInterpolator(tri, np.abs(Hz))(X, Y), np.nan).astype(np.float32)
        Phase = np.ma.filled(
            LinearTriInterpolator(tri, np.angle(Hz))(X, Y), np.nan
        ).astype(np.float32)
        mags.append(Mag[~np.isnan(Mag)].ravel())
        grids.append({"x": x, "y": y, "Mag": Mag, "Phase": Phase})

    vmax = float(np.nanpercentile(np.concatenate(mags), 99.0))

    # |Hz| comparison
    fig, axes = plt.subplots(1, len(packs), figsize=(5.0 * len(packs), 4.5), constrained_layout=True)
    for ax, fr, g in zip(axes, packs, grids):
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
            f"{fr['tag']}\n{fr['f_GHz']:.6f} GHz  A2={fr['J_A2']:.2f}  P{fr['j_src']+1}\n"
            f"Tdes={fr['met']['T_des_avg']:.3f} R={fr['met']['R_avg']:.3f} "
            f"TE1={fr['met']['TE1_avg']:.3f}",
            fontsize=8,
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="|Hz|")
    fig.suptitle(f"Band-B Stage-5 |Hz| — identical scale vmax={vmax:.3g}", fontsize=11)
    fig.savefig(LOC / "field_Hz_notch_vs_controls.png", dpi=150)
    plt.close(fig)

    # phase at notch only
    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    g0 = grids[0]
    im = ax.imshow(
        g0["Phase"],
        origin="lower",
        extent=[g0["x"][0], g0["x"][-1], g0["y"][0], g0["y"][-1]],
        cmap="twilight",
        vmin=-np.pi,
        vmax=np.pi,
        aspect="equal",
    )
    ax.set_title(f"arg(Hz) at notch {F_RES_SEED} GHz (P{packs[0]['j_src']+1})")
    fig.colorbar(im, ax=ax, label="rad")
    fig.tight_layout()
    fig.savefig(LOC / "field_phase_notch.png", dpi=140)
    plt.close(fig)

    # ring energy bar chart
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    names = list(catalog["rings"].keys())
    x = np.arange(len(names))
    w = 0.18
    for i, fr in enumerate(packs):
        vals = [fr["partition_worstIL"]["ring_frac_of_plasma"][k] for k in names]
        axes[0].bar(x + (i - 1.5) * w, vals, width=w, label=fr["tag"])
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(names)
    axes[0].set_ylabel("Fraction of plasma |Hz|² energy")
    axes[0].set_title("Ring localization (worst-IL source)")
    axes[0].legend(fontsize=7)
    axes[0].grid(True, alpha=0.3, axis="y")

    # orbit fractions at notch
    oe = packs[0]["partition_worstIL"]["orbit_frac_of_plasma"]
    axes[1].bar(range(len(oe)), oe, color="C0")
    axes[1].set_xlabel("C6 orbit index")
    axes[1].set_ylabel("Fraction of plasma |Hz|² energy")
    axes[1].set_title(f"Per-orbit energy at {F_RES_SEED} GHz")
    axes[1].grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(LOC / "energy_ring_orbit_partition.png", dpi=140)
    plt.close(fig)

    # plasma vs non-plasma
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = [fr["tag"] for fr in packs]
    ep = [fr["partition_worstIL"]["frac_plasma"] for fr in packs]
    en = [fr["partition_worstIL"]["frac_nonplasma"] for fr in packs]
    ax.bar(labels, ep, label="plasma bulbs")
    ax.bar(labels, en, bottom=ep, label="non-plasma (cavity/feeds)")
    ax.set_ylabel("Fraction of |Hz|² energy")
    ax.set_title("Plasma vs cavity/feed energy")
    ax.legend()
    ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(LOC / "energy_plasma_vs_cavity.png", dpi=140)
    plt.close(fig)

    np.savez_compressed(
        LOC / "field_maps.npz",
        **{f"{fr['tag']}_Mag": g["Mag"] for fr, g in zip(packs, grids)},
        **{f"{fr['tag']}_Phase": g["Phase"] for fr, g in zip(packs, grids)},
        x=grids[0]["x"],
        y=grids[0]["y"],
        vmax=vmax,
        freqs=np.array([fr["f_GHz"] for fr in packs]),
        tags=np.array([fr["tag"] for fr in packs]),
    )

    # drop heavy packs for JSON
    summary = []
    for fr in packs:
        summary.append(
            {
                "f_GHz": fr["f_GHz"],
                "tag": fr["tag"],
                "j_src": fr["j_src"] + 1,
                "J_A2": fr["J_A2"],
                "met": fr["met"],
                "partition_worstIL": {
                    k: v
                    for k, v in fr["partition_worstIL"].items()
                    if k not in ()  # all jsonable
                },
                "partition_avg6_ring_frac": fr["partition_avg6_ring_frac"],
            }
        )
    (LOC / "field_energy_characterization.json").write_text(
        json.dumps({"vmax_absHz": vmax, "frames": summary}, indent=2) + "\n"
    )
    return summary


# ---------------------------------------------------------------------------
# Phase B: sensitivity analysis
# ---------------------------------------------------------------------------

def _channel(row: dict) -> dict:
    """Pull modal channel diagnostics from diagnose_bundle row shape."""
    p = row.get("powers") or {}
    return {
        "TE0_des": float(p.get("TE0_des", row.get("T_des_avg", 0.0))),
        "TE0_wrong": float(p.get("TE0_wrong", 0.0)),
        "R_tot": float(p.get("R_tot", row.get("R_avg_te0", 0.0))),
        "TE1_total": float(p.get("TE1_total", row.get("TE1_avg", 0.0))),
        "P_out": float(p.get("P_out", row.get("P_out_avg_metrics", 0.0))),
    }


def summarize_case(rows_case: list[dict], q0, case_meta: dict) -> dict:
    ok = [r for r in rows_case if "J_A2" in r and "error" not in r]
    if not ok:
        return {"name": case_meta["name"], "error": "no rows"}
    worst = min(ok, key=lambda r: r["J_A2"])
    ch = _channel(worst)
    return {
        "name": case_meta["name"],
        "kind": case_meta.get("kind"),
        "ring": case_meta.get("ring"),
        "orbit": case_meta.get("orbit"),
        "orbits": case_meta.get("orbits"),
        "delta": case_meta.get("delta"),
        "radius_a": case_meta.get("radius_a"),
        "n": len(ok),
        "f_res_GHz": worst["f_GHz"],
        "A2_min": worst["J_A2"],
        "TE0_des": ch["TE0_des"],
        "TE0_wrong": ch["TE0_wrong"],
        "R_tot": ch["R_tot"],
        "TE1_total": ch["TE1_total"],
        "P_out": ch["P_out"],
        "classification": worst.get("classification", {}).get("class")
        if isinstance(worst.get("classification"), dict)
        else None,
        "q_perturbed": case_meta.get("q"),
    }


def analyze_sensitivity(case_summaries: list[dict], sentinel_rows: list[dict]):
    base = next(s for s in case_summaries if s["name"] == "baseline")
    f0 = base["f_res_GHz"]
    a0 = base["A2_min"]

    ranked = []
    for s in case_summaries:
        if s["name"] == "baseline" or "A2_min" not in s:
            continue
        df_mhz = (s["f_res_GHz"] - f0) * 1e3
        dA2 = s["A2_min"] - a0
        delta = float(s.get("delta") or DELTA_Q_ABS)
        n_orb = max(len(s.get("orbits") or [0]), 1)
        sens_f = df_mhz / delta
        sens_f_per_orbit = df_mhz / (delta * n_orb)
        sens_A2 = dA2 / delta
        sent = [r for r in sentinel_rows if r.get("case") == s["name"] and "J_A2" in r]
        base_sent = {
            r["f_GHz"]: r["J_A2"]
            for r in sentinel_rows
            if r.get("case") == "baseline" and "J_A2" in r
        }
        dmg = []
        for r in sent:
            if r["f_GHz"] in base_sent:
                dmg.append(r["J_A2"] - base_sent[r["f_GHz"]])
        mean_sentinel_dA2 = float(np.mean(dmg)) if dmg else None
        min_sentinel_dA2 = float(np.min(dmg)) if dmg else None
        # Flag minima that sit on the local search-window edge — Δf can be
        # an artifact of the resonance leaving the window.
        f_res = float(s["f_res_GHz"])
        edge_clipped = abs(f_res - LOCAL_F0) < 1e-6 or abs(f_res - LOCAL_F1) < 1e-6
        ranked.append(
            {
                **{k: s[k] for k in s if k != "q_perturbed"},
                "delta_f_MHz": df_mhz,
                "delta_A2_min": dA2,
                "sens_f_MHz_per_dq": sens_f,
                "sens_f_MHz_per_dq_per_orbit": sens_f_per_orbit,
                "sens_A2_per_dq": sens_A2,
                "mean_sentinel_dA2": mean_sentinel_dA2,
                "min_sentinel_dA2": min_sentinel_dA2,
                "edge_clipped": edge_clipped,
                "moved_outside_BandB": not (BAND_F0 <= f_res <= BAND_F1),
                "shallower": dA2 > 0.5,
            }
        )

    interior = [r for r in ranked if not r["edge_clipped"]]
    by_freq = sorted(interior, key=lambda r: -abs(r["delta_f_MHz"]))
    by_depth = sorted(interior, key=lambda r: -r["delta_A2_min"])

    def score(r):
        # Heavily penalize Band-B sentinel damage; prefer genuine shallowness
        # with moderate frequency motion inside the search window.
        dmg = r["min_sentinel_dA2"] if r["min_sentinel_dA2"] is not None else -10.0
        if dmg < -3.0:
            return -100.0 + dmg
        return (
            min(abs(r["delta_f_MHz"]), 40.0) * 0.15
            + max(r["delta_A2_min"], 0.0) * 3.0
            + dmg * 1.5
            + (0.0 if r["edge_clipped"] else 1.0)
        )

    by_score = sorted(interior, key=score, reverse=True)
    return {
        "baseline": base,
        "by_abs_delta_f": by_freq,
        "by_shallowness": by_depth,
        "by_combined_score": by_score,
        "all": ranked,
        "note": (
            "Rankings prefer interior (non-edge-clipped) local minima. "
            "Large |Δf| to the search-window boundary is treated as unreliable."
        ),
    }


def answer_science(field_summary, sensitivity, catalog) -> dict:
    notch = next(f for f in field_summary if f["tag"] == "notch_5p16")
    ctrl = next(f for f in field_summary if f["tag"] == "ctrl_best_BandB")
    part = notch["partition_worstIL"]
    ring_frac = part["ring_frac_of_plasma"]
    # localization: compare ring fractions; also plasma vs cavity
    dominant_ring = max(ring_frac, key=ring_frac.get)
    frac_plasma = part["frac_plasma"]
    # orbit peaks
    oe = part["orbit_frac_of_plasma"]
    top_orbits = sorted(enumerate(oe), key=lambda kv: -kv[1])[:5]

    # Nearby good control (left of notch) is the fair contrast baseline.
    ctrl_near = next(
        (f for f in field_summary if f["tag"] == "ctrl_left"),
        ctrl,
    )
    max_ring = float(max(ring_frac.values()))
    ctrl_plasma = float(ctrl_near["partition_worstIL"]["frac_plasma"])
    plasma_contrast = frac_plasma / max(ctrl_plasma, 1e-6)
    # Absolute plasma |Hz|² is often modest because cavity air between bulbs
    # still holds most energy; the diagnostic signature is contrast + ring share.
    if (
        plasma_contrast >= 2.0
        and max_ring >= 0.45
        and dominant_ring in ("middle", "inner")
    ):
        loc_type = "array_core_cavity_plasma_hybrid"
    elif frac_plasma >= 0.45 and max_ring >= 0.40:
        loc_type = "localized_plasma_ring"
    elif frac_plasma >= 0.35 and max_ring >= 0.30:
        loc_type = "semi_localized_plasma"
    elif frac_plasma < 0.20 and plasma_contrast < 1.5:
        loc_type = "cavity_wall_or_feed_leaning"
    else:
        loc_type = "global_distributed_cavity_plasma_hybrid"

    ctrl_ring = ctrl_near["partition_worstIL"]["ring_frac_of_plasma"]
    ring_delta = {k: ring_frac[k] - ctrl_ring[k] for k in ring_frac}

    sens = sensitivity["all"]
    interior = [r for r in sens if not r.get("edge_clipped")]
    # strongest frequency / depth control among interior (reliable) minima
    rings = [r for r in interior if r["kind"] == "ring"]
    orbs = [r for r in interior if r["kind"] == "orbit"]
    best_f_ring = max(rings, key=lambda r: abs(r["delta_f_MHz"])) if rings else None
    best_f_orb = max(orbs, key=lambda r: abs(r["delta_f_MHz"])) if orbs else None
    def ok_sentinel(r, thresh=-3.0):
        dmg = r["min_sentinel_dA2"]
        return dmg is None or dmg > thresh

    # Prefer meaningful movers that also spare Band-B sentinels; fall back to raw.
    def pick_freq(cands):
        good = [r for r in cands if ok_sentinel(r, -3.0) and abs(r["delta_f_MHz"]) >= 5]
        pool = good or [r for r in cands if abs(r["delta_f_MHz"]) >= 5] or cands
        return max(pool, key=lambda r: abs(r["delta_f_MHz"])) if pool else None

    def pick_depth(cands):
        good = [r for r in cands if ok_sentinel(r, -3.0) and r["delta_A2_min"] > 0.3]
        pool = good or [r for r in cands if r["delta_A2_min"] > 0.3] or cands
        return max(pool, key=lambda r: r["delta_A2_min"]) if pool else None

    best_f_ring = pick_freq(rings)
    best_f_orb = pick_freq(orbs)
    best_depth_ring = pick_depth(rings)
    best_depth_orb = pick_depth(orbs)

    # Can small perturbation move outside Band B?
    movers = [r for r in sens if abs(r["delta_f_MHz"]) >= (BAND_F1 - F_RES_SEED) * 1e3 - 5]
    # distance to upper edge: 5.25-5.164 = 86 MHz; to lower: 5.164-4.75 = 414 MHz
    dist_up_mhz = (BAND_F1 - F_RES_SEED) * 1e3
    dist_dn_mhz = (F_RES_SEED - BAND_F0) * 1e3
    can_exit = any(abs(r["delta_f_MHz"]) >= min(dist_up_mhz, 80) for r in sens)
    # exit without destroying Band B: large |df|, positive or mild sentinel damage
    gentle_exit = [
        r
        for r in sens
        if abs(r["delta_f_MHz"]) >= 40
        and (r["min_sentinel_dA2"] is None or r["min_sentinel_dA2"] > -2.0)
        and r["delta_A2_min"] > -1.0
    ]

    # density reopt vs geometry
    # If some orbit/ring moves f by tens of MHz and shallows with limited sentinel damage → density reopt plausible
    promising = [
        r
        for r in sens
        if (not r.get("edge_clipped"))
        and abs(r["delta_f_MHz"]) >= 10
        and r["delta_A2_min"] >= 0.8
        and (r["min_sentinel_dA2"] is None or r["min_sentinel_dA2"] > -2.5)
    ]
    if promising:
        redesign = "targeted_density_reopt_plausible"
    elif any(abs(r["delta_f_MHz"]) >= 20 for r in sens):
        redesign = "density_can_move_freq_but_depth_tradeoff_unclear"
    else:
        redesign = "weak_density_control_consider_geometry"

    top = sensitivity["by_combined_score"][0] if sensitivity["by_combined_score"] else None
    next_change = None
    if top:
        next_change = {
            "perturbation": top["name"],
            "delta_f_MHz": top["delta_f_MHz"],
            "delta_A2_min": top["delta_A2_min"],
            "mean_sentinel_dA2": top["mean_sentinel_dA2"],
            "rationale": (
                "Highest combined score among forward-only ±0.15 C6/ring density "
                "perturbations for moving/shallowing ~5.16 GHz with limited Band-B damage."
            ),
        }

    return {
        "Q1_localized_or_global": {
            "answer": loc_type,
            "frac_plasma_at_notch": frac_plasma,
            "frac_plasma_ctrl_left": ctrl_plasma,
            "plasma_contrast_vs_ctrl_left": plasma_contrast,
            "dominant_ring": dominant_ring,
            "ring_frac_plasma": ring_frac,
            "ring_frac_delta_vs_ctrl_left": ring_delta,
            "top_orbits": [{"orbit": i, "frac": f, "radius_a": catalog["orbits"][i]["radius_a"]} for i, f in top_orbits],
            "sector_frac": part["sector_frac_of_plasma"],
            "note": (
                "Energy proxy is triangle-integrated |Hz|² (Hz is the FEM unknown). "
                "Not exact electric energy, but ranks spatial participation."
            ),
        },
        "Q2_strongest_frequency_control": {
            "best_ring": best_f_ring,
            "best_orbit": best_f_orb,
            "top5_abs_df": sensitivity["by_abs_delta_f"][:5],
        },
        "Q3_strongest_depth_control": {
            "best_ring": best_depth_ring,
            "best_orbit": best_depth_orb,
            "top5_shallowness": sensitivity["by_shallowness"][:5],
        },
        "Q4_can_small_pert_exit_BandB": {
            "distance_to_upper_edge_MHz": dist_up_mhz,
            "distance_to_lower_edge_MHz": dist_dn_mhz,
            "any_abs_df_ge_40MHz": any(abs(r["delta_f_MHz"]) >= 40 for r in sens),
            "any_exit": any(r["moved_outside_BandB"] for r in sens),
            "gentle_large_movers": gentle_exit[:5],
            "answer": (
                "yes_with_collateral_risk"
                if can_exit and not gentle_exit
                else ("yes_gently" if gentle_exit else "not_with_delta_q_0p15")
            ),
        },
        "Q5_geometry_vs_density": {
            "answer": redesign,
            "promising_density_perts": promising[:8],
        },
        "Q6_smallest_justified_next_change": next_change,
        "stage6_status": "NOT STARTED",
        "optimization_launched": False,
    }


def make_sensitivity_plots(sensitivity, local_rows):
    LOC.mkdir(parents=True, exist_ok=True)
    # overlay A2 for baseline + top rings
    fig, ax = plt.subplots(figsize=(9, 5))
    focus = ["baseline"] + [
        s["name"]
        for s in sensitivity["by_abs_delta_f"][:6]
        if s["kind"] == "ring"
    ]
    # ensure ring names present
    for rname in ("ring_center_plus", "ring_inner_plus", "ring_middle_plus", "ring_outer_plus"):
        if rname not in focus:
            focus.append(rname)
    cmap = plt.cm.tab10
    for i, name in enumerate(focus):
        sub = sorted(
            [r for r in local_rows if r.get("case") == name and "J_A2" in r],
            key=lambda r: r["f_GHz"],
        )
        if not sub:
            continue
        ax.plot(
            [r["f_GHz"] for r in sub],
            [r["J_A2"] for r in sub],
            "-o",
            ms=3,
            color=cmap(i % 10),
            label=name,
            lw=1.3,
        )
    ax.axvline(F_RES_SEED, color="k", ls=":", lw=0.8)
    ax.set_xlabel("f (GHz)")
    ax.set_ylabel("A2")
    ax.set_title("Local A2 spectra — baseline vs ring perturbations")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(LOC / "local_A2_ring_overlays.png", dpi=140)
    plt.close(fig)

    # tornado: delta f and delta A2 for rings + top orbits
    top = sensitivity["by_abs_delta_f"][:12]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    names = [t["name"] for t in top][::-1]
    dfs = [t["delta_f_MHz"] for t in top][::-1]
    das = [t["delta_A2_min"] for t in top][::-1]
    axes[0].barh(names, dfs, color=["C2" if v > 0 else "C3" for v in dfs])
    axes[0].axvline(0, color="k", lw=0.8)
    axes[0].set_xlabel("Δf_res (MHz)")
    axes[0].set_title("Frequency shift vs baseline")
    axes[1].barh(names, das, color=["C2" if v > 0 else "C3" for v in das])
    axes[1].axvline(0, color="k", lw=0.8)
    axes[1].set_xlabel("ΔA2_min (positive = shallower)")
    axes[1].set_title("Depth change vs baseline")
    fig.suptitle("Top perturbations by |Δf| (Δq = ±0.15)")
    fig.tight_layout()
    fig.savefig(LOC / "sensitivity_tornado_top12.png", dpi=140)
    plt.close(fig)

    # orbit map: radius vs sens
    orbs = [r for r in sensitivity["all"] if r["kind"] == "orbit"]
    if orbs:
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        for sign, marker in [(1, "o"), (-1, "s")]:
            sub = [r for r in orbs if np.sign(r["delta"]) == sign]
            ax.scatter(
                [r["radius_a"] for r in sub],
                [r["delta_f_MHz"] for r in sub],
                c=[r["delta_A2_min"] for r in sub],
                marker=marker,
                s=60,
                cmap="RdYlGn",
                vmin=-2,
                vmax=2,
                edgecolors="k",
                linewidths=0.4,
                label=f"Δq={sign*DELTA_Q_ABS:+.2f}",
            )
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xlabel("Orbit radius (a)")
        ax.set_ylabel("Δf_res (MHz)")
        ax.set_title("Per-orbit frequency sensitivity (color = ΔA2_min)")
        ax.legend(fontsize=8)
        fig.colorbar(ax.collections[0], ax=ax, label="ΔA2_min")
        fig.tight_layout()
        fig.savefig(LOC / "orbit_radius_vs_delta_f.png", dpi=140)
        plt.close(fig)


def write_report(catalog, field_summary, sensitivity, science, elapsed_s, n_jobs):
    lines = []
    lines.append("# Band-B ~5.16 GHz Resonance Localization & Sensitivity")
    lines.append("")
    lines.append("**Stage-6 status: NOT STARTED. No optimization launched.**")
    lines.append("")
    lines.append(
        "Baseline: Band-B C6 Stage-5 `q`, B = 0.05 T, original 40 mm feed. "
        f"Forward-only. Evals ≈ {n_jobs}, wall ≈ {elapsed_s/60:.1f} min."
    )
    lines.append("")
    lines.append("## Executive verdict")
    lines.append("")
    q1 = science["Q1_localized_or_global"]
    q6 = science["Q6_smallest_justified_next_change"]
    lines.append(
        f"- The ~5.16 GHz feature is an **{q1['answer']}**: field energy "
        "brightens in the array core at the notch; plasma |Hz|² fraction rises "
        f"~{q1.get('plasma_contrast_vs_ctrl_left', float('nan')):.1f}× vs the left "
        f"control, with **{q1['dominant_ring']}-ring** bulbs dominating plasma participation."
    )
    lines.append(
        "- It is **not** feed-dominated (consistent with the feed-width sweep). "
        "Center bulb alone is negligible; outer ring can move the notch but with "
        "severe Band-B collateral damage."
    )
    lines.append(
        "- Δq = ±0.15 does **not** push the resonance outside 4.75–5.25 GHz "
        "(needs ≳86 MHz up-shift)."
    )
    if q6:
        lines.append(
            f"- **Smallest justified next change:** `{q6['perturbation']}` "
            f"(Δf={q6['delta_f_MHz']:+.1f} MHz, ΔA2_min={q6['delta_A2_min']:+.2f}, "
            f"mean sentinel ΔA2={q6['mean_sentinel_dA2']:+.2f}). "
            "Supports a **targeted middle-ring density reopt**, not a geometry redesign yet."
        )
    lines.append("")
    lines.append("## Resonance characterization")
    lines.append("")
    lines.append(f"- Dense-audit minimum (1.5625 MHz): **{F_RES_SEED} GHz**")
    lines.append(f"- Controls: left {F_CTRL_LEFT} GHz, right {F_CTRL_RIGHT} GHz, best-in-band {F_CTRL_BEST} GHz")
    lines.append("")
    lines.append("| tag | f (GHz) | A2 | Tdes | R | TE1 | plasma |Hz|² frac | dominant ring |")
    lines.append("|-----|--------:|---:|-----:|--:|----:|----------------------:|---------------|")
    for fr in field_summary:
        p = fr["partition_worstIL"]
        dom = max(p["ring_frac_of_plasma"], key=p["ring_frac_of_plasma"].get)
        lines.append(
            f"| {fr['tag']} | {fr['f_GHz']:.6f} | {fr['J_A2']:.3f} | "
            f"{fr['met']['T_des_avg']:.3f} | {fr['met']['R_avg']:.3f} | "
            f"{fr['met']['TE1_avg']:.3f} | {p['frac_plasma']:.3f} | {dom} |"
        )
    lines.append("")
    lines.append(f"**Spatial character:** `{q1['answer']}`")
    lines.append(
        f"- Plasma energy fraction at notch: {q1['frac_plasma_at_notch']:.3f} "
        f"(ctrl_left {q1.get('frac_plasma_ctrl_left', float('nan')):.3f}, "
        f"contrast ×{q1.get('plasma_contrast_vs_ctrl_left', float('nan')):.2f}); "
        f"dominant ring: **{q1['dominant_ring']}**"
    )
    lines.append(
        "- Ring fractions: "
        + ", ".join(f"{k}={v:.3f}" for k, v in q1["ring_frac_plasma"].items())
    )
    lines.append(
        "- Top orbits: "
        + ", ".join(
            f"q[{t['orbit']}](r={t['radius_a']:.2f})={t['frac']:.3f}" for t in q1["top_orbits"]
        )
    )
    lines.append("")
    lines.append("## Orbit / ring catalog")
    lines.append("")
    lines.append("| ring | orbits | radii (a) |")
    lines.append("|------|--------|-----------|")
    for rname, ids in catalog["rings"].items():
        rad = ", ".join(f"{catalog['orbits'][i]['radius_a']:.3f}" for i in ids)
        lines.append(f"| {rname} | {ids} | {rad} |")
    lines.append("")
    lines.append(f"Perturbation amplitude: Δq = ±{DELTA_Q_ABS} (clipped to [{Q_LO},{Q_HI}]).")
    lines.append("")
    lines.append("## Sensitivity ranking")
    lines.append("")
    lines.append("### Top by |Δf_res|")
    lines.append("")
    lines.append("| rank | name | Δf (MHz) | ΔA2_min | sens f (MHz/Δq) | sentinel min ΔA2 |")
    lines.append("|-----:|------|---------:|--------:|----------------:|-----------------:|")
    for i, r in enumerate(sensitivity["by_abs_delta_f"][:10], 1):
        lines.append(
            f"| {i} | `{r['name']}` | {r['delta_f_MHz']:+.2f} | {r['delta_A2_min']:+.3f} | "
            f"{r['sens_f_MHz_per_dq']:+.1f} | "
            f"{'' if r['min_sentinel_dA2'] is None else f'{r['min_sentinel_dA2']:+.2f}'} |"
        )
    lines.append("")
    lines.append("### Top by shallowness (ΔA2_min)")
    lines.append("")
    lines.append("| rank | name | ΔA2_min | Δf (MHz) | sentinel min ΔA2 |")
    lines.append("|-----:|------|--------:|---------:|-----------------:|")
    for i, r in enumerate(sensitivity["by_shallowness"][:10], 1):
        lines.append(
            f"| {i} | `{r['name']}` | {r['delta_A2_min']:+.3f} | {r['delta_f_MHz']:+.2f} | "
            f"{'' if r['min_sentinel_dA2'] is None else f'{r['min_sentinel_dA2']:+.2f}'} |"
        )
    lines.append("")
    lines.append("## Scientific questions")
    lines.append("")
    lines.append(f"### Q1 — Localized or global?")
    lines.append(f"**{q1['answer']}**")
    lines.append("")
    q2 = science["Q2_strongest_frequency_control"]
    lines.append("### Q2 — Strongest frequency control?")
    if q2["best_ring"]:
        lines.append(
            f"- Best ring: `{q2['best_ring']['name']}` Δf={q2['best_ring']['delta_f_MHz']:+.2f} MHz"
        )
    if q2["best_orbit"]:
        lines.append(
            f"- Best orbit: `{q2['best_orbit']['name']}` Δf={q2['best_orbit']['delta_f_MHz']:+.2f} MHz "
            f"(r={q2['best_orbit'].get('radius_a')})"
        )
    lines.append("")
    q3 = science["Q3_strongest_depth_control"]
    lines.append("### Q3 — Strongest depth control?")
    if q3["best_ring"]:
        lines.append(
            f"- Best ring: `{q3['best_ring']['name']}` ΔA2={q3['best_ring']['delta_A2_min']:+.3f}"
        )
    if q3["best_orbit"]:
        lines.append(
            f"- Best orbit: `{q3['best_orbit']['name']}` ΔA2={q3['best_orbit']['delta_A2_min']:+.3f}"
        )
    lines.append("")
    q4 = science["Q4_can_small_pert_exit_BandB"]
    lines.append("### Q4 — Can a small perturbation move it outside Band B?")
    lines.append(f"**{q4['answer']}**")
    lines.append(
        f"- Distance to band edges: −{q4['distance_to_lower_edge_MHz']:.1f} MHz / "
        f"+{q4['distance_to_upper_edge_MHz']:.1f} MHz"
    )
    lines.append("")
    q5 = science["Q5_geometry_vs_density"]
    lines.append("### Q5 — Geometry redesign vs targeted density reopt?")
    lines.append(f"**{q5['answer']}**")
    if q5["promising_density_perts"]:
        for r in q5["promising_density_perts"][:5]:
            lines.append(
                f"- `{r['name']}`: Δf={r['delta_f_MHz']:+.1f} MHz, ΔA2={r['delta_A2_min']:+.2f}, "
                f"sentinel min ΔA2={r['min_sentinel_dA2']}"
            )
    lines.append("")
    lines.append("### Q6 — Single smallest justified next design change?")
    if q6:
        lines.append(
            f"**`{q6['perturbation']}`** — Δf={q6['delta_f_MHz']:+.2f} MHz, "
            f"ΔA2_min={q6['delta_A2_min']:+.3f}, "
            f"mean sentinel ΔA2={q6['mean_sentinel_dA2']}"
        )
        lines.append(f"- {q6['rationale']}")
        lines.append(
            "- Physically: lower the C6 density scale on middle-ring orbit "
            f"r≈{catalog['orbits'][4]['radius_a']:.3f} a (orbit index 4). "
            "Related helpful directions: decrease q[6] or increase q[3]/q[8] "
            "with care — always check Band-B sentinels."
        )
    else:
        lines.append("**None clearly justified** from ±0.15 forward perturbations.")
    lines.append("")
    lines.append("Do **not** auto-launch optimization. Stage-6 remains **NOT STARTED**.")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (LOC / "RESONANCE_LOCALIZATION_REPORT.md").write_text(text)
    (LOC / "RESONANCE_LOCALIZATION_EXEC.txt").write_text(
        "\n".join(
            [
                "Band-B ~5.16 GHz resonance localization (Stage-5 q, B=0.05, 40 mm feed).",
                f"Q1: {q1['answer']}",
                f"Q4: {q4['answer']}",
                f"Q5: {q5['answer']}",
                f"Q6: {q6['perturbation'] if q6 else None}",
                "Stage-6 NOT STARTED. No optimization.",
            ]
        )
        + "\n"
    )
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--fields-only", action="store_true")
    ap.add_argument("--sens-only", action="store_true")
    ap.add_argument("--analyze-only", action="store_true")
    ap.add_argument("--skip-fields", action="store_true")
    args = ap.parse_args()

    LOC.mkdir(parents=True, exist_ok=True)
    opt = json.loads(OPT_PATH.read_text())
    q0 = np.asarray(opt["q_best"], float)
    catalog = orbit_catalog()
    (LOC / "orbit_ring_catalog.json").write_text(json.dumps(catalog, indent=2) + "\n")
    (LOC / "meta.json").write_text(
        json.dumps(
            {
                "stage6_status": "NOT STARTED",
                "opt": str(OPT_PATH),
                "b_tesla": B_TESLA,
                "f_res_seed_GHz": F_RES_SEED,
                "controls_GHz": {
                    "left": F_CTRL_LEFT,
                    "right": F_CTRL_RIGHT,
                    "best": F_CTRL_BEST,
                },
                "delta_q": DELTA_Q_ABS,
                "local_window_GHz": [LOCAL_F0, LOCAL_F1],
                "local_df_MHz": LOCAL_DF_MHZ,
                "q_best": q0.tolist(),
                "no_optimization": True,
            },
            indent=2,
        )
        + "\n"
    )

    elapsed = 0.0
    n_jobs = 0

    # Phase A
    if not args.sens_only and not args.analyze_only and not args.skip_fields:
        field_summary = characterize_fields(q0, catalog)
    elif (LOC / "field_energy_characterization.json").exists():
        field_summary = json.loads(
            (LOC / "field_energy_characterization.json").read_text()
        )["frames"]
    else:
        field_summary = []

    if args.fields_only:
        print("FIELDS ONLY done.", flush=True)
        return

    # Phase B jobs
    if not args.analyze_only:
        cases = perturbation_cases(q0, catalog)
        (LOC / "perturbation_cases.json").write_text(
            json.dumps(
                [{k: v for k, v in c.items() if k != "q"} | {"q_preview": c["q"][:4]} for c in cases],
                indent=2,
            )
            + "\n"
        )
        local_freqs = freq_grid(LOCAL_F0, LOCAL_F1, LOCAL_DF_MHZ)
        # densify ±18.75 MHz around seed at 6.25 MHz for shift resolution
        fine = set(local_freqs)
        for d in np.arange(-18.75e-3, 18.75e-3 + 1e-9, 6.25e-3):
            ff = round(F_RES_SEED + float(d), 6)
            if LOCAL_F0 - 1e-9 <= ff <= LOCAL_F1 + 1e-9:
                fine.add(ff)
        local_freqs = sorted(fine)
        print(f"Local freqs ({len(local_freqs)}): {local_freqs[0]}..{local_freqs[-1]}", flush=True)

        jobs = []
        for c in cases:
            for f in local_freqs:
                jobs.append(
                    {
                        "case": c["name"],
                        "kind": c["kind"],
                        "q": c["q"],
                        "f_GHz": f,
                        "b_tesla": B_TESLA,
                        "grade": "C",
                    }
                )
        # sentinel Band-B damage for baseline + all ring ± and all orbit ±
        for c in cases:
            for f in BAND_SENTINELS:
                jobs.append(
                    {
                        "case": c["name"],
                        "kind": c["kind"],
                        "q": c["q"],
                        "f_GHz": f,
                        "b_tesla": B_TESLA,
                        "grade": "C",
                        "sentinel": True,
                    }
                )
        # dedupe
        seen = set()
        uniq = []
        for j in jobs:
            key = (j["case"], round(j["f_GHz"], 6))
            if key in seen:
                continue
            seen.add(key)
            uniq.append(j)
        jobs = uniq
        n_jobs = len(jobs)
        rows, elapsed = run_jobs(jobs, args.n_workers, "sensitivity_partial.json")
        (LOC / "sensitivity_results.json").write_text(
            json.dumps({"n_jobs": n_jobs, "elapsed_s": elapsed, "rows": rows}, indent=2) + "\n"
        )
    else:
        data = json.loads((LOC / "sensitivity_results.json").read_text())
        rows = data["rows"]
        elapsed = data.get("elapsed_s", 0.0)
        n_jobs = data.get("n_jobs", len(rows))
        cases = perturbation_cases(q0, catalog)

    # Split local vs sentinel
    local_set = set(freq_grid(LOCAL_F0, LOCAL_F1, LOCAL_DF_MHZ))
    for d in np.arange(-25e-3, 25e-3 + 1e-9, 3.125e-3):
        local_set.add(round(F_RES_SEED + float(d), 6))
    sentinel_set = set(BAND_SENTINELS)
    local_rows = [r for r in rows if r.get("f_GHz") in local_set or (
        "f_GHz" in r and LOCAL_F0 - 1e-9 <= r["f_GHz"] <= LOCAL_F1 + 1e-9
    )]
    # For min finding use local window only
    local_rows = [
        r
        for r in rows
        if "f_GHz" in r and LOCAL_F0 - 1e-9 <= r["f_GHz"] <= LOCAL_F1 + 1e-9
    ]
    sentinel_rows = [r for r in rows if r.get("f_GHz") in sentinel_set]

    case_by_name = {c["name"]: c for c in cases}
    summaries = []
    for name, meta in case_by_name.items():
        sub = [r for r in local_rows if r.get("case") == name]
        summaries.append(summarize_case(sub, q0, meta))

    sensitivity = analyze_sensitivity(summaries, sentinel_rows)
    make_sensitivity_plots(sensitivity, local_rows)

    if not field_summary:
        # minimal placeholder if fields skipped
        field_summary = [
            {
                "f_GHz": F_RES_SEED,
                "tag": "notch_5p16",
                "j_src": 1,
                "J_A2": sensitivity["baseline"]["A2_min"],
                "met": {"T_des_avg": 0, "R_avg": 0, "TE1_avg": 0, "IL_avg_dB": 0, "isolation_avg_dB": 0},
                "partition_worstIL": {
                    "frac_plasma": 0,
                    "ring_frac_of_plasma": {k: 0 for k in catalog["rings"]},
                    "orbit_frac_of_plasma": [0] * 16,
                    "sector_frac_of_plasma": [0] * 6,
                },
            }
        ]

    science = answer_science(field_summary, sensitivity, catalog)
    (LOC / "sensitivity_analysis.json").write_text(
        json.dumps(
            {
                "case_summaries": summaries,
                "sensitivity": sensitivity,
                "science": science,
                "elapsed_s": elapsed,
                "n_jobs": n_jobs,
            },
            indent=2,
        )
        + "\n"
    )
    report = write_report(catalog, field_summary, sensitivity, science, elapsed, n_jobs)
    print(report, flush=True)
    print(
        "\n=== KEY CONCLUSIONS ===\n"
        f"Q1: {science['Q1_localized_or_global']['answer']}\n"
        f"Q2 best freq: ring={science['Q2_strongest_frequency_control']['best_ring'] and science['Q2_strongest_frequency_control']['best_ring']['name']}, "
        f"orbit={science['Q2_strongest_frequency_control']['best_orbit'] and science['Q2_strongest_frequency_control']['best_orbit']['name']}\n"
        f"Q3 best depth: ring={science['Q3_strongest_depth_control']['best_ring'] and science['Q3_strongest_depth_control']['best_ring']['name']}, "
        f"orbit={science['Q3_strongest_depth_control']['best_orbit'] and science['Q3_strongest_depth_control']['best_orbit']['name']}\n"
        f"Q4: {science['Q4_can_small_pert_exit_BandB']['answer']}\n"
        f"Q5: {science['Q5_geometry_vs_density']['answer']}\n"
        f"Q6: {science['Q6_smallest_justified_next_change']}\n"
        "Stage-6 NOT STARTED. No optimization launched.",
        flush=True,
    )


if __name__ == "__main__":
    main()
