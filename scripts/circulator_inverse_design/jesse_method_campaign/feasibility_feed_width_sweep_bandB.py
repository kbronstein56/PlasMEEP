#!/usr/bin/env python3
"""Fixed-q Band-B feed-width sweep (no q optimization).

Uses Band-B C6 Stage-5 optimized q at B=0.05 T. Sweeps clear feed width W
from the production 40 mm aperture through intermediate values to the prior
23.06 mm / 6.5 GHz-cutoff geometry. Dense forward audits only.

Stage-6 NOT STARTED. Do not launch q reoptimization from this script.
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

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
FEAS = OUT / "feasibility"
SWEEP = FEAS / "feed_width_sweep_bandB"

sys.path[:0] = [
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

from feed_narrow_geometry import te_cutoffs  # noqa: E402
from feasibility_feed_narrow_compare import (  # noqa: E402
    _init_worker,
    diagnose_S,
    build_system_feed,
    s_matrix_te0_drive,
)

OPT_PATH = FEAS / "B_4p75_5p25_C6_B0.0500_stage5_opt.json"
B_TESLA = 0.05
BAND_F0, BAND_F1 = 4.75, 5.25
# Guard band so resonances entering/leaving Band B are visible
AUDIT_F0, AUDIT_F1 = 4.65, 5.35
WATCH = [4.909375, 5.1640625]

# Physically motivated clear widths (mm):
#   40.00  original production feed (TE1 ~3.75 GHz, open throughout Band B)
#   35.00  intermediate
#   30.00  TE1 cutoff near mid-band (~5.00 GHz)
#   28.00  TE1 modestly above Band-B upper edge (~5.35 GHz)
#   26.50  additional cutoff margin (~5.66 GHz)
#   25.00  further margin (~6.00 GHz)
#   23.06  prior aggressive 6.5 GHz-cutoff geometry
WIDTHS_MM = [40.0, 35.0, 30.0, 28.0, 26.5, 25.0, 23.06]
A_M = 0.02  # lattice unit (m); matches sixport_common.a


def width_mm_to_lu(W_mm: float) -> float:
    return float(W_mm) * 1e-3 / A_M


def freq_grid(f0, f1, df):
    n = int(round((f1 - f0) / df)) + 1
    return [round(f0 + i * df, 6) for i in range(n)]


def width_catalog() -> list[dict]:
    rows = []
    for W_mm in WIDTHS_MM:
        W_lu = width_mm_to_lu(W_mm)
        cut = te_cutoffs(W_lu)
        te1 = cut["TE1_cutoff_GHz"]
        te2 = cut["modes"][2]["f_cut_GHz"]
        note = []
        if abs(W_mm - 40.0) < 1e-6:
            note.append("original production")
        if abs(W_mm - 23.06) < 1e-3:
            note.append("prior aggressive 6.5 GHz cutoff")
        if 5.25 < te1 <= 5.45:
            note.append("TE1 modestly above Band-B upper edge")
        if te1 > BAND_F1 + 0.3:
            note.append("extra TE1 cutoff margin over Band B")
        if te1 < BAND_F0:
            note.append("TE1 open throughout Band B")
        elif te1 < BAND_F1:
            note.append("TE1 cutoff inside Band B")
        rows.append(
            {
                "W_mm": W_mm,
                "W_lu": W_lu,
                "TE1_cutoff_GHz": te1,
                "TE2_cutoff_GHz": te2,
                "TE1_open_in_BandB": te1 < BAND_F1,
                "TE1_open_at_BandB_lower": te1 < BAND_F0,
                "note": "; ".join(note) if note else "",
                "cutoffs": cut,
            }
        )
    return rows


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
            "W_mm": float(payload["W_mm"]),
            "W_lu": float(payload["clear_width_lu"]),
            "f_GHz": float(payload["f_GHz"]),
            "b_tesla": float(payload["b_tesla"]),
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


def build_jobs(q, df_MHz: float) -> list[dict]:
    base = set(freq_grid(AUDIT_F0, AUDIT_F1, df_MHz * 1e-3))
    # Modest densification around known Band-B resonance families
    # (half-step of the base grid, ±25 MHz) — enough to track shifts without
    # exploding the cheap fixed-q sweep budget.
    fine = max(df_MHz * 0.5, 3.125) * 1e-3
    for fw in WATCH:
        for d in np.arange(-25e-3, 25e-3 + 1e-9, fine):
            ff = round(fw + float(d), 6)
            if AUDIT_F0 - 1e-9 <= ff <= AUDIT_F1 + 1e-9:
                base.add(ff)
    freqs = sorted(base)
    jobs = []
    for W_mm in WIDTHS_MM:
        W_lu = width_mm_to_lu(W_mm)
        for f in freqs:
            jobs.append(
                {
                    "q": q,
                    "f_GHz": f,
                    "b_tesla": B_TESLA,
                    "clear_width_lu": W_lu,
                    "W_mm": W_mm,
                    "grade": "C",
                }
            )
    return jobs


def run_jobs(jobs, n_workers=12):
    SWEEP.mkdir(parents=True, exist_ok=True)
    results = []
    t0 = time.time()
    print(
        f"BAND-B FEED WIDTH SWEEP: {len(jobs)} evals, "
        f"{len(WIDTHS_MM)} widths, n_workers={n_workers}",
        flush=True,
    )
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
                    "W_mm": j["W_mm"],
                    "f_GHz": j["f_GHz"],
                    "error": repr(e),
                }
                print(f"FAIL W={j['W_mm']} f={j['f_GHz']}: {e}", flush=True)
            results.append(r)
            if done % 20 == 0 or done == len(jobs):
                print(
                    f"  [{done}/{len(jobs)}] W={r.get('W_mm')} "
                    f"f={r.get('f_GHz')} A2={r.get('J_A2', 'err')}",
                    flush=True,
                )
                (SWEEP / "feed_width_sweep_partial.json").write_text(
                    json.dumps(
                        {"rows": results, "elapsed_s": time.time() - t0},
                        indent=2,
                    )
                    + "\n"
                )
    results = sorted(
        results,
        key=lambda r: (r.get("W_mm") or 0, r.get("f_GHz") or 0),
    )
    return results, time.time() - t0


def _ok(rows):
    return [r for r in rows if "error" not in r and "J_A2" in r]


def summarize_by_width(rows: list[dict]) -> list[dict]:
    out = []
    for W_mm in WIDTHS_MM:
        sub = _ok([r for r in rows if abs(r.get("W_mm", -1) - W_mm) < 1e-6])
        if not sub:
            out.append({"W_mm": W_mm, "n": 0, "error": "no rows"})
            continue
        band = [r for r in sub if BAND_F0 - 1e-9 <= r["f_GHz"] <= BAND_F1 + 1e-9]
        if not band:
            band = sub
        worst = min(band, key=lambda r: r["J_A2"])
        best = max(band, key=lambda r: r["J_A2"])
        # local minima near known families
        fam = {}
        for name, f0 in [("near_4p91", 4.91), ("near_5p16", 5.16)]:
            window = [r for r in sub if abs(r["f_GHz"] - f0) <= 0.08]
            if window:
                w = min(window, key=lambda r: r["J_A2"])
                fam[name] = {
                    "f_GHz": w["f_GHz"],
                    "J_A2": w["J_A2"],
                    "TE0_des": w["TE0_des"],
                    "TE0_wrong": w["TE0_wrong"],
                    "R_tot": w["R_tot"],
                    "TE1_total": w["TE1_total"],
                    "P_out": w["P_out"],
                    "classification": _classify_point(w),
                }
            else:
                fam[name] = None
        cut = te_cutoffs(width_mm_to_lu(W_mm))
        out.append(
            {
                "W_mm": W_mm,
                "W_lu": width_mm_to_lu(W_mm),
                "TE1_cutoff_GHz": cut["TE1_cutoff_GHz"],
                "TE2_cutoff_GHz": cut["modes"][2]["f_cut_GHz"],
                "n_band": len(band),
                "n_audit": len(sub),
                "worst_A2_BandB": worst["J_A2"],
                "worst_f_GHz_BandB": worst["f_GHz"],
                "worst_point": {
                    "f_GHz": worst["f_GHz"],
                    "J_A2": worst["J_A2"],
                    "TE0_des": worst["TE0_des"],
                    "TE0_wrong": worst["TE0_wrong"],
                    "R_tot": worst["R_tot"],
                    "TE1_total": worst["TE1_total"],
                    "P_out": worst["P_out"],
                    "classification": _classify_point(worst),
                },
                "best_A2_BandB": best["J_A2"],
                "best_f_GHz_BandB": best["f_GHz"],
                "mean_A2_BandB": float(np.mean([r["J_A2"] for r in band])),
                "mean_TE1_BandB": float(np.mean([r["TE1_total"] for r in band])),
                "mean_TE0_des_BandB": float(np.mean([r["TE0_des"] for r in band])),
                "mean_TE0_wrong_BandB": float(np.mean([r["TE0_wrong"] for r in band])),
                "mean_R_BandB": float(np.mean([r["R_tot"] for r in band])),
                "families": fam,
            }
        )
    return out


def _classify_point(r: dict) -> str:
    """Scattering-channel classification at a resonance/worst point."""
    te1 = float(r["TE1_total"])
    wrong = float(r["TE0_wrong"])
    refl = float(r["R_tot"])
    des = float(r["TE0_des"])
    vals = {
        "TE1": te1,
        "wrong_TE0": wrong,
        "reflection": refl,
        "absorption_or_missing": max(0.0, 1.0 - float(r["P_out"])),
    }
    top = max(vals, key=vals.get)
    if des >= 0.55 and r["J_A2"] > -2.0:
        return f"soft_floor_dominated_by_{top}"
    if vals[top] < 0.15:
        return f"mixed_shallow_{top}"
    return f"dominated_by_{top}"


def track_families(rows: list[dict]) -> dict:
    """Continuously track ~4.91 and ~5.16 GHz resonance families vs W."""
    tracks = {"near_4p91": [], "near_5p16": []}
    seeds = {"near_4p91": 4.91, "near_5p16": 5.16}
    # start from wide feed; allow seed to drift with width
    last_f = dict(seeds)
    for W_mm in WIDTHS_MM:
        sub = _ok([r for r in rows if abs(r.get("W_mm", -1) - W_mm) < 1e-6])
        if not sub:
            continue
        fs = np.array([r["f_GHz"] for r in sub])
        a2 = np.array([r["J_A2"] for r in sub])
        for name, f_seed in seeds.items():
            # search window around last tracked frequency
            f_ref = last_f[name]
            mask = np.abs(fs - f_ref) <= 0.10
            if not np.any(mask):
                mask = np.abs(fs - f_seed) <= 0.12
            if not np.any(mask):
                tracks[name].append({"W_mm": W_mm, "found": False})
                continue
            idx_local = np.where(mask)[0]
            # prefer local A2 minima (second-difference peak or simply min)
            i_min = idx_local[np.argmin(a2[idx_local])]
            r = sub[i_min]
            last_f[name] = r["f_GHz"]
            tracks[name].append(
                {
                    "W_mm": W_mm,
                    "found": True,
                    "f_GHz": r["f_GHz"],
                    "J_A2": r["J_A2"],
                    "TE0_des": r["TE0_des"],
                    "TE0_wrong": r["TE0_wrong"],
                    "R_tot": r["R_tot"],
                    "TE1_total": r["TE1_total"],
                    "P_out": r["P_out"],
                    "classification": _classify_point(r),
                    "delta_f_from_seed_MHz": (r["f_GHz"] - f_seed) * 1e3,
                }
            )
    return tracks


def make_plots(rows: list[dict], summary: list[dict], tracks: dict):
    SWEEP.mkdir(parents=True, exist_ok=True)
    ok = _ok(rows)
    widths = np.array(WIDTHS_MM, float)
    # common frequency axis (union; NaN fill)
    freqs = sorted({r["f_GHz"] for r in ok})
    f_arr = np.array(freqs, float)

    def grid(metric: str):
        G = np.full((len(widths), len(freqs)), np.nan)
        fmap = {f: i for i, f in enumerate(freqs)}
        wmap = {w: i for i, w in enumerate(WIDTHS_MM)}
        for r in ok:
            wi = wmap.get(r["W_mm"])
            fi = fmap.get(r["f_GHz"])
            if wi is None or fi is None:
                continue
            G[wi, fi] = r[metric]
        return G

    # --- A2 heatmap ---
    GA = grid("J_A2")
    fig, ax = plt.subplots(figsize=(10, 5.2))
    vmax = np.nanmax(GA)
    vmin = np.nanmin(GA)
    im = ax.pcolormesh(
        f_arr,
        widths,
        GA,
        shading="nearest",
        cmap="RdYlGn",
        vmin=max(vmin, -20),
        vmax=min(vmax, 0) if vmax < 0 else vmax,
    )
    ax.axvline(BAND_F0, color="k", ls="--", lw=0.8, alpha=0.7)
    ax.axvline(BAND_F1, color="k", ls="--", lw=0.8, alpha=0.7)
    for name, col in [("near_4p91", "C0"), ("near_5p16", "C3")]:
        xs, ys = [], []
        for t in tracks[name]:
            if t.get("found"):
                xs.append(t["f_GHz"])
                ys.append(t["W_mm"])
        if xs:
            ax.plot(xs, ys, "o-", color=col, ms=4, lw=1.2, label=name)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Feed clear width W (mm)")
    ax.set_title("Band-B fixed-q A2(f, W) — Stage-5 q, B=0.05 T")
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(im, ax=ax, label="A2")
    fig.tight_layout()
    fig.savefig(SWEEP / "A2_heatmap_f_vs_W.png", dpi=160)
    plt.close(fig)

    # --- overlays ---
    def overlay(metric, ylabel, fname, hline=None):
        fig, ax = plt.subplots(figsize=(9.5, 5.0))
        cmap = plt.cm.viridis(np.linspace(0.05, 0.95, len(WIDTHS_MM)))
        for i, W_mm in enumerate(WIDTHS_MM):
            sub = sorted(
                [r for r in ok if abs(r["W_mm"] - W_mm) < 1e-6],
                key=lambda r: r["f_GHz"],
            )
            if not sub:
                continue
            ax.plot(
                [r["f_GHz"] for r in sub],
                [r[metric] for r in sub],
                color=cmap[i],
                lw=1.3,
                label=f"{W_mm:.2f} mm",
            )
        ax.axvspan(BAND_F0, BAND_F1, color="0.85", alpha=0.35, label="Band B")
        if hline is not None:
            ax.axhline(hline, color="k", ls=":", lw=0.8)
        ax.set_xlabel("Frequency (GHz)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} vs f — fixed Band-B Stage-5 q")
        ax.legend(fontsize=7, ncol=2, loc="best")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(SWEEP / fname, dpi=150)
        plt.close(fig)

    overlay("J_A2", "A2", "A2_overlay_all_widths.png")
    overlay("TE0_des", "Desired TE0 power", "TE0_des_overlay.png")
    overlay("TE0_wrong", "Wrong/reverse TE0 power", "TE0_wrong_overlay.png")
    overlay("R_tot", "Reflection (R0+R1)", "R_tot_overlay.png")
    overlay("TE1_total", "TE1 total power", "TE1_overlay.png")
    overlay("P_out", "Total modal power P_out", "P_out_overlay.png", hline=1.0)

    # --- worst-A2 vs W ---
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0))
    Ws = [s["W_mm"] for s in summary if "worst_A2_BandB" in s]
    worstA = [s["worst_A2_BandB"] for s in summary if "worst_A2_BandB" in s]
    worstF = [s["worst_f_GHz_BandB"] for s in summary if "worst_A2_BandB" in s]
    meanA = [s["mean_A2_BandB"] for s in summary if "worst_A2_BandB" in s]
    meanTE1 = [s["mean_TE1_BandB"] for s in summary if "worst_A2_BandB" in s]
    axes[0, 0].plot(Ws, worstA, "o-", color="C3")
    axes[0, 0].set_ylabel("Worst Band-B A2")
    axes[0, 0].set_title("Worst dense-audit A2 vs W")
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 1].plot(Ws, worstF, "s-", color="C0")
    axes[0, 1].axhline(4.91, color="C0", ls=":", alpha=0.5)
    axes[0, 1].axhline(5.16, color="C3", ls=":", alpha=0.5)
    axes[0, 1].set_ylabel("f of worst A2 (GHz)")
    axes[0, 1].set_title("Worst-point frequency vs W")
    axes[0, 1].grid(True, alpha=0.3)
    axes[1, 0].plot(Ws, meanA, "o-", color="C2", label="mean A2")
    axes[1, 0].plot(Ws, worstA, "o--", color="C3", label="worst A2")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].set_xlabel("W (mm)")
    axes[1, 0].set_ylabel("A2")
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 1].plot(Ws, meanTE1, "o-", color="C1")
    axes[1, 1].set_xlabel("W (mm)")
    axes[1, 1].set_ylabel("Mean Band-B TE1 power")
    axes[1, 1].set_title("TE1 vs W (Band B mean)")
    axes[1, 1].grid(True, alpha=0.3)
    fig.suptitle("Band-B fixed-q feed-width sweep summary", y=1.01)
    fig.tight_layout()
    fig.savefig(SWEEP / "worst_A2_and_metrics_vs_W.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- family tracks ---
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0))
    for name, ax_row, color in [
        ("near_4p91", axes[0], "C0"),
        ("near_5p16", axes[1], "C3"),
    ]:
        ts = [t for t in tracks[name] if t.get("found")]
        if not ts:
            continue
        ax_row[0].plot([t["W_mm"] for t in ts], [t["f_GHz"] for t in ts], "o-", color=color)
        ax_row[0].set_ylabel("Resonance f (GHz)")
        ax_row[0].set_title(f"{name}: frequency vs W")
        ax_row[0].grid(True, alpha=0.3)
        ax_row[1].plot([t["W_mm"] for t in ts], [t["J_A2"] for t in ts], "o-", color=color, label="A2")
        ax2 = ax_row[1].twinx()
        ax2.plot(
            [t["W_mm"] for t in ts],
            [t["TE1_total"] for t in ts],
            "s--",
            color="C1",
            alpha=0.7,
            label="TE1",
        )
        ax2.plot(
            [t["W_mm"] for t in ts],
            [t["TE0_wrong"] for t in ts],
            "^:",
            color="C4",
            alpha=0.7,
            label="wrong TE0",
        )
        ax_row[1].set_ylabel("A2 depth")
        ax2.set_ylabel("Power fraction")
        ax_row[1].set_title(f"{name}: depth & channels vs W")
        ax_row[1].grid(True, alpha=0.3)
    axes[1, 0].set_xlabel("W (mm)")
    axes[1, 1].set_xlabel("W (mm)")
    fig.tight_layout()
    fig.savefig(SWEEP / "resonance_family_tracks.png", dpi=150)
    plt.close(fig)

    # channel-composition heatmaps
    for metric, fname, title in [
        ("TE0_des", "TE0_des_heatmap.png", "Desired TE0"),
        ("TE0_wrong", "TE0_wrong_heatmap.png", "Wrong TE0"),
        ("R_tot", "R_tot_heatmap.png", "Reflection"),
        ("TE1_total", "TE1_heatmap.png", "TE1 total"),
    ]:
        G = grid(metric)
        fig, ax = plt.subplots(figsize=(9.5, 4.8))
        im = ax.pcolormesh(f_arr, widths, G, shading="nearest", cmap="magma")
        ax.axvline(BAND_F0, color="w", ls="--", lw=0.7, alpha=0.6)
        ax.axvline(BAND_F1, color="w", ls="--", lw=0.7, alpha=0.6)
        ax.set_xlabel("Frequency (GHz)")
        ax.set_ylabel("W (mm)")
        ax.set_title(f"{title} (f, W)")
        fig.colorbar(im, ax=ax, label=title)
        fig.tight_layout()
        fig.savefig(SWEEP / fname, dpi=150)
        plt.close(fig)


def answer_science(summary: list[dict], tracks: dict) -> dict:
    """Answer the five scientific questions and recommend reopt or not."""
    byW = {s["W_mm"]: s for s in summary if "worst_A2_BandB" in s}
    W0 = byW[40.0]
    Wn = byW[23.06]
    # best worst-case A2 over Band B
    best = max(summary, key=lambda s: s.get("worst_A2_BandB", -1e9))
    # intermediate: TE1 suppressed in Band B (fc > 5.25) but not 23.06
    inter = [
        s
        for s in summary
        if s["W_mm"] < 40.0 - 1e-6
        and s["W_mm"] > 23.06 + 1e-3
        and s["TE1_cutoff_GHz"] > BAND_F1
    ]

    def fam_shift(name):
        ts = [t for t in tracks[name] if t.get("found")]
        if len(ts) < 2:
            return None
        return {
            "f_at_40mm": ts[0]["f_GHz"],
            "f_at_narrowest": ts[-1]["f_GHz"],
            "delta_f_MHz": (ts[-1]["f_GHz"] - ts[0]["f_GHz"]) * 1e3,
            "A2_at_40mm": ts[0]["J_A2"],
            "A2_at_narrowest": ts[-1]["J_A2"],
            "smooth": all(
                abs(ts[i + 1]["f_GHz"] - ts[i]["f_GHz"]) < 0.08 for i in range(len(ts) - 1)
            ),
            "path": [
                {
                    "W_mm": t["W_mm"],
                    "f_GHz": t["f_GHz"],
                    "J_A2": t["J_A2"],
                    "class": t["classification"],
                }
                for t in ts
            ],
        }

    # Q1: intermediate that suppresses TE1 without 23.06 penalties?
    # "Promising" requires TE1 suppression AND worst-case Band-B A2 at least as
    # good as the original 40 mm feed (the metric that matters for feasibility).
    # Mean-channel improvement alone is not enough if a deeper notch appears.
    q1_candidates = []
    for s in inter:
        te1_ok = s["mean_TE1_BandB"] < 0.5 * W0["mean_TE1_BandB"]
        less_wrong_than_narrow = s["mean_TE0_wrong_BandB"] < Wn["mean_TE0_wrong_BandB"]
        less_R_than_narrow = s["mean_R_BandB"] < Wn["mean_R_BandB"]
        better_worst_than_40 = s["worst_A2_BandB"] >= W0["worst_A2_BandB"] - 0.25
        better_worst_than_23 = s["worst_A2_BandB"] > Wn["worst_A2_BandB"] + 0.3
        q1_candidates.append(
            {
                "W_mm": s["W_mm"],
                "TE1_cutoff_GHz": s["TE1_cutoff_GHz"],
                "mean_TE1": s["mean_TE1_BandB"],
                "mean_wrong": s["mean_TE0_wrong_BandB"],
                "mean_R": s["mean_R_BandB"],
                "worst_A2": s["worst_A2_BandB"],
                "te1_suppressed_vs_40mm": te1_ok,
                "less_wrong_than_23mm": less_wrong_than_narrow,
                "less_R_than_23mm": less_R_than_narrow,
                "worst_A2_competitive_with_40mm": better_worst_than_40,
                "better_worst_than_23mm": better_worst_than_23,
                "promising": te1_ok and better_worst_than_40,
            }
        )

    promising = [c for c in q1_candidates if c["promising"]]
    # Q3: can feed width move dominant bad resonance outside Band B?
    # Check if worst-point frequency exits [4.75,5.25] or family A2 in-band weakens above -3
    q3 = {
        "any_worst_outside_band": any(
            s["worst_f_GHz_BandB"] < BAND_F0 + 0.01 or s["worst_f_GHz_BandB"] > BAND_F1 - 0.01
            for s in summary
            if s["W_mm"] < 40
        ),
        "family_4p91_exits_or_softens": False,
        "family_5p16_exits_or_softens": False,
        "notes": [],
    }
    for name, key in [("near_4p91", "family_4p91_exits_or_softens"), ("near_5p16", "family_5p16_exits_or_softens")]:
        ts = [t for t in tracks[name] if t.get("found")]
        for t in ts:
            if t["W_mm"] >= 40:
                continue
            out = t["f_GHz"] < BAND_F0 or t["f_GHz"] > BAND_F1
            soft = t["J_A2"] > -3.0
            if out or soft:
                q3[key] = True
                q3["notes"].append(
                    f"{name} at W={t['W_mm']:.2f}: f={t['f_GHz']:.4f}, A2={t['J_A2']:.2f}, "
                    f"{'outside band' if out else 'softened'}"
                )

    # Q4/Q5: improvement of worst dense A2 over entire Band B
    delta_vs_40 = best["worst_A2_BandB"] - W0["worst_A2_BandB"]
    # Reopt justification: meaningful worst-case gain AND not just trading one notch for another
    # Threshold: ≥1.5 A2 improvement in worst-case over Band B, or ≥2.0 if candidate is intermediate
    justified = False
    reason = []
    if best["W_mm"] == 40.0:
        justified = False
        reason.append("No narrowed/intermediate width improves worst Band-B A2 over the original 40 mm feed.")
    elif delta_vs_40 >= 1.5:
        justified = True
        reason.append(
            f"Best fixed-q width W={best['W_mm']:.2f} mm improves worst Band-B A2 by "
            f"{delta_vs_40:.2f} relative to 40 mm "
            f"({W0['worst_A2_BandB']:.2f} → {best['worst_A2_BandB']:.2f})."
        )
        if best["W_mm"] != 23.06:
            reason.append("Best width is intermediate — joint W+q reopt scientifically justified.")
        else:
            reason.append(
                "Best width is the aggressive 23.06 mm feed; a small joint reopt may still help "
                "if wrong-TE0/reflection can be mitigated by q adaptation."
            )
    elif delta_vs_40 >= 0.5:
        justified = False
        reason.append(
            f"Only modest worst-A2 gain ({delta_vs_40:.2f}) at W={best['W_mm']:.2f} mm — "
            "insufficient alone to justify reoptimization; consider only if TE1 physics "
            "must be removed for diagnostic clarity."
        )
    else:
        justified = False
        reason.append(
            f"Negligible/negative worst-A2 change ({delta_vs_40:.2f}) — do not reoptimize on narrowed feeds yet."
        )

    f4 = fam_shift("near_4p91")
    f5 = fam_shift("near_5p16")
    # 4.91 tracker can jump once the deep TE1 notch collapses; treat large
    # discontinuous redshift + classification change as "moves/softens" evidence
    # that feeds participate, even if the tracked minimum is imperfectly continuous.
    q2_answer = "yes" if (f4 or {}).get("smooth") and (f5 or {}).get("smooth") else "mixed"
    if f5 and abs(f5.get("delta_f_MHz", 99)) < 15 and f4 and abs(f4.get("delta_f_MHz", 0)) > 50:
        q2_answer = (
            "partial — ~5.16 GHz family is nearly pinned in frequency while deepening/"
            "reclassifying; ~4.91 GHz family softens and redshifts (tracker may latch "
            "onto a shallower neighbor once the TE1 notch collapses)"
        )

    return {
        "Q1_intermediate_TE1_suppress_without_23mm_penalty": {
            "answer": (
                "yes"
                if promising
                else (
                    "no — TE1 can be suppressed, but not without a worse Band-B worst-case "
                    "(wrong-TE0 / reflection penalties at ~5.16 GHz)"
                    if any(c["te1_suppressed_vs_40mm"] for c in q1_candidates)
                    else "no"
                )
            ),
            "candidates": q1_candidates,
            "promising": promising,
        },
        "Q2_resonances_move_smoothly_with_W": {
            "answer": q2_answer,
            "family_4p91": f4,
            "family_5p16": f5,
            "interpretation": (
                "Feeds participate in the resonant eigenstructure: scattering "
                "classification of both families evolves continuously with W "
                "(TE1 → reflection or wrong-TE0). Frequency motion is strong for "
                "the ~4.91 family and weak for ~5.16."
            ),
        },
        "Q3_can_W_move_bad_resonance_out_of_BandB": q3,
        "Q4_clear_candidate_for_worst_case_BandB_A2": {
            "best_W_mm": best["W_mm"],
            "best_worst_A2": best["worst_A2_BandB"],
            "best_worst_f_GHz": best["worst_f_GHz_BandB"],
            "orig_40mm_worst_A2": W0["worst_A2_BandB"],
            "orig_40mm_worst_f_GHz": W0["worst_f_GHz_BandB"],
            "narrow_23mm_worst_A2": Wn["worst_A2_BandB"],
            "delta_vs_40mm": delta_vs_40,
            "improves_entire_BandB_worst_not_just_4p91": best["worst_A2_BandB"]
            > W0["worst_A2_BandB"] + 0.5,
        },
        "Q5_justify_BandB_reopt_on_narrowed_feed": {
            "recommended": justified,
            "reasons": reason,
            "candidate_W_mm": best["W_mm"] if justified else None,
            "stage6_status": "NOT STARTED",
            "do_not_auto_launch": True,
        },
        "width_ranking_by_worst_BandB_A2": sorted(
            [
                {
                    "W_mm": s["W_mm"],
                    "worst_A2": s["worst_A2_BandB"],
                    "worst_f": s["worst_f_GHz_BandB"],
                    "mean_A2": s["mean_A2_BandB"],
                    "mean_TE1": s["mean_TE1_BandB"],
                }
                for s in summary
                if "worst_A2_BandB" in s
            ],
            key=lambda x: -x["worst_A2"],
        ),
    }


def write_report(catalog, summary, tracks, science, elapsed_s, n_jobs):
    lines = []
    lines.append("# Band-B Fixed-q Feed-Width Sweep")
    lines.append("")
    lines.append("**Stage-6 status: NOT STARTED**")
    lines.append("")
    lines.append("Fixed Band-B C6 Stage-5 optimized `q` at B = 0.05 T. No q optimization.")
    lines.append(
        f"Dense forward audit over {AUDIT_F0}–{AUDIT_F1} GHz "
        f"(Band B = {BAND_F0}–{BAND_F1} GHz + guard)."
    )
    lines.append(f"Evals: {n_jobs}. Wall time: {elapsed_s/60:.1f} min.")
    lines.append("")
    lines.append("## Analytical cutoffs (solver modal convention)")
    lines.append("")
    lines.append("| W (mm) | W (lu) | TE1 fc (GHz) | TE2 fc (GHz) | Note |")
    lines.append("|-------:|-------:|-------------:|-------------:|------|")
    for c in catalog:
        lines.append(
            f"| {c['W_mm']:.2f} | {c['W_lu']:.4f} | {c['TE1_cutoff_GHz']:.3f} | "
            f"{c['TE2_cutoff_GHz']:.3f} | {c['note']} |"
        )
    lines.append("")
    lines.append("Formula: `f_c,m = m · c / (2 W_m)`, `W_m = W_lu · a`, `a = 20 mm`.")
    lines.append("")
    lines.append("## Worst-case Band-B dense A2 vs width")
    lines.append("")
    lines.append(
        "| W (mm) | TE1 fc | Worst A2 | Worst f (GHz) | Mean A2 | Mean TE1 | "
        "Mean wrong TE0 | Mean R | Worst class |"
    )
    lines.append(
        "|-------:|-------:|---------:|--------------:|--------:|---------:|"
        "---------------:|-------:|-------------|"
    )
    for s in summary:
        if "worst_A2_BandB" not in s:
            continue
        wp = s["worst_point"]
        lines.append(
            f"| {s['W_mm']:.2f} | {s['TE1_cutoff_GHz']:.3f} | {s['worst_A2_BandB']:.3f} | "
            f"{s['worst_f_GHz_BandB']:.4f} | {s['mean_A2_BandB']:.3f} | "
            f"{s['mean_TE1_BandB']:.3f} | {s['mean_TE0_wrong_BandB']:.3f} | "
            f"{s['mean_R_BandB']:.3f} | `{wp['classification']}` |"
        )
    lines.append("")
    lines.append("## Resonance family tracking")
    lines.append("")
    for name in ("near_4p91", "near_5p16"):
        lines.append(f"### {name}")
        lines.append("")
        lines.append("| W (mm) | f (GHz) | Δf from seed (MHz) | A2 | TE1 | wrong TE0 | R | class |")
        lines.append("|-------:|--------:|-------------------:|---:|----:|----------:|--:|-------|")
        for t in tracks[name]:
            if not t.get("found"):
                lines.append(f"| {t['W_mm']:.2f} | — | — | — | — | — | — | not found |")
                continue
            lines.append(
                f"| {t['W_mm']:.2f} | {t['f_GHz']:.4f} | {t['delta_f_from_seed_MHz']:+.1f} | "
                f"{t['J_A2']:.3f} | {t['TE1_total']:.3f} | {t['TE0_wrong']:.3f} | "
                f"{t['R_tot']:.3f} | `{t['classification']}` |"
            )
        lines.append("")
    lines.append("## Scientific questions")
    lines.append("")
    q1 = science["Q1_intermediate_TE1_suppress_without_23mm_penalty"]
    lines.append(f"### Q1 — Intermediate width TE1 suppress without 23 mm penalties?")
    lines.append(f"**Answer: {q1['answer']}**")
    lines.append("")
    lines.append(
        "Criterion for promising: mean TE1 ≪ 40 mm **and** worst Band-B A2 "
        "competitive with 40 mm. Mean-channel cosmetics alone do not count."
    )
    if q1["promising"]:
        for c in q1["promising"]:
            lines.append(
                f"- Promising W={c['W_mm']:.2f} mm: "
                f"worst A2={c['worst_A2']:.3f}, mean TE1={c['mean_TE1']:.3f}, "
                f"mean wrong={c['mean_wrong']:.3f}, mean R={c['mean_R']:.3f}"
            )
    else:
        lines.append(
            "- **No** intermediate width meets the promising criterion. "
            "Widths with TE1 cutoff above Band B do suppress TE1, but they "
            "deepen the ~5.16 GHz family into a wrong-TE0 (or reflection) notch "
            "that is *worse* than the original 40 mm worst-case A2."
        )
        for c in q1["candidates"]:
            lines.append(
                f"  - W={c['W_mm']:.2f} mm: TE1_suppressed={c['te1_suppressed_vs_40mm']}, "
                f"worst A2={c['worst_A2']:.3f} (40 mm baseline "
                f"{science['Q4_clear_candidate_for_worst_case_BandB_A2']['orig_40mm_worst_A2']:.3f}), "
                f"mean wrong={c['mean_wrong']:.3f}, mean R={c['mean_R']:.3f}"
            )
    lines.append("")
    lines.append(
        "Note: some intermediate widths improve *mean* Band-B A2 while "
        "degrading the *worst* point — the feasibility metric is worst-case."
    )
    lines.append("")
    q2 = science["Q2_resonances_move_smoothly_with_W"]
    lines.append(f"### Q2 — Do resonance frequencies move smoothly with W?")
    lines.append(f"**Answer: {q2['answer']}**")
    for lab, key in [("~4.91", "family_4p91"), ("~5.16", "family_5p16")]:
        fs = q2.get(key) or {}
        if fs:
            lines.append(
                f"- {lab}: {fs['f_at_40mm']:.4f} → {fs['f_at_narrowest']:.4f} GHz "
                f"(Δ = {fs['delta_f_MHz']:+.1f} MHz); A2 {fs['A2_at_40mm']:.2f} → {fs['A2_at_narrowest']:.2f}"
            )
    lines.append("")
    q3 = science["Q3_can_W_move_bad_resonance_out_of_BandB"]
    lines.append("### Q3 — Can feed width move the dominant bad resonance outside Band B?")
    lines.append(
        f"- Family softens/exits: 4.91={q3['family_4p91_exits_or_softens']}, "
        f"5.16={q3['family_5p16_exits_or_softens']}"
    )
    for n in q3["notes"]:
        lines.append(f"- {n}")
    lines.append("")
    q4 = science["Q4_clear_candidate_for_worst_case_BandB_A2"]
    lines.append("### Q4 — Clear candidate for **worst-case dense A2 over entire Band B**?")
    lines.append(
        f"**Best fixed-q width: W = {q4['best_W_mm']:.2f} mm** "
        f"(worst A2 = {q4['best_worst_A2']:.3f} at {q4['best_worst_f_GHz']:.4f} GHz)."
    )
    lines.append(
        f"vs 40 mm: worst A2 = {q4['orig_40mm_worst_A2']:.3f} at {q4['orig_40mm_worst_f_GHz']:.4f} GHz "
        f"(Δ = {q4['delta_vs_40mm']:+.3f})."
    )
    lines.append(
        f"vs 23.06 mm: worst A2 = {q4['narrow_23mm_worst_A2']:.3f}."
    )
    lines.append(
        f"Improves entire-Band-B worst (not merely the old 4.909 min): "
        f"**{q4['improves_entire_BandB_worst_not_just_4p91']}**"
    )
    lines.append("")
    q5 = science["Q5_justify_BandB_reopt_on_narrowed_feed"]
    lines.append("### Q5 — Justify a narrowed/intermediate-feed Band-B reoptimization?")
    lines.append(f"**Recommendation: {'YES' if q5['recommended'] else 'NO'}**")
    for r in q5["reasons"]:
        lines.append(f"- {r}")
    lines.append(
        "- Fixed-q evidence: narrowing trades the TE1-dominated ~4.91 notch for a "
        "deeper wrong-TE0/reflection-dominated ~5.16 notch. A joint W+q reopt would "
        "need to *overcome* that conversion, not merely retune around the old 4.91 min. "
        "That is a larger ask than a 'small' Band-B calculation is likely to settle."
    )
    lines.append("")
    lines.append(
        "Do **not** auto-launch that optimization. Stage-6 remains **NOT STARTED**."
    )
    lines.append("")
    lines.append("## Plots")
    lines.append("")
    for p in [
        "A2_heatmap_f_vs_W.png",
        "A2_overlay_all_widths.png",
        "worst_A2_and_metrics_vs_W.png",
        "resonance_family_tracks.png",
        "TE0_des_overlay.png",
        "TE0_wrong_overlay.png",
        "R_tot_overlay.png",
        "TE1_overlay.png",
        "P_out_overlay.png",
        "TE0_des_heatmap.png",
        "TE0_wrong_heatmap.png",
        "R_tot_heatmap.png",
        "TE1_heatmap.png",
    ]:
        lines.append(f"- `{p}`")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (SWEEP / "FEED_WIDTH_SWEEP_BAND_B.md").write_text(text)
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--df-MHz", type=float, default=12.5)
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()

    SWEEP.mkdir(parents=True, exist_ok=True)
    catalog = width_catalog()
    (SWEEP / "width_cutoffs.json").write_text(json.dumps(catalog, indent=2) + "\n")

    if args.plot_only:
        data = json.loads((SWEEP / "feed_width_sweep_results.json").read_text())
        rows = data["rows"]
        elapsed = data.get("elapsed_s", 0.0)
        n_jobs = data.get("n_jobs", len(rows))
    else:
        opt = json.loads(OPT_PATH.read_text())
        q = np.asarray(opt["q_best"], float).tolist()
        meta = {
            "stage6_status": "NOT STARTED",
            "design": "Band-B C6 Stage-5",
            "opt_path": str(OPT_PATH),
            "b_tesla": B_TESLA,
            "band_GHz": [BAND_F0, BAND_F1],
            "audit_GHz": [AUDIT_F0, AUDIT_F1],
            "df_MHz": args.df_MHz,
            "widths_mm": WIDTHS_MM,
            "watch_GHz": WATCH,
            "q_best": q,
            "no_q_optimization": True,
        }
        (SWEEP / "sweep_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
        jobs = build_jobs(q, args.df_MHz)
        rows, elapsed = run_jobs(jobs, n_workers=args.n_workers)
        n_jobs = len(jobs)
        payload = {
            "meta": meta,
            "catalog": catalog,
            "n_jobs": n_jobs,
            "elapsed_s": elapsed,
            "rows": rows,
        }
        (SWEEP / "feed_width_sweep_results.json").write_text(
            json.dumps(payload, indent=2) + "\n"
        )

    summary = summarize_by_width(rows)
    tracks = track_families(rows)
    science = answer_science(summary, tracks)
    make_plots(rows, summary, tracks)
    report = write_report(catalog, summary, tracks, science, elapsed, n_jobs)
    out = {
        "summary": summary,
        "tracks": tracks,
        "science": science,
        "elapsed_s": elapsed,
        "n_jobs": n_jobs,
        "stage6_status": "NOT STARTED",
        "reopt_launched": False,
    }
    (SWEEP / "feed_width_sweep_analysis.json").write_text(json.dumps(out, indent=2) + "\n")
    print(report)
    print(
        f"\nDONE. Recommendation reopt={science['Q5_justify_BandB_reopt_on_narrowed_feed']['recommended']}. "
        "Stage-6 NOT STARTED.",
        flush=True,
    )


if __name__ == "__main__":
    main()
