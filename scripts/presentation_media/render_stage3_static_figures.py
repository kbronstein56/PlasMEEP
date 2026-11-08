#!/usr/bin/env python3
"""Stage-3 presentation STATIC figures (PNG) — no FEM solves.

Reads production audit / opt JSON from PlasMEEP campaign outputs.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from media_common import (  # noqa: E402
    C_ACCENT,
    C_BAD,
    C_GOOD,
    C_INK,
    C_MUTED,
    C_WARN,
    OUT_ROOT,
    PROD_CAMPAIGN,
    stamp_kind,
)

OUT = OUT_ROOT / "stage3_figures"


def load(name: str):
    return json.loads((PROD_CAMPAIGN / name).read_text())


def rows_of(d):
    return d["rows"] if isinstance(d, dict) else d


def xy(rows, key):
    f = np.array([r["f_GHz"] for r in rows], float)
    y = np.array([r.get(key) if r.get(key) is not None else np.nan for r in rows], float)
    return f, y


def mark_freqs(ax, freqs, color, label=None, first=True):
    for i, f in enumerate(freqs):
        ax.axvline(f, color=color, ls="--", lw=1.6, alpha=0.85, label=label if first and i == 0 else None)


def fig_s2_vs_s3_broadband():
    """Primary Stage-2 vs Stage-3 spectral comparison."""
    s2 = load("bb_stage2_audit_50MHz.json")
    # Prefer full dense 25 if present, else 50 + local 25
    dense = PROD_CAMPAIGN / "bb_stage3_dense_25MHz.json"
    if dense.exists():
        s3 = load("bb_stage3_dense_25MHz.json")
        s3_label = "Stage-3 dense 25 MHz"
    else:
        s3 = load("bb_stage3_audit_50MHz.json")
        s3_label = "Stage-3 audit 50 MHz"

    f2, a2 = xy(rows_of(s2), "J_A2")
    f3, a3 = xy(rows_of(s3), "J_A2")
    _, il2 = xy(rows_of(s2), "IL_dB")
    _, il3 = xy(rows_of(s3), "IL_dB")
    _, iso2 = xy(rows_of(s2), "isolation_dB")
    _, iso3 = xy(rows_of(s3), "isolation_dB")

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(
        "Stage-2 → Stage-3 broadband improvement (A2_hybrid)",
        fontsize=28,
        fontweight="bold",
        color=C_INK,
        y=0.97,
    )
    ax0 = fig.add_axes([0.07, 0.55, 0.88, 0.35])
    ax1 = fig.add_axes([0.07, 0.10, 0.42, 0.35])
    ax2 = fig.add_axes([0.55, 0.10, 0.40, 0.35])
    stamp_kind(ax0, "DATA", "OPTIMIZATION AUDIT")
    stamp_kind(ax1, "DATA")
    stamp_kind(ax2, "DATA")

    ax0.plot(f2, a2, color=C_MUTED, lw=2.2, label="Stage-2 (50 MHz)")
    ax0.plot(f3, a3, color=C_ACCENT, lw=2.6, label=s3_label)
    refine = PROD_CAMPAIGN / "bb_stage3_dense_12p5MHz_local.json"
    if refine.exists():
        fr, ar = xy(rows_of(load("bb_stage3_dense_12p5MHz_local.json")), "J_A2")
        # show only the sharpest refine failures as markers
        for fmark, col, lab in (
            (5.8125, C_BAD, "12.5 MHz hole 5.8125"),
            (4.775, C_WARN, "25 MHz hole 4.775"),
        ):
            j = int(np.argmin(np.abs(fr - fmark)))
            ax0.scatter(
                [fr[j]],
                [ar[j]],
                s=90,
                zorder=6,
                color=col,
                edgecolors="k",
                linewidths=0.6,
                label=lab,
            )
    mark_freqs(ax0, [5.75, 5.825], C_GOOD, "repaired nodes 5.75 / 5.825")
    mark_freqs(ax0, [5.8125], C_BAD, first=False)
    mark_freqs(ax0, [4.775], C_WARN, first=False)
    ax0.axhline(0, color="#94a3b8", lw=1)
    ax0.set_ylabel("A2 objective (higher better)")
    ax0.set_xlabel("Frequency (GHz)")
    ax0.legend(loc="lower right", fontsize=13, ncol=2)
    ax0.set_xlim(4.0, 6.0)

    ax1.plot(f2, il2, color=C_MUTED, lw=2.0, label="Stage-2")
    ax1.plot(f3, il3, color=C_ACCENT, lw=2.4, label="Stage-3")
    mark_freqs(ax1, [5.75, 5.825], C_GOOD)
    mark_freqs(ax1, [5.8125, 4.775], C_BAD)
    ax1.set_ylabel("IL (dB)")
    ax1.set_xlabel("Frequency (GHz)")
    ax1.legend(fontsize=13)
    ax1.set_xlim(4.0, 6.0)

    ax2.plot(f2, iso2, color=C_MUTED, lw=2.0, label="Stage-2")
    ax2.plot(f3, iso3, color=C_ACCENT, lw=2.4, label="Stage-3")
    mark_freqs(ax2, [5.75, 5.825], C_GOOD)
    mark_freqs(ax2, [5.8125, 4.775], C_BAD)
    ax2.axhline(0, color="#94a3b8", lw=1)
    ax2.set_ylabel("Isolation (dB)")
    ax2.set_xlabel("Frequency (GHz)")
    ax2.legend(fontsize=13)
    ax2.set_xlim(4.0, 6.0)

    out = OUT / "13_stage2_vs_stage3_broadband.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def _bar_pair(ax, labels, s2, s3, ylabel, title):
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, s2, w, color=C_MUTED, label="Stage-2")
    ax.bar(x + w / 2, s3, w, color=C_ACCENT, label="Stage-3")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=16)
    ax.legend(fontsize=12)
    stamp_kind(ax, "DATA")


def fig_hole_repair(f_ghz: float, stem: str, title: str):
    # Prefer exact 25 MHz audits
    s2 = load("bb_stage2_audit_25MHz_local.json")
    s3_path = PROD_CAMPAIGN / "bb_stage3_dense_25MHz.json"
    s3 = load("bb_stage3_dense_25MHz.json") if s3_path.exists() else load(
        "bb_stage3_audit_25MHz_local.json"
    )

    def at(rows, f):
        return min(rows, key=lambda r: abs(r["f_GHz"] - f))

    r2, r3 = at(rows_of(s2), f_ghz), at(rows_of(s3), f_ghz)
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(title, fontsize=26, fontweight="bold", color=C_INK, y=0.96)
    ax0 = fig.add_axes([0.08, 0.18, 0.38, 0.68])
    ax1 = fig.add_axes([0.55, 0.18, 0.38, 0.68])
    _bar_pair(
        ax0,
        ["A2", "iso (dB)", "IL (dB)"],
        [r2["J_A2"], r2["isolation_dB"], r2["IL_dB"]],
        [r3["J_A2"], r3["isolation_dB"], r3["IL_dB"]],
        "Metric value",
        f"f = {f_ghz:.3f} GHz  (exact audit nodes)",
    )
    ax0.axhline(0, color="#94a3b8", lw=1)
    # Annotation table
    ax1.axis("off")
    stamp_kind(ax1, "DATA", "QUANTITATIVE")
    txt = (
        f"Frequency: {f_ghz:.3f} GHz\n\n"
        f"Stage-2:\n"
        f"  A2 = {r2['J_A2']:.3f}\n"
        f"  isolation = {r2['isolation_dB']:.2f} dB\n"
        f"  IL = {r2['IL_dB']:.2f} dB\n"
        f"  T_des = {r2.get('T_des', float('nan')):.4f}\n\n"
        f"Stage-3:\n"
        f"  A2 = {r3['J_A2']:.3f}\n"
        f"  isolation = {r3['isolation_dB']:.2f} dB\n"
        f"  IL = {r3['IL_dB']:.2f} dB\n"
        f"  T_des = {r3.get('T_des', float('nan')):.4f}\n\n"
        f"ΔA2 = {r3['J_A2'] - r2['J_A2']:+.3f}\n"
        f"Δiso = {r3['isolation_dB'] - r2['isolation_dB']:+.2f} dB"
    )
    ax1.text(
        0.05,
        0.95,
        txt,
        transform=ax1.transAxes,
        va="top",
        ha="left",
        fontsize=18,
        family="monospace",
        color=C_INK,
        bbox=dict(boxstyle="round,pad=0.6", fc="#f8fafc", ec="#cbd5e1"),
    )
    out = OUT / f"{stem}.png"
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def fig_4p1_improvement():
    s2 = load("bb_stage2_audit_50MHz.json")
    s3 = load("bb_stage3_audit_50MHz.json")

    def at(rows, f):
        return min(rows, key=lambda r: abs(r["f_GHz"] - f))

    r2, r3 = at(rows_of(s2), 4.1), at(rows_of(s3), 4.1)
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(
        "4.1 GHz recovery preserved and improved (Stage-2 → Stage-3)",
        fontsize=26,
        fontweight="bold",
        color=C_INK,
        y=0.96,
    )
    ax = fig.add_axes([0.12, 0.18, 0.76, 0.68])
    labels = ["T_des", "IL (dB)", "iso (dB)", "A2"]
    # Scale T_des×10 for visibility on shared axis? Better twin or separate.
    # Use grouped normalized annotation bars for IL/iso/A2 and text for T_des.
    x = np.arange(3)
    w = 0.36
    ax.bar(x - w / 2, [r2["IL_dB"], r2["isolation_dB"], r2["J_A2"]], w, color=C_MUTED, label="Stage-2")
    ax.bar(x + w / 2, [r3["IL_dB"], r3["isolation_dB"], r3["J_A2"]], w, color=C_GOOD, label="Stage-3")
    ax.set_xticks(x)
    ax.set_xticklabels(["IL (dB)", "isolation (dB)", "A2"])
    ax.axhline(0, color="#94a3b8", lw=1)
    ax.legend(fontsize=14)
    stamp_kind(ax, "DATA")
    ax.text(
        0.02,
        0.98,
        f"T_des: {r2['T_des']:.3f} → {r3['T_des']:.3f}\n"
        f"IL: {r2['IL_dB']:.2f} → {r3['IL_dB']:.2f} dB",
        transform=ax.transAxes,
        va="top",
        fontsize=18,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cbd5e1"),
    )
    out = OUT / "16_stage3_4p1_improvement.png"
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def fig_device_summary():
    dense = PROD_CAMPAIGN / "bb_stage3_dense_25MHz.json"
    s3 = load("bb_stage3_dense_25MHz.json") if dense.exists() else load(
        "bb_stage3_audit_50MHz.json"
    )
    opt = load("opt_broadband_stage3_A2_C.json")
    f, a2 = xy(rows_of(s3), "J_A2")
    _, iso = xy(rows_of(s3), "isolation_dB")

    # Geometry scatter (same helpers as Stage-2 progress renderer)
    sys.path[:0] = [
        str(Path("/home/daq_user/PlasMEEP-viz/scripts/circulator_inverse_design")),
        str(Path("/home/daq_user/PlasMEEP-viz/scripts/validation/fem_full_device_validation")),
        str(Path("/home/daq_user/PlasMEEP-viz/scripts/validation")),
        str(Path("/home/daq_user/PlasMEEP-viz/scripts")),
    ]
    import sixport_common as sc
    from full91_geometry import bulb_centers_device
    from sixfold_orbits import build_orbits, expand_tied

    centers = np.asarray(bulb_centers_device(sc), float)
    orbits, _origin = build_orbits(centers)
    q = np.asarray(opt["q_best"], float)
    s = expand_tied(q, orbits)

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(
        "Stage-3 broadband plasma circulator — inverse-design summary",
        fontsize=26,
        fontweight="bold",
        color=C_INK,
        y=0.97,
    )
    ax0 = fig.add_axes([0.04, 0.12, 0.42, 0.78])
    ax1 = fig.add_axes([0.52, 0.55, 0.44, 0.35])
    ax2 = fig.add_axes([0.52, 0.12, 0.44, 0.35])
    stamp_kind(ax0, "DATA", "Stage-3 q_best")
    stamp_kind(ax1, "DATA")
    stamp_kind(ax2, "DATA")

    xy0 = centers - centers.mean(0)
    sca = ax0.scatter(
        xy0[:, 0],
        xy0[:, 1],
        c=s,
        s=70,
        cmap="viridis",
        vmin=0.4,
        vmax=1.45,
        edgecolors="k",
        linewidths=0.3,
    )
    # Port markers at guide angles 0°,60°,… (outside array)
    R = float(np.max(np.linalg.norm(xy0, axis=1))) * 1.28
    for i, ang_deg in enumerate(range(0, 360, 60)):
        ang = np.deg2rad(ang_deg)
        ax0.text(
            R * np.cos(ang),
            R * np.sin(ang),
            f"P{i+1}",
            ha="center",
            va="center",
            fontsize=15,
            fontweight="bold",
            color=C_ACCENT,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cbd5e1", alpha=0.92),
            zorder=6,
        )
    ax0.set_xlim(-1.45 * R / 1.28, 1.45 * R / 1.28)
    ax0.set_ylim(-1.45 * R / 1.28, 1.45 * R / 1.28)
    ax0.annotate(
        "",
        xy=(0.55 * R * np.cos(np.deg2rad(60)), 0.55 * R * np.sin(np.deg2rad(60))),
        xytext=(0.55 * R, 0),
        arrowprops=dict(arrowstyle="->", color=C_GOOD, lw=2.5),
    )
    ax0.text(
        0.02,
        0.02,
        "Desired circulation (B+): P1→P2→…",
        transform=ax0.transAxes,
        fontsize=13,
        color=C_GOOD,
    )
    ax0.set_aspect("equal")
    ax0.set_title("91-element C6 magnetized plasma array", fontsize=16)
    ax0.set_xlabel("x (a)")
    ax0.set_ylabel("y (a)")
    fig.colorbar(sca, ax=ax0, fraction=0.046, pad=0.04, label="relative density")

    # Prefer refine overlay markers even if plotted series is 25 MHz
    ax1.plot(f, a2, color=C_ACCENT, lw=2.4)
    mark_freqs(ax1, [5.75, 5.825], C_GOOD, "repaired nodes")
    mark_freqs(ax1, [5.8125], C_BAD, "hidden 5.8125", first=True)
    mark_freqs(ax1, [4.775], C_WARN, "hidden 4.775", first=True)
    ax1.set_ylabel("A2")
    ax1.set_xlim(4, 6)
    ax1.set_title(
        f"J_BB = {opt.get('J_BB_best'):.3f} · 30-freq soft-min A2",
        fontsize=15,
    )

    ax2.plot(f, iso, color=C_ACCENT, lw=2.4)
    ax2.axhline(0, color="#94a3b8", lw=1)
    mark_freqs(ax2, [5.75, 5.825], C_GOOD)
    mark_freqs(ax2, [5.8125, 4.775], C_BAD)
    ax2.set_ylabel("Isolation (dB)")
    ax2.set_xlabel("Frequency (GHz)")
    ax2.set_xlim(4, 6)

    out = OUT / "17_stage3_device_summary.png"
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def fig_optimization_evolution():
    """Stage-1 → Stage-2 → Stage-3 J_BB progression from histories."""
    s1 = load("opt_broadband_production_C.json")
    s2 = load("opt_broadband_stage2_A2_C.json")
    s3 = load("opt_broadband_stage3_A2_C.json")
    h2 = load("bb_stage2_history.json")["history"]
    h3 = load("bb_stage3_history.json")["history"]

    def accepted(h):
        xs, ys = [], []
        best = None
        for r in h:
            if r.get("kind") in ("start_eval", "accepted_iter"):
                j = float(r["J_BB"])
                best = j if best is None else max(best, j)
                xs.append(int(r["iteration"]))
                ys.append(best)
        return np.asarray(xs), np.asarray(ys)

    x2, y2 = accepted(h2)
    x3, y3 = accepted(h3)

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(
        "Broadband inverse-design evolution across stages",
        fontsize=26,
        fontweight="bold",
        color=C_INK,
        y=0.96,
    )
    ax = fig.add_axes([0.08, 0.12, 0.84, 0.75])
    stamp_kind(ax, "DATA", "OPTIMIZATION HISTORY")
    ax.plot(x2, y2, color=C_MUTED, lw=2.4, label="Stage-2 A2 (28 freq)")
    ax.plot(x3, y3, color=C_ACCENT, lw=2.8, label="Stage-3 A2 (30 freq)")
    ax.scatter(
        [0],
        [float(s1.get("J_BB_best", s1.get("J_BB", np.nan)))],
        s=120,
        color=C_WARN,
        zorder=5,
        label=f"Stage-1 final J≈{s1.get('J_BB_best', s1.get('J_BB'))}",
    )
    ax.set_xlabel("Accepted L-BFGS iteration")
    ax.set_ylabel("Best-so-far J_BB (soft-min A2)")
    ax.legend(fontsize=14, loc="lower right")
    ax.text(
        0.02,
        0.98,
        f"Stage-2 final {s2.get('J_BB_best'):.3f}\nStage-3 final {s3.get('J_BB_best'):.3f}",
        transform=ax.transAxes,
        va="top",
        fontsize=16,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cbd5e1"),
    )
    out = OUT / "18_stage3_optimization_evolution.png"
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig_s2_vs_s3_broadband()
    fig_hole_repair(5.75, "14_stage3_hole_repair_5p75", "5.75 GHz hole repair — Stage-2 vs Stage-3")
    fig_hole_repair(5.825, "15_stage3_hole_repair_5p825", "5.825 GHz hole repair — Stage-2 vs Stage-3")
    fig_4p1_improvement()
    fig_device_summary()
    fig_optimization_evolution()
    print("Stage-3 static figures complete", flush=True)


if __name__ == "__main__":
    main()
