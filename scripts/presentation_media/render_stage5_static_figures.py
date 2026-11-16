#!/usr/bin/env python3
"""Stage-5 presentation STATIC figures (PNG) — no FEM solves.

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

OUT = OUT_ROOT / "stage5_figures"
C_S3 = "#94a3b8"
C_S4 = "#64748b"
C_S5 = C_ACCENT


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
        ax.axvline(
            f, color=color, ls="--", lw=1.6, alpha=0.85, label=label if first and i == 0 else None
        )


def at(rows, f):
    return min(rows, key=lambda r: abs(r["f_GHz"] - f))


def fig_s3_s4_s5_broadband():
    s3 = load("bb_stage3_dense_25MHz.json")
    s4 = load("bb_stage4_dense_25MHz.json")
    s5 = load("bb_stage5_dense_6p25MHz.json")
    f3, a3 = xy(rows_of(s3), "J_A2")
    f4, a4 = xy(rows_of(s4), "J_A2")
    f5, a5 = xy(rows_of(s5), "J_A2")
    _, il3 = xy(rows_of(s3), "IL_dB")
    _, il4 = xy(rows_of(s4), "IL_dB")
    _, il5 = xy(rows_of(s5), "IL_dB")
    _, iso3 = xy(rows_of(s3), "isolation_dB")
    _, iso4 = xy(rows_of(s4), "isolation_dB")
    _, iso5 = xy(rows_of(s5), "isolation_dB")

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(
        "Stage-3 → Stage-4 → Stage-5 broadband comparison (A2_hybrid)",
        fontsize=26,
        fontweight="bold",
        color=C_INK,
        y=0.97,
    )
    ax0 = fig.add_axes([0.07, 0.55, 0.88, 0.35])
    ax1 = fig.add_axes([0.07, 0.10, 0.42, 0.35])
    ax2 = fig.add_axes([0.55, 0.10, 0.40, 0.35])
    for ax in (ax0, ax1, ax2):
        stamp_kind(ax, "DATA", "DENSE AUDIT")

    ax0.plot(f3, a3, color=C_S3, lw=1.8, label="Stage-3 25 MHz")
    ax0.plot(f4, a4, color=C_S4, lw=2.0, label="Stage-4 25 MHz")
    ax0.plot(f5, a5, color=C_S5, lw=2.4, label="Stage-5 6.25 MHz")
    fine = PROD_CAMPAIGN / "bb_stage5_dense_3p125MHz_local.json"
    if fine.exists():
        fr, ar = xy(rows_of(load("bb_stage5_dense_3p125MHz_local.json")), "J_A2")
        for fmark, col, lab in (
            (4.734375, C_BAD, "true worst 4.734375"),
            (5.08125, C_WARN, "5.08125 improved"),
            (5.8125, C_GOOD, "5.8125 watch"),
        ):
            j = int(np.argmin(np.abs(fr - fmark)))
            if abs(fr[j] - fmark) < 0.01:
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
    mark_freqs(ax0, [4.734375, 5.08125, 5.8125], C_BAD, first=False)
    ax0.axhline(0, color="#94a3b8", lw=1)
    ax0.set_ylabel("A2 objective (higher better)")
    ax0.set_xlabel("Frequency (GHz)")
    ax0.legend(loc="lower right", fontsize=11, ncol=2)
    ax0.set_xlim(4.0, 6.0)

    ax1.plot(f3, il3, color=C_S3, lw=1.8, label="Stage-3")
    ax1.plot(f4, il4, color=C_S4, lw=2.0, label="Stage-4")
    ax1.plot(f5, il5, color=C_S5, lw=2.4, label="Stage-5")
    mark_freqs(ax1, [4.734375, 5.08125, 5.8125], C_BAD)
    ax1.set_ylabel("IL (dB)")
    ax1.set_xlabel("Frequency (GHz)")
    ax1.legend(fontsize=12)
    ax1.set_xlim(4.0, 6.0)

    ax2.plot(f3, iso3, color=C_S3, lw=1.8, label="Stage-3")
    ax2.plot(f4, iso4, color=C_S4, lw=2.0, label="Stage-4")
    ax2.plot(f5, iso5, color=C_S5, lw=2.4, label="Stage-5")
    mark_freqs(ax2, [4.734375, 5.08125, 5.8125], C_BAD)
    ax2.axhline(0, color="#94a3b8", lw=1)
    ax2.set_ylabel("Isolation (dB)")
    ax2.set_xlabel("Frequency (GHz)")
    ax2.legend(fontsize=12)
    ax2.set_xlim(4.0, 6.0)

    out = OUT / "25_stage3_vs_stage4_vs_stage5_broadband.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def fig_pair_compare(f_ghz: float, stem: str, title: str, left_name: str, right_name: str,
                     left_label="Stage-4", right_label="Stage-5"):
    left = load(left_name)
    right = load(right_name)
    rL, rR = at(rows_of(left), f_ghz), at(rows_of(right), f_ghz)
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(title, fontsize=26, fontweight="bold", color=C_INK, y=0.96)
    ax0 = fig.add_axes([0.08, 0.18, 0.38, 0.68])
    ax1 = fig.add_axes([0.55, 0.18, 0.38, 0.68])
    x = np.arange(3)
    w = 0.38
    ax0.bar(x - w / 2, [rL["J_A2"], rL["isolation_dB"], rL["IL_dB"]], w, color=C_S4, label=left_label)
    ax0.bar(x + w / 2, [rR["J_A2"], rR["isolation_dB"], rR["IL_dB"]], w, color=C_S5, label=right_label)
    ax0.set_xticks(x)
    ax0.set_xticklabels(["A2", "iso (dB)", "IL (dB)"])
    ax0.set_ylabel("Metric value")
    ax0.set_title(f"nearest sample ≈ {f_ghz:.5g} GHz", fontsize=16)
    ax0.legend(fontsize=12)
    ax0.axhline(0, color="#94a3b8", lw=1)
    stamp_kind(ax0, "DATA")
    ax1.axis("off")
    stamp_kind(ax1, "DATA", "QUANTITATIVE")
    txt = (
        f"Requested: {f_ghz:.5g} GHz\n"
        f"{left_label} sample: {rL['f_GHz']:.5g} GHz\n"
        f"{right_label} sample: {rR['f_GHz']:.5g} GHz\n\n"
        f"{left_label}:\n"
        f"  A2 = {rL['J_A2']:.3f}\n"
        f"  isolation = {rL['isolation_dB']:.2f} dB\n"
        f"  IL = {rL['IL_dB']:.2f} dB\n"
        f"  T_des = {rL.get('T_des', float('nan')):.4e}\n\n"
        f"{right_label}:\n"
        f"  A2 = {rR['J_A2']:.3f}\n"
        f"  isolation = {rR['isolation_dB']:.2f} dB\n"
        f"  IL = {rR['IL_dB']:.2f} dB\n"
        f"  T_des = {rR.get('T_des', float('nan')):.4e}\n\n"
        f"ΔA2 = {rR['J_A2'] - rL['J_A2']:+.3f}\n"
        f"Δiso = {rR['isolation_dB'] - rL['isolation_dB']:+.2f} dB"
    )
    ax1.text(
        0.05,
        0.95,
        txt,
        transform=ax1.transAxes,
        va="top",
        ha="left",
        fontsize=16,
        family="monospace",
        color=C_INK,
        bbox=dict(boxstyle="round,pad=0.6", fc="#f8fafc", ec="#cbd5e1"),
    )
    out = OUT / f"{stem}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def fig_progression_5p8125():
    """Stage-3 → Stage-4 → Stage-5 at 5.8125 GHz."""
    specs = [
        ("Stage-3", "bb_stage3_localize_5p8_6p25MHz.json", C_S3),
        ("Stage-4", "bb_stage4_dense_6p25MHz_local.json", C_S4),
        ("Stage-5", "bb_stage5_dense_6p25MHz.json", C_S5),
    ]
    rows = []
    for lab, name, col in specs:
        r = at(rows_of(load(name)), 5.8125)
        rows.append((lab, r, col))
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(
        "5.8125 GHz progression · Stage-3 → Stage-4 → Stage-5",
        fontsize=26,
        fontweight="bold",
        color=C_INK,
        y=0.96,
    )
    ax0 = fig.add_axes([0.08, 0.18, 0.50, 0.68])
    ax1 = fig.add_axes([0.62, 0.18, 0.32, 0.68])
    x = np.arange(3)
    w = 0.25
    for i, (lab, r, col) in enumerate(rows):
        ax0.bar(
            x + (i - 1) * w,
            [r["J_A2"], r["isolation_dB"], r["IL_dB"]],
            w,
            color=col,
            label=f"{lab} @ {r['f_GHz']:.5g}",
        )
    ax0.set_xticks(x)
    ax0.set_xticklabels(["A2", "iso (dB)", "IL (dB)"])
    ax0.axhline(0, color="#94a3b8", lw=1)
    ax0.legend(fontsize=12)
    ax0.set_ylabel("Metric value")
    stamp_kind(ax0, "DATA")
    ax1.axis("off")
    stamp_kind(ax1, "DATA", "QUANTITATIVE")
    lines = []
    for lab, r, _ in rows:
        lines.append(
            f"{lab} ({r['f_GHz']:.5g} GHz):\n"
            f"  A2={r['J_A2']:.3f}  IL={r['IL_dB']:.2f}  iso={r['isolation_dB']:.2f}\n"
        )
    ax1.text(
        0.05,
        0.95,
        "\n".join(lines),
        transform=ax1.transAxes,
        va="top",
        fontsize=15,
        family="monospace",
        color=C_INK,
        bbox=dict(boxstyle="round,pad=0.6", fc="#f8fafc", ec="#cbd5e1"),
    )
    out = OUT / "28_stage5_progression_5p8125.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def fig_device_summary():
    s5 = load("bb_stage5_dense_6p25MHz.json")
    opt = load("opt_broadband_stage5_A2_C.json")
    summ = load("BB_STAGE5_DENSE_AUDIT_SUMMARY.json")
    f, a2 = xy(rows_of(s5), "J_A2")
    _, iso = xy(rows_of(s5), "isolation_dB")

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
        "Stage-5 broadband plasma circulator — inverse-design summary",
        fontsize=26,
        fontweight="bold",
        color=C_INK,
        y=0.97,
    )
    ax0 = fig.add_axes([0.04, 0.12, 0.42, 0.78])
    ax1 = fig.add_axes([0.52, 0.55, 0.44, 0.35])
    ax2 = fig.add_axes([0.52, 0.12, 0.44, 0.35])
    stamp_kind(ax0, "DATA", "Stage-5 q_best")
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
    R = float(np.max(np.linalg.norm(xy0, axis=1))) * 1.28
    for i, ang_deg in enumerate(range(0, 360, 60)):
        ang = np.deg2rad(ang_deg)
        ax0.text(
            R * np.cos(ang),
            R * np.sin(ang),
            f"P{i+1}",
            ha="center",
            va="center",
            fontsize=12,
            color=C_MUTED,
        )
    ax0.set_aspect("equal")
    ax0.set_title("Plasma bulb scale factors (C6)", fontsize=16)
    fig.colorbar(sca, ax=ax0, fraction=0.046, pad=0.04, label="s")

    ax1.plot(f, a2, color=C_S5, lw=2.4)
    ax1.axhline(0, color="#94a3b8", lw=1)
    tw = summ.get("true_worst") or {}
    if tw:
        ax1.axvline(tw["f_GHz"], color=C_BAD, ls="--", lw=1.6)
        ax1.scatter([tw["f_GHz"]], [tw["J_A2"]], color=C_BAD, s=80, zorder=5)
    ax1.set_ylabel("A2")
    ax1.set_xlabel("Frequency (GHz)")
    ax1.set_xlim(4.0, 6.0)
    ax1.set_title(
        f"J_BB={opt.get('J_BB_best'):.3f}  ·  true worst @ {tw.get('f_GHz', float('nan')):.5g} GHz",
        fontsize=14,
    )

    ax2.plot(f, iso, color=C_S5, lw=2.4)
    ax2.axhline(0, color="#94a3b8", lw=1)
    if tw:
        ax2.axvline(tw["f_GHz"], color=C_BAD, ls="--", lw=1.6)
    ax2.set_ylabel("Isolation (dB)")
    ax2.set_xlabel("Frequency (GHz)")
    ax2.set_xlim(4.0, 6.0)

    out = OUT / "30_stage5_device_summary.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def fig_true_worst():
    summ = load("BB_STAGE5_DENSE_AUDIT_SUMMARY.json")
    tw = summ["true_worst"]
    s5 = load("bb_stage5_dense_3p125MHz_local.json")
    f, a2 = xy(rows_of(s5), "J_A2")
    # focus neighborhood
    mask = (f >= tw["f_GHz"] - 0.05) & (f <= tw["f_GHz"] + 0.05)
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.suptitle(
        f"Stage-5 true dense-audit worst @ {tw['f_GHz']:.5g} GHz",
        fontsize=26,
        fontweight="bold",
        color=C_INK,
        y=0.96,
    )
    ax = fig.add_axes([0.08, 0.12, 0.84, 0.75])
    stamp_kind(ax, "DATA", "3.125 MHz LOCAL")
    ax.plot(f[mask], a2[mask], color=C_S5, lw=2.6)
    ax.scatter([tw["f_GHz"]], [tw["J_A2"]], s=120, color=C_BAD, zorder=5, label="true worst")
    # mark 6.25 MHz neighbors
    for fn in (4.73125, 4.7375):
        ax.axvline(fn, color=C_WARN, ls=":", lw=1.4, alpha=0.9)
    ax.axhline(0, color="#94a3b8", lw=1)
    ax.axvline(tw["f_GHz"], color=C_BAD, ls="--", lw=1.6)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("A2")
    ax.legend(fontsize=14)
    ax.text(
        0.02,
        0.05,
        f"A2={tw['J_A2']:.3f}\nIL={tw['IL_dB']:.2f} dB\niso={tw['isolation_dB']:.2f} dB\n"
        f"T_des={tw.get('T_des', float('nan')):.4e}\n"
        f"between 6.25 MHz nodes 4.73125 / 4.7375",
        transform=ax.transAxes,
        fontsize=16,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cbd5e1"),
    )
    out = OUT / "29_stage5_true_worst_4p734375.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    print("wrote", out, flush=True)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig_s3_s4_s5_broadband()
    fig_pair_compare(
        5.08125,
        "26_stage5_repair_5p08125",
        "Stage-4 5.08125 GHz failure → Stage-5",
        "bb_stage4_dense_6p25MHz_local.json",
        "bb_stage5_dense_6p25MHz.json",
    )
    fig_pair_compare(
        4.73125,
        "27_stage5_repair_4p73125",
        "Stage-4 ~4.73125 GHz failure → Stage-5",
        "bb_stage4_localize_4p73_6p25MHz.json",
        "bb_stage5_dense_6p25MHz.json",
    )
    # also true-worst vs Stage-4 nearest
    fig_pair_compare(
        4.734375,
        "27b_stage5_true_worst_vs_stage4",
        "True Stage-5 worst 4.734375 GHz vs Stage-4 neighborhood",
        "bb_stage4_localize_4p73_6p25MHz.json",
        "bb_stage5_dense_3p125MHz_local.json",
    )
    fig_progression_5p8125()
    fig_true_worst()
    fig_device_summary()
    print("Stage-5 static figures complete →", OUT, flush=True)


if __name__ == "__main__":
    main()
