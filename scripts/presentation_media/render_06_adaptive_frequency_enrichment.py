#!/usr/bin/env python3
"""06 — Adaptive frequency enrichment story from ACTUAL Stage-1 audit artifacts."""
from __future__ import annotations

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
    clear_frames,
    encode_mp4,
    frames_to_gif,
    load_json,
    save_frame,
    stamp_kind,
    write_poster,
)

OUT = OUT_ROOT / "06_adaptive_frequency_enrichment"
FRAMES = OUT / "frames"

# From BB_ENRICHMENT_RECOMMENDATION.md / production meta
STAGE1_OBJ = [4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0]
ENRICH = [4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95]
STAGE2_OBJ = sorted(set(STAGE1_OBJ + ENRICH))


def main():
    clear_frames(FRAMES)
    s100 = load_json(PROD_CAMPAIGN / "spectral_100MHz.json")
    s50 = load_json(PROD_CAMPAIGN / "spectral_50MHz.json")
    s25band = load_json(PROD_CAMPAIGN / "spectral_25MHz_5p50_6p00.json")
    cmp_ = load_json(PROD_CAMPAIGN / "spectral_grid_comparison.json")
    # between-grid hole evidence from enrichment doc numbers also in audit50 after stage1
    audit50 = load_json(PROD_CAMPAIGN / "bb_final_audit_50MHz.json")
    audit25 = load_json(PROD_CAMPAIGN / "bb_final_audit_25MHz_local.json")

    phases = [
        ("100 MHz objective / recon grid", 100, s100, STAGE1_OBJ),
        ("50 MHz dense audit", 50, s50, STAGE1_OBJ),
        ("25 MHz local band (5.50–6.00)", 25, s25band, STAGE1_OBJ),
        ("Enrichment → Stage-2 28-freq grid", 5.625, None, STAGE2_OBJ),
    ]

    paths = []
    fi = 0
    for pi, (title, tag, spec, obj_freqs) in enumerate(phases):
        for k in range(50):
            t = k / 49
            fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
            fig.suptitle(
                "Adaptive frequency enrichment (Stage-1 audits → Stage-2 grid)",
                fontsize=26, fontweight="bold", color=C_INK, y=0.97,
            )
            ax = fig.add_axes([0.07, 0.38, 0.86, 0.50])
            ax2 = fig.add_axes([0.07, 0.08, 0.86, 0.22])
            stamp_kind(ax, "DATA", "campaign spectral_* + enrichment rec")
            stamp_kind(ax2, "DATA")

            # progress ladder
            steps = ["100 MHz", "50 MHz", "25 MHz", "enrich +5.625…"]
            for si, lab in enumerate(steps):
                x0 = 0.05 + si * 0.24
                active = si <= pi
                ax2.add_patch(
                    plt.Rectangle(
                        (x0, 0.35), 0.20, 0.45,
                        transform=ax2.transAxes,
                        fc=C_ACCENT if active else "#dbe3f0",
                        ec=C_INK, lw=1.2, clip_on=False,
                    )
                )
                ax2.text(
                    x0 + 0.10, 0.57, lab, transform=ax2.transAxes,
                    ha="center", va="center", fontsize=14,
                    color="white" if active else C_MUTED, fontweight="bold",
                )
                if si < 3:
                    ax2.annotate(
                        "",
                        xy=(x0 + 0.22, 0.57),
                        xytext=(x0 + 0.20, 0.57),
                        xycoords=ax2.transAxes, textcoords=ax2.transAxes,
                        arrowprops=dict(arrowstyle="-|>", color=C_MUTED),
                    )
            ax2.axis("off")
            ax2.set_title(f"Phase: {title}", fontsize=16, loc="left")

            if spec is not None:
                rows = spec["rows"]
                f = np.array([r["f_GHz"] for r in rows], float)
                # prefer J_A if present else isolation
                if "J_A" in rows[0] and rows[0]["J_A"] is not None:
                    y = np.array([r["J_A"] for r in rows], float)
                    ylab = r"$J_A$"
                else:
                    y = np.array([r["isolation_dB"] for r in rows], float)
                    ylab = "isolation (dB)"
                ax.plot(f, y, "-o", ms=4, color=C_ACCENT, label=f"{tag} MHz samples" if tag != 5.625 else "")
                # mark objective nodes
                for fo in obj_freqs:
                    ax.axvline(fo, color="#94a3b8", lw=0.6, alpha=0.5)
                # highlight known holes
                for hole, col in ((5.65, C_WARN), (5.625, C_BAD), (4.1, C_MUTED)):
                    ax.axvline(hole, color=col, ls="--", lw=1.5)
                ax.set_ylabel(ylab)
                # annotate worst
                wi = int(np.nanargmin(y))
                ax.scatter([f[wi]], [y[wi]], s=120, color=C_BAD, zorder=5)
                ax.text(f[wi], y[wi], f"  worst @ {f[wi]:.3f}", color=C_BAD, fontsize=13, va="bottom")
            else:
                # show enrichment markers on a unit axis
                ax.set_ylim(-1.5, 1.5)
                ax.axhline(0, color="#cbd5e1")
                for fo in STAGE1_OBJ:
                    ax.plot([fo], [0], "o", color="#94a3b8", ms=7)
                for fo in ENRICH:
                    ax.plot([fo], [0], "D", color=C_BAD, ms=10)
                    ax.text(fo, 0.35, f"{fo:g}", rotation=90, ha="center", fontsize=11, color=C_BAD)
                ax.text(0.02, 0.92, "Gray circles: Stage-1 21×100 MHz objective nodes",
                        transform=ax.transAxes, fontsize=14, color=C_MUTED)
                ax.text(0.02, 0.82, "Red diamonds: enriched nodes (incl. 5.625 GHz between-grid hole)",
                        transform=ax.transAxes, fontsize=14, color=C_BAD)
                ax.set_ylabel("grid markers")

            ax.set_xlim(3.95, 6.05)
            ax.set_xlabel("f (GHz)")
            ax.set_title(title, fontsize=18)
            ax.grid(True, alpha=0.25)

            # footnote with documented numbers
            note = (
                "Documented between-grid failure (BB_ENRICHMENT_RECOMMENDATION.md): "
                "21×100 MHz worst J_A=−0.657 @ 5.7; 50 MHz audit −0.827 @ 5.65; "
                "25 MHz local −1.153 @ 5.625 → add enrichment set before Stage-2."
            )
            fig.text(0.07, 0.01, note, fontsize=12, color=C_MUTED)

            # fade-in highlight on enrichment phase
            if pi == 3:
                ax.set_facecolor((1.0, 0.97, 0.97, 0.35 + 0.2 * t))

            p = FRAMES / f"frame_{fi:04d}.png"
            save_frame(fig, p)
            paths.append(p)
            fi += 1

    # comparison callouts (static last frame duplicate with cmp numbers)
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    ax = fig.add_axes([0.08, 0.12, 0.84, 0.75])
    ax.axis("off")
    stamp_kind(ax, "DATA", "spectral_grid_comparison.json")
    fig.suptitle("Grid-sufficiency comparisons (saved)", fontsize=28, fontweight="bold", y=0.95)
    lines = [
        f"200 vs 100: Δiso={cmp_['200_vs_100']['d_worst_iso_dB']:.2f} dB  pass_iso={cmp_['200_vs_100']['pass_iso']}",
        f"100 vs 50:  Δiso={cmp_['100_vs_50']['d_worst_iso_dB']:.2f} dB  pass_iso={cmp_['100_vs_50']['pass_iso']}",
        f"50 vs 25 (5.50–6.00 band): Δiso={cmp_['50_vs_25_band_5p50_6p00']['d_worst_iso_dB']:.2f} dB  "
        f"pass_iso={cmp_['50_vs_25_band_5p50_6p00']['pass_iso']}",
        "",
        f"Stage-1 50 MHz audit rows: {len(audit50['rows'])}",
        f"Stage-1 25 MHz local audit rows: {len(audit25['rows'])}",
        "",
        "Stage-2 objective grid = 21 + 7 enriched = 28 frequencies (fixed during L-BFGS).",
    ]
    ax.text(0.05, 0.9, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=20,
            family="monospace", color=C_INK,
            bbox=dict(boxstyle="round,pad=0.6", fc="#f8fafc", ec=C_ACCENT))
    for _ in range(40):
        p = FRAMES / f"frame_{fi:04d}.png"
        if _ == 0:
            save_frame(fig, p)
            first = p
        else:
            import shutil
            shutil.copy(first, p)
        paths.append(p)
        fi += 1

    mp4 = OUT / "06_adaptive_frequency_enrichment.mp4"
    encode_mp4(str(FRAMES / "frame_*.png"), mp4)
    write_poster(paths[len(paths) // 2], OUT / "06_adaptive_frequency_enrichment_poster.png")
    frames_to_gif(paths, OUT / "06_adaptive_frequency_enrichment.gif", fps=10)
    print("wrote", mp4)


if __name__ == "__main__":
    main()
