#!/usr/bin/env python3
"""05 — Broadband frequency sweep from EXISTING spectral / Stage-2 records (no solves)."""
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

OUT = OUT_ROOT / "05_broadband_frequency_sweep"
FRAMES = OUT / "frames"


def rows_xy(rows, ykey):
    f = np.array([r["f_GHz"] for r in rows], float)
    y = np.array([r.get(ykey) if r.get(ykey) is not None else np.nan for r in rows], float)
    return f, y


def main():
    clear_frames(FRAMES)
    # Stage-1 dense audits / reconnaissance spectra (SF / production-era)
    spec100 = load_json(PROD_CAMPAIGN / "spectral_100MHz.json")
    spec50 = load_json(PROD_CAMPAIGN / "spectral_50MHz.json")
    # Stage-1 production final objective-grid metrics live in opt + audit
    prod = load_json(PROD_CAMPAIGN / "opt_broadband_production_C.json")
    audit50 = load_json(PROD_CAMPAIGN / "bb_final_audit_50MHz.json")
    stage2 = load_json(PROD_CAMPAIGN / "opt_broadband_stage2_A2_C.json")
    s2hist = load_json(PROD_CAMPAIGN / "bb_stage2_history.json")["history"]
    start = next(r for r in s2hist if r.get("kind") == "start_eval")
    end_acc = [r for r in s2hist if r.get("kind") == "accepted_iter"][-1]

    scenes = []

    # Scene A: stage-1 reconnaissance IL / iso on 100 & 50 MHz grids
    scenes.append(("stage1_recon", "Stage-1 spectral reconnaissance (existing SF scan data)"))
    # Scene B: stage-1 production design audit 50 MHz
    scenes.append(("stage1_prod", "Stage-1 production design — 50 MHz audit spectra"))
    # Scene C: stage-2 start vs final objective-grid J_A2 / IL proxy from history
    scenes.append(("stage2_compare", "Stage-2 A2: start vs final objective-grid metrics (history)"))

    paths = []
    fi = 0
    for scene_id, scene_title in scenes:
        for k in range(45):  # 1.5 s each
            fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
            fig.suptitle(
                f"Broadband frequency sweep — {scene_title}",
                fontsize=24, fontweight="bold", color=C_INK, y=0.97,
            )
            ax1 = fig.add_axes([0.07, 0.55, 0.40, 0.35])
            ax2 = fig.add_axes([0.55, 0.55, 0.40, 0.35])
            ax3 = fig.add_axes([0.07, 0.10, 0.88, 0.35])
            for ax in (ax1, ax2, ax3):
                stamp_kind(ax, "DATA")

            if scene_id == "stage1_recon":
                for spec, c, lab in (
                    (spec100, C_MUTED, "100 MHz grid"),
                    (spec50, C_ACCENT, "50 MHz grid"),
                ):
                    f, il = rows_xy(spec["rows"], "IL_dB")
                    f2, iso = rows_xy(spec["rows"], "isolation_dB")
                    ax1.plot(f, il, "-o", ms=3, color=c, label=lab)
                    ax2.plot(f2, iso, "-o", ms=3, color=c, label=lab)
                ax1.set_title("Insertion loss (dB)")
                ax2.set_title("Isolation (dB)")
                ax1.set_ylabel("IL (dB)")
                ax2.set_ylabel("iso (dB)")
                ax3.axis("off")
                ax3.text(
                    0.01, 0.8,
                    "Source: spectral_100MHz.json / spectral_50MHz.json (campaign reconnaissance).\n"
                    "These are measured spectra on saved designs — not a new FEM solve.\n"
                    f"50 MHz worst iso={spec50['worst_isolation_dB']['value']:.1f} dB @ "
                    f"{spec50['worst_isolation_dB']['f_GHz']} GHz; "
                    f"worst IL={spec50['worst_IL_dB']['value']:.1f} dB @ "
                    f"{spec50['worst_IL_dB']['f_GHz']} GHz.",
                    transform=ax3.transAxes, fontsize=16, va="top", color=C_INK,
                )
                ax1.legend(frameon=False, fontsize=12)
                ax2.legend(frameon=False, fontsize=12)

            elif scene_id == "stage1_prod":
                f, il = rows_xy(audit50["rows"], "IL_dB")
                f2, iso = rows_xy(audit50["rows"], "isolation_dB")
                f3, tdes = rows_xy(audit50["rows"], "T_des")
                ax1.plot(f, il, "-o", ms=3.5, color=C_WARN)
                ax2.plot(f2, iso, "-o", ms=3.5, color=C_ACCENT)
                ax3.plot(f3, tdes, "-o", ms=3.5, color=C_GOOD)
                ax1.axvline(4.1, color=C_BAD, ls="--", label="4.1 GHz hole")
                ax2.axvline(4.1, color=C_BAD, ls="--")
                ax3.axvline(4.1, color=C_BAD, ls="--")
                ax1.set_title("Stage-1 final design IL")
                ax2.set_title("Stage-1 final design isolation")
                ax3.set_title(r"Stage-1 final $T_{des}$ — motivates Stage-2 A2")
                ax1.set_ylabel("IL (dB)")
                ax2.set_ylabel("iso (dB)")
                ax3.set_ylabel(r"$T_{des}$")
                ax3.set_xlabel("f (GHz)")
                jbb = prod.get("J_BB_best", float("nan"))
                ax1.legend(frameon=False)
                ax3.text(
                    0.99, 0.05,
                    f"opt_broadband_production_C.json  J_BB_best={jbb:.3f}",
                    transform=ax3.transAxes, ha="right", fontsize=13, color=C_MUTED,
                )

            else:  # stage2_compare
                f0 = np.asarray(start["freqs_GHz"], float)
                j0 = np.asarray(start["J_per_freq"], float)
                f1 = np.asarray(end_acc["freqs_GHz"], float)
                j1 = np.asarray(end_acc["J_per_freq"], float)
                # IL from at_4p1 / worst stored — also plot IL reconstructed if present in worst only
                # Use history metrics: plot J_A2 and annotate IL@4.1
                ax1.plot(f0, j0, "-o", ms=3.5, color=C_MUTED, label="Stage-2 start")
                ax1.plot(f1, j1, "-o", ms=3.5, color=C_ACCENT, label="Stage-2 final (nit30)")
                ax1.axvline(4.1, color=C_GOOD, ls="--")
                ax1.axvline(float(end_acc["worst_f_GHz"]), color=C_BAD, ls=":")
                ax1.set_ylabel(r"$J_{A2}$")
                ax1.set_title("Objective-grid A2 spectra")
                ax1.legend(frameon=False, fontsize=12)

                # 4.1 recovery bars
                ax2.axis("off")
                a0, a1 = start["at_4p1"], end_acc["at_4p1"]
                txt = (
                    "4.1 GHz recovery (from Stage-2 history)\n"
                    f"  T_des: {a0['T_des']:.4e} → {a1['T_des']:.4f}\n"
                    f"  IL:    {a0['IL_dB']:.1f} dB → {a1['IL_dB']:.1f} dB\n"
                    f"  iso:   {a0['isolation_dB']:.1f} → {a1['isolation_dB']:.1f} dB\n\n"
                    f"Worst freq migrated: {start['worst_f_GHz']} → {end_acc['worst_f_GHz']} GHz\n"
                    f"best J_BB: {stage2.get('J_BB_best', float('nan')):.3f}"
                )
                ax2.text(0.05, 0.95, txt, transform=ax2.transAxes, va="top", fontsize=16,
                         family="monospace",
                         bbox=dict(boxstyle="round,pad=0.5", fc="#eef8f0", ec=C_GOOD))

                # T_des @4.1 and worst over accepted iters
                acc = [r for r in s2hist if r.get("kind") == "accepted_iter"]
                nits = [0] + [r["iteration"] for r in acc]
                tdes = [start["at_4p1"]["T_des"]] + [r["at_4p1"]["T_des"] for r in acc]
                wf = [start["worst_f_GHz"]] + [r["worst_f_GHz"] for r in acc]
                ax3.plot(nits, tdes, "-o", color=C_GOOD, label=r"$T_{des}$ @ 4.1 GHz")
                ax3b = ax3.twinx()
                ax3b.plot(nits, wf, "-s", color=C_BAD, ms=4, label="worst f (GHz)")
                ax3.set_xlabel("accepted iteration")
                ax3.set_ylabel(r"$T_{des}$ @ 4.1")
                ax3b.set_ylabel("worst f (GHz)")
                ax3.set_title("Stage-2 recovery trajectory (DATA)")
                h1, l1 = ax3.get_legend_handles_labels()
                h2, l2 = ax3b.get_legend_handles_labels()
                ax3.legend(h1 + h2, l1 + l2, loc="center right", frameon=False, fontsize=12)

            for ax in (ax1, ax2, ax3):
                if hasattr(ax, "grid"):
                    ax.grid(True, alpha=0.25)
            p = FRAMES / f"frame_{fi:04d}.png"
            save_frame(fig, p)
            paths.append(p)
            fi += 1

    mp4 = OUT / "05_broadband_frequency_sweep_stage1.mp4"
    encode_mp4(str(FRAMES / "frame_*.png"), mp4)
    # also alias comparison name
    import shutil
    shutil.copy(mp4, OUT / "05_broadband_frequency_sweep_stage2_comparison.mp4")
    write_poster(paths[len(paths) // 2], OUT / "05_broadband_frequency_sweep_poster.png")
    frames_to_gif(paths, OUT / "05_broadband_frequency_sweep.gif", fps=10)
    print("wrote", mp4)


if __name__ == "__main__":
    main()
