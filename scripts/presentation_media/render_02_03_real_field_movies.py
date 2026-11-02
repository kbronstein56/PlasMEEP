#!/usr/bin/env python3
"""02/03 — Real FEM complex-Hz field movies (no fabricated fields).

Reads production exports from:
  PlasMEEP/outputs/.../jesse_method_campaign/presentation_fields/

Animates: frame(t) = Re[ Hz * exp(-i ω t) ], ω = 2π f.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from media_common import (  # noqa: E402
    C_INK,
    C_MUTED,
    FPS,
    OUT_ROOT,
    PROD_CAMPAIGN,
    clear_frames,
    encode_mp4,
    frames_to_gif,
    save_frame,
    stamp_kind,
    write_poster,
)

FIELD_ROOT = PROD_CAMPAIGN / "presentation_fields"
N_PERIODS = 2
FRAMES_PER_PERIOD = 24  # 48 frames @ 2 periods → 1.6 s @ 30 fps


def load_export(stem: str):
    d = FIELD_ROOT / stem
    if not d.exists():
        # fuzzy match by prefix
        cands = sorted(FIELD_ROOT.glob(stem + "*")) if FIELD_ROOT.exists() else []
        if not cands:
            raise FileNotFoundError(f"missing field export: {d}")
        d = cands[0]
    meta = json.loads((d / "meta.json").read_text())
    x = np.load(d / "x.npy")
    y = np.load(d / "y.npy")
    Re = np.load(d / "Hz_real.npy")
    Im = np.load(d / "Hz_imag.npy")
    Hz = Re + 1j * Im
    return d, meta, x, y, Hz


def find_stem(design_tag: str, f_ghz: float) -> str:
    """Match export_presentation_fields naming."""
    # e.g. stage2_best_P1_4p100GHz_Bplus
    ftag = f"{f_ghz:.3f}".replace(".", "p")
    stem = f"{design_tag}_P1_{ftag}GHz_Bplus"
    exact = FIELD_ROOT / stem
    if exact.exists():
        return stem
    # tolerate float formatting variants
    for p in sorted(FIELD_ROOT.glob(f"{design_tag}_P1_*GHz_Bplus")):
        meta = json.loads((p / "meta.json").read_text())
        if abs(float(meta["f_GHz"]) - f_ghz) < 1e-6:
            return p.name
    raise FileNotFoundError(stem)


def field_extent(x, y):
    return [float(x.min()), float(x.max()), float(y.min()), float(y.max())]


def animate_pair(
    out_dir: Path,
    stem: str,
    title: str,
    left,
    right,
    left_title: str,
    right_title: str,
    subtitle: str,
    require_same_freq: bool = True,
):
    """Side-by-side Re[Hz e^{-iωt}] movies; shared color scale from both |Hz|.

    Phase advances in fractions of each panel's own RF period so unequal
    frequencies (e.g. 5.65 vs 5.75 GHz) still show correct oscillation.
    """
    clear_frames(out_dir / "frames")
    _, meta_l, x_l, y_l, Hz_l = left
    _, meta_r, x_r, y_r, Hz_r = right
    if require_same_freq:
        assert abs(meta_l["f_GHz"] - meta_r["f_GHz"]) < 1e-9
    vmax = float(
        np.nanpercentile(
            np.concatenate([np.abs(Hz_l).ravel(), np.abs(Hz_r).ravel()]), 99.5
        )
    )
    vmax = max(vmax, 1e-30)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    n_frames = N_PERIODS * FRAMES_PER_PERIOD
    frame_paths = []
    for i in range(n_frames):
        # Same frame index = same fraction of each panel's RF period
        phase = 2.0 * np.pi * (i / FRAMES_PER_PERIOD)
        Fl = np.real(Hz_l * np.exp(-1j * phase))
        Fr = np.real(Hz_r * np.exp(-1j * phase))
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(title, fontsize=26, fontweight="bold", color=C_INK, y=0.97)
        ax0 = fig.add_axes([0.04, 0.12, 0.44, 0.78])
        ax1 = fig.add_axes([0.52, 0.12, 0.44, 0.78])
        for ax, F, x, y, ttl, meta in (
            (ax0, Fl, x_l, y_l, left_title, meta_l),
            (ax1, Fr, x_r, y_r, right_title, meta_r),
        ):
            im = ax.imshow(
                F,
                origin="lower",
                extent=field_extent(x, y),
                cmap="RdBu_r",
                norm=norm,
                aspect="equal",
                interpolation="bilinear",
            )
            ax.set_title(ttl, fontsize=18)
            ax.set_xlabel("x (cm)")
            ax.set_ylabel("y (cm)")
            stamp_kind(ax, "SIMULATION", r"Re[Hz e$^{-i\omega t}$]")
            ax.text(
                0.98,
                0.02,
                f"P{meta['source_port']}  B={meta['B_tesla']:+.3f} T\n"
                f"f={meta['f_GHz']:.3f} GHz  φ={phase/np.pi:.2f}π",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=12,
                color=C_MUTED,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cbd5e1", alpha=0.9),
            )
        fig.colorbar(im, ax=[ax0, ax1], fraction=0.02, pad=0.02, label="Re[Hz] (rel.)")
        fig.text(0.5, 0.02, subtitle, ha="center", fontsize=13, color=C_MUTED)
        fp = out_dir / "frames" / f"frame_{i:04d}.png"
        save_frame(fig, fp)
        frame_paths.append(fp)

    mp4 = out_dir / f"{stem}.mp4"
    encode_mp4(str(out_dir / "frames" / "frame_*.png"), mp4, fps=FPS)
    write_poster(frame_paths[0], out_dir / f"{stem}_poster.png")
    frames_to_gif(frame_paths, out_dir / f"{stem}_preview.gif", fps=12)
    (out_dir / f"{stem}_PROVENANCE.txt").write_text(
        json.dumps(
            {
                "left": meta_l,
                "right": meta_r,
                "animation": "Re[Hz * exp(-1j*omega*t)]",
                "field_root": str(FIELD_ROOT),
                "n_frames": n_frames,
                "periods": N_PERIODS,
            },
            indent=2,
        )
        + "\n"
    )
    print("wrote", mp4)


def animate_single(
    out_dir: Path,
    stem: str,
    title: str,
    export,
    subtitle: str,
):
    clear_frames(out_dir / "frames")
    _, meta, x, y, Hz = export
    vmax = float(np.nanpercentile(np.abs(Hz), 99.5))
    vmax = max(vmax, 1e-30)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    n_frames = N_PERIODS * FRAMES_PER_PERIOD
    frame_paths = []
    for i in range(n_frames):
        phase = 2.0 * np.pi * (i / FRAMES_PER_PERIOD)
        F = np.real(Hz * np.exp(-1j * phase))
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(title, fontsize=26, fontweight="bold", color=C_INK, y=0.97)
        ax = fig.add_axes([0.08, 0.12, 0.78, 0.78])
        im = ax.imshow(
            F,
            origin="lower",
            extent=field_extent(x, y),
            cmap="RdBu_r",
            norm=norm,
            aspect="equal",
            interpolation="bilinear",
        )
        stamp_kind(ax, "SIMULATION", r"Re[Hz e$^{-i\omega t}$]")
        ax.set_xlabel("x (cm)")
        ax.set_ylabel("y (cm)")
        ax.text(
            0.98,
            0.02,
            f"P{meta['source_port']}  B={meta['B_tesla']:+.3f} T\n"
            f"f={meta['f_GHz']:.3f} GHz  "
            f"IL={meta['metrics_te0'].get('IL_avg_dB', float('nan')):.1f} dB  "
            f"iso={meta['metrics_te0'].get('isolation_avg_dB', float('nan')):.1f} dB",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=13,
            color=C_MUTED,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cbd5e1", alpha=0.9),
        )
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Re[Hz] (rel.)")
        fig.text(0.5, 0.02, subtitle, ha="center", fontsize=13, color=C_MUTED)
        fp = out_dir / "frames" / f"frame_{i:04d}.png"
        save_frame(fig, fp)
        frame_paths.append(fp)
    mp4 = out_dir / f"{stem}.mp4"
    encode_mp4(str(out_dir / "frames" / "frame_*.png"), mp4, fps=FPS)
    write_poster(frame_paths[0], out_dir / f"{stem}_poster.png")
    frames_to_gif(frame_paths, out_dir / f"{stem}_preview.gif", fps=12)
    (out_dir / f"{stem}_PROVENANCE.txt").write_text(
        json.dumps({"meta": meta, "animation": "Re[Hz * exp(-1j*omega*t)]"}, indent=2)
        + "\n"
    )
    print("wrote", mp4)


def main():
    if not FIELD_ROOT.exists() or not (FIELD_ROOT / "INDEX.json").exists():
        raise SystemExit(
            f"Field exports not ready at {FIELD_ROOT}. "
            "Run export_presentation_fields.py in production first."
        )

    # 02 — uniform vs Stage-2 at 4.1 GHz
    u41 = load_export(find_stem("uniform", 4.1))
    s41 = load_export(find_stem("stage2_best", 4.1))
    out02 = OUT_ROOT / "02_uniform_vs_optimized_field_movie"
    animate_pair(
        out02,
        "02_uniform_vs_optimized_field_movie",
        "02 — Uniform vs Stage-2 best · P1 · 4.1 GHz · B+",
        u41,
        s41,
        "Uniform (q=1)",
        "Stage-2 A2 best",
        "Production Direct-DtN FEM · grade C · real complex Hz exports · not schematic",
    )

    # 03a — successful circulation at 5.65 GHz (on-grid Stage-2 worst / representative)
    s565 = load_export(find_stem("stage2_best", 5.65))
    out03 = OUT_ROOT / "03_circulation_field_movie"
    animate_single(
        out03,
        "03_circulation_field_movie_5p65",
        "03 — Stage-2 circulation · P1 · 5.65 GHz · B+",
        s565,
        "On-objective Stage-2 frequency · SIMULATION Re[Hz e^{-iωt}]",
    )

    # 03b — off-grid failure at 5.75 GHz if present
    try:
        s575 = load_export(find_stem("stage2_best", 5.75))
        animate_pair(
            out03,
            "03_success_vs_offgrid_failure",
            "03 — Stage-2 · 5.65 GHz (on-grid) vs 5.75 GHz (audit hole) · P1 · B+",
            s565,
            s575,
            "5.65 GHz (on objective)",
            "5.75 GHz (off-grid failure)",
            "Same Stage-2 q_best · illustrates spectral hole motivating Stage-3 proposal",
            require_same_freq=False,
        )
    except FileNotFoundError as e:
        print("skip 5.75 comparison:", e)

    print("02/03 real field movies complete")


if __name__ == "__main__":
    main()
