#!/usr/bin/env python3
"""Stage-4 REAL FEM field movies — slow, loopable, with GIF + PNG poster.

Animates Re[Hz exp(-i ω t)] from production presentation_fields exports.
Target duration ~13 s @ 30 fps (2 visual periods + holds).
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
    FIELD_FRAMES_PER_PERIOD,
    FIELD_GIF_FPS,
    FIELD_HOLD_FRAMES,
    FIELD_N_PERIODS,
    FPS,
    OUT_ROOT,
    PROD_CAMPAIGN,
    clear_frames,
    encode_mp4,
    field_phase_schedule,
    frames_to_gif,
    save_frame,
    stamp_kind,
    write_poster,
)

FIELD_ROOT = PROD_CAMPAIGN / "presentation_fields"


def load_export(stem: str):
    d = FIELD_ROOT / stem
    if not d.exists():
        cands = sorted(FIELD_ROOT.glob(stem + "*")) if FIELD_ROOT.exists() else []
        if not cands:
            raise FileNotFoundError(f"missing field export: {d}")
        d = cands[0]
    meta = json.loads((d / "meta.json").read_text())
    x = np.load(d / "x.npy")
    y = np.load(d / "y.npy")
    Hz = np.load(d / "Hz_real.npy") + 1j * np.load(d / "Hz_imag.npy")
    return d, meta, x, y, Hz


def find_stem(design_tag: str, f_ghz: float) -> str:
    ftag = f"{f_ghz:.4f}".replace(".", "p")
    stem = f"{design_tag}_P1_{ftag}GHz_Bplus"
    if (FIELD_ROOT / stem).exists():
        return stem
    for p in sorted(FIELD_ROOT.glob(f"{design_tag}_P1_*GHz_Bplus")):
        meta = json.loads((p / "meta.json").read_text())
        if abs(float(meta["f_GHz"]) - f_ghz) < 1e-6:
            return p.name
    raise FileNotFoundError(stem)


def field_extent(x, y):
    return [float(x.min()), float(x.max()), float(y.min()), float(y.max())]


def metric_box(meta):
    m = meta.get("metrics_te0") or {}
    return (
        f"P{meta['source_port']}  B={meta['B_tesla']:+.3f} T\n"
        f"f={meta['f_GHz']:.5g} GHz\n"
        f"IL={m.get('IL_avg_dB', float('nan')):.1f} dB  "
        f"iso={m.get('isolation_avg_dB', float('nan')):.1f} dB"
    )


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
    shared_scale: bool = True,
):
    clear_frames(out_dir / "frames")
    _, meta_l, x_l, y_l, Hz_l = left
    _, meta_r, x_r, y_r, Hz_r = right
    if require_same_freq:
        assert abs(meta_l["f_GHz"] - meta_r["f_GHz"]) < 1e-9

    if shared_scale:
        vmax = float(
            np.nanpercentile(
                np.concatenate([np.abs(Hz_l).ravel(), np.abs(Hz_r).ravel()]), 99.5
            )
        )
        vmax = max(vmax, 1e-30)
        norms = (
            TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
            TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
        )
        scale_note = f"shared |Hz| p99.5 = {vmax:.3e}"
    else:
        vl = max(float(np.nanpercentile(np.abs(Hz_l), 99.5)), 1e-30)
        vr = max(float(np.nanpercentile(np.abs(Hz_r), 99.5)), 1e-30)
        norms = (
            TwoSlopeNorm(vmin=-vl, vcenter=0.0, vmax=vl),
            TwoSlopeNorm(vmin=-vr, vcenter=0.0, vmax=vr),
        )
        scale_note = f"independent scales L={vl:.3e} R={vr:.3e}"

    phases = field_phase_schedule(
        FIELD_N_PERIODS, FIELD_FRAMES_PER_PERIOD, FIELD_HOLD_FRAMES
    )
    frame_paths = []
    for i, phase in enumerate(phases):
        Fl = np.real(Hz_l * np.exp(-1j * phase))
        Fr = np.real(Hz_r * np.exp(-1j * phase))
        fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
        fig.suptitle(title, fontsize=24, fontweight="bold", color=C_INK, y=0.97)
        ax0 = fig.add_axes([0.04, 0.12, 0.44, 0.78])
        ax1 = fig.add_axes([0.52, 0.12, 0.44, 0.78])
        im = None
        for ax, F, x, y, ttl, meta, norm in (
            (ax0, Fl, x_l, y_l, left_title, meta_l, norms[0]),
            (ax1, Fr, x_r, y_r, right_title, meta_r, norms[1]),
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
                metric_box(meta) + f"\nφ={phase / np.pi:.2f}π",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=11,
                color=C_MUTED,
                bbox=dict(
                    boxstyle="round,pad=0.25", fc="white", ec="#cbd5e1", alpha=0.9
                ),
            )
        fig.colorbar(im, ax=[ax0, ax1], fraction=0.02, pad=0.02, label="Re[Hz] (rel.)")
        fig.text(
            0.5,
            0.02,
            f"{subtitle}  ·  {scale_note}  ·  ~{len(phases)/FPS:.1f}s @ {FPS}fps",
            ha="center",
            fontsize=12,
            color=C_MUTED,
        )
        fp = out_dir / "frames" / f"frame_{i:04d}.png"
        save_frame(fig, fp)
        frame_paths.append(fp)

    mp4 = out_dir / f"{stem}.mp4"
    encode_mp4(str(out_dir / "frames" / "frame_*.png"), mp4, fps=FPS)
    write_poster(frame_paths[int(FIELD_HOLD_FRAMES)], out_dir / f"{stem}_poster.png")
    frames_to_gif(frame_paths, out_dir / f"{stem}_preview.gif", fps=FIELD_GIF_FPS)
    (out_dir / f"{stem}_PROVENANCE.txt").write_text(
        json.dumps(
            {
                "left": meta_l,
                "right": meta_r,
                "animation": "Re[Hz * exp(-1j*omega*t)]",
                "n_frames": len(phases),
                "duration_s": len(phases) / FPS,
                "periods": FIELD_N_PERIODS,
                "hold_frames": FIELD_HOLD_FRAMES,
                "shared_color_scale": shared_scale,
                "field_root": str(FIELD_ROOT),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {mp4} ({len(phases)/FPS:.1f}s)", flush=True)
    return mp4


def main():
    if not FIELD_ROOT.exists():
        raise SystemExit(f"missing {FIELD_ROOT}")

    jobs = []
    specs = [
        (
            "20_stage4_notch_5p8125",
            "20 — 5.8125 GHz notch · Stage-3 vs Stage-4 · P1 · B+",
            "stage3_best",
            "stage4_best",
            5.8125,
            "Stage-3",
            "Stage-4",
            "REAL FEM · shared color scale · A2 −12.02 → −4.48",
        ),
        (
            "21_stage4_notch_4p775",
            "21 — 4.775 GHz · Stage-3 vs Stage-4 · P1 · B+",
            "stage3_best",
            "stage4_best",
            4.775,
            "Stage-3",
            "Stage-4",
            "REAL FEM · shared color scale · iso −8.4 → +9.4 dB",
        ),
        (
            "22_stage4_new_hole_5p08125",
            "22 — NEW dense-audit worst @ 5.08125 GHz · Stage-3 vs Stage-4 · P1 · B+",
            "stage3_best",
            "stage4_best",
            5.08125,
            "Stage-3",
            "Stage-4 (true worst)",
            "REAL FEM · shared color scale · Stage-4 A2≈−11.04",
        ),
        (
            "24_stage4_new_hole_4p725",
            "24 — NEW off-grid hole @ 4.725 GHz · Stage-3 vs Stage-4 · P1 · B+",
            "stage3_best",
            "stage4_best",
            4.725,
            "Stage-3",
            "Stage-4",
            "REAL FEM · shared color scale · Stage-4 A2≈−9.16",
        ),
    ]
    for stem, title, left_tag, right_tag, f, lt, rt, sub in specs:
        try:
            jobs.append(
                (
                    stem,
                    title,
                    find_stem(left_tag, f),
                    find_stem(right_tag, f),
                    lt,
                    rt,
                    sub,
                    True,
                )
            )
        except FileNotFoundError as e:
            print(f"skip {stem}: {e}", flush=True)

    if not jobs:
        raise SystemExit("No Stage-4 field jobs ready — wait for exports.")

    for stem, title, left_stem, right_stem, lt, rt, sub, shared in jobs:
        out = OUT_ROOT / stem
        animate_pair(
            out,
            stem,
            title,
            load_export(left_stem),
            load_export(right_stem),
            lt,
            rt,
            sub,
            require_same_freq=True,
            shared_scale=shared,
        )

    print("Stage-4 field movies complete", flush=True)


if __name__ == "__main__":
    main()
