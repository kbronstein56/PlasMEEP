#!/usr/bin/env python3
"""Shared helpers for PlasMEEP presentation media (no FEM/Meep solves)."""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Paths: viz worktree writes locally; production outputs are read-only.
VIZ_ROOT = Path(__file__).resolve().parents[2]
PROD_ROOT = Path("/home/daq_user/PlasMEEP")
PROD_CAMPAIGN = (
    PROD_ROOT / "outputs/circulator_inverse_design/jesse_method_campaign"
)
PROD_PRES_FIGS = (
    PROD_ROOT / "outputs/validation/fem_meep_validation/presentation_figures"
)
OUT_ROOT = VIZ_ROOT / "outputs/presentation_media"
SCRATCH = OUT_ROOT / "_scratch"

W, H = 1920, 1080
FPS = 30
DPI = 100  # 1920x1080 at figsize=(19.2, 10.8)

# Presentation palette (light, high-contrast labels)
C_BG = "#f7f8fb"
C_INK = "#1a1a1a"
C_ACCENT = "#0b5fff"
C_GOOD = "#1b7f3a"
C_WARN = "#c27803"
C_BAD = "#b42318"
C_MUTED = "#5c6778"
C_SCHEMATIC = "#7a3e00"
C_DATA = "#0b5fff"
C_SIM = "#5b2c6f"

mpl.rcParams.update(
    {
        "figure.facecolor": C_BG,
        "axes.facecolor": "#ffffff",
        "savefig.facecolor": C_BG,
        "font.size": 18,
        "axes.titlesize": 28,
        "axes.labelsize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 16,
        "axes.linewidth": 1.4,
        "lines.linewidth": 2.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    }
)


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def load_json(path: Path):
    return json.loads(path.read_text())


def stamp_kind(ax, kind: str, extra: str = "") -> None:
    """Large corner badge: SCHEMATIC / DATA / ANALYTIC / SIMULATION."""
    colors = {
        "SCHEMATIC": C_SCHEMATIC,
        "DATA": C_DATA,
        "ANALYTIC": C_SIM,
        "SIMULATION": C_SIM,
        "MIXED": C_WARN,
    }
    color = colors.get(kind.upper(), C_MUTED)
    label = kind.upper() if not extra else f"{kind.upper()} · {extra}"
    ax.text(
        0.01,
        0.985,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.35", fc=color, ec="none", alpha=0.95),
        zorder=20,
    )


def new_frame(title: str | None = None):
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    ax = fig.add_axes([0.07, 0.10, 0.86, 0.80])
    if title:
        fig.suptitle(title, fontsize=30, fontweight="bold", color=C_INK, y=0.965)
    return fig, ax


def save_frame(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def encode_mp4(frame_glob: str, mp4_path: Path, fps: int = FPS) -> None:
    """Encode H.264 MP4, no audio. Requires ffmpeg."""
    mp4_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(fps),
        "-pattern_type",
        "glob",
        "-i",
        frame_glob,
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "20",
        "-movflags",
        "+faststart",
        "-an",
        str(mp4_path),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def frames_to_gif(frame_paths: list[Path], gif_path: Path, fps: int = 10) -> None:
    """Small GIF preview via ffmpeg (palette)."""
    if not frame_paths:
        return
    gif_path.parent.mkdir(parents=True, exist_ok=True)
    lst = gif_path.with_suffix(".txt")
    # subsample for light GIF
    step = max(1, len(frame_paths) // 40)
    chosen = frame_paths[::step]
    if chosen[-1] != frame_paths[-1]:
        chosen.append(frame_paths[-1])
    tmp = gif_path.parent / "_gif_frames"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    for i, src in enumerate(chosen):
        shutil.copy(src, tmp / f"g_{i:04d}.png")
    pal = tmp / "palette.png"
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            str(tmp / "g_%04d.png"),
            "-vf",
            "scale=960:-1:flags=lanczos,palettegen",
            str(pal),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            str(tmp / "g_%04d.png"),
            "-i",
            str(pal),
            "-lavfi",
            "scale=960:-1:flags=lanczos[x];[x][1:v]paletteuse",
            "-loop",
            "0",
            str(gif_path),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    shutil.rmtree(tmp, ignore_errors=True)
    lst.unlink(missing_ok=True)


def write_poster(frame_path: Path, poster_path: Path) -> None:
    shutil.copy(frame_path, poster_path)


def clear_frames(frame_dir: Path) -> None:
    if frame_dir.exists():
        shutil.rmtree(frame_dir)
    frame_dir.mkdir(parents=True, exist_ok=True)


def accepted_stage2_records(history: list[dict]) -> list[dict]:
    """Real Stage-2 checkpoints: start_eval + accepted_iter (and keep last new_best per nit)."""
    out = []
    for rec in history:
        kind = rec.get("kind")
        if kind in ("start_eval", "accepted_iter"):
            out.append(rec)
    return out


def best_so_far_curve(accepted: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    nits = []
    best = []
    jbest = None
    for rec in accepted:
        j = float(rec["J_BB"])
        jbest = j if jbest is None else max(jbest, j)
        nits.append(int(rec["iteration"]))
        best.append(jbest)
    return np.asarray(nits), np.asarray(best)
