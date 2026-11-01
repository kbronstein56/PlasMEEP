#!/usr/bin/env python3
"""Stubs for 02/03 real field movies.

Intentional non-invention: these scripts document the required complex-Hz
exports and the correct animation formula Re[Hz exp(-iωt)].

They do NOT synthesize fake device fields. Running this writes a placeholder
poster + QC note only.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from media_common import C_INK, C_SCHEMATIC, OUT_ROOT, save_frame, stamp_kind  # noqa: E402

NEEDED = """
REQUIRED EXPORTS (not present yet — do not invent):

02_uniform_vs_optimized_field_movie
  - Complex Hz on a common XY grid for:
      * uniform start design, source P1, selected f (e.g. 4.1 and/or final worst)
      * Stage-2 best design (q_best from opt_broadband_stage2_A2_C.json), same f/src
  - Files ideally: Hz_real.npy, Hz_imag.npy, x.npy, y.npy, meta.json
    meta.json: {f_GHz, B_tesla, source_port, design_tag, units, solver, mesh_grade}

03_circulation_field_movie (B+ sense)
  - Same complex Hz export for Stage-2 best at a well-behaved frequency
  - Optional companion: B reversed (B_MINUS) for Onsager/circulation sense check

ANIMATION RULE (when exports arrive):
  frame(t) = Re[ (Hz_real + 1j*Hz_imag) * exp(-1j * omega * t) ]
  DO NOT animate |Hz| as oscillating.
  Label every frame: SIMULATION · Re[Hz e^{-iωt}] · f=… · B=± …
  Keep port numbering P1…P6 and desired=(j+1)%6 for B_PLUS consistent.
"""


def placeholder(stem: str, title: str, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    ax = fig.add_axes([0.08, 0.1, 0.84, 0.8])
    ax.axis("off")
    fig.suptitle(title, fontsize=28, fontweight="bold", color=C_INK, y=0.95)
    stamp_kind(ax, "SCHEMATIC", "AWAITING FIELD EXPORTS")
    ax.text(
        0.5, 0.55,
        "No fabricated EM fields.\nWaiting for complex Hz exports from Stage-2 / production solvers.",
        ha="center", va="center", fontsize=22, color=C_SCHEMATIC,
        transform=ax.transAxes,
        bbox=dict(boxstyle="round,pad=0.6", fc="#fff7ed", ec=C_SCHEMATIC),
    )
    ax.text(0.02, 0.02, NEEDED, fontsize=11, family="monospace", color="#334155",
            transform=ax.transAxes, va="bottom")
    save_frame(fig, out_dir / f"{stem}_poster.png")
    (out_dir / f"{stem}_NEEDS_EXPORTS.txt").write_text(NEEDED)


def main():
    placeholder(
        "02_uniform_vs_optimized_field_movie",
        "02 — Uniform vs optimized field movie (STUB)",
        OUT_ROOT / "02_field_movies_stubs",
    )
    placeholder(
        "03_circulation_field_movie",
        "03 — Circulation field movie (STUB)",
        OUT_ROOT / "03_field_movies_stubs",
    )
    print("wrote field-movie stubs (no fabricated fields)")


if __name__ == "__main__":
    main()
