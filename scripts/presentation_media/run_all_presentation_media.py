#!/usr/bin/env python3
"""Build all presentation media with light CPU (no FEM/Meep)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "render_01_device_circulation_overview.py",
    "render_stage2_progress.py",
    "render_05_broadband_frequency_sweep.py",
    "render_06_adaptive_frequency_enrichment.py",
    "render_07_TE0_TE1_modes.py",
    "render_08_validation_ladder.py",
    "render_02_03_field_movie_stubs.py",
]


def main():
    py = sys.executable
    for name in SCRIPTS:
        print("=" * 60, flush=True)
        print("RUN", name, flush=True)
        subprocess.run([py, str(HERE / name)], check=True)
    print("ALL DONE")


if __name__ == "__main__":
    main()
