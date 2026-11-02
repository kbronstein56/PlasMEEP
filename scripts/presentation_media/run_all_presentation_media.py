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
]


def main():
    py = sys.executable
    for name in SCRIPTS:
        print("=" * 60, flush=True)
        print("RUN", name, flush=True)
        subprocess.run([py, str(HERE / name)], check=True)

    # 02/03: real field movies if production exports exist; else stubs only
    from pathlib import Path

    field_index = Path(
        "/home/daq_user/PlasMEEP/outputs/circulator_inverse_design/"
        "jesse_method_campaign/presentation_fields/INDEX.json"
    )
    if field_index.exists():
        name = "render_02_03_real_field_movies.py"
    else:
        name = "render_02_03_field_movie_stubs.py"
    print("=" * 60, flush=True)
    print("RUN", name, flush=True)
    subprocess.run([py, str(HERE / name)], check=True)
    print("ALL DONE")


if __name__ == "__main__":
    main()
