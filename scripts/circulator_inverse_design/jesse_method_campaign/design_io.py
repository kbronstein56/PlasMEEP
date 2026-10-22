#!/usr/bin/env python3
"""Load design vectors from multimode optimization JSON artifacts."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OPT = ROOT / "outputs" / "circulator_inverse_design" / "multimode_optimization"


def load_design_q(design: str) -> np.ndarray:
    """Prefer best['q'], then q_best / q / q0."""
    p = OPT / f"{design}.json"
    d = json.loads(p.read_text())
    if isinstance(d.get("best"), dict) and "q" in d["best"]:
        return np.asarray(d["best"]["q"], float)
    for key in ("q_best", "q", "q0"):
        if key in d:
            return np.asarray(d[key], float)
    raise KeyError(f"no design vector in {p}")
