#!/usr/bin/env python3
"""Overlay Meep prism solids against the old FEM edge-tube mask."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "horn_localization"
doc = json.loads((OUT / "meep_horn_polygons.json").read_text())
a = float(doc["a_m"])

fig, ax = plt.subplots(figsize=(8.2, 7.2))
area = 0.0
for horn in doc["horns"]:
    for wall in horn["walls"]:
        q = np.asarray(wall["vertices_m"], dtype=float)
        ax.add_patch(Polygon(q, closed=True, facecolor="C0", edgecolor="k", alpha=0.85, lw=0.4))
        x, y = q[:, 0], q[:, 1]
        area += 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    src = np.asarray(horn["source_center_a"], dtype=float) * a
    mon = np.asarray(horn["monitor_center_a"], dtype=float) * a
    throat = np.asarray(horn["throat_a"], dtype=float) * a
    t = np.asarray(horn["tangent"], dtype=float)
    span = 0.5 * float(horn["span_a"]) * a
    for c, color in ((src, "C3"), (mon, "C2"), (throat, "C1")):
        ax.plot([c[0] - span * t[0], c[0] + span * t[0]], [c[1] - span * t[1], c[1] + span * t[1]], color=color, lw=1.2)

ax.set_aspect("equal")
ax.autoscale()
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_title("Meep Add_Prism PEC walls (exact vertices)\nred source, green monitor, orange throat")
fig.tight_layout()
fig.savefig(OUT / "meep_horn_polygons.png", dpi=140)
print("shoelace area m^2", area)
print("wrote", OUT / "meep_horn_polygons.png")
