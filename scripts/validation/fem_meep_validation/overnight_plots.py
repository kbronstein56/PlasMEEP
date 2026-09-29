#!/usr/bin/env python3
"""Port-power and Meep-FEM error plots for the overnight B=0 series."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parents[3] / "outputs" / "validation" / "fem_meep_validation" / "overnight_b0"


def db(a, b):
    return 10.0 * np.log10(np.asarray(a) / np.asarray(b))


def load_meep(tag):
    rows = []
    for path in sorted(OUT.glob(f"meep_{tag}_ppc*.json")):
        data = json.loads(path.read_text())
        rows.append(data)
    rows.sort(key=lambda r: r["points_per_cm"])
    return rows


def fit(dx, y):
    # linear in dx and linear in dx^2, intercept at dx=0
    x = np.asarray(dx, dtype=float)
    y = np.asarray(y, dtype=float)
    a1 = np.polyfit(x, y, 1)
    a2 = np.polyfit(x**2, y, 1)
    pred1 = np.polyval(a1, x)
    pred2 = np.polyval(a2, x**2)
    return {
        "linear_dx_intercept": float(a1[1]),
        "linear_dx_slope": float(a1[0]),
        "linear_dx_rms": float(np.sqrt(np.mean((pred1 - y) ** 2))),
        "linear_dx2_intercept": float(a2[1]),
        "linear_dx2_slope": float(a2[0]),
        "linear_dx2_rms": float(np.sqrt(np.mean((pred2 - y) ** 2))),
    }


def series(tag, fem_level):
    fem_rows = json.loads((OUT / f"fem_{tag}.json").read_text())
    fem = next(r for r in fem_rows if r["level"] == fem_level)
    meeps = load_meep(tag)
    return fem, meeps


def plot(tag, fem_level, ports=(1, 2, 3)):
    fem, meeps = series(tag, fem_level)
    if len(meeps) < 2:
        print(tag, "not enough Meep points")
        return
    dx = np.array([10.0 / m["points_per_cm"] for m in meeps])
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    names = [f"P{i+1}" for i in ports]
    for i, name in zip(ports, names):
        y = np.array([m["normalized_power"][i] for m in meeps])
        axes[0].plot(dx, y, "o-", label=f"Meep {name}")
        axes[0].axhline(fem["normalized"][i], ls="--", label=f"FEM {name}")
        err = db(y, fem["normalized"][i])
        axes[1].plot(dx, err, "o-", label=name)
        if len(meeps) >= 3:
            fits = fit(dx, err)
            print(tag, name, "err_dB", [round(float(v), 4) for v in err], fits)
    axes[0].set_xlabel("dx (mm)")
    axes[0].set_ylabel("normalized port power")
    axes[0].set_title(f"{tag}: power vs dx")
    axes[1].set_xlabel("dx (mm)")
    axes[1].set_ylabel("Meep − FEM (dB)")
    axes[1].axhline(0.0, color="k", lw=0.6)
    axes[1].set_title(f"{tag} vs {fem_level}")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / f"{tag}_vs_dx.png", dpi=140)
    print("wrote", OUT / f"{tag}_vs_dx.png")


def main():
    # The Meep JSON key is normalized_powers or powers. Detect it.
    sample = next(OUT.glob("meep_quartz_ppc*.json"), None)
    if sample:
        keys = json.loads(sample.read_text()).keys()
        print("meep keys", sorted(keys))
    plot("quartz", "FEM-VB")
    if (OUT / "fem_plasma.json").exists() and list(OUT.glob("meep_plasma_ppc*.json")):
        plot("plasma", "FEM-VC")


if __name__ == "__main__":
    main()
