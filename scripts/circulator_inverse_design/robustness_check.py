#!/usr/bin/env python3
"""Small robustness perturbations around the best design (characterization)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from circulator_core import B_PLUS, F_HZ, S_HI, S_LO, load_device, solve_six  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"


def main():
    opt = json.loads(Path(sys.argv[1]).read_text())
    grade = sys.argv[2] if len(sys.argv) > 2 else "M"
    s0 = np.asarray(opt["best"]["s"], float)
    P_ref = float(opt["P_ref"])
    dev = load_device(grade)
    base = solve_six(dev, s0, F_HZ, B_PLUS, P_ref=P_ref)
    rows = [{"name": "nominal", "J": base.J, "desired_avg": base.metrics["desired_avg"],
             "isolation_avg_dB": base.metrics["isolation_avg_dB"],
             "desired_min": base.metrics["desired_min"]}]
    # density scale ±5%
    for fac, name in [(0.95, "s*0.95"), (1.05, "s*1.05")]:
        s = np.clip(s0 * fac, S_LO, S_HI)
        bun = solve_six(dev, s, F_HZ, B_PLUS, P_ref=P_ref)
        rows.append({"name": name, "J": bun.J, "desired_avg": bun.metrics["desired_avg"],
                     "isolation_avg_dB": bun.metrics["isolation_avg_dB"],
                     "desired_min": bun.metrics["desired_min"]})
    # frequency ±50 MHz
    for f, name in [(3.80e9, "f=3.80"), (3.90e9, "f=3.90")]:
        bun = solve_six(dev, s0, f, B_PLUS, P_ref=P_ref)
        rows.append({"name": name, "J": bun.J, "desired_avg": bun.metrics["desired_avg"],
                     "isolation_avg_dB": bun.metrics["isolation_avg_dB"],
                     "desired_min": bun.metrics["desired_min"]})
    # B ±10%
    for b, name in [(0.045, "B=0.045"), (0.055, "B=0.055")]:
        bun = solve_six(dev, s0, F_HZ, b, P_ref=P_ref)
        rows.append({"name": name, "J": bun.J, "desired_avg": bun.metrics["desired_avg"],
                     "isolation_avg_dB": bun.metrics["isolation_avg_dB"],
                     "desired_min": bun.metrics["desired_min"]})
    # a few random per-rod 2% noise (break symmetry slightly)
    rng = np.random.default_rng(1)
    for k in range(3):
        s = np.clip(s0 * (1.0 + 0.02 * rng.standard_normal(len(s0))), S_LO, S_HI)
        bun = solve_six(dev, s, F_HZ, B_PLUS, P_ref=P_ref)
        rows.append({"name": f"noise2pct_{k}", "J": bun.J, "desired_avg": bun.metrics["desired_avg"],
                     "isolation_avg_dB": bun.metrics["isolation_avg_dB"],
                     "desired_min": bun.metrics["desired_min"]})
    out = {"grade": grade, "nominal_J": base.J, "rows": rows}
    path = OUT / f"robustness_{Path(sys.argv[1]).stem}_{grade}.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE", path.name, flush=True)
    for r in rows:
        print(r, flush=True)


if __name__ == "__main__":
    main()
