#!/usr/bin/env python3
"""Evaluate a design vector on a chosen mesh: 6x6 power, B-reversal, frequency sweep."""
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
from circulator_core import B_PLUS, F_HZ, load_device, solve_six  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"


def evaluate(s, grade, b_tesla, f_hz, P_ref=None, tag="eval"):
    dev = load_device(grade)
    bun = solve_six(dev, s, f_hz, b_tesla, P_ref=P_ref)
    out = {
        "tag": tag,
        "grade": grade,
        "B_T": b_tesla,
        "f_Hz": f_hz,
        "P_ref": bun.P_ref,
        "J": bun.J,
        "metrics": bun.metrics,
        "power": bun.power.tolist(),
        "factor_s": bun.factor_s,
        "solve_s": bun.solve_s,
        "s_min": float(np.min(s)),
        "s_max": float(np.max(s)),
        "s_mean": float(np.mean(s)),
    }
    path = OUT / f"{tag}.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE", path.name, "J", bun.J, "des", bun.metrics["desired_avg"], flush=True)
    return out


def main():
    # usage: evaluate_design.py <opt_json> [grade] [mode]
    # mode: point | breeversal | sweep | all
    opt_path = Path(sys.argv[1])
    grade = sys.argv[2] if len(sys.argv) > 2 else "F"
    mode = sys.argv[3] if len(sys.argv) > 3 else "all"
    opt = json.loads(opt_path.read_text())
    s = np.asarray(opt["best"]["s"], float)
    P_ref = float(opt["P_ref"])
    stem = opt_path.stem
    results = {}
    if mode in ("point", "all"):
        results["plus"] = evaluate(s, grade, B_PLUS, F_HZ, P_ref=P_ref, tag=f"{stem}_eval_{grade}_Bp0.05")
    if mode in ("breeversal", "all"):
        results["minus"] = evaluate(s, grade, -B_PLUS, F_HZ, P_ref=P_ref, tag=f"{stem}_eval_{grade}_Bm0.05")
        results["zero"] = evaluate(s, grade, 0.0, F_HZ, P_ref=P_ref, tag=f"{stem}_eval_{grade}_B0")
        # check reversal: at -B, reverse port should exceed desired for each source
        Pm = np.array(results["minus"]["power"], float)
        Pp = np.array(results["plus"]["power"], float)
        ok = True
        for j in range(6):
            des = (j + 1) % 6
            rev = (j - 1) % 6
            if not (Pm[rev, j] > Pm[des, j] and Pp[des, j] > Pp[rev, j]):
                ok = False
        results["reversal_pass"] = bool(ok)
        print("B_REVERSAL", ok, flush=True)
    if mode in ("sweep", "all"):
        sweep = []
        for f in [3.70e9, 3.75e9, 3.80e9, 3.85e9, 3.90e9, 3.95e9, 4.00e9]:
            r = evaluate(s, grade, B_PLUS, f, P_ref=P_ref, tag=f"{stem}_sweep_{grade}_f{f/1e9:.2f}")
            sweep.append(r)
        results["sweep"] = [
            {
                "f_Hz": r["f_Hz"],
                "J": r["J"],
                "desired_avg": r["metrics"]["desired_avg"],
                "reverse_avg": r["metrics"]["reverse_avg"],
                "through_avg": r["metrics"]["through_avg"],
                "isolation_avg_dB": r["metrics"]["isolation_avg_dB"],
            }
            for r in sweep
        ]
    (OUT / f"{stem}_eval_{grade}_summary.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
