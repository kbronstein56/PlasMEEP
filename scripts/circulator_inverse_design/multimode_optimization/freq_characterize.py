#!/usr/bin/env python3
"""Phase 9: frequency characterization of best multimode design.

TE1 cutoff ≈ c/(2W) with W=2a. Below cutoff TE1 is evanescent and must not
be counted as outgoing propagating power in engineering FoMs.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from multimode_opt_core import B_PLUS, OUT, dump_json, eval_design, te0_metrics_from_S

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from modes import mode_constants  # noqa: E402

C0 = 299792458.0
FREQS_GHZ = [3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00]


def te1_cutoff_hz():
    W = float(sc.clear_width) * sc.a  # meters
    return C0 / (2.0 * W)


def classify_modes(f_hz: float):
    k0 = 2.0 * np.pi * (f_hz * sc.a / C0)
    modes = mode_constants(k0, float(sc.clear_width), n_modes=4)
    out = []
    for m, md in enumerate(modes):
        beta = complex(md["beta"])
        prop = abs(beta.imag) < 1e-8 and beta.real > 1e-8
        out.append(
            {
                "m": m,
                "beta": [float(beta.real), float(beta.imag)],
                "propagating": bool(prop),
                "kc": float(md.get("kc", m * np.pi / float(sc.clear_width))),
            }
        )
    n_prop = sum(1 for x in out if x["propagating"])
    return out, n_prop


def metrics_propagating_only(S, n_prop_per_port: int):
    """TE0 engineering metrics; TE1 counted only if propagating."""
    met = te0_metrics_from_S(S)
    if n_prop_per_port < 2:
        # zero out TE1 contribution for reporting
        for r in met["per_source"]:
            r["TE1_total_propagating"] = 0.0
            r["TE1_total_raw"] = r["TE1_total"]
            r["TE1_total"] = 0.0
        met["TE1_avg"] = 0.0
        met["TE1_note"] = "TE1 below cutoff — not counted as propagating power"
    else:
        for r in met["per_source"]:
            r["TE1_total_propagating"] = r["TE1_total"]
        met["TE1_note"] = "TE1 propagating — included"
    return met


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--opt-json", required=True)
    ap.add_argument("--grade", default="C")
    ap.add_argument("--label", default="best")
    ap.add_argument("--freqs", default=",".join(str(f) for f in FREQS_GHZ))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    opt = json.loads(Path(args.opt_json).read_text())
    s = np.asarray(opt["best"]["s"], float)
    fc = te1_cutoff_hz()
    freqs = [float(x) * 1e9 for x in args.freqs.split(",")]
    t0 = time.time()
    rows = []
    for f in freqs:
        modes, n_prop = classify_modes(f)
        near = abs(f - fc) / fc < 0.01
        print(
            f"f={f/1e9:.3f} GHz n_prop/port≈{n_prop} near_cutoff={near} ...",
            flush=True,
        )
        try:
            bun = eval_design(args.grade, s, f, B_PLUS, with_absorption=True)
            met = metrics_propagating_only(bun["S"], n_prop)
            err = None
            S = bun["S"]
        except Exception as e:
            met = None
            err = repr(e)
            S = None
            print(f"  FAILED: {e}", flush=True)
        rows.append(
            {
                "f_Hz": f,
                "f_GHz": f / 1e9,
                "TE1_cutoff_Hz": fc,
                "TE1_cutoff_GHz": fc / 1e9,
                "modes": modes,
                "n_prop_modes_in_basis": n_prop,
                "near_cutoff": near,
                "beta_TE1": modes[1]["beta"] if len(modes) > 1 else None,
                "metrics": met,
                "error": err,
            }
        )
        if S is not None:
            np.savez_compressed(
                OUT / f"S12_freq_{args.label}_{args.grade}_f{f/1e9:.2f}.npz",
                S=S,
                f_Hz=f,
                s=s,
            )

    report = {
        "opt_json": str(args.opt_json),
        "label": args.label,
        "grade": args.grade,
        "TE1_cutoff_Hz": fc,
        "TE1_cutoff_GHz": fc / 1e9,
        "rows": rows,
        "elapsed_s": time.time() - t0,
        "broadband_recommendation": None,
    }
    # recommendation
    above = [r for r in rows if r["f_Hz"] > fc * 1.01 and r["metrics"]]
    below = [r for r in rows if r["f_Hz"] < fc * 0.99 and r["metrics"]]
    near_rows = [r for r in rows if r["near_cutoff"]]
    rec = {
        "option_A": (
            "Optimize entirely above TE1 cutoff "
            f"(safely > {fc/1e9:.4f} GHz), e.g. 3.80–4.00 GHz if TE1 stays propagating "
            "and the objective always uses 12 channels."
        ),
        "option_B": (
            "Allow crossing TE1 cutoff with a frequency-dependent propagating-channel "
            "objective (count TE1 only when β_TE1 is real). Requires a new frozen "
            "derivation before broadband opt."
        ),
        "observed": {
            "n_above": len(above),
            "n_below": len(below),
            "n_near": len(near_rows),
            "note": (
                "At 3.85 GHz TE1 is barely above cutoff; broadband that dips to 3.75 "
                "needs cutoff-aware treatment."
            ),
        },
        "recommendation": "NEEDS_CUTOFF_DERIVATION",
        "rationale": (
            "Single-frequency design lives just above TE1 cutoff. Any band that "
            "approaches or crosses ~3.747 GHz changes the number of propagating "
            "channels; do not reuse a fixed 12-channel objective without an explicit "
            "frozen cutoff treatment."
        ),
    }
    report["broadband_recommendation"] = rec
    dump_json(OUT / f"freq_char_{args.label}_{args.grade}.json", report)

    lines = [
        "# MULTIMODE_FREQ_CHARACTERIZATION",
        "",
        f"Design `{args.opt_json}`, grade={args.grade}.",
        f"TE1 cutoff (analytic parallel-plate): **{fc/1e9:.4f} GHz** "
        f"(W = clear_width × a = {sc.clear_width} × {sc.a} m).",
        "",
        "| f (GHz) | n_prop | β_TE1 (re, im) | T_des | IL | iso | RL | TE1 | note |",
        "|--------:|-------:|----------------|------:|---:|----:|---:|----:|------|",
    ]
    for r in rows:
        b1 = r["beta_TE1"]
        bstr = f"({b1[0]:.4g}, {b1[1]:.4g})" if b1 else "—"
        if r["metrics"] is None:
            lines.append(
                f"| {r['f_GHz']:.2f} | {r['n_prop_modes_in_basis']} | {bstr} | — | — | — | — | — | FAIL |"
            )
            continue
        m = r["metrics"]
        note = "near cutoff" if r["near_cutoff"] else ("TE1 evan." if r["n_prop_modes_in_basis"] < 2 else "ok")
        lines.append(
            f"| {r['f_GHz']:.2f} | {r['n_prop_modes_in_basis']} | {bstr} | "
            f"{m['T_des_avg']:.4f} | {m['IL_avg_dB']:.2f} | {m['isolation_avg_dB']:.2f} | "
            f"{m['RL_avg_dB']:.2f} | {m['TE1_avg']:.4f} | {note} |"
        )
    lines += [
        "",
        "## Broadband decision (Phase 10)",
        "",
        f"**BROADBAND_READY: {rec['recommendation']}**",
        "",
        rec["rationale"],
        "",
        "- Option A: " + rec["option_A"],
        "- Option B: " + rec["option_B"],
        "",
        "Broadband optimization itself is deferred to the next campaign.",
        "",
    ]
    (OUT / "MULTIMODE_FREQ_CHARACTERIZATION.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"TE1_cutoff_GHz": fc / 1e9, "recommendation": rec["recommendation"], "elapsed_s": report["elapsed_s"]}, indent=2))


if __name__ == "__main__":
    main()
