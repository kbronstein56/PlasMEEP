#!/usr/bin/env python3
"""Sensitivity study of soft-min τ / α on fixed spectra (no FEM solves).

τ controls how closely soft-min ≈ worst-frequency J.
α controls mean-performance vs worst-frequency weighting.

Does NOT choose τ/α merely because FD differentiation passed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "broadband_method"
sys.path.insert(0, str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"))

from broadband_objective import softmin_bb  # noqa: E402

# Synthetic spectra representative of circulator FoM shapes over 4–6 GHz
SPECTRA = {
    "flat_good": np.full(11, 0.55),
    "flat_poor": np.full(11, 0.10),
    "one_sharp_dip": np.array([0.55, 0.54, 0.53, 0.10, 0.52, 0.55, 0.56, 0.54, 0.53, 0.55, 0.54]),
    "two_dips": np.array([0.50, 0.15, 0.48, 0.47, 0.12, 0.49, 0.50, 0.48, 0.47, 0.49, 0.50]),
    "slow_roll": np.linspace(0.60, 0.25, 11),
    "edge_dip": np.array([0.05, 0.50, 0.52, 0.53, 0.54, 0.55, 0.54, 0.53, 0.52, 0.51, 0.50]),
}

TAUS = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50]
ALPHAS = [0.0, 0.10, 0.25, 0.50, 1.0]


def analyze(Jf: np.ndarray, tau: float, alpha: float) -> dict:
    J_bb, w, J_soft, J_mean, p = softmin_bb(Jf, tau=tau, alpha=alpha)
    worst = float(np.min(Jf))
    # How close softmin is to worst (0 = exact; large = soft)
    soft_gap = J_soft - worst
    # Weight mass on the worst sample
    i_worst = int(np.argmin(Jf))
    return {
        "J_BB": J_bb,
        "J_softmin": J_soft,
        "J_mean": J_mean,
        "J_worst": worst,
        "softmin_minus_worst": soft_gap,
        "w_on_worst": float(w[i_worst]),
        "p_on_worst": float(p[i_worst]),
        "w_entropy": float(-np.sum(w * np.log(np.clip(w, 1e-30, None)))),
        "effective_n": float(np.exp(-np.sum(p * np.log(np.clip(p, 1e-30, None))))),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, Jf in SPECTRA.items():
        for tau in TAUS:
            for alpha in ALPHAS:
                r = analyze(Jf, tau, alpha)
                r.update({"spectrum": name, "tau": tau, "alpha": alpha})
                rows.append(r)

    # Recommendations from fixed spectra (not from FD pass)
    # Prefer: softmin tracks worst within ~0.05 on one_sharp_dip, while
    # retaining some mean signal (alpha in 0.1–0.25) and not collapsing to
    # single-point gradient (effective_n >~2 on flat spectra).
    rec = []
    for tau in TAUS:
        for alpha in ALPHAS:
            dip = analyze(SPECTRA["one_sharp_dip"], tau, alpha)
            flat = analyze(SPECTRA["flat_good"], tau, alpha)
            two = analyze(SPECTRA["two_dips"], tau, alpha)
            ok_track = dip["softmin_minus_worst"] <= 0.08
            ok_flat_spread = flat["effective_n"] >= 5.0  # near-uniform on flat
            ok_not_pure_mean = alpha < 0.9
            ok_not_pure_min = tau > 0.005 and dip["w_on_worst"] < 0.95
            score = (
                ok_track
                and ok_flat_spread
                and ok_not_pure_mean
                and (0.05 <= alpha <= 0.35)
                and (0.02 <= tau <= 0.10)
            )
            rec.append(
                {
                    "tau": tau,
                    "alpha": alpha,
                    "dip_soft_gap": dip["softmin_minus_worst"],
                    "dip_w_worst": dip["w_on_worst"],
                    "flat_eff_n": flat["effective_n"],
                    "two_w_worst": two["w_on_worst"],
                    "candidate": bool(score),
                }
            )

    candidates = [r for r in rec if r["candidate"]]
    # Prefer current prototype if it qualifies; else best by smallest dip gap
    preferred = None
    for r in candidates:
        if abs(r["tau"] - 0.05) < 1e-12 and abs(r["alpha"] - 0.25) < 1e-12:
            preferred = r
            break
    if preferred is None and candidates:
        preferred = min(candidates, key=lambda r: r["dip_soft_gap"])

    report = {
        "spectra": {k: v.tolist() for k, v in SPECTRA.items()},
        "tau_grid": TAUS,
        "alpha_grid": ALPHAS,
        "rows": rows,
        "candidate_summary": candidates,
        "preferred": preferred,
        "notes": {
            "tau": "smaller τ → soft-min closer to min_i J(f_i); larger τ → smoother average",
            "alpha": "α=0 → pure soft-min; α=1 → pure mean; α≈0.25 mixes robustness + average IL",
            "fd_check": "FD soft-min weight derivative rel~5.7e-11 is necessary but not sufficient",
        },
    }
    (OUT / "softmin_tau_alpha_sensitivity.json").write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# Soft-min τ / α sensitivity (fixed spectra)",
        "",
        "No FEM solves. Synthetic J(f) shapes over 11 samples (proxy for 4–6 GHz).",
        "",
        "## Meanings",
        "",
        "- **τ**: how closely soft-min approximates worst-frequency J. "
        "Small τ → near hard min; large τ → near mean.",
        "- **α**: weight on mean vs soft-min. α=0 pure soft-min; α=1 pure mean.",
        "",
        "FD differentiation PASS (rel ≈ 5.7e-11) is **not** a reason to freeze values.",
        "",
        "## Preferred from this study",
        "",
    ]
    if preferred:
        lines += [
            f"- τ = **{preferred['tau']}**, α = **{preferred['alpha']}**",
            f"- one-sharp-dip: softmin−worst = {preferred['dip_soft_gap']:.4f}, "
            f"w_worst = {preferred['dip_w_worst']:.3f}",
            f"- flat spectrum effective n ≈ {preferred['flat_eff_n']:.2f}",
            "",
            "Current prototype (τ=0.05, α=0.25) "
            + ("**qualifies**." if preferred["tau"] == 0.05 and preferred["alpha"] == 0.25 else "is overridden by this study."),
            "",
        ]
    else:
        lines += ["No candidate in the scored box; inspect JSON.", ""]

    lines += [
        "## Candidate table (scored box)",
        "",
        "| τ | α | dip soft−worst | dip w_worst | flat eff_n |",
        "|--:|--:|---------------:|------------:|-----------:|",
    ]
    for r in candidates:
        lines.append(
            f"| {r['tau']} | {r['alpha']} | {r['dip_soft_gap']:.4f} | "
            f"{r['dip_w_worst']:.3f} | {r['flat_eff_n']:.2f} |"
        )
    lines += [
        "",
        "Artifact: `softmin_tau_alpha_sensitivity.json`.",
    ]
    (OUT / "SOFTMIN_TAU_ALPHA_SENSITIVITY.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"preferred": preferred, "n_candidates": len(candidates)}, indent=2))


if __name__ == "__main__":
    main()
