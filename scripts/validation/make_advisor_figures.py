#!/usr/bin/env python3
"""Generate advisor-update figures from existing validation JSON (no Meep runs)."""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "validation" / "advisor_update"
PORTS = ROOT / "outputs" / "validation" / "ports"
RECIP = ROOT / "outputs" / "validation" / "reciprocity"
FARADAY32 = ROOT / "outputs" / "validation" / "faraday_res32.json"
FARADAY64 = ROOT / "outputs" / "validation" / "faraday_res64.json"


def _style():
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 200,
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "axes.grid": True,
            "grid.alpha": 0.3,
        }
    )


def fig_faraday():
    f32 = json.loads(FARADAY32.read_text())
    f64 = json.loads(FARADAY64.read_text())
    theory = f32["theory"]["+B"]["kappa_rad_per_a"]

    cases = ["B=0", "+B", "−B"]
    x = np.arange(len(cases))
    width = 0.28

    def vals(bundle):
        return [
            bundle["results"]["B=0"]["kappa_rad_per_a"],
            bundle["results"]["+B"]["kappa_rad_per_a"],
            bundle["results"]["-B"]["kappa_rad_per_a"],
        ]

    theory_vals = [0.0, theory, -theory]
    m32 = vals(f32)
    m64 = vals(f64)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    ax.axhline(0, color="k", lw=0.8)
    ax.bar(x - width, theory_vals, width, label="Theory κ", color="#4C78A8")
    ax.bar(x, m32, width, label="Measured res=32", color="#F58518")
    ax.bar(x + width, m64, width, label="Measured res=64", color="#54A24B")
    ax.set_xticks(x)
    ax.set_xticklabels(cases)
    ax.set_ylabel("κ (rad / a)")
    ax.set_title("Faraday rotation rate vs bias")
    ax.legend(loc="best")
    ax.text(
        0.02,
        0.02,
        "Note: measured κ sign is opposite theory convention;\n"
        "PASS gates use |κ| and B-reversal only.",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
    )

    ax = axes[1]
    err32 = abs(abs(m32[1]) - abs(theory)) / abs(theory) * 100
    err64 = abs(abs(m64[1]) - abs(theory)) / abs(theory) * 100
    ax.bar([0, 1], [err32, err64], color=["#F58518", "#54A24B"], width=0.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["res=32", "res=64"])
    ax.set_ylabel("|κ| relative error vs theory (%)")
    ax.set_title("Magnitude accuracy")
    for i, e in enumerate([err32, err64]):
        ax.text(i, e + 0.01, f"{e:.3f}%", ha="center", va="bottom", fontsize=10)
    ax.set_ylim(0, max(err32, err64) * 1.45 + 0.05)

    fig.suptitle(
        "Homogeneous magnetized-plasma Faraday benchmark (PlasMEEP GyrotropicDrude)",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    path = OUT / "fig1_faraday_benchmark.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_reciprocity():
    # Baseline overnight
    base = {
        32: json.loads((RECIP / "pair_P1P2_res32_rt40.json").read_text())[
            "reciprocity"
        ]["max_abs_diff_dB"],
        48: json.loads((RECIP / "pair_P1P2_res48_rt80.json").read_text())[
            "reciprocity"
        ]["max_abs_diff_dB"],
        64: json.loads((RECIP / "pair_P1P2_res64_rt80.json").read_text())[
            "reciprocity"
        ]["max_abs_diff_dB"],
        96: json.loads((RECIP / "pair_P1P2_res96_rt80.json").read_text())[
            "reciprocity"
        ]["max_abs_diff_dB"],
    }
    p2p3_base = json.loads((RECIP / "pair_P2P3_res32_rt40.json").read_text())[
        "reciprocity"
    ]["max_abs_diff_dB"]

    te1 = {}
    for res, name in [
        (32, "pair_P1P2_te1_hz_line_res32_rt40.json"),
        (48, "pair_P1P2_te1_hz_line_res48_rt80.json"),
        (64, "pair_P1P2_te1_hz_line_res64_rt80.json"),
    ]:
        te1[res] = json.loads((PORTS / name).read_text())["reciprocity"][
            "max_abs_diff_dB"
        ]
    p2p3_te1 = json.loads((PORTS / "pair_P2P3_te1_hz_line_res32_rt40.json").read_text())[
        "reciprocity"
    ]["max_abs_diff_dB"]

    horns_only = None
    diag_path = ROOT / "outputs" / "validation" / "diagnostics" / "pair_horns_te1hz_P1P2_res32_rt40.json"
    if diag_path.exists():
        horns_only = json.loads(diag_path.read_text())["reciprocity"]["max_abs_diff_dB"]

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    res_b = sorted(base)
    res_t = sorted(te1)
    ax.plot(
        res_b,
        [base[r] for r in res_b],
        "o-",
        color="#4C78A8",
        lw=2,
        ms=7,
        label="Baseline (Hz-line + axis flux) P1↔P2",
    )
    ax.plot(
        res_t,
        [te1[r] for r in res_t],
        "s-",
        color="#E45756",
        lw=2,
        ms=7,
        label="TE1 launch + axis flux P1↔P2",
    )
    ax.axhline(
        0.2,
        color="0.35",
        ls="--",
        lw=1,
        label="Gate: axis↔diag ≲ 0.2 dB",
    )
    # Orientation control markers at res=32
    ax.scatter(
        [32],
        [p2p3_base],
        marker="^",
        s=80,
        color="#4C78A8",
        zorder=5,
        label=f"Baseline P2↔P3 @32 ({p2p3_base:.4f} dB)",
    )
    ax.scatter(
        [32],
        [p2p3_te1],
        marker="D",
        s=70,
        color="#E45756",
        zorder=5,
        label=f"TE1 P2↔P3 @32 ({p2p3_te1:.4f} dB)",
    )
    if horns_only is not None:
        ax.scatter(
            [32],
            [horns_only],
            marker="x",
            s=90,
            color="#000000",
            lw=2,
            zorder=6,
            label=f"TE1 horns_only P1↔P2 @32 ({horns_only:.2f} dB)",
        )
    ax.set_xlabel("Resolution (pixels / a)")
    ax.set_ylabel(r"max $|P_{ij}-P_{ji}|$ (dB)")
    ax.set_title("B=0 reciprocity residual vs resolution")
    ax.set_xticks([32, 48, 64, 96])
    ax.legend(loc="upper right", framealpha=0.95)
    ax.set_ylim(bottom=-0.05)
    fig.text(
        0.5,
        -0.02,
        "P1↔P2 = axis↔diagonal (worst class). P2↔P3 = diagonal↔diagonal control. "
        "Uniform ρ, fp=8 GHz, B=0.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    path = OUT / "fig2_reciprocity_vs_resolution.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_incident():
    base_rows = json.loads((RECIP / "incident_power_vs_res.json").read_text())["rows"]
    base = {int(r["res"]): r["axis_over_diag"] for r in base_rows}

    te1 = {}
    for path in [
        PORTS / "incident_te1_hz_line_res32_48.json",
        PORTS / "incident_te1_hz_line_res64.json",
        PORTS / "smoke_incident_te1_hz_line_res32.json",
    ]:
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        for r in data["rows"]:
            te1[int(r["res"])] = r["axis_over_diag"]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    res_b = sorted(base)
    res_t = sorted(te1)
    ax.axhline(1.0, color="0.2", ls="-", lw=1.2, label="Ideal axis/diag = 1")
    ax.axhspan(0.97, 1.03, color="#54A24B", alpha=0.15, label="±3% band")
    ax.plot(
        res_b,
        [base[r] for r in res_b],
        "o-",
        color="#4C78A8",
        lw=2,
        ms=7,
        label="Baseline launcher",
    )
    ax.plot(
        res_t,
        [te1[r] for r in res_t],
        "s-",
        color="#E45756",
        lw=2,
        ms=7,
        label="TE1 launcher",
    )
    ax.set_xlabel("Resolution (pixels / a)")
    ax.set_ylabel("Incident power axis / diagonal")
    ax.set_title("Port launch-power orientation uniformity")
    ax.set_xticks(sorted(set(res_b) | set(res_t)))
    ax.legend(loc="best")
    ax.set_ylim(0.9, 1.2)
    fig.text(
        0.5,
        -0.02,
        "Axis = mean(P1,P4); diagonal = mean(P2,P3,P5,P6). "
        "TE1 removes the ~10% baseline orientation bias.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    path = OUT / "fig3_incident_power_uniformity.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


LORENTZ = ROOT / "outputs" / "validation" / "lorentz_direct"
MODEPROF = ROOT / "outputs" / "validation" / "mode_profiles"


def fig_direct_lorentz():
    lorentz_path = LORENTZ / "lorentz_direct_res32.json"
    if not lorentz_path.exists():
        return None
    data = json.loads(lorentz_path.read_text())
    port_te1 = None
    diag_path = ROOT / "outputs" / "validation" / "diagnostics" / "pair_horns_te1hz_P1P2_res32_rt40.json"
    if diag_path.exists():
        port_te1 = json.loads(diag_path.read_text())["reciprocity"]["max_abs_diff_dB"]

    pairs = data["pairs"]
    labels = [f"P{p['port_a']+1}↔P{p['port_b']+1}" for p in pairs]
    lorentz_db = [abs(p["amp_err_dB"]) for p in pairs]

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, lorentz_db, width, label="Direct Lorentz (Hz point)", color="#54A24B")
    if port_te1 is not None:
        ax.bar(
            [0],
            [port_te1],
            width,
            label="te1_hz_line + flux (P1↔P2 only)",
            color="#E45756",
        )
    ax.set_yscale("log")
    ax.set_ylabel("|amplitude error| (dB)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title("Direct reciprocity vs port-normalized metric (horns_only, res32)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.text(
        0.5,
        -0.02,
        "Direct test: localized Hz sources, complex Hz at monitor — no flux normalization.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout()
    path = OUT / "fig4_direct_lorentz_vs_port.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    _style()
    OUT.mkdir(parents=True, exist_ok=True)
    paths = [fig_faraday(), fig_reciprocity(), fig_incident()]
    p4 = fig_direct_lorentz()
    if p4 is not None:
        paths.append(p4)
    caption = OUT / "FIGURE_CAPTIONS.md"
    caption.write_text(
        """# Advisor-update figure captions

## fig1_faraday_benchmark.png
Homogeneous Faraday rotation benchmark using the PlasMEEP `GyrotropicDrudeSusceptibility`
path (fs=5 GHz, fp=2 GHz, |B|=0.05 T). Left: measured κ at B=0 / +B / −B versus independent
cold-plasma theory. Right: |κ| relative error at res=32 (~0.09%) and res=64 (~0.02%).
PASS criteria use magnitude and B-reversal; absolute Stokes/theory sign convention is opposite.

## fig2_reciprocity_vs_resolution.png
B=0 six-port reciprocity residual for the worst orientation class (P1↔P2, axis↔diagonal)
versus Meep resolution. Baseline notebook-faithful Hz-line sources peak at **1.64 dB** at
res=64; TE1 launch reduces that to **1.03 dB** but remains above the 0.2 dB gate.
Diagonal↔diagonal control pair P2↔P3 stays ~0 at res=32 for both formulations.

## fig3_incident_power_uniformity.png
Ratio of launched incident power for axis-aligned ports to diagonal (60°) ports.
Baseline retains a ~10–14% orientation bias at all resolutions. TE1 cos-profile launch
brings the ratio within a few percent of unity (inside the shaded ±3% band).
""",
        encoding="utf-8",
    )
    for p in paths:
        print("Wrote", p)
    print("Wrote", caption)


if __name__ == "__main__":
    main()
