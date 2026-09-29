#!/usr/bin/env python3
"""Compare complex fields on shared lines. Powers already disagree, so this localizes where."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parents[3] / "outputs" / "validation" / "fem_meep_validation" / "horn_localization"
meep = np.load(OUT / "meep_horns_only_fields.npz")
fem = np.load(OUT / "fem_horns_only_fields_FEM-M.npz")


def align(ref, other):
    """Remove one global complex scale so shapes can be compared."""
    den = np.vdot(other, other)
    if abs(den) < 1e-30:
        return other
    scale = np.vdot(other, ref) / den
    return other * scale


def stats(name, span_half):
    fHz = fem[f"{name}_Hz"]
    mHz = meep[f"{name}_Hz"]
    s_f = fem[f"{name}_s"]
    s_m = np.linspace(s_f[0], s_f[-1], len(mHz))
    # Meep volume is the same physical segment; resample Meep onto FEM abscissa.
    m_on_f = np.interp(s_f, s_m, np.real(mHz)) + 1j * np.interp(s_f, s_m, np.imag(mHz))
    fitted = align(m_on_f, fHz)
    num = np.linalg.norm(m_on_f - fitted)
    den = np.linalg.norm(m_on_f) + 1e-30
    phase = np.angle(m_on_f * np.conj(fitted))
    # weight phase by amplitude
    w = np.abs(m_on_f)
    w = w / (w.sum() + 1e-30)
    # Sx for port-0 lines (outward = +x): 0.5 Re(Ey conj(Hz))
    def sdot(Hz, Ex, Ey):
        return 0.5 * np.real(Ey * np.conj(Hz))

    mEx = np.interp(s_f, s_m, np.real(meep[f"{name}_Ex"])) + 1j * np.interp(s_f, s_m, np.imag(meep[f"{name}_Ex"]))
    mEy = np.interp(s_f, s_m, np.real(meep[f"{name}_Ey"])) + 1j * np.interp(s_f, s_m, np.imag(meep[f"{name}_Ey"]))
    return {
        "name": name,
        "rel_L2_after_complex_scale": float(num / den),
        "abs_corr": float(np.corrcoef(np.abs(m_on_f), np.abs(fHz))[0, 1]),
        "phase_std_rad": float(np.sqrt(np.sum(w * (phase - np.sum(w * phase)) ** 2))),
        "meep_Hz": m_on_f,
        "fem_Hz_aligned": fitted,
        "s": s_f,
        "meep_S": sdot(m_on_f, mEx, mEy),
        "fem_S": sdot(fitted, fem[f"{name}_Ex"] * (fitted[np.argmax(np.abs(fitted))] / (fHz[np.argmax(np.abs(fHz))] + 1e-30)), fem[f"{name}_Ey"]),
    }


# Recompute S with the same complex scale applied to all components.
rows = []
fig, axes = plt.subplots(4, 2, figsize=(9.5, 11), sharex=False)
for i, name in enumerate(("source", "throat", "monitor", "center")):
    fHz = fem[f"{name}_Hz"]
    fEx = fem[f"{name}_Ex"]
    fEy = fem[f"{name}_Ey"]
    mHz = meep[f"{name}_Hz"]
    mEx = meep[f"{name}_Ex"]
    mEy = meep[f"{name}_Ey"]
    s = fem[f"{name}_s"]
    s_m = np.linspace(s[0], s[-1], len(mHz))

    def resamp(z):
        return np.interp(s, s_m, np.real(z)) + 1j * np.interp(s, s_m, np.imag(z))

    mHz, mEx, mEy = resamp(mHz), resamp(mEx), resamp(mEy)
    den = np.vdot(fHz, fHz)
    scale = np.vdot(fHz, mHz) / den if abs(den) > 1e-30 else 1
    fHz, fEx, fEy = fHz * scale, fEx * scale, fEy * scale
    rel = float(np.linalg.norm(mHz - fHz) / (np.linalg.norm(mHz) + 1e-30))
    corr = float(np.corrcoef(np.abs(mHz), np.abs(fHz))[0, 1]) if np.std(np.abs(mHz)) > 0 else float("nan")
    phase = np.angle(mHz * np.conj(fHz))
    w = np.abs(mHz)
    w = w / (w.sum() + 1e-30)
    mean_ph = float(np.sum(w * phase))
    phase_std = float(np.sqrt(np.sum(w * (phase - mean_ph) ** 2)))
    meep_S = 0.5 * np.real(mEy * np.conj(mHz))
    fem_S = 0.5 * np.real(fEy * np.conj(fHz))
    rows.append(
        {
            "line": name,
            "rel_L2_Hz": rel,
            "abs_Hz_corr": corr,
            "phase_std_rad": phase_std,
            "phase_std_deg": float(np.degrees(phase_std)),
        }
    )
    ax_a, ax_p = axes[i]
    ax_a.plot(s, np.abs(mHz), label="Meep")
    ax_a.plot(s, np.abs(fHz), label="FEM-M aligned")
    ax_a.set_ylabel("|Hz|")
    ax_a.set_title(f"{name}  rel L2={rel:.3f}  corr={corr:.4f}")
    ax_p.plot(s, np.degrees(np.unwrap(np.angle(mHz))), label="Meep")
    ax_p.plot(s, np.degrees(np.unwrap(np.angle(fHz))), label="FEM-M")
    ax_p.set_ylabel("phase deg")
    if i == 0:
        ax_a.legend(fontsize=8)
        ax_p.legend(fontsize=8)
    ax_p.set_xlabel("position along line (a)")

fig.tight_layout()
fig.savefig(OUT / "field_lines_horns_only.png", dpi=140)
(OUT / "field_line_compare.json").write_text(json.dumps(rows, indent=2) + "\n")
print(json.dumps(rows, indent=2))
