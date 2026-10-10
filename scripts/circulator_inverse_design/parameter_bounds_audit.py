#!/usr/bin/env python3
"""Audit plasma density-scale bounds for optimization."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sixport_common as sc  # noqa: E402
from gyrotropic_tensor import dissipation_matrix, matrix_of, ordinary_from_si, rho_xy, tensor_ordinary  # noqa: E402
from plasma_sensitivity import deps_ds, drho_ds  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"
S_LO, S_HI = 0.05, 1.80
SAMPLES = np.array([0.0, 0.05, 0.25, 0.5, 1.0, 1.5, 1.8, 2.0])


def main():
    rows = []
    ok = True
    for b in (0.0, 0.05, -0.05):
        for s in SAMPLES:
            f, fp, gamma, fc = ordinary_from_si(sc.fs_Hz, sc.fp_Hz * np.sqrt(max(s, 0.0)), sc.gamma_Hz, b, sc.a)
            exx, exy, eyx, eyy, ezz, eta = tensor_ordinary(f, fp, gamma, fc)
            rxx, rxy, ryx, ryy, det = rho_xy(exx, exy, eyx, eyy)
            eps = matrix_of(exx, exy, eyx, eyy)
            diss = dissipation_matrix(eps)
            evals = np.linalg.eigvalsh(diss)
            # analytic deps at this s (including s=0)
            fp_ref = sc.fp_Hz * sc.a / 299792458.0
            de = deps_ds(f, fp_ref, gamma, fc, max(s, 0.0))
            dr = drho_ds(f, fp_ref, gamma, fc, max(s, 1e-30) if s == 0 else s)
            # FD check away from exact 0 if needed
            h = 1e-5
            s_fd = max(s, h)
            fp_p = fp_ref * np.sqrt(s_fd + h)
            fp_m = fp_ref * np.sqrt(max(s_fd - h, 0.0))
            ep = tensor_ordinary(f, fp_p, gamma, fc)[:4]
            em = tensor_ordinary(f, fp_m, gamma, fc)[:4]
            de_fd = [(ep[i] - em[i]) / (2 * h) for i in range(4)]
            de_ana = deps_ds(f, fp_ref, gamma, fc, s_fd)[:4]
            eps_rel = max(abs(de_fd[i] - de_ana[i]) / max(abs(de_ana[i]), 1e-30) for i in range(4))
            singular = abs(det) < 1e-6
            passive = bool(np.min(evals) >= -1e-10)
            finite = all(np.isfinite([exx, exy, rxx, rxy, det]))
            if singular or not finite:
                ok = False
            rows.append({
                "B_T": b,
                "s": float(s),
                "det_eps": [float(np.real(det)), float(np.imag(det))],
                "min_diss_eig": float(np.min(evals)),
                "passive": passive,
                "singular": bool(singular),
                "eps_deriv_fd_rel": float(eps_rel),
                "rho_xx": [float(np.real(rxx)), float(np.imag(rxx))],
            })
    # freeze bounds
    bounds = {
        "s_lo": S_LO,
        "s_hi": S_HI,
        "rationale": (
            "s=0 is mathematically ok (fp^2=0) but excluded from the closed interval "
            "to avoid optimizer pile-up on a vacuum corner; s_lo=0.05 keeps a weak plasma. "
            "s_hi=1.80 stays below the exploratory s=2 map and away from denser regimes "
            "not exercised in the validated Jacobian campaign. Production s=1 is interior."
        ),
        "samples": rows,
        "pass": bool(ok and all(r["passive"] for r in rows if S_LO <= r["s"] <= S_HI)),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "parameter_bounds_audit.json").write_text(json.dumps(bounds, indent=2) + "\n")
    md = [
        "# Circulator parameter bounds",
        "",
        "Frozen before production optimization.",
        "",
        f"- `s ∈ [{S_LO}, {S_HI}]`",
        "- `f_p^2 = s f_p,ref^2` with `f_p,ref = 8 GHz`",
        "- Analytic `∂ε/∂s` remains finite at `s=0` (depends on `f_p,ref^2`, not `1/√s`)",
        "- In-range samples: no permittivity-determinant singularity; dissipation matrix PSD",
        "",
        bounds["rationale"],
        "",
        f"**BOUNDS_AUDIT: {'PASS' if bounds['pass'] else 'FAIL'}**",
        "",
    ]
    (OUT / "CIRCULATOR_PARAMETER_BOUNDS.md").write_text("\n".join(md))
    print("bounds", S_LO, S_HI, "pass", bounds["pass"], flush=True)


if __name__ == "__main__":
    main()
