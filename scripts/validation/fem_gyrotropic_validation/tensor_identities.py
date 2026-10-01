#!/usr/bin/env python3
"""Tensor identities for the independently derived gyrotropic permittivity."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from gyrotropic_tensor import (  # noqa: E402
    circular_split,
    dissipation_matrix,
    faraday_kappa,
    matrix_of,
    ordinary_from_si,
    rho_xy,
    tensor_ordinary,
    voigt_eps_eff,
)
from faraday_benchmark import gyrotropic_drude_eps_eta  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_gyrotropic_validation"


def scalar_drude(f, fp, gamma):
    return 1.0 - fp * fp / (f * f + 1j * f * gamma)


def one_case(f, fp, gamma, fc) -> dict:
    exx, exy, eyx, eyy, ezz, eta = tensor_ordinary(f, fp, gamma, fc)
    exx0, exy0, eyx0, eyy0, ezz0, _ = tensor_ordinary(f, fp, gamma, 0.0)
    exxm, exym, eyxm, eyym, _, etam = tensor_ordinary(f, fp, gamma, -fc)
    rxx, rxy, ryx, ryy, det = rho_xy(exx, exy, eyx, eyy)
    eps = matrix_of(exx, exy, eyx, eyy)
    rho = matrix_of(rxx, rxy, ryx, ryy)
    ident = rho @ eps - np.eye(2)
    inv = np.linalg.inv(eps)
    ons = eps - matrix_of(exxm, eyxm, exym, eyym)
    # epsilon(B) should equal epsilon(-B).T, which swaps xy and yx.
    ons_t = eps - matrix_of(exxm, eyxm, exym, eyym).T
    # Build epsilon(-B).T explicitly.
    eps_mb_t = matrix_of(exxm, exym, eyxm, eyym).T
    diss = dissipation_matrix(eps)
    eig = np.linalg.eigvalsh(0.5 * (diss + diss.conj().T))
    # Compare diagonal and |η| with the existing helper. Off-diagonal sign is reported, not fitted.
    helper_e, helper_eta = gyrotropic_drude_eps_eta(f, fp, gamma, fc)
    stored = -3.317759870618328 + 0.001121496070290476j
    return {
        "fc": fc,
        "exx_eq_eyy": abs(exx - eyy),
        "exy_eq_minus_eyx": abs(exy + eyx),
        "diag_even_in_B": max(abs(exx - exxm), abs(eyy - eyym)),
        "off_odd_in_B": abs(exy + exym),
        "b0_xx": abs(exx0 - scalar_drude(f, fp, gamma)),
        "b0_off": abs(exy0) + abs(eyx0),
        "rho_eps_resid": float(np.linalg.norm(ident)),
        "inv_resid": float(np.linalg.norm(rho - inv)),
        "onsager_resid": float(np.linalg.norm(eps - eps_mb_t)),
        "det_abs": abs(det),
        "diss_min_eig": float(np.min(np.real(eig))),
        "helper_diag_diff": abs(exx - helper_e),
        "helper_eta_magnitude_diff": abs(abs(eta) - abs(helper_eta)),
        "helper_eta_sum": abs(eta + helper_eta),
        "voigt_n2": [voigt_eps_eff(exx, eta).real, voigt_eps_eff(exx, eta).imag],
        "ezz_b0_like": abs(ezz0 - scalar_drude(f, fp, gamma)),
        "stored_production_unused": abs(stored),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f_si, fp_si, g_si, fc_si = ordinary_from_si(3.85e9, 8.00e9, 1.00e6, 0.05)
    rows = []
    grid = []
    for f in (0.20, 0.2568443533025771, 0.45):
        for fp in (0.15, 0.5337025523170433):
            for gamma in (0.0, 6.671281903963041e-05, 0.01):
                for fc in (0.0, 0.02, 0.08, -0.08):
                    grid.append((f, fp, gamma, fc))
    # Production SI point at 0.05 T, plus the validated B=0 numbers.
    grid.append((f_si, fp_si, g_si, 0.0))
    grid.append((f_si, fp_si, g_si, fc_si))
    grid.append((f_si, fp_si, g_si, -fc_si))
    worst = {}
    for f, fp, gamma, fc in grid:
        rec = one_case(f, fp, gamma, fc)
        rec.update({"f": f, "fp": fp, "gamma": gamma})
        rows.append(rec)
    keys = [
        "exx_eq_eyy",
        "exy_eq_minus_eyx",
        "diag_even_in_B",
        "off_odd_in_B",
        "b0_xx",
        "b0_off",
        "rho_eps_resid",
        "inv_resid",
        "onsager_resid",
        "helper_diag_diff",
        "helper_eta_magnitude_diff",
    ]
    for key in keys:
        worst[key] = max(abs(r[key]) for r in rows)
    # Dissipation must not go negative beyond roundoff. Lossless rows sit at 0.
    diss_min = min(r["diss_min_eig"] for r in rows)
    # Faraday sign: κ(+B) = -κ(-B), and κ(0) = 0, on a propagating dielectric-like point.
    f, fp, gamma, fc = 0.45, 0.15, 1e-4, 0.05
    ep, eta = tensor_ordinary(f, fp, gamma, fc)[0], tensor_ordinary(f, fp, gamma, fc)[5]
    em, etam = tensor_ordinary(f, fp, gamma, -fc)[0], tensor_ordinary(f, fp, gamma, -fc)[5]
    e0, eta0 = tensor_ordinary(f, fp, gamma, 0.0)[0], tensor_ordinary(f, fp, gamma, 0.0)[5]
    kp = circular_split(ep, eta)
    km = circular_split(em, etam)
    k0 = circular_split(e0, eta0)
    kap = faraday_kappa(f, kp[0], kp[1])
    kam = faraday_kappa(f, km[0], km[1])
    ka0 = faraday_kappa(f, k0[0], k0[1])
    # Production B=0 epsilon versus the frozen value.
    fs, fps, gs, _ = ordinary_from_si(3.85e9, 8.0e9, 1.0e6, 0.0)
    e_prod = tensor_ordinary(fs, fps, gs, 0.0)[0]
    stored = -3.317759870618328 + 0.001121496070290476j
    summary = {
        "n_cases": len(rows),
        "worst": worst,
        "dissipation_min_eigenvalue": diss_min,
        "faraday_kappa_plus": [kap.real, kap.imag],
        "faraday_kappa_minus": [kam.real, kam.imag],
        "faraday_kappa_zero": [ka0.real, ka0.imag],
        "faraday_reversal_resid": abs(kap + kam),
        "production_b0_eps_diff": abs(e_prod - stored),
        "production_fc_0.05T": fc_si,
        "note": "helper_eta_sum near 2|η| means the existing helper stores the opposite off-diagonal sign.",
    }
    # helper eta sum on one magnetized row
    sample = one_case(f, fp, gamma, fc)
    summary["sample_helper_eta_sum"] = sample["helper_eta_sum"]
    summary["sample_helper_diag_diff"] = sample["helper_diag_diff"]
    (OUT / "tensor_identities.json").write_text(json.dumps({"summary": summary, "rows": rows}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
