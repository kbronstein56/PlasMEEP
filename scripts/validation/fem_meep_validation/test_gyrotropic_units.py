#!/usr/bin/env python3
"""
Gyrotropic unit tests: B=0 reduction, B→−B off-diagonal reversal, Faraday |κ| check.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
]

import sixport_common as sc  # noqa: E402
from faraday_benchmark import (  # noqa: E402
    circular_eigenpermittivities,
    gyrotropic_drude_eps_eta,
    propagation_constants,
    signed_faraday_kappa,
)
from fem_validated_solver import eps_tensor_at_bias, rho_from_eps  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_meep_validation" / "phase9"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    checks = []

    # 1) B=0 → η=0 → isotropic
    exx, exy, eyx, eyy = eps_tensor_at_bias(0.0)
    eps_p, eta = gyrotropic_drude_eps_eta(sc.fs_a, sc.fp_a, sc.gamma_a, 0.0)
    ok1 = abs(eta) < 1e-14 and abs(exy) < 1e-14 and abs(eyx) < 1e-14 and abs(exx - eyy) < 1e-14
    rxx, rxy, ryx, ryy = rho_from_eps(exx, exy, eyx, eyy)
    ok1b = abs(rxy) < 1e-12 and abs(ryx) < 1e-12 and abs(rxx - 1 / exx) < 1e-10
    checks.append({"name": "B0_isotropic", "pass": bool(ok1 and ok1b), "eps_perp": [exx.real, exx.imag], "eta": [eta.real, eta.imag]})

    # 2) B → −B reverses off-diagonal
    bias = 0.05
    e1 = eps_tensor_at_bias(bias)
    e2 = eps_tensor_at_bias(-bias)
    ok2 = (
        abs(e1[0] - e2[0]) < 1e-12
        and abs(e1[3] - e2[3]) < 1e-12
        and abs(e1[1] + e2[1]) < 1e-12
        and abs(e1[2] + e2[2]) < 1e-12
    )
    r1 = rho_from_eps(*e1)
    r2 = rho_from_eps(*e2)
    ok2b = abs(r1[1] + r2[1]) < 1e-10 and abs(r1[2] + r2[2]) < 1e-10
    checks.append(
        {
            "name": "B_reversal_offdiag",
            "pass": bool(ok2 and ok2b),
            "eps_xy_plus": [e1[1].real, e1[1].imag],
            "eps_xy_minus": [e2[1].real, e2[1].imag],
            "rho_xy_plus": [r1[1].real, r1[1].imag],
            "rho_xy_minus": [r2[1].real, r2[1].imag],
        }
    )

    # 3) Faraday |κ| vs analytic at bias=0.05 (homogeneous theory, not FEM solve)
    eps_perp, eta = gyrotropic_drude_eps_eta(sc.fs_a, sc.fp_a, sc.gamma_a, bias)
    ep, em = circular_eigenpermittivities(eps_perp, eta)
    kp, km = propagation_constants(sc.fs_a, ep, em)
    kappa = signed_faraday_kappa(kp, km)
    # Reverse B
    eps_perp2, eta2 = gyrotropic_drude_eps_eta(sc.fs_a, sc.fp_a, sc.gamma_a, -bias)
    ep2, em2 = circular_eigenpermittivities(eps_perp2, eta2)
    kp2, km2 = propagation_constants(sc.fs_a, ep2, em2)
    kappa2 = signed_faraday_kappa(kp2, km2)
    ok3 = abs(abs(kappa) - abs(kappa2)) < 1e-10 and abs(kappa + kappa2) < 1e-8
    # Load trusted Faraday res64 if present
    faraday_path = ROOT / "outputs" / "validation" / "faraday_res64.json"
    faraday_cmp = None
    if faraday_path.is_file():
        fj = json.loads(faraday_path.read_text())
        # find theory kappa magnitude if present
        faraday_cmp = {
            "file": str(faraday_path),
            "keys": list(fj.keys())[:20],
        }
        for key in ("kappa_theory", "kappa_meep", "theory_kappa", "kappa_signed_theory", "results"):
            if key in fj:
                faraday_cmp[key] = fj[key]
    checks.append(
        {
            "name": "faraday_kappa_B_reversal",
            "pass": bool(ok3),
            "kappa_plusB": [kappa.real, kappa.imag],
            "kappa_minusB": [kappa2.real, kappa2.imag],
            "abs_kappa": abs(kappa),
            "faraday_reference": faraday_cmp,
        }
    )

    # 4) Tensor form matches documented Meep convention ε=[[ε⊥,-iη],[iη,ε⊥]]
    ok4 = abs(e1[1] - (-1j * eta)) < 1e-12 and abs(e1[2] - (1j * eta)) < 1e-12
    checks.append({"name": "meep_tensor_signs", "pass": bool(ok4), "eta": [eta.real, eta.imag]})

    summary = {
        "all_pass": all(c["pass"] for c in checks),
        "checks": checks,
        "fs_a": sc.fs_a,
        "fp_a": sc.fp_a,
        "gamma_a": sc.gamma_a,
        "bias_test": bias,
    }
    path = OUT / "gyrotropic_unit_tests.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if summary["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
