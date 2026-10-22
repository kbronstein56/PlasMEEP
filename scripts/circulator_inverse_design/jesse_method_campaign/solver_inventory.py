#!/usr/bin/env python3
"""Inventory sparse direct solvers available without installing new packages."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"


def try_import(name):
    try:
        __import__(name)
        return True, None
    except Exception as e:
        return False, str(e)


def tiny_factor_bench(n=800, density=0.01):
    rng = np.random.default_rng(0)
    A = sparse.random(n, n, density=density, format="csc", dtype=np.complex128, random_state=rng)
    A = A + sparse.eye(n) * (2.0 + 1j)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    results = {}
    from scipy.sparse.linalg import splu, spsolve

    t0 = time.perf_counter()
    lu = splu(A)
    x = lu.solve(b)
    results["scipy_superlu"] = {
        "factor_solve_s": time.perf_counter() - t0,
        "residual": float(np.linalg.norm(A @ x - b)),
        "available": True,
    }
    # spsolve may use UMFPACK if scikit-umfpack present
    t0 = time.perf_counter()
    x2 = spsolve(A, b)
    results["scipy_spsolve"] = {
        "solve_s": time.perf_counter() - t0,
        "residual": float(np.linalg.norm(A @ x2 - b)),
        "available": True,
    }
    return results


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inv = {
        "scipy_sparse_linalg_splu": {"available": True, "backend": "SuperLU (default)"},
        "packages": {},
    }
    for pkg in [
        "pypardiso",
        "pardiso",
        "sklearn.utils.fixes",
        "sksparse",
        "sksparse.cholmod",
        "petsc4py",
        "mumps",
        "pymumps",
        "jax",
        "autograd",
    ]:
        ok, err = try_import(pkg.split(".")[0] if pkg != "sksparse.cholmod" else "sksparse")
        if pkg == "sksparse.cholmod":
            try:
                from sksparse.cholmod import cholesky  # noqa: F401

                ok, err = True, None
            except Exception as e:
                ok, err = False, str(e)
        inv["packages"][pkg] = {"available": ok, "error": err}
    inv["tiny_complex_bench"] = tiny_factor_bench()
    # symbolic reuse note
    inv["symbolic_factor_reuse"] = {
        "scipy_superlu": (
            "SciPy SuperLU factorizes numerically each call; no public API to "
            "reuse symbolic analysis across matrices with identical sparsity "
            "but different values (nearby frequencies). Would need a solver "
            "with separate symbolic/numeric phases (PARDISO, MUMPS, PETSc)."
        ),
        "recommendation": (
            "Do not change production SuperLU until an installed alternative "
            "beats FULL91-C LU (~1008 s) with validated residual/gradient."
        ),
    }
    (OUT / "solver_inventory.json").write_text(json.dumps(inv, indent=2) + "\n")
    print(json.dumps(inv, indent=2))


if __name__ == "__main__":
    main()
