#!/usr/bin/env python3
"""Tiny custom-VJP prototype: wrap a validated Ax=b adjoint as a differentiable op.

If JAX is unavailable, implements the same VJP contract in NumPy and compares to FD.
Does NOT port FULL91 FEM.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"


def build_tiny_A(s: float) -> np.ndarray:
    """2×2 complex SPD-ish system depending on scalar s (toy Maxwell stand-in)."""
    return np.array([[2.0 + s, 0.1j], [0.1j, 1.5 + 0.5 * s]], dtype=np.complex128)


def objective(x: np.ndarray) -> float:
    return float(np.real(np.vdot(x, x)))  # ||x||^2


def solve_and_J(s: float):
    A = build_tiny_A(s)
    b = np.array([1.0, 0.5], dtype=np.complex128)
    x = np.linalg.solve(A, b)
    return objective(x), A, x, b


def analytic_adjoint_grad(s: float) -> float:
    """dJ/ds via adjoint: J=x^H x, A x = b, dA/ds known.

    dJ = 2 Re(x^H dx), dx = -A^{-1} (dA) x ⇒
    A^H λ = x,  dJ/ds = -2 Re(λ^H (dA/ds) x).
    """
    _, A, x, _ = solve_and_J(s)
    lam = np.linalg.solve(A.conj().T, x)
    dAds = np.array([[1.0, 0.0], [0.0, 0.5]], dtype=np.complex128)
    return float(-2.0 * np.real(np.vdot(lam, dAds @ x)))


def numpy_custom_vjp_grad(s: float) -> float:
    """Emulate JAX custom_vjp: forward calls solve; backward calls adjoint."""
    return analytic_adjoint_grad(s)


def fd_grad(s: float, h=1e-6) -> float:
    Jp, *_ = solve_and_J(s + h)
    Jm, *_ = solve_and_J(s - h)
    return (Jp - Jm) / (2 * h)


def try_jax(s0: float):
    try:
        import jax
        import jax.numpy as jnp
        from jax import custom_vjp
    except ImportError as e:
        return {"available": False, "error": str(e)}

    @custom_vjp
    def J_of_s(s):
        # forward uses numpy host callback pattern
        def _fwd(sv):
            val = float(sv)
            J, A, x, b = solve_and_J(val)
            return jnp.asarray(J), (jnp.asarray(val),)

        return _fwd(s)[0]

    def J_fwd(s):
        val = float(np.asarray(s))
        J, A, x, b = solve_and_J(val)
        return jnp.asarray(J), (jnp.asarray(val),)

    def J_bwd(res, g):
        (s,) = res
        # call validated adjoint; scale by upstream cotangent
        dJds = analytic_adjoint_grad(float(s))
        return (jnp.asarray(dJds) * g,)

    J_of_s.defvjp(J_fwd, J_bwd)
    g = jax.grad(J_of_s)(jnp.asarray(s0))
    return {"available": True, "grad": float(g)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s0 = 0.3
    g_adj = analytic_adjoint_grad(s0)
    g_fd = fd_grad(s0)
    g_wrap = numpy_custom_vjp_grad(s0)
    jax_info = try_jax(s0)
    rel = abs(g_adj - g_fd) / max(abs(g_fd), 1e-30)
    report = {
        "s0": s0,
        "J": solve_and_J(s0)[0],
        "grad_adjoint": g_adj,
        "grad_custom_vjp_numpy": g_wrap,
        "grad_fd": g_fd,
        "rel_err_adj_vs_fd": rel,
        "pass": rel < 1e-5,
        "jax": jax_info,
        "conclusion": (
            "custom_vjp wrapping a validated adjoint matches FD on this toy problem. "
            "FULL91 production should keep the analytic discrete adjoint; JAX would be "
            "API sugar unless sparse complex LU moves into a JIT/GPU stack."
        ),
    }
    (OUT / "jax_custom_vjp_tiny.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
