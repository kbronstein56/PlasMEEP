#!/usr/bin/env python3
"""Cold-plasma tensor from the Lorentz electron equation, e^{-iωt}, B along +z.

Ordinary frequencies are in one common unit. ε_xy = +i η for f_c > 0.
"""
from __future__ import annotations

import numpy as np

E_CHARGE = 1.602176634e-19
E_MASS = 9.1093837015e-31
C_LIGHT = 299792458.0


def cyclotron_hz(b_tesla: float) -> float:
    """f_c = e B / (2π m), signed with B_z."""
    return E_CHARGE * b_tesla / (2.0 * np.pi * E_MASS)


def ordinary_from_si(f_hz: float, fp_hz: float, gamma_hz: float, b_tesla: float, a_m: float = 0.020):
    scale = a_m / C_LIGHT
    return f_hz * scale, fp_hz * scale, gamma_hz * scale, cyclotron_hz(b_tesla) * scale


def tensor_ordinary(f: float, fp: float, gamma: float, fc: float):
    """Return ε_xx, ε_xy, ε_yx, ε_yy, ε_zz and the scalar η."""
    u = f + 1j * gamma
    denom = u * u - fc * fc
    e_perp = 1.0 - (fp * fp) * u / (f * denom)
    eta = (fp * fp) * fc / (f * denom)
    e_zz = 1.0 - (fp * fp) / (f * u)
    return e_perp, 1j * eta, -1j * eta, e_perp, e_zz, eta


def rho_xy(e_xx: complex, e_xy: complex, e_yx: complex, e_yy: complex):
    det = e_xx * e_yy - e_xy * e_yx
    return (
        e_yy / det,
        -e_xy / det,
        -e_yx / det,
        e_xx / det,
        det,
    )


def matrix_of(exx, exy, eyx, eyy) -> np.ndarray:
    return np.array([[exx, exy], [eyx, eyy]], dtype=np.complex128)


def voigt_eps_eff(e_perp: complex, eta: complex) -> complex:
    """n² for k ⟂ B in the Hz polarization."""
    return (e_perp * e_perp - eta * eta) / e_perp


def circular_split(e_perp: complex, eta: complex):
    """Eigenvalues of the xy tensor. (1, i) has ε_⊥ - η."""
    return e_perp - eta, e_perp + eta


def faraday_kappa(f: float, e_plus: complex, e_minus: complex) -> complex:
    """κ = (k_+ - k_-)/2 for propagation along +z. k = 2π f √ε, Re(k) >= 0."""
    def k_of(eps: complex) -> complex:
        k = 2.0 * np.pi * f * np.sqrt(eps)
        if np.real(k) < 0.0:
            k = -k
        return complex(k)

    return 0.5 * (k_of(e_plus) - k_of(e_minus))


def dissipation_matrix(eps: np.ndarray) -> np.ndarray:
    """(ε - ε†) / (2i). Positive semidefinite for a passive e^{-iωt} medium."""
    return (eps - eps.conj().T) / (2j)
