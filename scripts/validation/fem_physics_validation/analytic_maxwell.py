#!/usr/bin/env python3
"""Independent analytic Maxwell references for the Hz polarization.

Time convention e^{-iωt}. Outgoing waves use H^{(1)}.
The 2D unknown is Hz. Isotropic media satisfy

    div((1/ε) grad Hz) + k0^2 Hz = 0,

with Hz and (1/ε) ∂n Hz continuous. Tangential E is

    (Ex, Ey) = (i/ω) (1/ε) (∂y Hz, -∂x Hz).

Nothing in this module imports the FEM assembler, a mesh, or Meep.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import hankel1, h1vp, jv, jvp

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def kz_of(k0: float, eps: complex, ky: float = 0.0) -> complex:
    """Normal wave number with Im(kx) >= 0, so exp(i kx x) is outgoing or decaying as x → +∞."""
    val = complex(k0**2 * eps - ky**2)
    kx = np.sqrt(val)
    if np.real(kx) < 0 and abs(np.imag(kx)) < 1e-14:
        kx = -kx
    if np.imag(kx) < 0:
        kx = -kx
    return complex(kx)


def fresnel_ht(k0: float, eps1: complex, eps2: complex, ky: float = 0.0) -> dict:
    """Reflection and transmission of Hz at a planar interface x=0.

    Incident wave exp(i kx1 x + i ky y) in medium 1.
    r = (kx1/ε1 - kx2/ε2) / (kx1/ε1 + kx2/ε2), t = 1 + r.
    Power uses Re(kx/ε) |Hz|^2, the normal Poynting factor for this polarization.
    """
    kx1 = kz_of(k0, eps1, ky)
    kx2 = kz_of(k0, eps2, ky)
    a1 = kx1 / eps1
    a2 = kx2 / eps2
    r = (a1 - a2) / (a1 + a2)
    t = 1.0 + r
    p1 = np.real(a1)
    p2 = np.real(a2)
    R = abs(r) ** 2
    T = (p2 / p1) * abs(t) ** 2 if abs(p1) > 1e-14 else float("nan")
    return {
        "r": complex(r),
        "t": complex(t),
        "R_power": float(R),
        "T_power": float(T),
        "power_sum": float(R + T) if np.isfinite(T) else float("nan"),
        "kx1": kx1,
        "kx2": kx2,
    }


def stack_response(k0: float, stack: list[tuple[float, complex]]) -> dict:
    """stack is a list of (thickness, epsilon), left to right. Both exteriors are air.

    Phase reference: incident wave e^{i k0 x} with unit amplitude at x=0, and the
    front face is at x=0. Downstream field for x greater than the stack thickness
    is t exp(i k0 x). Reflected field for x<0 is r exp(-i k0 x).
    """
    # State just to the right of the back face: outgoing, Hz=t_raw at the back face,
    # then scale so the incident amplitude is 1.
    x = 0.0
    faces = []
    for thickness, eps in stack:
        faces.append((x, x + thickness, complex(eps)))
        x += thickness
    L = x
    # Integrate the state from the back face toward the front.
    # At x=L+, Hz = 1, (1/ε) dHz/dx = i k0 / 1 * 1, for a unit outgoing wave.
    # Propagate this state backward through the layers. Then match to A, B on the left.
    hz = 1.0 + 0j
    slope = 1j * k0  # (1/ε) dHz/dx in air, ε=1, for exp(i k0 (x-L)) renormalized at x=L
    # Step backward through reversed layers. Within a layer of constant k,
    # Hz(x) = P cos(k (x-x_right)) + Q sin(k (x-x_right)) is inconvenient going backward.
    # Use the transfer over a distance -d:
    # [Hz, s] at left = M [Hz, s] at right, s=(1/ε) Hz_x.
    for x_left, x_right, eps in reversed(faces):
        d = x_right - x_left
        k = kz_of(k0, eps, 0.0)
        # Forward transfer by +d: 
        # Hz(x+d) = Hz cos(kd) + s (ε/k) sin(kd)
        # s(x+d) = s cos(kd) - Hz (k/ε) sin(kd)
        # so backward, d → -d, sin flips sign.
        kd = k * d
        c = np.cos(kd)
        s = np.sin(kd)
        # Going from right to left by distance d:
        # Hz_left = c * Hz_right - (eps/k) s * slope_right
        # slope_left = (k/eps) s * Hz_right + c * slope_right
        hz_left = c * hz - (eps / k) * s * slope
        slope_left = (k / eps) * s * hz + c * slope
        hz, slope = hz_left, slope_left
    # At x=0-, air: Hz = A + B, s = i k0 (A - B)
    # A + B = hz, i k0 (A - B) = slope
    A = 0.5 * (hz + slope / (1j * k0))
    B = 0.5 * (hz - slope / (1j * k0))
    # The outgoing wave was scaled to Hz(L)=1, so t_at_back = 1/A if incident A is normalized to 1,
    # and the field for x>L is exp(i k0 (x-L)) / A * (incident becomes 1).
    # Downstream ratio to a vacuum wave exp(i k0 x) is [exp(i k0 (x-L))/A] / exp(i k0 x) = exp(-i k0 L)/A.
    r = B / A
    t_downstream = np.exp(-1j * k0 * L) / A
    return {
        "r": complex(r),
        "t_over_vacuum": complex(t_downstream),
        "thickness": float(L),
        "A_left": complex(A),
        "B_left": complex(B),
    }


def layer_field(x: np.ndarray, k0: float, stack: list[tuple[float, complex]], xs: float) -> np.ndarray:
    """Field of a unit-incident stack, including layer interiors, incident phase exp(-i k0 xs) at x=0."""
    # Reconstruct coefficients by a forward march from the known left state.
    resp = stack_response(k0, stack)
    # Incident amplitude at the global exp(i k0 x) basis is exp(-i k0 xs), so A_left_phys = exp(-i k0 xs).
    scale = np.exp(-1j * k0 * xs) / resp["A_left"]
    # State at x=0+ equals state at x=0- : Hz = A+B, s = ik (A-B), then times scale*A_left
    A = resp["A_left"] * scale
    B = resp["B_left"] * scale
    hz0 = A + B
    s0 = 1j * k0 * (A - B)
    # Walk layers, storing (x_left, x_right, eps, hz_left, s_left)
    pieces = []
    x_cursor = 0.0
    hz, slope = hz0, s0
    for thickness, eps in stack:
        pieces.append((x_cursor, x_cursor + thickness, complex(eps), hz, slope, kz_of(k0, complex(eps))))
        d = thickness
        k = kz_of(k0, complex(eps))
        kd = k * d
        c = np.cos(kd)
        sn = np.sin(kd)
        hz_right = c * hz + (eps / k) * sn * slope
        slope_right = c * slope - (k / eps) * sn * hz
        hz, slope = hz_right, slope_right
        x_cursor += thickness
    L = x_cursor
    out = np.empty(len(x), dtype=np.complex128)
    for i, xv in enumerate(np.asarray(x, float)):
        if xv <= 0.0:
            out[i] = A * np.exp(1j * k0 * xv) + B * np.exp(-1j * k0 * xv)
        elif xv >= L:
            # hz at back face is the outgoing value; continue as exp(i k0 (x-L))
            out[i] = hz * np.exp(1j * k0 * (xv - L))
        else:
            placed = False
            for x0, x1, eps, hz_l, s_l, k in pieces:
                if x0 <= xv <= x1 + 1e-12:
                    d = xv - x0
                    kd = k * d
                    out[i] = np.cos(kd) * hz_l + (eps / k) * np.sin(kd) * s_l
                    placed = True
                    break
            if not placed:
                out[i] = np.nan
    return out


def mie_T(k0: float, eps: complex, radius: float, nmax: int) -> np.ndarray:
    """Diagonal T_n, n = -nmax..nmax, for a bare cylinder in air. b_n = T_n a_n for incident J_n."""
    ns = np.arange(-nmax, nmax + 1)
    T = np.zeros(len(ns), dtype=np.complex128)
    kp = k0 * np.sqrt(complex(eps))
    for i, n in enumerate(ns):
        n = int(n)
        j = jv(n, k0 * radius)
        jp = jvp(n, k0 * radius, 1)
        h = hankel1(n, k0 * radius)
        hp = h1vp(n, k0 * radius, 1)
        jin = jv(n, kp * radius)
        jnp = jvp(n, kp * radius, 1)
        right = (kp / eps) * (jnp / jin)
        T[i] = (right * j - k0 * jp) / (k0 * hp - right * h)
    return T


def coated_T(k0: float, radii: list[float], eps: list[complex], nmax: int) -> np.ndarray:
    """T_n for concentric layers. radii and eps are core → shell, exterior air is implicit.

    radii: [R1, R2, ..., Rm] interfaces. eps: [eps_core, eps_layer2, ..., eps_outer_shell], same length.
    """
    if len(radii) != len(eps):
        raise ValueError("radii and eps must match")
    ns = np.arange(-nmax, nmax + 1)
    T = np.zeros(len(ns), dtype=np.complex128)
    for i, n in enumerate(ns):
        T[i] = _coated_one(k0, radii, eps, int(n))
    return T


def _coated_one(k0: float, radii: list[float], eps: list[complex], n: int) -> complex:
    """Solve multilayer matching for one azimuthal order. Returns the exterior T_n."""
    # Regions: 0..M-1 shells, plus exterior air. Unknowns: core A, then (B,C) per shell after core, plus S.
    # Region 0: A J only. Regions 1..M-1: B J + C H. Exterior: J + S H.
    m = len(radii)
    # unknowns: A, then 2 per additional interior region, then S. Interior regions after core: m-1.
    n_unknown = 1 + 2 * (m - 1) + 1
    mat = np.zeros((n_unknown, n_unknown), dtype=np.complex128)
    rhs = np.zeros(n_unknown, dtype=np.complex128)
    # Map region j=0..m-1. For j>=1, coeff index of J is 1+2*(j-1), H is that+1. S is last.
    def idx_J(region: int) -> int:
        if region == 0:
            return 0
        return 1 + 2 * (region - 1)

    row = 0
    for interface, R in enumerate(radii):
        left = interface
        right_is_exterior = interface == m - 1
        eps_l = complex(eps[left])
        k_l = k0 * np.sqrt(eps_l)
        if right_is_exterior:
            eps_r = 1.0 + 0j
            k_r = k0 * 1.0
        else:
            eps_r = complex(eps[left + 1])
            k_r = k0 * np.sqrt(eps_r)
        # Hz continuity and (k/ε) derivative continuity.
        for kind in ("value", "deriv"):
            if left == 0:
                z = k_l * R
                basis = jv(n, z) if kind == "value" else (k_l / eps_l) * jvp(n, z, 1)
                mat[row, 0] += basis
            else:
                j_i = idx_J(left)
                zl = k_l * R
                if kind == "value":
                    mat[row, j_i] += jv(n, zl)
                    mat[row, j_i + 1] += hankel1(n, zl)
                else:
                    mat[row, j_i] += (k_l / eps_l) * jvp(n, zl, 1)
                    mat[row, j_i + 1] += (k_l / eps_l) * h1vp(n, zl, 1)
            if right_is_exterior:
                zr = k_r * R
                if kind == "value":
                    mat[row, -1] -= hankel1(n, zr)
                    rhs[row] += jv(n, zr)
                else:
                    mat[row, -1] -= k_r * h1vp(n, zr, 1)
                    rhs[row] += k_r * jvp(n, zr, 1)
            else:
                j_i = idx_J(left + 1)
                zr = k_r * R
                if kind == "value":
                    mat[row, j_i] -= jv(n, zr)
                    mat[row, j_i + 1] -= hankel1(n, zr)
                else:
                    mat[row, j_i] -= (k_r / eps_r) * jvp(n, zr, 1)
                    mat[row, j_i + 1] -= (k_r / eps_r) * h1vp(n, zr, 1)
            row += 1
    sol = np.linalg.solve(mat, rhs)
    return complex(sol[-1])


def addition_check(k0: float = 1.3, nmax: int = 12) -> float:
    """Graf check: H0(k|r-c|) = sum_m J_m(kρ) H_m(kR) e^{im(θ-φ)} for ρ<R."""
    c = np.array([2.0, 0.4])
    p = np.array([0.2, -0.1])
    R = np.linalg.norm(c)
    phi = np.arctan2(c[1], c[0])
    rho = np.linalg.norm(p)
    theta = np.arctan2(p[1], p[0])
    exact = hankel1(0, k0 * np.linalg.norm(p - c))
    acc = 0j
    for m in range(-nmax, nmax + 1):
        acc += jv(m, k0 * rho) * hankel1(m, k0 * R) * np.exp(1j * m * (theta - phi))
    err0 = abs(acc - exact)
    # Higher-order check. φ is the angle of (wave center − expansion center).
    worst = err0
    for n in (1, -2, 3):
        dvec = p - c
        exact_n = hankel1(n, k0 * np.linalg.norm(dvec)) * np.exp(1j * n * np.arctan2(dvec[1], dvec[0]))
        acc_n = 0j
        for m in range(-nmax, nmax + 1):
            acc_n += jv(m, k0 * rho) * np.exp(1j * m * theta) * hankel1(m - n, k0 * R) * np.exp(-1j * (m - n) * phi)
        worst = max(worst, abs(acc_n - exact_n))
    return float(worst)


def solve_clusters(centers: np.ndarray, k0: float, radius: float, eps: complex, source: np.ndarray, nmax: int,
                   coated: tuple[list[float], list[complex]] | None = None) -> np.ndarray:
    """Scattered coefficients b[l, n_index] for identical cylinders and a unit H0 point source.

    Incident field is H0(k|r-source|), expanded with the Graf formula validated above.
    """
    centers = np.asarray(centers, float)
    n_obj = len(centers)
    ns = np.arange(-nmax, nmax + 1)
    n_mode = len(ns)
    if coated is None:
        T = mie_T(k0, eps, radius, nmax)
    else:
        T = coated_T(k0, coated[0], coated[1], nmax)
    # System: b_m^l - T_m sum_{j≠l} sum_n G_{mn}^{lj} b_n^j = T_m a_m^l
    dim = n_obj * n_mode
    M = np.eye(dim, dtype=np.complex128)
    rhs = np.zeros(dim, dtype=np.complex128)
    ext = np.zeros((n_obj, n_mode), dtype=np.complex128)
    for l in range(n_obj):
        # Angle of the source as seen from the cylinder. The Graf check uses this direction.
        d = source - centers[l]
        R = np.linalg.norm(d)
        phi = np.arctan2(d[1], d[0])
        for im, m in enumerate(ns):
            ext[l, im] = hankel1(int(m), k0 * R) * np.exp(-1j * int(m) * phi)
            rhs[l * n_mode + im] = T[im] * ext[l, im]
    for l in range(n_obj):
        for j in range(n_obj):
            if j == l:
                continue
            d = centers[j] - centers[l]
            R = np.linalg.norm(d)
            phi = np.arctan2(d[1], d[0])
            for im, m in enumerate(ns):
                row = l * n_mode + im
                for inn, n in enumerate(ns):
                    # Verified against a direct Hankel evaluation: coefficient of
                    # J_m(k ρ_l) e^{im θ_l} in H_n(k ρ_j) e^{in θ_j}.
                    G = hankel1(int(m - n), k0 * R) * np.exp(-1j * int(m - n) * phi)
                    M[row, j * n_mode + inn] -= T[im] * G
    solve_clusters.last_cond = float(np.linalg.cond(M))
    b = np.linalg.solve(M, rhs).reshape(n_obj, n_mode)
    return b


def cluster_field(xy: np.ndarray, centers: np.ndarray, k0: float, radius: float, eps: complex,
                  source: np.ndarray, b: np.ndarray, nmax: int) -> np.ndarray:
    """Exterior field. Points inside a cylinder are left nan; use a near-field routine there."""
    ns = np.arange(-nmax, nmax + 1)
    xy = np.asarray(xy, float)
    out = hankel1(0, k0 * np.linalg.norm(xy - source, axis=1))
    for l, c in enumerate(np.asarray(centers, float)):
        d = xy - c
        rho = np.linalg.norm(d, axis=1)
        theta = np.arctan2(d[:, 1], d[:, 0])
        inside = rho <= radius
        for im, n in enumerate(ns):
            hn = np.zeros(len(xy), dtype=np.complex128)
            hn[~inside] = hankel1(int(n), k0 * rho[~inside])
            out = out + b[l, im] * hn * np.exp(1j * int(n) * theta)
        out = np.where(inside, np.nan, out)
    return out


def self_checks() -> dict:
    # Air Fresnel.
    fr = fresnel_ht(1.7, 1.0 + 0j, 1.0 + 0j, 0.0)
    # Dielectric normal, power conservation.
    fd = fresnel_ht(1.7, 1.0 + 0j, 3.8 + 0j, 0.0)
    # Negative epsilon.
    fn = fresnel_ht(1.7, 1.0 + 0j, -3.3 + 0.001j, 0.3)
    # Slab of air must be invisible.
    air = stack_response(1.7, [(0.4, 1.0 + 0j)])
    # Bare Mie air cylinder.
    Tair = mie_T(1.7, 1.0 + 0j, 0.23, 8)
    # Coated with everything air matches a trivial scatterer.
    Tcoat_air = coated_T(1.7, [0.23, 0.325, 0.375], [1, 1, 1], 6)
    # Coated with shell and gap equal to the core epsilon matches a bare cylinder of the outer radius.
    eps = -3.3 + 0.001j
    Tbare = mie_T(1.7, eps, 0.375, 8)
    Tcoat = coated_T(1.7, [0.23, 0.325, 0.375], [eps, eps, eps], 8)
    # Quartz removed (eps=1) matches the bare plasma core. Gap is already air.
    Tcore = mie_T(1.7, eps, 0.23, 8)
    Tgap = coated_T(1.7, [0.23, 0.325, 0.375], [eps, 1.0, 1.0], 8)
    graf = addition_check()
    # One-cylinder cluster solver matches Mie at a test point.
    k0 = 1.7
    b = solve_clusters(np.array([[0.0, 0.0]]), k0, 0.23, eps, np.array([-2.0, 0.0]), 8)
    Tm = mie_T(k0, eps, 0.23, 8)
    d = np.array([-2.0, 0.0]) - np.array([0.0, 0.0])
    R = np.linalg.norm(d)
    phi = np.arctan2(d[1], d[0])
    ns = np.arange(-8, 9)
    ext = np.array([hankel1(int(m), k0 * R) * np.exp(-1j * int(m) * phi) for m in ns])
    mie_vs_cluster = float(np.max(np.abs(b[0] - Tm * ext)))
    return {
        "fresnel_air_r": abs(fr["r"]),
        "fresnel_dielectric_power_sum": fd["power_sum"],
        "fresnel_plasma_r_abs": abs(fn["r"]),
        "air_slab_r": abs(air["r"]),
        "air_slab_t": abs(air["t_over_vacuum"] - 1),
        "bare_air_T": float(np.max(np.abs(Tair))),
        "coated_air_T": float(np.max(np.abs(Tcoat_air))),
        "coated_equals_bare_outer": float(np.max(np.abs(Tcoat - Tbare))),
        "coated_air_shell_equals_core": float(np.max(np.abs(Tgap - Tcore))),
        "graf_error": graf,
        "cluster_reduces_to_mie": mie_vs_cluster,
    }


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    rec = self_checks()
    print(json.dumps(rec, indent=2))
    (OUT / "analytic_self_checks.json").write_text(json.dumps(rec, indent=2) + "\n")
