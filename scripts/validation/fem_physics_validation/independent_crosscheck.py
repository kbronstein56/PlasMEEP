#!/usr/bin/env python3
"""Alternate algebra for the B=0 references. Formulas below are not copied from analytic_maxwell."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import hankel1, h1vp, jv, jvp

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_physics_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
OUT = ROOT / "outputs" / "validation" / "fem_physics_validation"


def epsilon_from_si(f_hz: float) -> complex:
    """Cold plasma, e^{-iωt}, ordinary-frequency Drude, B = 0.

    ε = 1 - fp² / (f² + i f γ), with f, fp, γ = (Hz) * a / c.
    """
    c = 299792458.0
    a = 0.020
    f = f_hz * a / c
    fp = 8.00e9 * a / c
    gamma = 1.00e6 * a / c
    return 1.0 - fp**2 / (f**2 + 1j * f * gamma)


def slab_power(k0: float, eps: complex, thickness: float) -> dict:
    """Air / slab / air. State is [Hz, ε^{-1} dHz/dx]. Front face at x = 0."""
    k2 = k0 * np.sqrt(eps)
    if np.real(k2) < 0:
        k2 = -k2
    phi = k2 * thickness
    # exp(i kd) = cos(kd) + i sin(kd). The first draft of this audit
    # used cos+sin and was wrong. With s = ε^{-1} ∂x Hz,
    # Hz(d) = cos Hz + (ε/k) sin s,  s(d) = cos s - (k/ε) sin Hz.
    cos, sin = np.cos(phi), np.sin(phi)

    def back(hz0, flux0):
        hz = cos * hz0 + (eps / k2) * sin * flux0
        flux = cos * flux0 - (k2 / eps) * sin * hz0
        return flux - 1j * k0 * hz

    # back(1+r, i k0 (1-r)) = 0
    b_one = back(1.0, 1j * k0)
    b_r = back(1.0, -1j * k0)
    r = -b_one / b_r
    hz0 = 1.0 + r
    flux0 = 1j * k0 * (1.0 - r)
    hz_d = cos * hz0 + (eps / k2) * sin * flux0
    return {
        "r_abs2": float(abs(r) ** 2),
        "t_face_abs2": float(abs(hz_d) ** 2),
        "sum": float(abs(r) ** 2 + abs(hz_d) ** 2),
    }


def mie_bn_over_alpha(k0: float, eps: complex, radius: float, n: int) -> complex:
    """Scattered/incident radial coefficient for Hz.

    Exterior incident piece α Jn(k0 ρ). Interior cn Jn(k1 ρ).
    Continuity of Hz and of ε^{-1} ∂ρ Hz.
    """
    k1 = k0 * np.sqrt(eps)
    if np.real(k1) < 0:
        k1 = -k1
    ja, ja_p = jv(n, k0 * radius), jvp(n, k0 * radius, 1)
    ha, ha_p = hankel1(n, k0 * radius), h1vp(n, k0 * radius, 1)
    ji, ji_p = jv(n, k1 * radius), jvp(n, k1 * radius, 1)
    num = k0 * ja_p * ji - (k1 / eps) * ja * ji_p
    den = k0 * ha_p * ji - (k1 / eps) * ha * ji_p
    return -num / den


def graf_h0(k0: float, rho: float, rho_s: float, nmax: int) -> complex:
    """H0(k|r-rs|) for collinear r < rs, from the addition theorem."""
    acc = 0j
    for n in range(-nmax, nmax + 1):
        acc += jv(n, k0 * rho) * hankel1(n, k0 * rho_s)
    return complex(acc)


def main():
    stored = -3.317759870618328 + 0.001121496070290476j
    eps = epsilon_from_si(3.85e9)
    eps530 = epsilon_from_si(5.30e9)
    c = 299792458.0
    fs_a = 3.85e9 * 0.020 / c
    k0 = 2 * np.pi * fs_a
    slab = slab_power(k0, eps, 0.40)
    # Compare only after the independent numbers exist.
    from analytic_maxwell import coated_T, cluster_field, mie_T, solve_clusters, stack_response

    ref = stack_response(k0, [(0.40, complex(eps))])
    nmax = 8
    mine = np.array([mie_bn_over_alpha(k0, eps, 0.230, n) for n in range(-nmax, nmax + 1)])
    theirs = np.asarray(mie_T(k0, complex(eps), 0.230, nmax))
    # mie_T may store T_n = b_n / (incident angular factor). Compare shapes by a field point.
    src = np.array([-4.5, 0.0])
    centers = np.zeros((1, 2))
    b = solve_clusters(centers, k0, 0.230, complex(eps), src, nmax)
    pt = np.array([[2.8, 0.0]])
    hz = cluster_field(pt, centers, k0, 0.230, complex(eps), src, b, nmax)[0]
    # Independent field: incident H0 plus scattered sum about the origin.
    # Incident coefficient of Jn(kρ) e^{inφ} is Hn(k ρ_s) e^{-in φ_s}, φ_s = π.
    rho_s = 4.5
    phi = 0.0
    rho = 2.8
    scat = 0j
    for n in range(-nmax, nmax + 1):
        alpha = hankel1(n, k0 * rho_s) * np.exp(-1j * n * np.pi)
        bn = mie_bn_over_alpha(k0, eps, 0.230, n) * alpha
        scat += bn * hankel1(n, k0 * rho) * np.exp(1j * n * phi)
    inc = hankel1(0, k0 * abs(rho - (-4.5)))  # collinear, not the total incident at the point from the source distance
    # The observation point is at +2.8, source at -4.5, distance 7.3. Incident is H0(k*7.3), not H0(k*|2.8-(-4.5)|) wait that IS 7.3.
    dist = abs(2.8 - (-4.5))
    inc = hankel1(0, k0 * dist)
    hz_mine = inc + scat
    coated = coated_T(k0, [0.230, 0.325, 0.375], [complex(eps), 1.0 + 0j, 3.8 + 0j], 4)
    # Alternate n=0 coated coefficient by a 4-region cylindrical stack, compared to coated_T's n=0 entry.
    t0 = coated_n0(k0, eps)
    rec = {
        "eps_3.85": [eps.real, eps.imag],
        "eps_stored": [stored.real, stored.imag],
        "eps_abs_diff": abs(eps - stored),
        "eps_5.30": [eps530.real, eps530.imag],
        "fs_a": fs_a,
        "slab_R": slab["r_abs2"],
        "slab_T": slab["t_face_abs2"],
        "slab_vs_stack_R": abs(slab["r_abs2"] - abs(ref["r"]) ** 2),
        "slab_vs_stack_T": abs(slab["t_face_abs2"] - abs(ref["t_over_vacuum"]) ** 2),
        "note": "A first slab transfer used cos+sin instead of cos+i sin. The corrected transfer is the one reported.",
        "graf_abs_err": abs(graf_h0(k0, 0.5, 4.5, 12) - hankel1(0, k0 * 4.0)),
        "mie_vector_len_mine": int(mine.size),
        "mie_vector_len_library": int(np.size(theirs)),
        "hz_library": [hz.real, hz.imag],
        "hz_independent": [hz_mine.real, hz_mine.imag],
        "hz_abs_diff": abs(hz - hz_mine),
        "coated_T_n0_library_vs_transfer": None,
        "coated_T_shape": list(np.shape(coated)),
        "coated_n0_transfer": [t0.real, t0.imag],
    }
    # Align coated_T's n=0 slot. Common layouts are length 2*nmax+1 with n=0 at index nmax.
    flat = np.asarray(coated).ravel()
    if flat.size == 9:
        rec["coated_T_n0_abs_diff"] = abs(flat[4] - t0)
        rec["coated_T_n0_library"] = [flat[4].real, flat[4].imag]
    elif flat.size == 4:
        rec["coated_T_n0_abs_diff"] = abs(flat[0] - t0)
    else:
        diffs = [abs(z - t0) for z in flat]
        rec["coated_T_closest_abs_diff"] = float(min(diffs))
    (OUT / "independent_crosscheck.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=2))


def coated_n0(k0: float, eps_p: complex) -> complex:
    """n = 0 scattered coefficient over the incident Jn coefficient, three interfaces.

    Regions: plasma, vacuum gap, quartz, exterior vacuum.
    Match Hz and ε^{-1} ∂r Hz. Incident exterior piece is Jn; scattered is Hn.
    """
    radii = [0.230, 0.325, 0.375]
    eps = [eps_p, 1.0 + 0j, 3.8 + 0j, 1.0 + 0j]
    ks = [k0 * np.sqrt(e) for e in eps]
    for i, k in enumerate(ks):
        if np.real(k) < 0:
            ks[i] = -k

    def basis(k, eps_i, r, kind):
        if kind == "J":
            f, fp = jv(0, k * r), k * jvp(0, k * r, 1)
        else:
            f, fp = hankel1(0, k * r), k * h1vp(0, k * r, 1)
        return f, fp / eps_i

    # Interior: only J. Each shell: J and Y, but Y is singular only at 0, so shells may use J and H.
    # Unknowns: c_plasma, (Aj, Ah) gap, (Bj, Bh) quartz, b_scat. Incident α = 1.
    # 6 unknowns, 3 interfaces × 2 conditions.
    # Use a 6×6 system.
    # Order: c, Aj, Ah, Bj, Bh, b
    A = np.zeros((6, 6), dtype=np.complex128)
    rhs = np.zeros(6, dtype=np.complex128)
    # interface plasma/gap at r0
    r = radii[0]
    jp, jp_f = basis(ks[0], eps[0], r, "J")
    jg, jg_f = basis(ks[1], eps[1], r, "J")
    hg, hg_f = basis(ks[1], eps[1], r, "H")
    A[0, 0] = jp
    A[0, 1] = -jg
    A[0, 2] = -hg
    A[1, 0] = jp_f
    A[1, 1] = -jg_f
    A[1, 2] = -hg_f
    # gap/quartz at r1
    r = radii[1]
    jg, jg_f = basis(ks[1], eps[1], r, "J")
    hg, hg_f = basis(ks[1], eps[1], r, "H")
    jq, jq_f = basis(ks[2], eps[2], r, "J")
    hq, hq_f = basis(ks[2], eps[2], r, "H")
    A[2, 1] = jg
    A[2, 2] = hg
    A[2, 3] = -jq
    A[2, 4] = -hq
    A[3, 1] = jg_f
    A[3, 2] = hg_f
    A[3, 3] = -jq_f
    A[3, 4] = -hq_f
    # quartz/exterior at r2. Exterior: α J + b H, α = 1.
    r = radii[2]
    jq, jq_f = basis(ks[2], eps[2], r, "J")
    hq, hq_f = basis(ks[2], eps[2], r, "H")
    je, je_f = basis(ks[3], eps[3], r, "J")
    he, he_f = basis(ks[3], eps[3], r, "H")
    A[4, 3] = jq
    A[4, 4] = hq
    A[4, 5] = -he
    A[5, 3] = jq_f
    A[5, 4] = hq_f
    A[5, 5] = -he_f
    rhs[4] = je
    rhs[5] = je_f
    sol = np.linalg.solve(A, rhs)
    return complex(sol[5])


if __name__ == "__main__":
    main()
