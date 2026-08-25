#!/usr/bin/env python3
"""
Homogeneous magnetized-plasma Faraday rotation benchmark.

Uses the PlasMEEP Get_Med() -> GyrotropicDrudeSusceptibility material path in a
Meep tutorial-style geometry: linear Ex polarization, propagation along +z,
bias along +z.  This isolates gyrotropic material physics from the six-port
circulator geometry.

Default frequencies are chosen so both circular eigenmodes propagate.  The
circulator notebook point (fs=3.85 GHz, fp=8 GHz, B=0.05 T) has Re(eps+/-)<0
and is rejected by the pre-flight eigenpermittivity check.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import meep as mp
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from plasmeep.lib import Plasmeep as pm  # noqa: E402

C = 299_792_458.0
E_CHARGE = 1.602_176_634e-19
ME = 9.109_383_7015e-31

# Keep PlasMEEP length unit consistent with the circulator notebook.
DEFAULT_A = 0.028

# Propagating Faraday validation point (NOT the notebook 3.85/8 GHz case).
# At |B|=0.05 T, fc ~ 1.40 GHz; both R/L eigenpermittivities are > 0.
DEFAULT_FS_HZ = 5.0e9
DEFAULT_FP_HZ = 2.0e9
DEFAULT_GAMMA_HZ = 1.0e6
DEFAULT_B0_T = 0.05

# Require Re(eps+/-) above this for a meaningful propagating rotation test.
MIN_RE_EPS = 0.05

# Homogeneous-material PASS gates.
B0_ABS_TOL = 0.002
B0_REL_TOL = 0.03
ASYMMETRY_TOL = 0.10
THEORY_MAG_TOL = 0.15
FIT_RMSE_FRAC_OF_SPAN = 0.05
FIT_R2_MIN = 0.98


def git_head() -> str:
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=ROOT,
                stderr=subprocess.DEVNULL,
                text=True,
            ).strip()
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def nondimensionalize(a: float, hz_values: np.ndarray) -> np.ndarray:
    helper = pm(a, 32, 1.0, 10.0, 10.0, B=np.zeros(3))
    return helper.Nondimensionalize_Freq(hz_values)


def build_material(a: float, fp_a: float, gamma_a: float, bz_t: float):
    helper = pm(
        a,
        32,
        1.0,
        10.0,
        10.0,
        B=np.array([0.0, 0.0, bz_t]),
    )
    return helper.Get_Med(eps=1.0, wp=fp_a, gamma=gamma_a)


def gyrotropic_drude_eps_eta(
    fs_a: float,
    fp_a: float,
    gamma_a: float,
    bias_a: float,
    eps_inf: float = 1.0,
    sigma: float = 1.0,
) -> tuple[complex, complex]:
    """
    Meep gyrotropic Drude tensor components in Meep ordinary-frequency units.

    From Meep Materials (Gyrotropic Drude-Lorentz model), the polarization ODE is

        d²P/dt² + γ dP/dt - (dP/dt)×b + ω_n² P = σ ω_n² E          (Lorentz)

    For the Drude form the restoring term ω_n² P on the LHS is omitted; the RHS
    still uses the Meep `frequency` parameter as ω_n in σ ω_n² (plasma scale).

    Frequency-domain Lorentz susceptibilities (bias along z):

        Δ_n = ω_n² - ω² - i ω γ
        χ_⊥ = ω_n² Δ_n σ / (Δ_n² - ω² b²)
        η   = ω_n² ω b σ / (Δ_n² - ω² b²)

    Drude limit (omit ω_n² on LHS only):

        Δ = -ω² - i ω γ
        χ_⊥ = ω_p² Δ σ / (Δ² - ω² b²)
        η   = ω_p² ω b σ / (Δ² - ω² b²)

    Meep's Python/tutorial algebra uses ordinary frequencies f = ω/2π (and the
    same units for `gamma` and bias magnitude |b|/2π).  Equivalent form:

        Δ_f = -f² - i f γ
        ε_⊥ = ε_∞ + σ f_p² Δ_f / (Δ_f² - (f b)²)
        η   = σ f_p² f b / (Δ_f² - (f b)²)

    Algebraically this is the cold-plasma result with the 1/ω factors:

        ε_⊥ = ε_∞ - ω_p² (ω + iγ) / [ω ((ω + iγ)² - ω_c²)]
        η   =         ω_p² ω_c     / [ω ((ω + iγ)² - ω_c²)]

    with ω_c = 2π b.  The previous benchmark omitted those 1/ω denominators.
    """
    delta_f = -fs_a**2 - 1j * fs_a * gamma_a
    denom = delta_f**2 - (fs_a * bias_a) ** 2
    eps_perp = eps_inf + sigma * fp_a**2 * delta_f / denom
    eta = sigma * fp_a**2 * fs_a * bias_a / denom
    return eps_perp, eta


def circular_eigenpermittivities(
    eps_perp: complex,
    eta: complex,
) -> tuple[complex, complex]:
    """
    Circular eigenpermittivities for k || B.

    For ε = [[ε_⊥, -iη], [iη, ε_⊥]]:
      Ey = +i Ex  ->  ε_+ = ε_⊥ + η
      Ey = -i Ex  ->  ε_- = ε_⊥ - η
    """
    return eps_perp + eta, eps_perp - eta


def propagation_constants(
    fs_a: float,
    eps_plus: complex,
    eps_minus: complex,
) -> tuple[complex, complex]:
    """
    Complex k± in Meep inverse-length units: k = 2π f √ε.

    Branch: Re(k) >= 0.
    """
    def k_of(eps: complex) -> complex:
        root = np.sqrt(eps)
        k = 2.0 * np.pi * fs_a * root
        if np.real(k) < 0.0:
            k = -k
        return complex(k)

    return k_of(eps_plus), k_of(eps_minus)


def signed_faraday_kappa(k_plus: complex, k_minus: complex) -> complex:
    """
    Linear-polarization rotation rate κ from circular-mode beat:

        κ = (k_+ - k_-) / 2

    Sign reverses when B -> -B because η -> -η swaps ε_+ and ε_-.
    """
    return 0.5 * (k_plus - k_minus)


def modes_are_propagating(
    eps_plus: complex,
    eps_minus: complex,
    k_plus: complex,
    k_minus: complex,
    min_re_eps: float = MIN_RE_EPS,
) -> tuple[bool, str]:
    """Abort criterion: both circular modes must support propagating waves."""
    if np.real(eps_plus) < min_re_eps or np.real(eps_minus) < min_re_eps:
        return (
            False,
            (
                f"evanescent/non-propagating circular eigenpermittivity: "
                f"Re(eps+)={np.real(eps_plus):.4f}, Re(eps-)={np.real(eps_minus):.4f} "
                f"(need both > {min_re_eps})"
            ),
        )

    for label, k in (("+", k_plus), ("-", k_minus)):
        re_k = float(np.real(k))
        im_k = float(np.imag(k))
        if re_k <= 0.0:
            return False, f"Re(k{label})={re_k:.4f} is not positive"
        if abs(im_k) > 0.5 * re_k:
            return (
                False,
                (
                    f"|Im(k{label})|={abs(im_k):.4f} is comparable to "
                    f"Re(k{label})={re_k:.4f}; rotation fit would be unreliable"
                ),
            )

    return True, "both circular modes propagate"


def theory_bundle(fs_a: float, fp_a: float, gamma_a: float, bias_a: float) -> dict:
    eps_perp, eta = gyrotropic_drude_eps_eta(fs_a, fp_a, gamma_a, bias_a)
    eps_plus, eps_minus = circular_eigenpermittivities(eps_perp, eta)
    k_plus, k_minus = propagation_constants(fs_a, eps_plus, eps_minus)
    kappa = signed_faraday_kappa(k_plus, k_minus)
    ok, reason = modes_are_propagating(eps_plus, eps_minus, k_plus, k_minus)
    return {
        "eps_perp": complex(eps_perp),
        "eta": complex(eta),
        "eps_plus": complex(eps_plus),
        "eps_minus": complex(eps_minus),
        "k_plus": complex(k_plus),
        "k_minus": complex(k_minus),
        "kappa_complex": complex(kappa),
        "kappa_rad_per_a": float(np.real(kappa)),
        "propagating": ok,
        "propagating_reason": reason,
    }


def stokes_orientation_unwrapped(ex: np.ndarray, ey: np.ndarray) -> np.ndarray:
    """
    Linear-polarization orientation from complex Ex, Ey via Stokes parameters.

    Orientation is π-periodic.  Compute

        2ψ = atan2(S2, S1)

    with S1 = |Ex|² - |Ey|² and S2 = 2 Re(Ex Ey*), unwrap that 2π-periodic
    quantity, then divide by 2.
    """
    s1 = np.abs(ex) ** 2 - np.abs(ey) ** 2
    s2 = 2.0 * np.real(ex * np.conj(ey))
    two_psi = np.arctan2(s2, s1)
    return 0.5 * np.unwrap(two_psi)


def fit_orientation_vs_z(
    z: np.ndarray,
    ex: np.ndarray,
    ey: np.ndarray,
    dpml_a: float,
    src_z: float,
) -> dict:
    """Unwrap Stokes orientation vs z and fit κ = dψ/dz in the interior."""
    psi = stokes_orientation_unwrapped(ex, ey)

    z_lo = -0.5 * (z.max() - z.min()) + dpml_a + 1.5
    z_hi = 0.5 * (z.max() - z.min()) - dpml_a - 0.5
    mask = (z > max(z_lo, src_z + 1.0)) & (z < z_hi)

    if np.count_nonzero(mask) < 8:
        raise RuntimeError("Not enough interior DFT samples for orientation fit.")

    z_fit = z[mask]
    psi_fit = psi[mask]
    slope, intercept = np.polyfit(z_fit, psi_fit, 1)
    resid = psi_fit - (slope * z_fit + intercept)
    rmse = float(np.sqrt(np.mean(resid**2)))
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((psi_fit - np.mean(psi_fit)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0.0 else float("nan")
    psi_span = float(np.max(psi_fit) - np.min(psi_fit))

    return {
        "kappa_rad_per_a": float(slope),
        "kappa_deg_per_a": float(np.degrees(slope)),
        "fit_intercept_rad": float(intercept),
        "fit_rmse_rad": rmse,
        "fit_r2": r2,
        "psi_span_rad": psi_span,
        "z_fit_span_a": [float(z_fit.min()), float(z_fit.max())],
        "n_fit_points": int(z_fit.size),
    }


def run_faraday_case(
    *,
    a: float,
    res: int,
    fs_a: float,
    fp_a: float,
    gamma_a: float,
    bz_t: float,
    length_a: float,
    dpml_a: float,
    run_time: float,
    src_offset_a: float,
) -> dict:
    material = build_material(a, fp_a, gamma_a, bz_t)

    cell = mp.Vector3(0.0, 0.0, length_a)
    pml = [mp.PML(thickness=dpml_a, direction=mp.Z)]

    src_z = -0.5 * length_a + dpml_a + src_offset_a
    sources = [
        mp.Source(
            mp.ContinuousSource(frequency=fs_a, is_integrated=True),
            component=mp.Ex,
            center=mp.Vector3(0.0, 0.0, src_z),
        )
    ]

    sim = mp.Simulation(
        cell_size=cell,
        geometry=[],
        sources=sources,
        boundary_layers=pml,
        default_material=material,
        resolution=res,
    )

    # ContinuousSource does not turn off: stop with finite until=, not
    # until_after_sources.
    dft_vol = mp.Volume(
        center=mp.Vector3(0.0, 0.0, 0.0),
        size=mp.Vector3(0.0, 0.0, length_a),
    )
    # Meep 1.30 accepts either (fcen, df, nfreq) or a freq list/array in *args.
    dft_fields = sim.add_dft_fields([mp.Ex, mp.Ey], [fs_a], where=dft_vol)

    sim.run(until=run_time)

    ex = np.asarray(sim.get_dft_array(dft_fields, mp.Ex, 0)).astype(np.complex128).flatten()
    ey = np.asarray(sim.get_dft_array(dft_fields, mp.Ey, 0)).astype(np.complex128).flatten()
    _x, _y, z, _w = sim.get_array_metadata(dft_cell=dft_fields)
    z = np.asarray(z, dtype=float).flatten()

    if ex.size != z.size or ey.size != z.size:
        raise RuntimeError(
            f"DFT field/geometry size mismatch: Ex={ex.size}, Ey={ey.size}, z={z.size}"
        )
    if z.size < 2 or not np.all(np.diff(z) > 0.0):
        raise RuntimeError(
            "DFT z metadata is not strictly monotonic increasing; "
            f"z[:5]={z[:5]}, z[-5:]={z[-5:]}"
        )

    fit = fit_orientation_vs_z(z, ex, ey, dpml_a, src_z)
    fit["Bz_T"] = bz_t
    return fit


def assess_results(
    results: dict,
    kappa_theory_plus: float,
    kappa_theory_minus: float,
) -> tuple[str, dict]:
    """
    PASS/FAIL uses convention-independent magnitude and B-reversal checks.

    Absolute measured κ vs theory sign is reported as a global mapping, not a gate.
    """
    k0 = results["B=0"]["kappa_rad_per_a"]
    kp = results["+B"]["kappa_rad_per_a"]
    km = results["-B"]["kappa_rad_per_a"]

    theory_mag = max(abs(kappa_theory_plus), abs(kappa_theory_minus), 1e-12)
    meas_mag = 0.5 * (abs(kp) + abs(km))
    scale = max(theory_mag, meas_mag, 1e-12)

    diagnostics = {
        "scale_rad_per_a": scale,
        "b0_residual": abs(k0),
        "b0_tol": max(B0_ABS_TOL, B0_REL_TOL * scale),
        "asymmetry": abs(abs(kp) - abs(km)) / scale,
        "theory_mag_error": abs(meas_mag - theory_mag) / theory_mag,
        "reversal_ok": bool(np.sign(kp) == -np.sign(km) and abs(kp) > 0.0 and abs(km) > 0.0),
        "global_sign_map": None,
    }

    # Global sign mapping: measured κ(+B) vs theory κ(+B). Not a PASS gate.
    if abs(kp) > 0.0 and abs(kappa_theory_plus) > 0.0:
        if np.sign(kp) == np.sign(kappa_theory_plus):
            diagnostics["global_sign_map"] = "measured_matches_theory_kappa_sign"
        else:
            diagnostics["global_sign_map"] = "measured_opposite_theory_kappa_sign"
    else:
        diagnostics["global_sign_map"] = "undefined"

    failures = []

    if abs(k0) > diagnostics["b0_tol"]:
        failures.append(
            f"B=0 residual {abs(k0):.6g} > tol {diagnostics['b0_tol']:.6g}"
        )

    if not diagnostics["reversal_ok"]:
        failures.append("κ(+B) and κ(-B) do not reverse sign")

    if diagnostics["asymmetry"] > ASYMMETRY_TOL:
        failures.append(
            f"±B |κ| asymmetry {diagnostics['asymmetry']:.3f} > {ASYMMETRY_TOL}"
        )

    if diagnostics["theory_mag_error"] > THEORY_MAG_TOL:
        failures.append(
            f"|κ| vs theory error {diagnostics['theory_mag_error']:.3f} > {THEORY_MAG_TOL}"
        )

    for label in ("B=0", "+B", "-B"):
        fit = results[label]
        rmse = fit["fit_rmse_rad"]
        r2 = fit["fit_r2"]
        span = max(abs(fit["psi_span_rad"]), 1e-12)
        # B=0 has near-zero span; require small absolute RMSE instead of R².
        if label == "B=0":
            if rmse > max(B0_ABS_TOL, 0.01):
                failures.append(f"{label} fit RMSE {rmse:.6g} too large for flat ψ(z)")
        else:
            if rmse > FIT_RMSE_FRAC_OF_SPAN * span:
                failures.append(
                    f"{label} fit RMSE {rmse:.6g} > "
                    f"{FIT_RMSE_FRAC_OF_SPAN} * ψ-span {span:.6g}"
                )
            if not np.isfinite(r2) or r2 < FIT_R2_MIN:
                failures.append(f"{label} fit R² {r2:.6g} < {FIT_R2_MIN}")

    diagnostics["failures"] = failures
    status = "PASS" if not failures else "FAIL/INVESTIGATE"
    return status, diagnostics


def _cfmt(z: complex) -> str:
    return f"{np.real(z):+.6f}{np.imag(z):+.6f}j"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=float, default=DEFAULT_A)
    parser.add_argument("--fs-hz", type=float, default=DEFAULT_FS_HZ)
    parser.add_argument("--fp-hz", type=float, default=DEFAULT_FP_HZ)
    parser.add_argument("--gamma-hz", type=float, default=DEFAULT_GAMMA_HZ)
    parser.add_argument("--b0-t", type=float, default=DEFAULT_B0_T)
    parser.add_argument("--res", type=int, default=64)
    parser.add_argument("--length-a", type=float, default=20.0)
    parser.add_argument("--dpml-a", type=float, default=1.0)
    parser.add_argument("--run-time", type=float, default=100.0)
    parser.add_argument("--json-out", type=str, default="")
    args = parser.parse_args()

    freqs_hz = np.array([args.fs_hz, args.fp_hz, args.gamma_hz])
    fs_a, fp_a, gamma_a = nondimensionalize(args.a, freqs_hz)
    fc_hz = E_CHARGE * args.b0_t / (2.0 * np.pi * ME)
    fc_a = float(nondimensionalize(args.a, np.array([fc_hz]))[0])

    # PlasMEEP stores Meep bias in ordinary-frequency units (fc_a), matching
    # GyrotropicDrudeSusceptibility / Faraday tutorial conventions.
    theory_plus = theory_bundle(fs_a, fp_a, gamma_a, +fc_a)
    theory_minus = theory_bundle(fs_a, fp_a, gamma_a, -fc_a)
    theory_zero = theory_bundle(fs_a, fp_a, gamma_a, 0.0)

    print("# Faraday Benchmark")
    print(f"Commit: {git_head()}")
    print(f"Meep version: {mp.__version__}")
    print(f"Resolution: {args.res}")
    print(f"a = {args.a} m")
    print(f"fs = {args.fs_hz/1e9:.4f} GHz -> {fs_a:.6f} c/a")
    print(f"fp = {args.fp_hz/1e9:.4f} GHz -> {fp_a:.6f} c/a")
    print(f"gamma = {args.gamma_hz/1e6:.4f} MHz -> {gamma_a:.6e} c/a")
    print(f"|B0| = {args.b0_t:.4f} T -> fc = {fc_a:.6f} c/a")
    print()
    print("## Circular eigenpermittivity pre-flight")
    for label, bundle in (
        ("B=0", theory_zero),
        ("+B", theory_plus),
        ("-B", theory_minus),
    ):
        print(f"  [{label}]")
        print(f"    eps+ = {_cfmt(bundle['eps_plus'])}")
        print(f"    eps- = {_cfmt(bundle['eps_minus'])}")
        print(f"    k+   = {_cfmt(bundle['k_plus'])}")
        print(f"    k-   = {_cfmt(bundle['k_minus'])}")
        print(
            f"    theory kappa = {bundle['kappa_rad_per_a']:.6f} rad/a "
            f"({np.degrees(bundle['kappa_rad_per_a']):.4f} deg/a)"
        )
        print(f"    propagating: {bundle['propagating_reason']}")

    abort_labels = [
        label
        for label, bundle in (
            ("B=0", theory_zero),
            ("+B", theory_plus),
            ("-B", theory_minus),
        )
        if not bundle["propagating"]
    ]
    if abort_labels:
        print()
        print("Status: ABORT")
        print(
            "Both circular eigenmodes must propagate for B=0, +B, and -B. "
            f"Failed cases: {', '.join(abort_labels)}"
        )
        print(
            "Note: the circulator notebook point fs=3.85 GHz, fp=8 GHz, "
            "B=0.05 T is evanescent for homogeneous plasma and is not used."
        )
        return 2

    print()
    print(
        f"  theory |kappa|(+/-B) = {abs(theory_plus['kappa_rad_per_a']):.6f} rad/a"
    )
    print()

    shared = dict(
        a=args.a,
        res=args.res,
        fs_a=fs_a,
        fp_a=fp_a,
        gamma_a=gamma_a,
        length_a=args.length_a,
        dpml_a=args.dpml_a,
        run_time=args.run_time,
        src_offset_a=0.5,
    )

    results = {}
    for label, bz in [("B=0", 0.0), ("+B", +args.b0_t), ("-B", -args.b0_t)]:
        print(f"Running {label} (Bz = {bz:+.4f} T)...")
        fit = run_faraday_case(bz_t=bz, **shared)
        results[label] = fit
        print(
            f"  kappa = {fit['kappa_rad_per_a']:.6f} rad/a "
            f"({fit['kappa_deg_per_a']:.4f} deg/a), "
            f"RMSE = {fit['fit_rmse_rad']:.6g} rad, "
            f"R² = {fit['fit_r2']:.6f}"
        )

    status, diagnostics = assess_results(
        results,
        theory_plus["kappa_rad_per_a"],
        theory_minus["kappa_rad_per_a"],
    )
    k0 = results["B=0"]["kappa_rad_per_a"]
    kp = results["+B"]["kappa_rad_per_a"]
    km = results["-B"]["kappa_rad_per_a"]
    theory_mag = abs(theory_plus["kappa_rad_per_a"])
    meas_mag = 0.5 * (abs(kp) + abs(km))

    print()
    print("## Summary")
    print(f"B=0 kappa: {k0:.6f} rad/a  (theory {theory_zero['kappa_rad_per_a']:.6f})")
    print(
        f"+B kappa:  {kp:.6f} rad/a  "
        f"(theory signed {theory_plus['kappa_rad_per_a']:.6f})"
    )
    print(
        f"-B kappa:  {km:.6f} rad/a  "
        f"(theory signed {theory_minus['kappa_rad_per_a']:.6f})"
    )
    print(f"Sign reversal (+B vs -B): {diagnostics['reversal_ok']}")
    print(f"Global sign mapping (not a PASS gate): {diagnostics['global_sign_map']}")
    print(f"|+B| vs |-B| relative asymmetry: {diagnostics['asymmetry']:.3f}")
    print(f"Measured |κ| = {meas_mag:.6f} rad/a")
    print(f"Theory   |κ| = {theory_mag:.6f} rad/a")
    print(f"|κ| vs theory relative error: {diagnostics['theory_mag_error']:.3f}")
    print(
        f"B=0 residual / tol: "
        f"{diagnostics['b0_residual']:.6g} / {diagnostics['b0_tol']:.6g}"
    )
    if diagnostics["failures"]:
        print("Gate failures:")
        for item in diagnostics["failures"]:
            print(f"  - {item}")
    print(f"Status: {status}")

    def cjson(z: complex) -> dict:
        return {"re": float(np.real(z)), "im": float(np.imag(z))}

    payload = {
        "commit": git_head(),
        "meep_version": mp.__version__,
        "parameters": {
            "a_m": args.a,
            "fs_hz": args.fs_hz,
            "fp_hz": args.fp_hz,
            "gamma_hz": args.gamma_hz,
            "b0_t": args.b0_t,
            "fc_hz": fc_hz,
            "res": args.res,
            "length_a": args.length_a,
            "run_time": args.run_time,
            "note": (
                "Defaults are a propagating Faraday point; "
                "notebook 3.85/8 GHz homogeneous plasma is evanescent."
            ),
        },
        "gates": {
            "b0_abs_tol": B0_ABS_TOL,
            "b0_rel_tol": B0_REL_TOL,
            "asymmetry_tol": ASYMMETRY_TOL,
            "theory_mag_tol": THEORY_MAG_TOL,
            "fit_rmse_frac_of_span": FIT_RMSE_FRAC_OF_SPAN,
            "fit_r2_min": FIT_R2_MIN,
        },
        "theory": {
            "+B": {
                "eps_perp": cjson(theory_plus["eps_perp"]),
                "eta": cjson(theory_plus["eta"]),
                "eps_plus": cjson(theory_plus["eps_plus"]),
                "eps_minus": cjson(theory_plus["eps_minus"]),
                "k_plus": cjson(theory_plus["k_plus"]),
                "k_minus": cjson(theory_plus["k_minus"]),
                "kappa_rad_per_a": theory_plus["kappa_rad_per_a"],
            },
            "-B": {
                "kappa_rad_per_a": theory_minus["kappa_rad_per_a"],
            },
            "B=0": {
                "kappa_rad_per_a": theory_zero["kappa_rad_per_a"],
            },
        },
        "results": results,
        "diagnostics": diagnostics,
        "status": status,
    }

    if args.json_out:
        out_dir = os.path.dirname(os.path.abspath(args.json_out))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.json_out, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
        print(f"Wrote {args.json_out}")

    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
