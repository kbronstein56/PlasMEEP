#!/usr/bin/env python3
"""Six-port circulator objective, port-power adjoint, and multi-source FEM solve.

Uses the validated FULL91 assembler and vacuum-feed guide-normal power.
Design variables: per-rod density scales s_i (material-only).
"""
from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

import sixport_common as sc  # noqa: E402
from diff_fem_core import assemble_dA_rod, per_rod_plasma_masks, rho_from_scales  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from full91_geometry import bulb_centers_device, production_constants  # noqa: E402
from full91_solve import inject, load_mesh, material_masks, port_lines  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"

# Frozen before optimization (see CIRCULATOR_OBJECTIVE_DERIVATION.md)
W_DES = 1.0
W_REV = 2.0
W_THR = 1.5
W_LEAK = 1.0
W_ACC = 0.25

F_HZ = 3.85e9
B_PLUS = 0.05
S_LO = 0.05
S_HI = 1.80


def port_roles(j: int):
    return {
        "desired": (j + 1) % 6,
        "reverse": (j - 1) % 6,
        "through": (j + 3) % 6,
        "leak_a": (j + 2) % 6,
        "leak_b": (j + 4) % 6,
        "source": j,
    }


@dataclass
class DeviceData:
    grade: str
    points: np.ndarray
    tris: np.ndarray
    rod_masks: list
    quartz: np.ndarray
    walls: np.ndarray
    centers: np.ndarray
    sampler: ElementSampler
    lines: list
    rhs: list
    r_plasma: float


@dataclass
class SolveBundle:
    rho: list
    A: object
    lu: object
    fields: list  # x_j for each source
    k0: float
    f_ord: float
    gamma: float
    fc: float
    fp_ref_ord: float
    s_vec: np.ndarray
    power: np.ndarray  # (6,6) P[i,j]
    P_ref: float
    metrics: dict
    J: float
    factor_s: float = 0.0
    solve_s: float = 0.0


def load_device(grade: str) -> DeviceData:
    _sc, const = production_constants()
    points, tris = load_mesh(grade)
    plasma, quartz, walls, centers = material_masks(points, tris, sc)
    centers = np.asarray(centers, float)
    r_p = const["r_plasma_material_a"]
    rod_masks = per_rod_plasma_masks(points, tris, centers, r_p)
    sampler = ElementSampler(points, tris)
    lines = port_lines(n_points=41)
    rhs = [inject(points, tris, sampler, line["source_xy"], line["weights"]) for line in lines]
    return DeviceData(
        grade=grade,
        points=points,
        tris=tris,
        rod_masks=rod_masks,
        quartz=quartz,
        walls=walls,
        centers=centers,
        sampler=sampler,
        lines=lines,
        rhs=rhs,
        r_plasma=r_p,
    )


def port_power_one(sampler, points, tris, x, rho, k0, line) -> float:
    """Outward guide-normal flux on one monitor line (vacuum formula)."""
    xy = line["monitor_xy"]
    n_hat = np.asarray(line["outward"], float)
    n_hat /= np.linalg.norm(n_hat)
    weights = np.asarray(line["weights"], float)
    hz, ex, ey = sampler.fields(x, *rho, k0, xy)
    sx = 0.5 * np.real(ey * np.conj(hz))
    sy = -0.5 * np.real(ex * np.conj(hz))
    return float(np.sum((n_hat[0] * sx + n_hat[1] * sy) * weights))


def power_matrix(dev: DeviceData, fields, rho, k0) -> np.ndarray:
    P = np.zeros((6, 6), float)
    for j, x in enumerate(fields):
        for i, line in enumerate(dev.lines):
            P[i, j] = port_power_one(dev.sampler, dev.points, dev.tris, x, rho, k0, line)
    return P


def power_adjoint_q(sampler, tris, x, rho, k0, line) -> np.ndarray:
    """Complex vector q with dP = Re(q† dx) for vacuum-feed port power."""
    n = len(sampler.points)
    q = np.zeros(n, np.complex128)
    xy = np.asarray(line["monitor_xy"], float)
    n_hat = np.asarray(line["outward"], float)
    n_hat /= np.linalg.norm(n_hat)
    weights = np.asarray(line["weights"], float)
    alpha = -1j / k0
    tid = sampler.locate(xy)
    for s in range(len(xy)):
        t = int(tid[s])
        if t < 0:
            continue
        nodes = tris[t]
        wbar = sampler._bary(t, xy[s])
        ell = np.asarray(wbar, float)
        bx = sampler.bx[t]
        by = sampler.by[t]
        gn = n_hat[0] * bx + n_hat[1] * by  # length-3
        un = x[nodes]
        hz = complex(np.dot(ell, un))
        dn = complex(np.dot(gn, un))
        ws = float(weights[s])
        # dP = (ws/2) Re[ α (gn·dx) Hz* ] + (ws/2) Re[ α dn (ℓ·dx)* ]
        # => q_k += (ws/2) conj(α) Hz gn_k  +  (ws/2) α dn ℓ_k
        for k in range(3):
            nd = int(nodes[k])
            q[nd] += 0.5 * ws * np.conj(alpha) * hz * gn[k]
            q[nd] += 0.5 * ws * alpha * dn * ell[k]
    return q


def metrics_from_power(P: np.ndarray, P_ref: float) -> dict:
    t = P / P_ref
    rows = []
    acc = []
    for j in range(6):
        r = port_roles(j)
        des = float(t[r["desired"], j])
        rev = float(t[r["reverse"], j])
        thr = float(t[r["through"], j])
        leak = float(t[r["leak_a"], j] + t[r["leak_b"], j])
        a = float(-t[r["source"], j])
        acc.append(a)
        rows.append({"desired": des, "reverse": rev, "through": thr, "leak": leak, "accepted": a})
    des_avg = float(np.mean([r["desired"] for r in rows]))
    rev_avg = float(np.mean([r["reverse"] for r in rows]))
    thr_avg = float(np.mean([r["through"] for r in rows]))
    leak_avg = float(np.mean([r["leak"] for r in rows]))
    acc_avg = float(np.mean(acc))
    J = (
        W_DES * des_avg
        - W_REV * rev_avg
        - W_THR * thr_avg
        - W_LEAK * leak_avg
        + W_ACC * acc_avg
    )
    eps = 1e-30
    return {
        "per_source": rows,
        "desired_avg": des_avg,
        "reverse_avg": rev_avg,
        "through_avg": thr_avg,
        "leak_avg": leak_avg,
        "accepted_avg": acc_avg,
        "desired_min": float(min(r["desired"] for r in rows)),
        "isolation_avg_dB": float(
            np.mean([10 * np.log10(max(r["desired"], eps) / max(r["reverse"], eps)) for r in rows])
        ),
        "insertion_avg_dB": float(np.mean([10 * np.log10(max(r["desired"], eps)) for r in rows])),
        "J": float(J),
        "P_ref": float(P_ref),
        "weights": {"des": W_DES, "rev": W_REV, "thr": W_THR, "leak": W_LEAK, "acc": W_ACC},
    }


def assemble_and_factor(dev: DeviceData, s_vec, f_hz: float, b_tesla: float):
    k0 = 2.0 * np.pi * (f_hz * sc.a / 299792458.0)
    rho, f, gamma, fc = rho_from_scales(
        len(dev.tris),
        dev.rod_masks,
        s_vec,
        dev.quartz,
        dev.walls,
        f_hz,
        b_tesla,
        sc.a,
        sc.fp_Hz,
        sc.gamma_Hz,
    )
    t0 = time.perf_counter()
    A = assemble_anisotropic(dev.points, dev.tris, *rho, k0, np.zeros(len(dev.points), dtype=bool))
    lu = splu(A.tocsc())
    t_fac = time.perf_counter() - t0
    fp_ref_ord = sc.fp_Hz * sc.a / 299792458.0
    return rho, A, lu, k0, f, gamma, fc, fp_ref_ord, t_fac


def solve_six(dev: DeviceData, s_vec, f_hz: float, b_tesla: float, P_ref: Optional[float] = None) -> SolveBundle:
    s_vec = np.asarray(s_vec, float)
    rho, A, lu, k0, f, gamma, fc, fp_ref_ord, t_fac = assemble_and_factor(dev, s_vec, f_hz, b_tesla)
    fields = []
    t_sol = 0.0
    for b in dev.rhs:
        t1 = time.perf_counter()
        fields.append(lu.solve(b))
        t_sol += time.perf_counter() - t1
    P = power_matrix(dev, fields, rho, k0)
    if P_ref is None:
        P_ref = float(np.mean([-P[j, j] for j in range(6)]))
        P_ref = max(P_ref, 1e-30)
    metrics = metrics_from_power(P, P_ref)
    return SolveBundle(
        rho=rho,
        A=A,
        lu=lu,
        fields=fields,
        k0=k0,
        f_ord=f,
        gamma=gamma,
        fc=fc,
        fp_ref_ord=fp_ref_ord,
        s_vec=s_vec.copy(),
        power=P,
        P_ref=P_ref,
        metrics=metrics,
        J=metrics["J"],
        factor_s=t_fac,
        solve_s=t_sol,
    )


def objective_q_for_source(dev: DeviceData, x, rho, k0, j: int, P_ref: float) -> np.ndarray:
    """Adjoint RHS q_j for source j contribution to J (before 1/6 is in metrics; include 1/6 here)."""
    r = port_roles(j)
    scale = 1.0 / (6.0 * P_ref)
    q = np.zeros(len(dev.points), np.complex128)
    q += scale * W_DES * power_adjoint_q(dev.sampler, dev.tris, x, rho, k0, dev.lines[r["desired"]])
    q -= scale * W_REV * power_adjoint_q(dev.sampler, dev.tris, x, rho, k0, dev.lines[r["reverse"]])
    q -= scale * W_THR * power_adjoint_q(dev.sampler, dev.tris, x, rho, k0, dev.lines[r["through"]])
    q -= scale * W_LEAK * power_adjoint_q(dev.sampler, dev.tris, x, rho, k0, dev.lines[r["leak_a"]])
    q -= scale * W_LEAK * power_adjoint_q(dev.sampler, dev.tris, x, rho, k0, dev.lines[r["leak_b"]])
    # a_j = -P_jj/P_ref => +w_acc a contributes -w_acc/P_ref * dP_jj
    q -= scale * W_ACC * power_adjoint_q(dev.sampler, dev.tris, x, rho, k0, dev.lines[r["source"]])
    return q


def gradient_s(dev: DeviceData, bundle: SolveBundle, rod_indices: Optional[Sequence[int]] = None) -> np.ndarray:
    """∂J/∂s_k for selected rods (default all). Reuses bundle.lu."""
    n = len(dev.rod_masks)
    idx = list(range(n)) if rod_indices is None else list(rod_indices)
    # adjoint fields
    lams = []
    for j, x in enumerate(bundle.fields):
        q = objective_q_for_source(dev, x, bundle.rho, bundle.k0, j, bundle.P_ref)
        try:
            lam = bundle.lu.solve(q, trans="H")
        except TypeError:
            lam = splu(bundle.A.conj().T.tocsc()).solve(q)
        lams.append(lam)
    g = np.zeros(n, float)
    for k in idx:
        dA = assemble_dA_rod(
            dev.points,
            dev.tris,
            dev.rod_masks[k],
            bundle.f_ord,
            bundle.fp_ref_ord,
            bundle.gamma,
            bundle.fc,
            bundle.k0,
            float(bundle.s_vec[k]),
        )
        acc = 0.0
        for j in range(6):
            acc += -np.real(np.vdot(lams[j], dA @ bundle.fields[j]))
        g[k] = acc
    return g


def forward_power_deriv(dev, x, dx, rho, k0, line) -> float:
    """dP from product rule on vacuum power (for FD cross-check via dx)."""
    xy = line["monitor_xy"]
    n_hat = np.asarray(line["outward"], float)
    n_hat /= np.linalg.norm(n_hat)
    weights = np.asarray(line["weights"], float)
    hz, ex, ey = dev.sampler.fields(x, *rho, k0, xy)
    dhz, dex, dey = dev.sampler.fields(dx, *rho, k0, xy)
    dsx = 0.5 * np.real(dey * np.conj(hz) + ey * np.conj(dhz))
    dsy = -0.5 * np.real(dex * np.conj(hz) + ex * np.conj(dhz))
    return float(np.sum((n_hat[0] * dsx + n_hat[1] * dsy) * weights))


def gradient_s_forward(dev: DeviceData, bundle: SolveBundle, rod_indices: Sequence[int]) -> np.ndarray:
    """Forward-mode chain rule for selected rods (expensive validation path)."""
    g = np.zeros(len(dev.rod_masks), float)
    P_ref = bundle.P_ref
    for k in rod_indices:
        dA = assemble_dA_rod(
            dev.points,
            dev.tris,
            dev.rod_masks[k],
            bundle.f_ord,
            bundle.fp_ref_ord,
            bundle.gamma,
            bundle.fc,
            bundle.k0,
            float(bundle.s_vec[k]),
        )
        acc = 0.0
        for j, x in enumerate(bundle.fields):
            dx = bundle.lu.solve(-(dA @ x))
            r = port_roles(j)
            dP = {}
            for name, ip in r.items():
                dP[name] = forward_power_deriv(dev, x, dx, bundle.rho, bundle.k0, dev.lines[ip])
            dJ_j = (
                W_DES * dP["desired"]
                - W_REV * dP["reverse"]
                - W_THR * dP["through"]
                - W_LEAK * (dP["leak_a"] + dP["leak_b"])
                + W_ACC * (-dP["source"])
            ) / (6.0 * P_ref)
            acc += dJ_j
        g[k] = acc
    return g
