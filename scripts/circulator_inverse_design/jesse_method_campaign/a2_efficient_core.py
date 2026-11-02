#!/usr/bin/env python3
"""Shared Stage-2 A2 efficient eval/grad with optional detailed timing.

Used by production Stage-2 opt, performance profiling, and field exports.
Does not change A2 mathematics.
"""
from __future__ import annotations

import os
import time
from typing import Any

import numpy as np
from scipy.sparse.linalg import splu

from adjoint_direct_dtn import build_system
from broadband_objective import softmin_bb
from diff_fem_core import assemble_dA_rod
from direct_dtn_s import build_Q_columns
from efficient_adjoint import S_from_U, _pack_C, _solve_multi
from multimode_opt_core import B_PLUS, expand_tied, orbits_and_masks, te0_metrics_from_S
from objective_A2 import EPS, objective_A2_from_S, scalar_adjoint_rhs_A2


def rss_gib() -> float:
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / (1024.0 * 1024.0)
    except Exception:
        return float("nan")
    return float("nan")


def eval_and_grad_A2(
    grade,
    q,
    orbits,
    f_hz,
    b_tesla=B_PLUS,
    eps=EPS,
    profile: bool = False,
    return_fields: bool = False,
):
    """Single-frequency A2 objective + C6 gradient (1 LU + 6+6 RHS)."""
    tim: dict[str, Any] = {}
    rss0 = rss_gib() if profile else None
    t_all = time.perf_counter()

    t0 = time.perf_counter()
    s = expand_tied(np.asarray(q, float), orbits)
    if profile:
        tim["expand_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    st = build_system(grade, s, f_hz, b_tesla, "projected")
    if profile:
        tim["assemble_build_system"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    _, _, o_masks = orbits_and_masks(rod_masks=st["rod_masks"])
    points, A = st["points"], st["A"]
    cols = build_Q_columns(points, st["ports"], st["modes"], n_prop=2)
    C, betas = _pack_C(points, cols)
    te0_cols = [2 * j for j in range(6)]
    if profile:
        tim["pack_ports"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    lu = splu(A.tocsc())
    if profile:
        tim["factor_lu"] = time.perf_counter() - t0
        tim["rss_GiB_after_factor"] = rss_gib()

    t0 = time.perf_counter()
    U = _solve_multi(lu, C[:, te0_cols], trans="N")
    if profile:
        tim["forward_te0_rhs"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    S = S_from_U(C, U, betas, te0_cols)
    J = objective_A2_from_S(S, eps=eps)
    met = te0_metrics_from_S(S)
    if profile:
        tim["S_and_objective"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    R = np.column_stack(
        [scalar_adjoint_rhs_A2(C, betas, S, j, eps=eps) for j in range(6)]
    )
    Lam = _solve_multi(lu, R, trans="T")
    if profile:
        tim["adjoint_te0_rhs"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    g = np.zeros(len(orbits), float)
    for m, o in enumerate(orbits):
        s_m = float(s[o["rods"][0]])
        dA = assemble_dA_rod(
            points,
            st["tris"],
            o_masks[m],
            st["f_ord"],
            st["fp_ref_ord"],
            st["gamma"],
            st["fc"],
            st["k0"],
            s=s_m,
        )
        acc = 0.0
        for j_loc in range(6):
            dAu = dA @ U[:, j_loc]
            acc += 2.0 * np.real(np.dot(Lam[:, j_loc], dAu))
        g[m] = float(acc / 6.0)
    if profile:
        tim["dA_contraction"] = time.perf_counter() - t0
        tim["total"] = time.perf_counter() - t_all
        tim["rss_GiB_start"] = rss0
        tim["rss_GiB_end"] = rss_gib()
        tim["n_lu"] = 1
        tim["n_forward_rhs"] = 6
        tim["n_adjoint_rhs"] = 6
        tim["n_dof"] = int(A.shape[0])
        tim["omp_num_threads"] = os.environ.get("OMP_NUM_THREADS")
        tim["openblas_num_threads"] = os.environ.get("OPENBLAS_NUM_THREADS")
        tim["mkl_num_threads"] = os.environ.get("MKL_NUM_THREADS")

    bun = {"J": J, "S": S, "metrics": met, "timings": tim if profile else None}
    if return_fields:
        bun["points"] = np.asarray(points, float)
        bun["tris"] = np.asarray(st["tris"], int)
        bun["U_te0"] = U  # columns = TE0 drives for ports 0..5
        bun["ports"] = st["ports"]
        bun["betas"] = betas
        bun["k0"] = st["k0"]
        bun["b_tesla"] = float(b_tesla)
        bun["f_hz"] = float(f_hz)
    return bun, g, s


def _top_weights(freqs_hz, w, Js, k=5):
    w = np.asarray(w, float)
    order = np.argsort(-w)[:k]
    return [
        {
            "f_GHz": float(freqs_hz[i] / 1e9),
            "weight": float(w[i]),
            "J_A2": float(Js[i]),
        }
        for i in order
    ]


def _pack_bb(Js, gs, mets, freqs_hz, tau, alpha):
    Js = np.asarray(Js, float)
    J_bb, w, J_soft, J_mean, p = softmin_bb(Js, tau=tau, alpha=alpha)
    g_bb = sum(wi * gi for wi, gi in zip(w, gs))
    worst_i = int(np.argmin(Js))
    idx_41 = int(np.argmin(np.abs(np.asarray(freqs_hz) / 1e9 - 4.1)))
    return {
        "J_BB": float(J_bb),
        "J_softmin": float(J_soft),
        "J_mean": float(J_mean),
        "J_per_freq": Js.tolist(),
        "weights": np.asarray(w, float).tolist(),
        "softmax_p": np.asarray(p, float).tolist(),
        "freqs_Hz": [float(f) for f in freqs_hz],
        "freqs_GHz": [float(f) / 1e9 for f in freqs_hz],
        "metrics_per_freq": mets,
        "grad_q": np.asarray(g_bb, float),
        "grad_norm": float(np.linalg.norm(g_bb)),
        "tau": tau,
        "alpha": alpha,
        "worst_J": float(Js[worst_i]),
        "worst_f_Hz": float(freqs_hz[worst_i]),
        "worst_f_GHz": float(freqs_hz[worst_i] / 1e9),
        "worst_metrics": mets[worst_i],
        "metrics_4p1": mets[idx_41],
        "J_4p1": float(Js[idx_41]),
        "max_weight": float(np.max(w)),
        "N_eff": float(1.0 / np.sum(np.asarray(w, float) ** 2)),
        "top5_weights": _top_weights(freqs_hz, w, Js, 5),
    }


def eval_bb_A2_serial(
    grade,
    q,
    orbits,
    freqs_hz,
    tau=0.05,
    alpha=0.25,
    profile: bool = False,
):
    """Serial frequency loop (production default)."""
    t0 = time.perf_counter()
    Js, gs, mets, per = [], [], [], []
    for f in freqs_hz:
        bun, g, _ = eval_and_grad_A2(
            grade, q, orbits, f_hz=float(f), profile=profile
        )
        Js.append(bun["J"])
        gs.append(g)
        mets.append(bun["metrics"])
        if profile and bun.get("timings"):
            per.append({"f_Hz": float(f), **bun["timings"]})
    out = _pack_bb(Js, gs, mets, freqs_hz, tau, alpha)
    out["wall_s"] = time.perf_counter() - t0
    out["n_workers"] = 1
    out["mode"] = "serial"
    if profile:
        out["per_freq_timings"] = per
    return out


def _ensure_sys_path():
    """Spawned workers need the same script path layout as production runners."""
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    paths = [
        str(Path(__file__).resolve().parent),
        str(root / "scripts" / "circulator_inverse_design" / "broadband_method"),
        str(root / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
        str(root / "scripts" / "circulator_inverse_design" / "multimode_ports"),
        str(root / "scripts" / "circulator_inverse_design"),
        str(root / "scripts" / "validation" / "fem_full_device_validation"),
        str(root / "scripts" / "validation"),
        str(root / "scripts"),
    ]
    for p in reversed(paths):
        if p not in sys.path:
            sys.path.insert(0, p)


def _init_parallel_worker(omp_threads: int):
    """Runs in each spawned worker before tasks; env must be set before BLAS use."""
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(int(omp_threads))
    _ensure_sys_path()


def _worker_payload(payload: dict):
    """Process-pool worker: one frequency A2 eval+grad."""
    _ensure_sys_path()
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(payload.get("omp_threads", 1))
    from multimode_opt_core import orbits_and_masks as _orb

    orbits, _, _ = _orb()
    bun, g, _ = eval_and_grad_A2(
        payload["grade"],
        np.asarray(payload["q"], float),
        orbits,
        f_hz=float(payload["f_Hz"]),
        profile=bool(payload.get("profile", False)),
    )
    return {
        "f_Hz": float(payload["f_Hz"]),
        "J": float(bun["J"]),
        "g": np.asarray(g, float),
        "metrics": bun["metrics"],
        "timings": bun.get("timings"),
    }


def eval_bb_A2_process_parallel(
    grade,
    q,
    orbits,
    freqs_hz,
    tau=0.05,
    alpha=0.25,
    n_workers: int = 2,
    omp_threads_per_worker: int = 1,
    profile: bool = False,
):
    """Process-parallel over frequencies. Default-off production option.

    Each frequency remains scientifically independent (own assemble + LU + 6+6).
    Numerical result must match serial to floating-point aggregation order.
    """
    # spawn so each worker gets a clean BLAS/OpenMP init after OMP=1 is set
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed

    n_workers = max(1, int(n_workers))
    if n_workers == 1:
        return eval_bb_A2_serial(
            grade, q, orbits, freqs_hz, tau=tau, alpha=alpha, profile=profile
        )

    t0 = time.perf_counter()
    payloads = [
        {
            "grade": grade,
            "q": np.asarray(q, float).tolist(),
            "f_Hz": float(f),
            "omp_threads": int(omp_threads_per_worker),
            "profile": profile,
        }
        for f in freqs_hz
    ]
    by_f = {}
    ctx = mp.get_context("spawn")
    # Inherit OMP=omp_threads_per_worker into spawned children at process start
    # (before they import SciPy/OpenBLAS), then restore parent env afterward.
    thread_keys = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    saved_env = {k: os.environ.get(k) for k in thread_keys}
    for k in thread_keys:
        os.environ[k] = str(int(omp_threads_per_worker))
    try:
        with ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=ctx,
            initializer=_init_parallel_worker,
            initargs=(int(omp_threads_per_worker),),
        ) as ex:
            futs = {ex.submit(_worker_payload, p): p["f_Hz"] for p in payloads}
            for fut in as_completed(futs):
                r = fut.result()
                by_f[r["f_Hz"]] = r
    finally:
        for k, v in saved_env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    # Preserve frequency order (deterministic)
    Js, gs, mets, per = [], [], [], []
    for f in freqs_hz:
        r = by_f[float(f)]
        Js.append(r["J"])
        gs.append(r["g"])
        mets.append(r["metrics"])
        if profile and r.get("timings"):
            per.append({"f_Hz": float(f), **r["timings"]})
    out = _pack_bb(Js, gs, mets, freqs_hz, tau, alpha)
    out["wall_s"] = time.perf_counter() - t0
    out["n_workers"] = n_workers
    out["omp_threads_per_worker"] = omp_threads_per_worker
    out["mode"] = "process_parallel"
    if profile:
        out["per_freq_timings"] = per
    return out


def eval_bb_A2(
    grade,
    q,
    orbits,
    freqs_hz,
    tau=0.05,
    alpha=0.25,
    n_workers: int = 1,
    omp_threads_per_worker: int = 1,
    profile: bool = False,
):
    """Production entry: serial unless n_workers>1."""
    if int(n_workers) <= 1:
        return eval_bb_A2_serial(
            grade, q, orbits, freqs_hz, tau=tau, alpha=alpha, profile=profile
        )
    return eval_bb_A2_process_parallel(
        grade,
        q,
        orbits,
        freqs_hz,
        tau=tau,
        alpha=alpha,
        n_workers=n_workers,
        omp_threads_per_worker=omp_threads_per_worker,
        profile=profile,
    )
