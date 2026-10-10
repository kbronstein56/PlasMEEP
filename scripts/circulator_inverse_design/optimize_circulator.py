#!/usr/bin/env python3
"""Sixfold-tied single-frequency (and optional broadband) circulator optimization."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from circulator_core import (  # noqa: E402
    B_PLUS,
    F_HZ,
    S_HI,
    S_LO,
    DeviceData,
    assemble_and_factor,
    load_device,
    metrics_from_power,
    objective_q_for_source,
    power_matrix,
    solve_six,
)
from diff_fem_core import assemble_dA_rod  # noqa: E402
from sixfold_orbits import build_orbits, expand_tied  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"


def orbit_masks(dev: DeviceData, orbits):
    masks = []
    for o in orbits:
        m = np.zeros(len(dev.tris), dtype=bool)
        for i in o["rods"]:
            m |= dev.rod_masks[i]
        masks.append(m)
    return masks


def gradient_q(dev, bundle, orbits, o_masks) -> np.ndarray:
    """∂J/∂q_m using one dA per orbit (all rods in orbit share s)."""
    lams = []
    for j, x in enumerate(bundle.fields):
        q = objective_q_for_source(dev, x, bundle.rho, bundle.k0, j, bundle.P_ref)
        try:
            lam = bundle.lu.solve(q, trans="H")
        except TypeError:
            lam = splu(bundle.A.conj().T.tocsc()).solve(q)
        lams.append(lam)
    g = np.zeros(len(orbits), float)
    # representative s for orbit (all equal under tying)
    for m, o in enumerate(orbits):
        s_m = float(bundle.s_vec[o["rods"][0]])
        dA = assemble_dA_rod(
            dev.points,
            dev.tris,
            o_masks[m],
            bundle.f_ord,
            bundle.fp_ref_ord,
            bundle.gamma,
            bundle.fc,
            bundle.k0,
            s_m,
        )
        acc = 0.0
        for j in range(6):
            acc += -np.real(np.vdot(lams[j], dA @ bundle.fields[j]))
        g[m] = acc
    return g


class OptState:
    def __init__(self, dev, orbits, o_masks, f_hz, b_tesla, P_ref, freqs=None):
        self.dev = dev
        self.orbits = orbits
        self.o_masks = o_masks
        self.f_hz = f_hz
        self.b_tesla = b_tesla
        self.P_ref = P_ref
        self.freqs = freqs  # None => single frequency
        self.history = []
        self.best = {"J": -1e300, "q": None, "metrics": None, "s": None}
        self.n_fev = 0
        self.n_gev = 0

    def s_from_q(self, q):
        return expand_tied(np.asarray(q, float), self.orbits)

    def eval_single(self, q):
        s = self.s_from_q(q)
        bun = solve_six(self.dev, s, self.f_hz, self.b_tesla, P_ref=self.P_ref)
        gq = gradient_q(self.dev, bun, self.orbits, self.o_masks)
        return bun, gq

    def eval_broadband(self, q, freqs):
        s = self.s_from_q(q)
        Js = []
        g_acc = np.zeros(len(self.orbits), float)
        metrics_band = []
        for f in freqs:
            bun = solve_six(self.dev, s, f, self.b_tesla, P_ref=self.P_ref)
            gq = gradient_q(self.dev, bun, self.orbits, self.o_masks)
            Js.append(bun.J)
            g_acc += gq
            metrics_band.append({"f_Hz": f, "J": bun.J, "desired_avg": bun.metrics["desired_avg"],
                                 "isolation_avg_dB": bun.metrics["isolation_avg_dB"]})
        # robust: mean + 0.5 * min (differentiable soft via exact min for now)
        Jmean = float(np.mean(Js))
        Jmin = float(np.min(Js))
        # For gradient of min: use the active frequency only (subgradient)
        imin = int(np.argmin(Js))
        # Recompute gradient properly: d(mean)/dq + 0.5 d(min)/dq
        # We need per-freq gradients stored — redo cheaply
        g_mean = g_acc / len(freqs)
        # get g at imin
        bun_min = solve_six(self.dev, s, freqs[imin], self.b_tesla, P_ref=self.P_ref)
        g_min = gradient_q(self.dev, bun_min, self.orbits, self.o_masks)
        J = Jmean + 0.5 * Jmin
        g = g_mean + 0.5 * g_min
        bun_min.J = J  # stash
        bun_min.metrics = {
            "J": J,
            "J_mean": Jmean,
            "J_min": Jmin,
            "band": metrics_band,
            "desired_avg": float(np.mean([m["desired_avg"] for m in metrics_band])),
            "isolation_avg_dB": float(np.mean([m["isolation_avg_dB"] for m in metrics_band])),
            "per_source": bun_min.metrics["per_source"],
            "reverse_avg": bun_min.metrics["reverse_avg"],
            "through_avg": bun_min.metrics["through_avg"],
            "leak_avg": bun_min.metrics["leak_avg"],
            "accepted_avg": bun_min.metrics["accepted_avg"],
            "desired_min": bun_min.metrics["desired_min"],
            "insertion_avg_dB": bun_min.metrics["insertion_avg_dB"],
            "P_ref": self.P_ref,
            "weights": bun_min.metrics["weights"],
        }
        return bun_min, g


def run_opt(grade="M", b_tesla=B_PLUS, f_hz=F_HZ, maxiter=40, q0=None, tag="single", freqs=None):
    OUT.mkdir(parents=True, exist_ok=True)
    dev = load_device(grade)
    orbits, _ = build_orbits(dev.centers)
    o_masks = orbit_masks(dev, orbits)
    nq = len(orbits)
    # freeze P_ref at uniform s=1, production frequency (even for broadband, use center f)
    base = solve_six(dev, np.ones(len(dev.rod_masks)), f_hz if freqs is None else F_HZ, b_tesla, P_ref=None)
    P_ref = base.P_ref
    print(f"P_ref frozen={P_ref:.6e} J0={base.J:.6e} orbits={nq}", flush=True)

    if q0 is None:
        q0 = np.ones(nq)
    else:
        q0 = np.asarray(q0, float)

    st = OptState(dev, orbits, o_masks, f_hz, b_tesla, P_ref, freqs=freqs)
    bounds = [(S_LO, S_HI)] * nq

    def fun(q):
        st.n_fev += 1
        if freqs is None:
            bun, gq = st.eval_single(q)
        else:
            bun, gq = st.eval_broadband(q, freqs)
        st._last_g = gq
        st._last_bun = bun
        rec = {
            "fev": st.n_fev,
            "J": bun.J,
            "desired_avg": bun.metrics["desired_avg"],
            "reverse_avg": bun.metrics["reverse_avg"],
            "through_avg": bun.metrics["through_avg"],
            "leak_avg": bun.metrics["leak_avg"],
            "accepted_avg": bun.metrics["accepted_avg"],
            "isolation_avg_dB": bun.metrics.get("isolation_avg_dB"),
            "q_min": float(np.min(q)),
            "q_max": float(np.max(q)),
            "grad_norm": float(np.linalg.norm(gq)),
        }
        st.history.append(rec)
        if bun.J > st.best["J"]:
            st.best = {
                "J": bun.J,
                "q": np.asarray(q, float).tolist(),
                "s": st.s_from_q(q).tolist(),
                "metrics": bun.metrics,
                "power": bun.power.tolist(),
            }
        print(
            f"[{tag}] fev={st.n_fev} J={bun.J:.6e} des={bun.metrics['desired_avg']:.4f} "
            f"rev={bun.metrics['reverse_avg']:.4f} thr={bun.metrics['through_avg']:.4f} "
            f"|g|={rec['grad_norm']:.3e}",
            flush=True,
        )
        return -bun.J  # maximize J

    def jac(q):
        st.n_gev += 1
        # scipy may call jac after fun with same q; reuse
        if hasattr(st, "_last_g") and hasattr(st, "_last_q") and np.allclose(q, st._last_q):
            return -st._last_g
        if freqs is None:
            bun, gq = st.eval_single(q)
        else:
            bun, gq = st.eval_broadband(q, freqs)
        st._last_g = gq
        st._last_bun = bun
        st._last_q = np.asarray(q, float).copy()
        return -gq

    def fun_and_cache(q):
        st._last_q = np.asarray(q, float).copy()
        return fun(q)

    t0 = time.perf_counter()
    res = minimize(
        fun_and_cache,
        q0,
        method="L-BFGS-B",
        jac=jac,
        bounds=bounds,
        options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-8, "disp": True},
    )
    wall = time.perf_counter() - t0
    out = {
        "tag": tag,
        "grade": grade,
        "B_T": b_tesla,
        "f_Hz": f_hz,
        "freqs": freqs,
        "P_ref": P_ref,
        "J0": base.J,
        "metrics0": base.metrics,
        "power0": base.power.tolist(),
        "success": bool(res.success),
        "message": str(res.message),
        "nit": int(res.nit),
        "nfev": int(res.nfev),
        "njev": int(getattr(res, "njev", st.n_gev)),
        "wall_s": wall,
        "history": st.history,
        "best": st.best,
        "final_q": np.asarray(res.x, float).tolist(),
        "bounds": [S_LO, S_HI],
        "orbits": [{"rods": o["rods"], "radius_a": o["radius_a"]} for o in orbits],
    }
    # evaluate final/best on same mesh for complete metrics
    s_best = np.asarray(st.best["s"], float)
    bun_best = solve_six(dev, s_best, f_hz if freqs is None else F_HZ, b_tesla, P_ref=P_ref)
    out["best_reeval"] = {"J": bun_best.J, "metrics": bun_best.metrics, "power": bun_best.power.tolist()}
    path = OUT / f"opt_{tag}_{grade}_B{b_tesla:+.4f}.json".replace("+", "p").replace("-", "m")
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE", path.name, "best J", st.best["J"], flush=True)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "single"
    grade = sys.argv[2] if len(sys.argv) > 2 else "M"
    if mode == "single":
        run_opt(grade=grade, tag="single_uniform", maxiter=35)
    elif mode == "pert":
        # small symmetry-preserving perturbation start
        rng = np.random.default_rng(0)
        # need n orbits
        from full91_geometry import bulb_centers_device
        import sixport_common as sc
        orbits, _ = build_orbits(bulb_centers_device(sc))
        q0 = np.clip(1.0 + 0.05 * rng.standard_normal(len(orbits)), S_LO, S_HI)
        run_opt(grade=grade, tag="single_pert", maxiter=25, q0=q0)
    elif mode == "broadband":
        # warm start from best single if present
        cands = sorted(OUT.glob("opt_single_*_Bp0.0500.json"))
        q0 = None
        if cands:
            best = max((json.loads(p.read_text()) for p in cands), key=lambda d: d["best"]["J"])
            q0 = best["best"]["q"]
            print("warm start from J", best["best"]["J"], flush=True)
        freqs = [3.75e9, 3.80e9, 3.85e9, 3.90e9, 3.95e9]
        run_opt(grade=grade, tag="broadband", maxiter=20, q0=q0, freqs=freqs)
    elif mode == "free91":
        # release tying — separate script path using per-rod LBFGS would be heavy;
        # use tied warm start then a few free steps via finite orbits expansion
        print("use refine_free91.py", flush=True)
    else:
        raise SystemExit(mode)


if __name__ == "__main__":
    main()
