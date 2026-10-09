#!/usr/bin/env python3
"""Tiny in-silico inverse-design proof: one density scale, maximize |Hz_probe|^2.

Not a circulator optimization. Only shows grad → optimizer → objective moves
in the predicted direction on the validated one-cylinder FEM.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation" / "fem_gyrotropic_validation"),
    str(ROOT / "scripts" / "validation" / "fem_meep_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
from adjoint_check import probe_row  # noqa: E402
from fem_validated_solver import ElementSampler, assemble_anisotropic  # noqa: E402
from jacobian_one_cylinder import build_mesh, solve_state  # noqa: E402
from plasma_sensitivity import element_drho  # noqa: E402
import sixport_common as sc  # noqa: E402

OUT = ROOT / "outputs" / "validation" / "fem_full_device_validation"
C = 299792458.0


def main():
    pts, tris, inside, center = build_mesh()
    sc.fs_a = sc.fs_Hz * sc.a / C
    k0 = 2 * np.pi * sc.fs_a
    sampler = ElementSampler(pts, tris)
    src = center + np.array([-1.5, 0.2])
    bvec = np.zeros(len(pts), np.complex128)
    t_src = int(sampler.locate(src[None, :])[0])
    w = sampler._bary(t_src, src)
    for k in range(3):
        bvec[int(tris[t_src, k])] += complex(w[k])
    probe = center + np.array([1.3, 0.4])
    L = probe_row(sampler, tris, probe, len(pts))
    b_tesla = 0.05
    s = 1.0
    s_lo, s_hi = 0.2, 2.0
    history = []
    lr = 0.15
    for it in range(12):
        x, A, ords, rho = solve_state(pts, tris, inside, s, b_tesla, k0, bvec)
        f, fp, gamma, fc = ords
        y = complex((L @ x)[0])
        J = abs(y) ** 2
        fp_ref = sc.fp_Hz * sc.a / C
        dr = element_drho(inside, f, fp_ref, gamma, fc, s)
        dA = assemble_anisotropic(
            pts, tris, *dr, k0, np.zeros(len(pts), dtype=bool), mass_scale=0.0, pin_empty=False
        )
        lu = splu(A.tocsc())
        rhs = np.asarray(L.conj().T @ np.array([y])).ravel()
        try:
            lam = lu.solve(rhs, trans="H")
        except TypeError:
            lam = splu(A.conj().T.tocsc()).solve(rhs)
        g = -2.0 * np.real(np.vdot(lam, dA @ x))
        history.append({"iter": it, "s": float(s), "J": float(J), "grad": float(g)})
        print(f"it={it} s={s:.6f} J={J:.6e} g={g:.6e}", flush=True)
        # ascent on J
        s = float(np.clip(s + lr * np.sign(g) * min(abs(g) * 5.0, 0.25), s_lo, s_hi))
    # predicted direction: J should not decrease on the first successful step with same sign
    improved = history[-1]["J"] > history[0]["J"]
    out = {
        "objective": "|Hz_probe|^2 ascent",
        "B_T": b_tesla,
        "history": history,
        "J0": history[0]["J"],
        "J_final": history[-1]["J"],
        "improved": bool(improved),
        "pass": bool(improved),
        "note": "Toy proof only. Not a circulator design result.",
    }
    (OUT / "opt_proof_one_cylinder.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE opt_proof_one_cylinder.json", out["pass"], flush=True)


if __name__ == "__main__":
    main()
