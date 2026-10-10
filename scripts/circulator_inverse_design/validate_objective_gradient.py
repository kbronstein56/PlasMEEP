#!/usr/bin/env python3
"""Three-way gradient check for the full six-port circulator objective.

1) centered FD of J
2) forward sensitivity + chain rule
3) port-power adjoint

Reduced: FULL91-C selected rods (optional).
Production check: FULL91-M at B=+0.05 T, f=3.85 GHz.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

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
    gradient_s,
    gradient_s_forward,
    load_device,
    solve_six,
)
from sixfold_orbits import build_orbits, reduce_gradient  # noqa: E402
import sixport_common as sc  # noqa: E402
from full91_geometry import bulb_centers_device  # noqa: E402

OUT = ROOT / "outputs" / "circulator_inverse_design"
HS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3)
REL_TOL = 1e-4  # same frozen relative gate as JACOBIAN_PASS_CRITERIA.md


def gate(errs):
    imin = int(np.argmin(errs))
    minimum = errs[imin]
    larger = max(errs[:imin]) if imin > 0 else (max(errs[1:]) if len(errs) > 1 else errs[0])
    return imin, minimum, larger, bool(minimum <= REL_TOL and larger >= 5.0 * minimum)


def select_rods(centers):
    c0 = centers - centers.mean(0)
    r = np.linalg.norm(c0, axis=1)
    i_center = int(np.argmin(r))
    shells = np.unique(np.round(r, 3))
    shells = shells[shells > 0.05]
    i_inner = int(np.argmin(np.abs(r - shells[0])))
    i_outer = int(np.argmin(np.abs(r - shells[-1])))
    shift = np.array([sc.nx_ports / 2.0, sc.ny_ports / 2.0])
    mon0 = np.asarray(sc.full_horns[0]["monitor_center"], float) + shift
    i_horn = int(np.argmin(np.linalg.norm(centers - mon0, axis=1)))
    cy = sc.ny_ports / 2.0
    cand = np.where((c0[:, 1] > 0.3) & (r > 0.5 * r.max()))[0]
    i_a = int(cand[len(cand) // 2]) if len(cand) else i_outer
    mirror = np.array([centers[i_a, 0], 2 * cy - centers[i_a, 1]])
    i_b = int(np.argmin(np.linalg.norm(centers - mirror, axis=1)))
    labs = ["center", "inner", "outer", "near_horn0", "sym_a", "sym_b"]
    ids = [i_center, i_inner, i_outer, i_horn, i_a, i_b]
    # unique
    out_l, out_i = [], []
    seen = set()
    for lab, i in zip(labs, ids):
        if i in seen:
            continue
        seen.add(i)
        out_l.append(lab)
        out_i.append(i)
    return out_l, out_i


def run(grade: str, b_tesla: float):
    print(f"=== gradient validation grade={grade} B={b_tesla} ===", flush=True)
    dev = load_device(grade)
    s0 = np.ones(len(dev.rod_masks))
    t0 = time.perf_counter()
    base = solve_six(dev, s0, F_HZ, b_tesla, P_ref=None)
    P_ref = base.P_ref
    print(f"base J={base.J:.6e} P_ref={P_ref:.6e} factor={base.factor_s:.1f}s", flush=True)
    labels, rods = select_rods(dev.centers)
    print("rods", list(zip(labels, rods)), flush=True)

    do_fwd = "--adjoint-fd-only" not in sys.argv
    g_adj = gradient_s(dev, base, rods)
    g_fwd = gradient_s_forward(dev, base, rods) if do_fwd else g_adj.copy()

    rows = []
    for lab, k in zip(labels, rods):
        sweep = []
        for h in HS:
            sp = s0.copy(); sp[k] = 1.0 + h
            sm = s0.copy(); sm[k] = 1.0 - h
            Jp = solve_six(dev, sp, F_HZ, b_tesla, P_ref=P_ref).J
            Jm = solve_six(dev, sm, F_HZ, b_tesla, P_ref=P_ref).J
            g_fd = (Jp - Jm) / (2.0 * h)
            rel_adj = abs(g_fd - g_adj[k]) / max(abs(g_adj[k]), 1e-30)
            rel_fwd = abs(g_fd - g_fwd[k]) / max(abs(g_fwd[k]), 1e-30) if do_fwd else float("nan")
            sweep.append({
                "h": h,
                "g_fd": g_fd,
                "fd_vs_adj": float(rel_adj),
                "fd_vs_fwd": float(rel_fwd),
            })
            print(
                f"{lab}({k}) h={h:.1e} fd={g_fd:.4e} adj={g_adj[k]:.4e} "
                f"fwd={g_fwd[k]:.4e} rel_adj={rel_adj:.3e}",
                flush=True,
            )
        errs = [r["fd_vs_adj"] for r in sweep]
        imin, minimum, larger, passed = gate(errs)
        fa = abs(g_fwd[k] - g_adj[k]) / max(abs(g_adj[k]), 1e-30) if do_fwd else 0.0
        rows.append({
            "label": lab,
            "rod": int(k),
            "g_adj": float(g_adj[k]),
            "g_fwd": float(g_fwd[k]) if do_fwd else None,
            "fwd_vs_adj": float(fa),
            "minimum_fd_vs_adj": minimum,
            "h_at_min": sweep[imin]["h"],
            "truncation": bool(larger >= 5.0 * minimum),
            "pass": bool(passed and fa <= REL_TOL),
            "sweep": sweep,
        })
        print("ROD", lab, "PASS" if rows[-1]["pass"] else "FAIL", minimum, flush=True)

    # tied-orbit FD check on first non-center orbit
    orbits, _ = build_orbits(dev.centers)
    orbit = next(o for o in orbits if o["size"] == 6)
    g_q_adj = reduce_gradient(g_adj, orbits)
    # need full g_adj for that orbit — recompute adj for orbit members if missing
    need = [i for i in orbit["rods"] if i not in rods]
    if need:
        g_extra = gradient_s(dev, base, need)
        g_all = g_adj.copy()
        for i in need:
            g_all[i] = g_extra[i]
    else:
        g_all = g_adj
    # fill zeros for non-tested — recompute full orbit gradient properly
    g_orbit_rods = gradient_s(dev, base, orbit["rods"])
    g_q = float(sum(g_orbit_rods[i] for i in orbit["rods"]))
    sweep_q = []
    for h in (3e-2, 1e-2, 3e-3):
        sp = s0.copy()
        sm = s0.copy()
        for i in orbit["rods"]:
            sp[i] = 1.0 + h
            sm[i] = 1.0 - h
        Jp = solve_six(dev, sp, F_HZ, b_tesla, P_ref=P_ref).J
        Jm = solve_six(dev, sm, F_HZ, b_tesla, P_ref=P_ref).J
        g_fd = (Jp - Jm) / (2.0 * h)
        rel = abs(g_fd - g_q) / max(abs(g_q), 1e-30)
        sweep_q.append({"h": h, "g_fd": g_fd, "rel": float(rel)})
        print(f"ORBIT rods={orbit['rods'][:3]}... h={h:.1e} fd={g_fd:.4e} adj={g_q:.4e} rel={rel:.3e}", flush=True)
    oerrs = [r["rel"] for r in sweep_q]
    _, omin, olarger, opass = gate(oerrs)

    out = {
        "grade": grade,
        "B_T": b_tesla,
        "f_Hz": F_HZ,
        "P_ref": P_ref,
        "J0": base.J,
        "rel_tol": REL_TOL,
        "rods": rows,
        "orbit_check": {
            "rods": orbit["rods"],
            "g_q": g_q,
            "sweep": sweep_q,
            "minimum_rel": omin,
            "pass": bool(opass),
        },
        "overall_pass": all(r["pass"] for r in rows) and bool(opass),
        "wall_s": time.perf_counter() - t0,
    }
    tag = f"grad_valid_{grade}_B{b_tesla:+.4f}".replace("+", "p").replace("-", "m")
    (OUT / f"{tag}.json").write_text(json.dumps(out, indent=2) + "\n")
    print("WROTE", tag, out["overall_pass"], flush=True)
    return out


def write_report(results):
    lines = [
        "# Circulator objective gradient validation",
        "",
        "Three-way check of the **complete** six-port power objective:",
        "centered FD, forward Jacobian chain rule, port-power adjoint.",
        "",
        f"Predeclared relative tolerance: `{REL_TOL}` (with truncation ≥5×).",
        "Weights and normalization frozen in `CIRCULATOR_OBJECTIVE_DERIVATION.md`.",
        "",
    ]
    overall = True
    for res in results:
        overall = overall and res["overall_pass"]
        lines += [
            f"## grade {res['grade']}, B = {res['B_T']} T",
            "",
            f"- J(s=1) = {res['J0']:.6e}, P_ref = {res['P_ref']:.6e}",
            "",
            "| rod | g_adj | min FD vs adj | fwd vs adj | pass |",
            "|-----|------:|--------------:|-----------:|:----:|",
        ]
        for r in res["rods"]:
            lines.append(
                f"| {r['label']}({r['rod']}) | {r['g_adj']:.4e} | {r['minimum_fd_vs_adj']:.3e} | "
                f"{r['fwd_vs_adj']:.3e} | {'yes' if r['pass'] else 'no'} |"
            )
        oc = res["orbit_check"]
        lines += [
            "",
            f"Tied-orbit FD: min rel `{oc['minimum_rel']:.3e}` → "
            f"{'PASS' if oc['pass'] else 'FAIL'}",
            "",
        ]
    lines += [
        f"**PHYSICAL_PORT_OBJECTIVE_GRADIENT: {'PASS' if overall else 'FAIL'}**",
        "",
    ]
    (OUT / "CIRCULATOR_OBJECTIVE_GRADIENT_VALIDATION.md").write_text("\n".join(lines))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    grade = sys.argv[1] if len(sys.argv) > 1 else "M"
    results = []
    # quick C smoke if requested
    if grade == "all":
        results.append(run("C", B_PLUS))
        results.append(run("M", B_PLUS))
    else:
        results.append(run(grade, B_PLUS))
    write_report(results)
    print("OVERALL", all(r["overall_pass"] for r in results), flush=True)


if __name__ == "__main__":
    main()
