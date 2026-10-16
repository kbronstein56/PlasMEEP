#!/usr/bin/env python3
"""Profile baseline vs efficient adjoint; require numerical agreement.

Writes ADJOINT_EFFICIENCY_BENCHMARK.md. Does not modify frozen SF code paths
used by the running single-frequency campaign.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "broadband_method"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from efficient_adjoint import baseline_eval_and_grad, eval_and_grad_efficient  # noqa: E402
from multimode_opt_core import load_start, orbits_and_masks  # noqa: E402

REL_BAR = 1e-3
ABS_BAR = 1e-5


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--grade", default="C")
    ap.add_argument(
        "--design",
        default="outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json",
    )
    ap.add_argument("--fd-orbit", type=int, default=1)
    ap.add_argument("--fd-h", type=float, default=3e-4)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orbits, _, _ = orbits_and_masks()
    # load q from opt json or Pref tag
    p = Path(args.design)
    if p.suffix == ".json" and p.exists():
        d = json.loads(p.read_text())
        q = np.asarray(d["best"]["q"], float)
        tag = p.stem
    else:
        q, _, tag = load_start(args.design, orbits)

    print("baseline eval+grad ...", flush=True)
    t0 = time.perf_counter()
    bun_b, g_b, s = baseline_eval_and_grad(args.grade, q, orbits)
    t_base = time.perf_counter() - t0
    print(f"  J={bun_b['J']:.6f} |g|={np.linalg.norm(g_b):.4e} dt={t_base:.1f}s", flush=True)

    print("efficient eval+grad ...", flush=True)
    t0 = time.perf_counter()
    bun_e, g_e, _ = eval_and_grad_efficient(
        args.grade, q, orbits, profile=True
    )
    t_eff = time.perf_counter() - t0
    tim = bun_e["timings"]
    print(f"  J={bun_e['J']:.6f} |g|={np.linalg.norm(g_e):.4e} dt={t_eff:.1f}s", flush=True)
    print("  timings:", json.dumps({k: round(v, 3) for k, v in tim.items()}), flush=True)

    dJ = abs(bun_e["J"] - bun_b["J"])
    # Compare TE0 columns of S
    S_b, S_e = bun_b["S"], bun_e["S"]
    te0 = [2 * j for j in range(6)]
    s_err = float(np.max(np.abs(S_e[:, te0] - S_b[:, te0])))
    g_rel = float(np.linalg.norm(g_e - g_b) / max(np.linalg.norm(g_b), 1e-30))
    g_abs = float(np.linalg.norm(g_e - g_b))

    # FD spot on one orbit using efficient J
    from multimode_opt_core import S_HI as SH, S_LO as SL

    oid = args.fd_orbit
    h = args.fd_h
    qp, qm = q.copy(), q.copy()
    qp[oid] = min(SH, qp[oid] + h)
    qm[oid] = max(SL, qm[oid] - h)
    hp = float(qp[oid] - q[oid])
    hm = float(q[oid] - qm[oid])
    print(f"FD orbit {oid} h={h} ...", flush=True)
    bun_p, _, _ = eval_and_grad_efficient(args.grade, qp, orbits)
    bun_m, _, _ = eval_and_grad_efficient(args.grade, qm, orbits)
    g_fd = float((bun_p["J"] - bun_m["J"]) / (hp + hm))
    g_adj = float(g_e[oid])
    fd_rel = abs(g_fd - g_adj) / max(abs(g_fd), abs(g_adj), 1e-30)

    ok = (
        dJ < 1e-10
        and s_err < 1e-8
        and g_rel <= REL_BAR
        and fd_rel <= REL_BAR
    )
    report = {
        "grade": args.grade,
        "design": tag,
        "baseline_s": t_base,
        "efficient_s": t_eff,
        "speedup": t_base / max(t_eff, 1e-9),
        "timings_efficient": tim,
        "J_baseline": bun_b["J"],
        "J_efficient": bun_e["J"],
        "dJ": dJ,
        "S_te0_col_max_abs_err": s_err,
        "grad_rel_err": g_rel,
        "grad_abs_err": g_abs,
        "fd_orbit": oid,
        "fd_h": h,
        "fd_grad": g_fd,
        "adj_grad": g_adj,
        "fd_rel_err": fd_rel,
        "n_forward_rhs_eff": 6,
        "n_adjoint_rhs_eff": 6,
        "n_lu_eff": 1,
        "PASS": ok,
        "rel_bar": REL_BAR,
    }
    (OUT / "adjoint_efficiency_benchmark.json").write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# ADJOINT_EFFICIENCY_BENCHMARK",
        "",
        f"Grade `{args.grade}`, design `{tag}`.",
        "",
        "## Wall time (one eval+grad)",
        "",
        f"| path | seconds | RHS / LU |",
        f"|------|--------:|----------|",
        f"| baseline (`_channel_basis` + `direct_S_from_factor`) | {t_base:.1f} | ~12+12 fwd/adj + duplicate S solves |",
        f"| efficient (TE0-only + scalar adjoint) | {t_eff:.1f} | **6 fwd + 6 adj, 1 LU** |",
        f"| speedup | **{report['speedup']:.2f}×** | |",
        "",
        "## Efficient path breakdown",
        "",
        "| stage | s |",
        "|-------|--:|",
    ]
    for k, v in tim.items():
        lines.append(f"| {k} | {v:.3f} |")
    lines += [
        "",
        "## Accuracy vs frozen baseline",
        "",
        f"| check | value | bar |",
        f"|-------|------:|----:|",
        f"| ΔJ | {dJ:.3e} | ~0 |",
        f"| max ‖ΔS[:,TE0]‖ | {s_err:.3e} | 1e-8 |",
        f"| ‖Δg‖/‖g‖ | {g_rel:.3e} | {REL_BAR} |",
        f"| FD orbit {oid} rel | {fd_rel:.3e} | {REL_BAR} |",
        "",
        f"**EFFICIENT_ADJOINT: {'PASS' if ok else 'FAIL'}**",
        "",
        "Redundant solves removed: TE1 incident columns not factored into opt;",
        "`direct_S_from_factor` not called (S from U); channel-wise V replaced by",
        "six scalar-objective adjoint RHS.",
        "",
    ]
    (OUT / "ADJOINT_EFFICIENCY_BENCHMARK.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: report[k] for k in ("PASS", "speedup", "grad_rel_err", "fd_rel_err", "efficient_s", "baseline_s")}, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
