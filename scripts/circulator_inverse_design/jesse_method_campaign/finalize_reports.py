#!/usr/bin/env python3
"""Update campaign markdown reports from completed JSON artifacts."""
from __future__ import annotations

import json
from pathlib import Path

OUT = Path(__file__).resolve().parents[3] / "outputs/circulator_inverse_design/jesse_method_campaign"


def load(name):
    p = OUT / name
    if not p.exists():
        return None
    return json.loads(p.read_text())


def main():
    lines = ["# Campaign results snapshot", "", f"Generated from artifacts in `{OUT.name}`.", ""]

    fm = load("forward_mode_bench.json")
    if fm:
        lines += ["## Forward-mode vs adjoint", ""]
        adj = fm["adjoint"]
        lines.append(f"- Adjoint wall: **{adj['wall_s']:.1f} s**, LU={adj['n_LU']}, RHS={adj['n_forward_rhs']}+{adj['n_adjoint_rhs']}")
        if "forward" in fm:
            f = fm["forward"]
            lines.append(
                f"- Forward full est: **{f['wall_full_est_s']:.1f} s**, sens RHS={f['n_sensitivity_rhs']}, "
                f"rel err vs adj={f['rel_err_vs_adjoint']:.3e}, pass={f['pass']}"
            )
            for row in fm.get("scaling_estimate", []):
                lines.append(
                    f"- Scaling {row['label']} (n={row['n_param']}): adj~{row['est_adjoint_s']:.0f}s "
                    f"fwd~{row['est_forward_s']:.0f}s → prefer **{row['preferred']}**"
                )
        lines.append("")

    mb = load("mode_baselines.json")
    if mb:
        lines += ["## Mode A vs B baselines", "", "| design | J_A | J_B | des0 | des1 | iso TE0 | iso tot |",
                  "|--------|----:|----:|-----:|-----:|--------:|--------:|"]
        for r in mb["rows"]:
            lines.append(
                f"| {r['design']} | {r['J_A']:.4f} | {r['J_B']:.4f} | {r['des0']:.4f} | "
                f"{r['des1']:.4f} | {r['iso_te0_dB']:.2f} | {r['iso_tot_dB']:.2f} |"
            )
        lines.append("")

    for tag in ("200MHz", "100MHz", "50MHz"):
        sp = load(f"spectral_{tag}.json")
        if not sp:
            continue
        w = sp.get("worst_isolation_dB", {})
        lines += [
            f"## Spectral {tag}",
            f"- n={sp.get('n_freq')} worst iso={w.get('value')} dB @ {w.get('f_GHz')} GHz",
            "",
        ]

    comp = load("spectral_grid_comparison.json")
    if comp:
        lines += ["## Grid comparison", f"```json\n{json.dumps(comp, indent=2)}\n```", ""]

    rot = load("rotation_source_check.json")
    if rot:
        lines += [
            "## Rotation source check",
            f"- C6 max rel err={rot.get('C6_max_rel_err')}, pass={rot.get('C6_pass_1e-6')}",
            "",
        ]

    sm = load("softmin_e2e_fd.json")
    if sm:
        lines += [
            "## Soft-min E2E FD",
            f"- J_BB={sm.get('J_BB')}, weight FD pass={sm.get('softmin_weight_fd_pass')}, "
            f"overall pass={sm.get('pass')}",
            "",
        ]

    ab = load("mode_ab_short_opt.json")
    if ab:
        lines += ["## Mode A/B short opt", ""]
        for k, v in ab.items():
            lines.append(f"- {k}: J_best={v.get('J_best')}, nit={v.get('nit')}")
        lines.append("")

    c36 = load("c3_c6_short_opt.json")
    if c36:
        lines += ["## C3 vs C6 short opt", f"- map |Δs|={c36.get('map_max_abs_ds')}"]
        for k in ("C6", "C3"):
            if k in c36:
                lines.append(f"- {k}: n_var={c36[k].get('n_var')} J_best={c36[k].get('J_best')}")
        lines.append("")

    (OUT / "CAMPAIGN_RESULTS_SNAPSHOT.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUT / "CAMPAIGN_RESULTS_SNAPSHOT.md")


if __name__ == "__main__":
    main()
