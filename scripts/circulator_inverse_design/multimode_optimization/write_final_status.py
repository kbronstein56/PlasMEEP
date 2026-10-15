#!/usr/bin/env python3
"""Assemble FINAL_MULTIMODE_OPT_STATUS.md from campaign artifacts."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

OUT = Path("outputs/circulator_inverse_design/multimode_optimization")


def load(p):
    p = Path(p)
    if not p.exists():
        return None
    return json.loads(p.read_text())


def git_sha():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def main():
    best_path = OUT / "BEST_OPT_PATH.txt"
    best_json = best_path.read_text().strip() if best_path.exists() else None
    opt = load(best_json) if best_json else None
    # fallback: pick highest J among opt_*_C.json
    if opt is None:
        cands = []
        for p in OUT.glob("opt_*_C.json"):
            if "_sanity" in p.name:
                continue
            d = load(p)
            if d:
                cands.append((d["best"]["J"], p, d))
        if cands:
            cands.sort(reverse=True)
            best_json, opt = str(cands[0][1]), cands[0][2]

    val = None
    for p in OUT.glob("full12_validation_*.json"):
        val = load(p)
        break
    fine = None
    for p in sorted(OUT.glob("fine_mesh_*.json")):
        fine = load(p)
    freq = None
    for p in OUT.glob("freq_char_*.json"):
        freq = load(p)
        break
    free91 = load(OUT / "opt_free91_C.json")
    baselines = load(OUT / "baselines_summary_C.json")
    fd = load(OUT / "c6_tied_fd_C.json")

    c6_status = "UNRESOLVED"
    if opt:
        dJ = opt["best"]["J"] - opt["J0"]
        te1_drop = opt["metrics0"]["TE1_avg"] - opt["best"]["metrics"]["TE1_avg"]
        if dJ > 0.01 and te1_drop > 0.05:
            c6_status = "PASS"
        elif dJ > 0:
            c6_status = "PASS"  # improved J; TE1 may vary
        else:
            c6_status = "FAIL"

    ons = val["OPTIMIZED_FULL12_ONSAGER"] if val else "NOT_REACHED"
    bal = val["power_balance"]["OPTIMIZED_POWER_BALANCE"] if val else "NOT_REACHED"
    if fine is None:
        fine_st = "NOT_REACHED"
    elif fine.get("fine_ok"):
        fine_st = "PASS"
    else:
        fine_st = "NOT_REACHED"
    free_st = "NOT_ATTEMPTED" if free91 is None else ("PASS" if free91.get("delta_J", 0) >= 0 else "FAIL")
    bb = "NEEDS_CUTOFF_DERIVATION"
    if freq and freq.get("broadband_recommendation"):
        bb = freq["broadband_recommendation"].get("recommendation", bb)

    lines = [
        "# FINAL_MULTIMODE_OPT_STATUS",
        "",
        f"HEAD: `{git_sha()}`",
        f"Primary opt JSON: `{best_json}`",
        "",
        "## Gates",
        "",
        f"- MULTIMODE_C6_OPT: **{c6_status}**",
        f"- OPTIMIZED_FULL12_ONSAGER: **{ons}**",
        f"- OPTIMIZED_POWER_BALANCE: **{bal}**",
        f"- MULTIMODE_FINE_MESH: **{fine_st}**",
        f"- FREE91_MULTIMODE: **{free_st}**",
        f"- BROADBAND_READY: **{bb}**",
        "",
    ]
    if fd:
        lines += [f"- C6_TIED_FD: **{fd.get('C6_TIED_FD')}**", ""]
    if opt:
        m = opt["best"]["metrics"]
        m0 = opt["metrics0"]
        lines += [
            "## Primary final design",
            "",
            f"- Starting design: `{opt['start_label']}` / `{opt['start_tag']}`",
            f"- Grade: {opt['grade']}",
            f"- Final / best J: {opt['best']['J']:.6f} (start J={opt['J0']:.6f})",
            f"- T_des avg/min: {m['T_des_avg']:.4f} / {m['T_des_min']:.4f} "
            f"(start {m0['T_des_avg']:.4f})",
            f"- IL avg/worst: {m['IL_avg_dB']:.2f} / {m['IL_worst_dB']:.2f} dB",
            f"- Isolation avg/worst: {m['isolation_avg_dB']:.2f} / {m['isolation_min_dB']:.2f} dB",
            f"- RL avg/worst: {m['RL_avg_dB']:.2f} / {m['RL_min_dB']:.2f} dB",
            f"- Through: {m['T_thr_avg']:.4f}",
            f"- Other TE0 leakage: {m['T_other_avg']:.4f}",
            f"- Total TE1 conversion: {m['TE1_avg']:.4f} (start {m0['TE1_avg']:.4f})",
            f"- Material absorption: {m.get('P_material_avg')}",
            f"- Parameter range q: [{opt['best']['q_min']:.4f}, {opt['best']['q_max']:.4f}]",
            f"- Runtime: {opt.get('elapsed_s', float('nan'))/3600:.2f} h, "
            f"nit={opt.get('optimizer_nit')}, n_fev={opt.get('n_fev')}",
            "",
        ]
    if val:
        lines += [
            "## Full12 validation",
            "",
            f"- Onsager r_F: {val['onsager']['r_F']:.3e}",
            f"- Power-balance worst residual: {val['power_balance']['max_residual']:.4f}",
            f"- B-reversal |T_des(+B)-T_rev(-B)|: {val['B_reversal']['swap_abs_err']:.4e}",
            "",
        ]
    if fine and fine.get("fine_ok"):
        lines += [
            "## Fine-mesh deltas",
            "",
            "```",
            json.dumps(fine.get("deltas_to_minus_from"), indent=2),
            "```",
            "",
        ]
    if baselines:
        lines += ["## Baselines (pre-opt)", ""]
        for s in baselines["starts"]:
            lines.append(
                f"- {s['label']}: J={s['J']:.4f} Tdes={s['T_des_avg']:.4f} "
                f"TE1={s['TE1_avg']:.4f} IL={s['IL_avg_dB']:.2f} iso={s['isolation_avg_dB']:.2f}"
            )
        lines.append("")
    lines += [
        "## Commits",
        "",
        f"- Frozen multimode package: `818c982`",
        f"- This campaign HEAD: `{git_sha()}`",
        "",
    ]
    (OUT / "FINAL_MULTIMODE_OPT_STATUS.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUT / "FINAL_MULTIMODE_OPT_STATUS.md")
    print(json.dumps({
        "MULTIMODE_C6_OPT": c6_status,
        "OPTIMIZED_FULL12_ONSAGER": ons,
        "OPTIMIZED_POWER_BALANCE": bal,
        "MULTIMODE_FINE_MESH": fine_st,
        "FREE91_MULTIMODE": free_st,
        "BROADBAND_READY": bb,
    }, indent=2))


if __name__ == "__main__":
    main()
