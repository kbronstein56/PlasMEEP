#!/usr/bin/env bash
# Finish single-frequency campaign after the three warm starts completed.
# Does NOT alter the frozen objective.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
export PYTHONPATH="scripts/circulator_inverse_design/multimode_optimization:scripts/circulator_inverse_design/multimode_ports:scripts/circulator_inverse_design:scripts/validation/fem_full_device_validation:scripts/validation/fem_gyrotropic_validation:scripts/validation/fem_meep_validation:scripts/validation:scripts:${PYTHONPATH:-}"
PY=/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python
OUT=outputs/circulator_inverse_design/multimode_optimization
LOG="$OUT/campaign_finish.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== FINISH SF CAMPAIGN $(date -Is) ====="
BEST_JSON=$(cat "$OUT/BEST_OPT_PATH.txt")
echo "BEST_JSON=$BEST_JSON"

echo "===== C6 opt report ====="
$PY - <<'PY'
import json
from pathlib import Path
OUT=Path("outputs/circulator_inverse_design/multimode_optimization")
rows=[]
for p in sorted(OUT.glob("opt_*_C.json")):
    if "_sanity" in p.name: continue
    d=json.loads(p.read_text())
    if "best" not in d: continue
    m=d["best"]["metrics"]; m0=d["metrics0"]
    rows.append((d["best"]["J"], p.name, d, m, m0))
rows.sort(reverse=True)
lines=["# MULTIMODE_C6_OPT_REPORT","",
"Frozen multimode direct-DtN C6 L-BFGS-B on FULL91-C.","",
"| start | J0 | J_best | T_des | IL | iso | RL | TE1_0→TE1 | h |",
"|-------|---:|-------:|------:|---:|----:|---:|----------:|--:|"]
for J,name,d,m,m0 in rows:
    lines.append(
        f"| {d['start_label']} | {d['J0']:.4f} | {J:.4f} | {m['T_des_avg']:.4f} | "
        f"{m['IL_avg_dB']:.2f} | {m['isolation_avg_dB']:.2f} | {m['RL_avg_dB']:.2f} | "
        f"{m0['TE1_avg']:.3f}→{m['TE1_avg']:.3f} | {d['elapsed_s']/3600:.2f} |"
    )
lines += ["","Primary = highest J among warm starts (uniform unexpectedly best).",
"Pref designs remain warm starts only under the multimode S.",""]
(OUT/"MULTIMODE_C6_OPT_REPORT.md").write_text("\n".join(lines)+"\n")
print("wrote C6 report", len(rows), "runs")
PY

echo "===== PHASE 6: full12 validation ====="
$PY scripts/circulator_inverse_design/multimode_optimization/validate_opt_full12.py \
  --opt-json "$BEST_JSON" --grade C --label primary

echo "===== PHASE 7: fine-mesh C→M then M→F ====="
$PY scripts/circulator_inverse_design/multimode_optimization/eval_fine_mesh.py \
  --opt-json "$BEST_JSON" --from-grade C --to-grade M --label primary
$PY scripts/circulator_inverse_design/multimode_optimization/eval_fine_mesh.py \
  --opt-json "$BEST_JSON" --from-grade M --to-grade F --label primary || true

echo "===== PHASE 9: frequency characterization ====="
$PY scripts/circulator_inverse_design/multimode_optimization/freq_characterize.py \
  --opt-json "$BEST_JSON" --grade C --label primary

echo "===== PHASE 8: free-91 short ====="
$PY scripts/circulator_inverse_design/multimode_optimization/refine_free91_multimode.py \
  --opt-json "$BEST_JSON" --grade C --maxiter 5 || echo "free91 skipped/failed"

echo "===== Final status ====="
$PY scripts/circulator_inverse_design/multimode_optimization/write_final_status.py
echo "===== SF CAMPAIGN COMPLETE $(date -Is) ====="
