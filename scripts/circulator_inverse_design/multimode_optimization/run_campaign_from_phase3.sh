#!/usr/bin/env bash
# Resume campaign after baselines + C6 FD PASS.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
export PYTHONPATH="scripts/circulator_inverse_design/multimode_optimization:scripts/circulator_inverse_design/multimode_ports:scripts/circulator_inverse_design:scripts/validation/fem_full_device_validation:scripts/validation/fem_gyrotropic_validation:scripts/validation/fem_meep_validation:scripts/validation:scripts:${PYTHONPATH:-}"
PY=/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python
OUT=outputs/circulator_inverse_design/multimode_optimization
mkdir -p "$OUT"
LOG="$OUT/campaign_phase3on.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== RESUME FROM PHASE 3 $(date -Is) ====="
echo "branch=$(git branch --show-current) HEAD=$(git rev-parse --short HEAD)"

# Require FD pass artifact
if [[ ! -f "$OUT/C6_TIED_FD_CHECK.md" ]] || ! grep -q 'C6_TIED_FD: PASS' "$OUT/C6_TIED_FD_CHECK.md"; then
  echo "ERROR: C6_TIED_FD must PASS before optimization. Re-run c6_tied_fd_check.py"
  exit 1
fi

echo "===== PHASE 3: sanity opt on C (maxiter=4) ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start pref_peak --maxiter 4 --no-absorption --tag-suffix _sanity

SANITY_JSON="$OUT/opt_pref_peak_C_sanity.json"
$PY - <<'PY'
import json
from pathlib import Path
p=Path("outputs/circulator_inverse_design/multimode_optimization/opt_pref_peak_C_sanity.json")
d=json.loads(p.read_text())
print("sanity J0", d["J0"], "best", d["best"]["J"], "delta", d["best"]["J"]-d["J0"])
te1_0=d["metrics0"]["TE1_avg"]; te1=d["best"]["metrics"]["TE1_avg"]
print("TE1", te1_0, "->", te1)
if d["best"]["J"] <= d["J0"] - 1e-6:
    raise SystemExit("sanity did not improve J")
print("SANITY_OK")
PY

echo "===== PHASE 4: production C6 on C — prioritize pref_peak ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start pref_peak --maxiter 30 --no-absorption

echo "===== PHASE 4b: broadband warm start (maxiter=20) ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start pref_broadband --maxiter 20 --no-absorption

echo "===== PHASE 4c: uniform warm start (maxiter=20) ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start uniform --maxiter 20 --no-absorption

echo "===== Select best ====="
BEST=$($PY - <<'PY'
import json
from pathlib import Path
out=Path("outputs/circulator_inverse_design/multimode_optimization")
cands=[]
for p in out.glob("opt_*_C.json"):
    if "_sanity" in p.name: continue
    d=json.loads(p.read_text())
    cands.append((d["best"]["J"], str(p)))
cands.sort(reverse=True)
for j,p in cands: print(f"  {p}: J={j:.6f}")
print(cands[0][1])
PY
)
BEST_JSON=$(echo "$BEST" | tail -n1)
echo "BEST_JSON=$BEST_JSON"
echo "$BEST_JSON" > "$OUT/BEST_OPT_PATH.txt"

echo "===== PHASE 6: full12 validation ====="
$PY scripts/circulator_inverse_design/multimode_optimization/validate_opt_full12.py \
  --opt-json "$BEST_JSON" --grade C --label primary

echo "===== PHASE 7: fine-mesh C→M then M→F ====="
$PY scripts/circulator_inverse_design/multimode_optimization/eval_fine_mesh.py \
  --opt-json "$BEST_JSON" --from-grade C --to-grade M --label primary
$PY scripts/circulator_inverse_design/multimode_optimization/eval_fine_mesh.py \
  --opt-json "$BEST_JSON" --from-grade M --to-grade F --label primary || true

echo "===== PHASE 4d: short M refine from best C q ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade M --start custom --q-json "$BEST_JSON" --start-label from_best_C \
  --maxiter 6 --no-absorption --tag-suffix _Mrefine || true

echo "===== PHASE 9: frequency characterization ====="
$PY scripts/circulator_inverse_design/multimode_optimization/freq_characterize.py \
  --opt-json "$BEST_JSON" --grade C --label primary

echo "===== Write reports ====="
$PY - <<'PY'
import json
from pathlib import Path
OUT=Path("outputs/circulator_inverse_design/multimode_optimization")
rows=[]
for p in sorted(OUT.glob("opt_*_C.json")):
    if "_sanity" in p.name: continue
    d=json.loads(p.read_text()); m=d["best"]["metrics"]; m0=d["metrics0"]
    rows.append((d["best"]["J"], p.name, d, m, m0))
rows.sort(reverse=True)
lines=["# MULTIMODE_C6_OPT_REPORT","","Frozen multimode direct-DtN C6 L-BFGS-B.","",
"| start | grade | J0 | J_best | T_des | IL | iso | RL | TE1_0→TE1 | h |",
"|-------|-------|---:|-------:|------:|---:|----:|---:|----------:|--:|"]
for J,name,d,m,m0 in rows:
    lines.append(f"| {d['start_label']} | {d['grade']} | {d['J0']:.4f} | {J:.4f} | {m['T_des_avg']:.4f} | {m['IL_avg_dB']:.2f} | {m['isolation_avg_dB']:.2f} | {m['RL_avg_dB']:.2f} | {m0['TE1_avg']:.3f}→{m['TE1_avg']:.3f} | {d['elapsed_s']/3600:.2f} |")
for p in sorted(OUT.glob("opt_*_M*.json")):
    d=json.loads(p.read_text()); m=d["best"]["metrics"]; m0=d["metrics0"]
    lines.append(f"| {d['start_label']} | {d['grade']} | {d['J0']:.4f} | {d['best']['J']:.4f} | {m['T_des_avg']:.4f} | {m['IL_avg_dB']:.2f} | {m['isolation_avg_dB']:.2f} | {m['RL_avg_dB']:.2f} | {m0['TE1_avg']:.3f}→{m['TE1_avg']:.3f} | {d['elapsed_s']/3600:.2f} |")
lines += ["","Production mesh: FULL91-C (~7 min/eval+grad). M used for transfer + short refine.",""]
(OUT/"MULTIMODE_C6_OPT_REPORT.md").write_text("\n".join(lines)+"\n")
print("wrote C6 report")
PY
$PY scripts/circulator_inverse_design/multimode_optimization/write_final_status.py

echo "===== PHASE 8: free-91 short ====="
$PY scripts/circulator_inverse_design/multimode_optimization/refine_free91_multimode.py \
  --opt-json "$BEST_JSON" --grade C --maxiter 5 || echo "free91 failed/skipped"
$PY scripts/circulator_inverse_design/multimode_optimization/write_final_status.py || true

echo "===== CAMPAIGN END $(date -Is) ====="
