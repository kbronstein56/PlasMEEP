#!/usr/bin/env bash
# Autonomous multimode C6 optimization campaign.
# Logs to outputs/circulator_inverse_design/multimode_optimization/campaign.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
export PYTHONPATH="scripts/circulator_inverse_design/multimode_optimization:scripts/circulator_inverse_design/multimode_ports:scripts/circulator_inverse_design:scripts/validation/fem_full_device_validation:scripts/validation/fem_gyrotropic_validation:scripts/validation/fem_meep_validation:scripts/validation:scripts:${PYTHONPATH:-}"
PY=/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python
OUT=outputs/circulator_inverse_design/multimode_optimization
mkdir -p "$OUT"
LOG="$OUT/campaign.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== CAMPAIGN START $(date -Is) ====="
echo "branch=$(git branch --show-current) HEAD=$(git rev-parse --short HEAD)"

echo "===== PHASE 2: baselines on C ====="
$PY scripts/circulator_inverse_design/multimode_optimization/run_baselines.py --grade C

echo "===== PHASE 1: C6 tied FD spot check ====="
$PY scripts/circulator_inverse_design/multimode_optimization/c6_tied_fd_check.py \
  --grade C --orbits 0,1 --hs 1e-3 --design opt_single_uniform_M_Bp0.0500

echo "===== PHASE 3: sanity opt on C (maxiter=4) ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start pref_peak --maxiter 4 --no-absorption --tag-suffix _sanity

# If sanity improved J, continue production
SANITY_JSON="$OUT/opt_pref_peak_C_sanity.json"
if [[ ! -f "$SANITY_JSON" ]]; then
  echo "SANITY FAILED: missing $SANITY_JSON"; exit 1
fi
$PY - <<'PY'
import json
from pathlib import Path
p=Path("outputs/circulator_inverse_design/multimode_optimization/opt_pref_peak_C_sanity.json")
d=json.loads(p.read_text())
print("sanity J0", d["J0"], "best", d["best"]["J"], "delta", d["best"]["J"]-d["J0"])
if d["best"]["J"] <= d["J0"] - 1e-6:
    raise SystemExit("sanity did not improve J")
print("SANITY_OK")
PY

echo "===== PHASE 4: production C6 on C (pref_peak, maxiter=30) ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start pref_peak --maxiter 30 --no-absorption

echo "===== PHASE 4b: production warm-start broadband ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start pref_broadband --maxiter 25 --no-absorption

echo "===== PHASE 4c: production warm-start uniform ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade C --start uniform --maxiter 25 --no-absorption

echo "===== Select best among three starts ====="
BEST=$($PY - <<'PY'
import json
from pathlib import Path
out=Path("outputs/circulator_inverse_design/multimode_optimization")
cands=[]
for p in out.glob("opt_*_C.json"):
    if "_sanity" in p.name: continue
    d=json.loads(p.read_text())
    cands.append((d["best"]["J"], p.name, str(p)))
cands.sort(reverse=True)
for j,n,p in cands:
    print(f"  {n}: J={j:.6f}")
print(cands[0][2])
PY
)
BEST_JSON=$(echo "$BEST" | tail -n1)
echo "BEST_JSON=$BEST_JSON"
echo "$BEST_JSON" > "$OUT/BEST_OPT_PATH.txt"

echo "===== PHASE 6: full12 validation of best ====="
$PY scripts/circulator_inverse_design/multimode_optimization/validate_opt_full12.py \
  --opt-json "$BEST_JSON" --grade C --label primary

echo "===== PHASE 7: fine-mesh transfer C→M (affordable) then try F ====="
$PY scripts/circulator_inverse_design/multimode_optimization/eval_fine_mesh.py \
  --opt-json "$BEST_JSON" --from-grade C --to-grade M --label primary
$PY scripts/circulator_inverse_design/multimode_optimization/eval_fine_mesh.py \
  --opt-json "$BEST_JSON" --from-grade M --to-grade F --label primary || true

echo "===== PHASE 4d: optional short M refinement from best C ====="
$PY scripts/circulator_inverse_design/multimode_optimization/optimize_multimode_c6.py \
  --grade M --start pref_peak --maxiter 8 --no-absorption --tag-suffix _Mrefine || true
# If we have a better C best, warm-start M from that q via a small helper:
$PY - <<'PY' || true
import json, numpy as np
from pathlib import Path
import sys
sys.path[:0]=["scripts/circulator_inverse_design/multimode_optimization"]
from multimode_opt_core import OUT, dump_json
best_path=Path(Path("outputs/circulator_inverse_design/multimode_optimization/BEST_OPT_PATH.txt").read_text().strip())
d=json.loads(best_path.read_text())
# write a synthetic start json compatible with load_s via opt dump
start={
  "best": {"s": d["best"]["s"], "q": d["best"]["q"], "J": d["best"]["J"]},
}
dump_json(OUT/"_warmstart_from_best_C.json", start)
print("wrote warmstart")
PY

echo "===== PHASE 9: frequency characterization on C ====="
$PY scripts/circulator_inverse_design/multimode_optimization/freq_characterize.py \
  --opt-json "$BEST_JSON" --grade C --label primary

echo "===== Write C6 opt report + final status ====="
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
lines=["# MULTIMODE_C6_OPT_REPORT","","Frozen multimode direct-DtN C6 L-BFGS-B on FULL91-C.","",
"| start | J0 | J_best | T_des | IL | iso | RL | TE1_0→TE1 | elapsed_h |","|-------|---:|-------:|------:|---:|----:|---:|----------:|---------:|"]
for J,name,d,m,m0 in rows:
    lines.append(f"| {d['start_label']} | {d['J0']:.4f} | {J:.4f} | {m['T_des_avg']:.4f} | {m['IL_avg_dB']:.2f} | {m['isolation_avg_dB']:.2f} | {m['RL_avg_dB']:.2f} | {m0['TE1_avg']:.3f}→{m['TE1_avg']:.3f} | {d['elapsed_s']/3600:.2f} |")
lines += ["","Best start retained as primary; Pref designs are warm starts only.",""]
(OUT/"MULTIMODE_C6_OPT_REPORT.md").write_text("\n".join(lines)+"\n")
print("wrote C6 report", len(rows), "runs")
PY
$PY scripts/circulator_inverse_design/multimode_optimization/write_final_status.py

echo "===== PHASE 8: optional free-91 (short, if primary exists) ====="
if [[ -f "$OUT/BEST_OPT_PATH.txt" ]]; then
  BEST_JSON=$(cat "$OUT/BEST_OPT_PATH.txt")
  $PY scripts/circulator_inverse_design/multimode_optimization/refine_free91_multimode.py \
    --opt-json "$BEST_JSON" --grade C --maxiter 5 || echo "free91 skipped/failed"
  $PY scripts/circulator_inverse_design/multimode_optimization/write_final_status.py || true
fi

echo "===== CAMPAIGN END $(date -Is) ====="
