#!/usr/bin/env bash
# Feasibility bandwidth ladder — NOT Stage-6.
# Runs sequentially to avoid oversubscribing 12-worker FEM jobs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
# shellcheck disable=SC1091
source /home/daq_user/miniconda3/etc/profile.d/conda.sh
conda activate plasmeep_mpi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export PLASMEEP_A2_FAST_KERNEL=1

FEAS=outputs/circulator_inverse_design/jesse_method_campaign/feasibility
mkdir -p "$FEAS"
LOG="$FEAS/ladder_chain_nohup.out"

echo "FEASIBILITY LADDER START $(date -Is) — Stage-6 NOT STARTED" | tee -a "$LOG"

run_band () {
  local band="$1"
  local start="$2"
  local maxiter="$3"
  local tag="${band}_C6_B0.0500_${start}"
  local opt="$FEAS/${tag}_opt.json"
  echo "==== OPT ${tag} maxiter=${maxiter} $(date -Is) ====" | tee -a "$LOG"
  if [[ -f "$opt" ]] && python -c "import json;print(json.load(open('$opt')).get('status'))" | grep -q complete; then
    echo "skip complete $opt" | tee -a "$LOG"
  else
    python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_opt.py \
      --band "$band" --start "$start" --maxiter "$maxiter" 2>&1 | tee -a "$LOG"
  fi
  echo "==== AUDIT ${tag} $(date -Is) ====" | tee -a "$LOG"
  python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_audit.py \
    --opt-json "$opt" --df-MHz 3.125 2>&1 | tee -a "$LOG" || true
}

# Narrow → wider. Reuse Stage-5 warm start; one alternate start on A only.
run_band A_4p9_5p1 stage5 20
run_band A_4p9_5p1 uniform 15
run_band B_4p75_5p25 stage5 20
run_band C_4p5_5p5 stage5 15

# Optional DOF probe on band B (manageable)
echo "==== C3 DOF probe on B $(date -Is) ====" | tee -a "$LOG"
python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_opt.py \
  --band B_4p75_5p25 --start stage5 --symmetry C3 --maxiter 15 \
  --tag B_4p75_5p25_C3_B0.0500_stage5 2>&1 | tee -a "$LOG" || true

# B-field probe on band C: production 0.05 T already done; literature-motivated 0.10 T labeled
echo "==== B=0.10 T probe on C $(date -Is) ====" | tee -a "$LOG"
python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_opt.py \
  --band C_4p5_5p5 --start stage5 --b-tesla 0.10 --maxiter 15 \
  --tag C_4p5_5p5_C6_B0.1000_stage5 2>&1 | tee -a "$LOG" || true

echo "FEASIBILITY LADDER COMPLETE $(date -Is) — Stage-6 NOT STARTED" | tee -a "$LOG"
