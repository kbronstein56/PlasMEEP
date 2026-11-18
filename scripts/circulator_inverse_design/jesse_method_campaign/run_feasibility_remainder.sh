#!/usr/bin/env bash
# Feasibility remainder: C3 DOF + B=0.10 probes (+ optional Band D).
# Assumes A/B/C C6@0.05 already complete. Stage-6 NOT STARTED.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
# shellcheck disable=SC1091
source /home/daq_user/miniconda3/etc/profile.d/conda.sh
conda activate plasmeep_mpi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1 PLASMEEP_A2_FAST_KERNEL=1

FEAS=outputs/circulator_inverse_design/jesse_method_campaign/feasibility
LOG="$FEAS/ladder_remainder_nohup.out"
B_S5="$FEAS/B_4p75_5p25_C6_B0.0500_stage5_opt.json"
mkdir -p "$FEAS"
echo "REMAINDER START $(date -Is) STAGE-6 NOT STARTED" | tee -a "$LOG"

is_complete () {
  [[ -f "$1" ]] || return 1
  python -c "import json,sys; d=json.load(open(sys.argv[1])); sys.exit(0 if d.get('status')=='complete' else 1)" "$1"
}

run_case () {
  local tag="$1"; local band="$2"; local start="$3"; local maxiter="$4"; shift 4
  local opt="$FEAS/${tag}_opt.json"
  echo "==== OPT ${tag} $(date -Is) ====" | tee -a "$LOG"
  if is_complete "$opt"; then
    echo "skip complete $opt" | tee -a "$LOG"
  else
    python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_opt.py \
      --band "$band" --start "$start" --maxiter "$maxiter" --tag "$tag" "$@" 2>&1 | tee -a "$LOG"
  fi
  local sum="$FEAS/audit_${tag}_opt_3.125MHz_SUMMARY.json"
  if [[ -f "$sum" ]]; then
    echo "skip audit $sum" | tee -a "$LOG"
  elif [[ -f "$opt" ]]; then
    python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_audit.py \
      --opt-json "$opt" --df-MHz 3.125 2>&1 | tee -a "$LOG" || true
  fi
}

# Phase 4: C3 vs C6 on Band B
run_case B_4p75_5p25_C3_B0.0500_stage5 B_4p75_5p25 stage5 12 --symmetry C3
if [[ -f "$B_S5" ]]; then
  run_case B_4p75_5p25_C3_B0.0500_Bs5 B_4p75_5p25 "file:${B_S5}" 10 --symmetry C3
fi

# Phase 3: B=0.10 T on Band C (literature-motivated; labeled extrapolation)
run_case C_4p5_5p5_C6_B0.1000_stage5 C_4p5_5p5 stage5 12 --b-tesla 0.10
if [[ -f "$B_S5" ]]; then
  run_case C_4p5_5p5_C6_B0.1000_Bs5 C_4p5_5p5 "file:${B_S5}" 8 --b-tesla 0.10
fi

# Band D from best 0.5 GHz design (B stage5)
if [[ -f "$B_S5" ]]; then
  run_case D_4p25_5p75_C6_B0.0500_Bs5 D_4p25_5p75 "file:${B_S5}" 8
fi

echo "REMAINDER COMPLETE $(date -Is) STAGE-6 NOT STARTED" | tee -a "$LOG"
