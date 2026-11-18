#!/usr/bin/env bash
# Feasibility bandwidth ladder — NOT Stage-6.
# Sequential only (each job uses n_workers=12). Skips already-complete opts.
#
# Important: Stage-5 q_best is a poor warm start on narrow bands (local basin with
# residual notches). Band A uniform finds a much better 200 MHz solution; wider
# bands therefore also test a warm start from that A-uniform q.
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
A_UNI_OPT="$FEAS/A_4p9_5p1_C6_B0.0500_uniform_opt.json"

echo "FEASIBILITY LADDER START $(date -Is) — Stage-6 NOT STARTED" | tee -a "$LOG"

is_complete () {
  local opt="$1"
  [[ -f "$opt" ]] || return 1
  python -c "import json,sys; d=json.load(open(sys.argv[1])); sys.exit(0 if d.get('status')=='complete' else 1)" "$opt"
}

# Usage: run_case <tag> <band> <start> <maxiter> [extra args]
run_case () {
  local tag="$1"
  local band="$2"
  local start="$3"
  local maxiter="$4"
  shift 4
  local opt="$FEAS/${tag}_opt.json"
  echo "==== OPT ${tag} maxiter=${maxiter} $(date -Is) ====" | tee -a "$LOG"
  if is_complete "$opt"; then
    echo "skip complete $opt" | tee -a "$LOG"
  else
    python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_opt.py \
      --band "$band" --start "$start" --maxiter "$maxiter" --tag "$tag" "$@" 2>&1 | tee -a "$LOG"
  fi
  echo "==== AUDIT ${tag} $(date -Is) ====" | tee -a "$LOG"
  local sum="$FEAS/audit_${tag}_opt_3.125MHz_SUMMARY.json"
  if [[ -f "$sum" ]]; then
    echo "skip audit summary exists $sum" | tee -a "$LOG"
  elif [[ -f "$opt" ]]; then
    python -u scripts/circulator_inverse_design/jesse_method_campaign/feasibility_band_audit.py \
      --opt-json "$opt" --df-MHz 3.125 2>&1 | tee -a "$LOG" || true
  else
    echo "missing opt $opt — skip audit" | tee -a "$LOG"
  fi
}

# Narrow → wider
run_case A_4p9_5p1_C6_B0.0500_stage5 A_4p9_5p1 stage5 20
run_case A_4p9_5p1_C6_B0.0500_uniform A_4p9_5p1 uniform 15

# Band B: Stage-5 basin + A-uniform warm start (materially different)
run_case B_4p75_5p25_C6_B0.0500_stage5 B_4p75_5p25 stage5 20
if [[ -f "$A_UNI_OPT" ]]; then
  run_case B_4p75_5p25_C6_B0.0500_Auni B_4p75_5p25 "file:${A_UNI_OPT}" 15
else
  echo "WARN: missing $A_UNI_OPT — skip B Auni start" | tee -a "$LOG"
fi

# Band C: Stage-5 + A-uniform warm start
run_case C_4p5_5p5_C6_B0.0500_stage5 C_4p5_5p5 stage5 15
if [[ -f "$A_UNI_OPT" ]]; then
  run_case C_4p5_5p5_C6_B0.0500_Auni C_4p5_5p5 "file:${A_UNI_OPT}" 12
fi

# DOF probe on band B (warm from Stage-5 mapped C6→C3; also from Auni if available)
echo "==== C3 DOF probe on B $(date -Is) ====" | tee -a "$LOG"
run_case B_4p75_5p25_C3_B0.0500_stage5 B_4p75_5p25 stage5 15 --symmetry C3
if [[ -f "$A_UNI_OPT" ]]; then
  run_case B_4p75_5p25_C3_B0.0500_Auni B_4p75_5p25 "file:${A_UNI_OPT}" 12 --symmetry C3
fi

# B-field probe on band C (0.05 covered; 0.10 = Li et al. literature-motivated)
# Prefer A-uniform warm start so B-sensitivity is not confounded by Stage-5 basin.
echo "==== B=0.10 T probe on C $(date -Is) ====" | tee -a "$LOG"
if [[ -f "$A_UNI_OPT" ]]; then
  run_case C_4p5_5p5_C6_B0.1000_Auni C_4p5_5p5 "file:${A_UNI_OPT}" 12 --b-tesla 0.10
else
  run_case C_4p5_5p5_C6_B0.1000_stage5 C_4p5_5p5 stage5 15 --b-tesla 0.10
fi

# Band D only after intermediate results (still same physics)
echo "==== Band D 4.25–5.75 $(date -Is) ====" | tee -a "$LOG"
if [[ -f "$A_UNI_OPT" ]]; then
  run_case D_4p25_5p75_C6_B0.0500_Auni D_4p25_5p75 "file:${A_UNI_OPT}" 10
else
  run_case D_4p25_5p75_C6_B0.0500_stage5 D_4p25_5p75 stage5 12
fi

echo "FEASIBILITY LADDER COMPLETE $(date -Is) — Stage-6 NOT STARTED" | tee -a "$LOG"
echo "STAGE-6: NOT STARTED / AWAITING AUTHORIZATION" | tee -a "$LOG"
