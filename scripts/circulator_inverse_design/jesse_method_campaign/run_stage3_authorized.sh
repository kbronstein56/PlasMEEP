#!/usr/bin/env bash
# Authorized Stage-3 optimization + post dense audit.
# Does not change physics. Resume-safe via bb_stage3_checkpoint.json.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"  # .../PlasMEEP
cd "$ROOT"
echo "ROOT=$ROOT"
# shellcheck disable=SC1091
source /home/daq_user/miniconda3/etc/profile.d/conda.sh
conda activate plasmeep_mpi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1

OUT=outputs/circulator_inverse_design/jesse_method_campaign
mkdir -p "$OUT"

echo "========================================================================"
echo "STAGE-3 AUTHORIZED PIPELINE"
echo "cwd=$ROOT"
echo "OMP=$OMP_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS MKL=$MKL_NUM_THREADS"
echo "========================================================================"

python -u scripts/circulator_inverse_design/jesse_method_campaign/production_broadband_opt_stage3.py \
  --n-workers 12 --omp-threads-per-worker 1 --maxiter 30

echo "STAGE-3 OPT COMPLETE — starting dense audit"
python -u scripts/circulator_inverse_design/jesse_method_campaign/audit_bb_stage3_design.py --stage all

echo "STAGE-3 PIPELINE COMPLETE"
