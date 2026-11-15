#!/usr/bin/env bash
# Authorized Stage-5 optimization + post dense audit.
# Resume-safe via bb_stage5_checkpoint.json. Does not launch Stage-6.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
echo "ROOT=$ROOT"
# shellcheck disable=SC1091
source /home/daq_user/miniconda3/etc/profile.d/conda.sh
conda activate plasmeep_mpi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export PLASMEEP_A2_FAST_KERNEL=1

OUT=outputs/circulator_inverse_design/jesse_method_campaign
mkdir -p "$OUT"

echo "========================================================================"
echo "STAGE-5 AUTHORIZED PIPELINE"
echo "cwd=$ROOT"
echo "commit=$(git rev-parse HEAD)"
echo "branch=$(git branch --show-current)"
echo "OMP=$OMP_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS MKL=$MKL_NUM_THREADS"
echo "PLASMEEP_A2_FAST_KERNEL=$PLASMEEP_A2_FAST_KERNEL"
echo "========================================================================"

python -u scripts/circulator_inverse_design/jesse_method_campaign/production_broadband_opt_stage5.py \
  --n-workers 12 --omp-threads-per-worker 1 --maxiter 30

echo "STAGE-5 OPT COMPLETE — starting dense audit (6.25 + 3.125 MHz refine)"
python -u scripts/circulator_inverse_design/jesse_method_campaign/audit_bb_stage5_dense.py --stage all

echo "STAGE-5 PIPELINE COMPLETE — Stage-6 NOT STARTED"
