#!/usr/bin/env bash
# Jesse method campaign launcher. Long FEM jobs run sequentially by default
# to avoid thrashing memory on FULL91-C.
set -euo pipefail
ROOT=/home/daq_user/PlasMEEP
OUT=$ROOT/outputs/circulator_inverse_design/jesse_method_campaign
LOG=$OUT/campaign.log
PY=/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python
export PYTHONPATH="\
$ROOT/scripts/circulator_inverse_design/jesse_method_campaign:\
$ROOT/scripts/circulator_inverse_design/broadband_method:\
$ROOT/scripts/circulator_inverse_design/multimode_optimization:\
$ROOT/scripts/circulator_inverse_design/multimode_ports:\
$ROOT/scripts/circulator_inverse_design:\
$ROOT/scripts/validation/fem_full_device_validation:\
$ROOT/scripts/validation/fem_gyrotropic_validation:\
$ROOT/scripts/validation/fem_meep_validation:\
$ROOT/scripts/validation:\
$ROOT/scripts:${PYTHONPATH:-}"
cd "$ROOT"
mkdir -p "$OUT"
exec >>"$LOG" 2>&1
echo "==== $(date -Is) campaign start ===="

run() { echo "---- $(date -Is) $*"; "$@"; }

run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/threefold_orbits.py
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/solver_inventory.py
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/jax_custom_vjp_tiny.py

run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/forward_mode_bench.py --grade C
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/eval_mode_baselines.py --grade C
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/rotation_source_check.py --grade C
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/spectral_recon.py --grade C --stage 200
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/spectral_recon.py --grade C --stage 100
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/spectral_recon.py --grade C --stage compare

run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/short_opt_mode_ab.py --grade C --maxiter 5 --which both
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/short_opt_c3_c6.py --grade C --maxiter 5 --which both
run "$PY" scripts/circulator_inverse_design/jesse_method_campaign/softmin_e2e_fd.py --grade C --orbits 1 --h 3e-4

echo "==== $(date -Is) campaign end ===="
