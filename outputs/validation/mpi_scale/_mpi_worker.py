
import json, os, sys, time
sys.path.insert(0, os.environ["SIXPORT_VAL"])
sys.path.insert(0, os.environ["SIXPORT_SCRIPTS"])
import meep as mp
from sixport_common import (
    ensure_normalizations,
    default_uniform_rho,
    build_circulator_device,
    simulate_circulator,
    physical_resolution_report,
)
from port_formulations import get_formulation

res = int(os.environ["BENCH_RES"])
run_time = float(os.environ["BENCH_RUNTIME"])
port = int(os.environ["BENCH_PORT"])
formul = os.environ.get("BENCH_FORMULATION", "te1_hz_line")

t0 = time.time()
cache = ensure_normalizations(
    res,
    run_time=run_time,
    force=True,
    ports=[port],
    verbose=False,
    formulation=formul,
)
t_norm = time.time() - t0

rho = default_uniform_rho()
t1 = time.time()
# Single-port device excitation via the same path as reciprocity studies
result = simulate_circulator(
    rho,
    [0.0, 0.0, 0.0],
    res=res,
    run_time=run_time,
    verbose=False,
    incident_cache=cache,
    ports=[port],
    formulation=formul,
)
t_dev = time.time() - t1
flux = float(abs(result["power_matrix"][0, 0]))

out = {
    "rank_env": int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1"))),
    "t_norm_s": t_norm,
    "t_device_s": t_dev,
    "flux": flux,
    "res": res,
    "run_time": run_time,
    "port": port,
    "formulation": formul,
    "physical_resolution": physical_resolution_report(res),
}
path = os.environ["BENCH_OUT"]
rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
if rank == 0:
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out))
