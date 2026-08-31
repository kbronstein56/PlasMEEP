
import json, os, sys, time
sys.path.insert(0, os.environ["VAL"])
sys.path.insert(0, os.environ["SCRIPTS"])
from plasmeep.ports.lorentz_probe import run_hz_transfer
from sixport_common import (
    build_circulator_device, default_uniform_rho, fs_a, source_df,
    horn_for_port, monitor_center_for_port, set_geometry_context,
)

res = int(os.environ["BENCH_RES"])
run_time = float(os.environ["BENCH_RT"])
port_a = int(os.environ["PORT_A"])
port_b = int(os.environ["PORT_B"])

set_geometry_context(res=res, horn_walls="prism")
rho = default_uniform_rho()
B = [0.0, 0.0, 0.0]
t0 = time.perf_counter()
_pmm, dev, _ = build_circulator_device(rho, B, res=res, device_mode="horns_only")
t_build = time.perf_counter() - t0

src = horn_for_port(port_a, res)["source_center"]
mon = monitor_center_for_port(port_b, res)
t1 = time.perf_counter()
hz = run_hz_transfer(
    dev, res=res, source_xy=src, monitor_xy=mon,
    frequency=fs_a, fwidth=source_df, run_time=run_time,
)
t_xfer = time.perf_counter() - t1

out = {
    "ranks": int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1"))),
    "t_build_s": t_build,
    "t_transfer_s": t_xfer,
    "t_total_s": t_build + t_xfer,
    "hz_abs": abs(hz),
}
rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
if rank == 0:
    with open(os.environ["OUT"], "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
