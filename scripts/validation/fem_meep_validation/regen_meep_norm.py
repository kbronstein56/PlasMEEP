#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "scripts/validation"), str(ROOT / "scripts")]
from sixport_common import set_geometry_context, ensure_normalizations
import sixport_common as sc
from physical_units import meep_resolution_from_points_per_cm

def rank():
    try:
        from mpi4py import MPI
        return int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        return 0

def main():
    res = int(meep_resolution_from_points_per_cm(25.0, a_m=sc.a))
    set_geometry_context(res=res, horn_walls="prism", grid_offset_cells=(0,0),
                         monitor_offset_cells=(0,0), coord_rotation_deg=0.0)
    out = ROOT / "outputs/validation/fem_meep_validation/phase2/norm_p0_fresh.pkl"
    if rank()==0:
        print("Regenerating norm force=True ->", out, flush=True)
    cache = ensure_normalizations(
        res=res, run_time=20.0, force=True, ports=[0],
        cache_path=str(out), verbose=True, formulation="num_mode_guide_normal",
    )
    if rank()==0:
        fd = cache["incident_flux_data_by_port"][0]
        E = np.asarray(fd.E); H = np.asarray(fd.H)
        print("P_inc", cache["incident_power_by_port"][0],
              "E_nnz", np.count_nonzero(E), "H_nnz", np.count_nonzero(H),
              "E_max", float(np.max(np.abs(E))), "H_max", float(np.max(np.abs(H))), flush=True)
        print("DONE", flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
