# Validation cache policy

## Safe to cache (invalidate on listed key changes)

| Artifact | Cache key fields | Reuse scope |
|---|---|---|
| Source normalization pickle | res, run_time, formulation, grid_offset, monitor_offset, horn_walls, coord_rotation, ports | Same geometry + source type |
| Numerical port mode JSON | res, frequency, horn_walls, grid_offset, port symmetry class | Launch/measure for matching horns |
| Lorentz direct reciprocity JSON | res, run_time, horn_walls, grid_offset | Diagnostic replay / plotting |
| Mode profile JSON | res, run_time, horn_walls, ports | Overlap analysis without re-sim |
| MPI timing summary | res, case type | Rank selection only |

## Must recalculate when changing

- `res` (grid spacing)
- `fs_a` / source frequency or `source_df`
- `horn_walls`, `grid_offset_cells`, `coord_rotation_deg`
- `device_mode` (horns_only vs full)
- plasma `rho`, `B`, or bulb geometry (full device)
- `run_time` below convergence-safe minimum

## Symmetry reuse

- Ports P2/P3 and P5/P6 (+60°/−60° pairs): numerical mode profiles may share one reference per |angle|.
- P2↔P3 reciprocity control: run once; do not duplicate for symmetric copies.

## Do not cache

- Full field HDF5 arrays (unless explicitly needed)
- Meep `Simulation` objects across processes
- Normalizations computed with `force=True` mid-sweep without updating cache key

## Overnight cheap defaults (after profiling)

- horns_only res32: use MPI ranks from `02_mpi_summary_res32.json`
- run_time from `03_runtime_convergence_res32.json` `recommended_run_time`
- Skip completed JSON phases in `overnight_campaign.py`
