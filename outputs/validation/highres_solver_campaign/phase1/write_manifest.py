#!/usr/bin/env python3
"""Document Phase1 environment status for the campaign."""
import json
from pathlib import Path
from datetime import datetime, timezone

out = {
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "production_env": "plasmeep_mpi — NOT MODIFIED",
    "campaign_venv": "/home/daq_user/PlasMEEP/.conda_envs/fdfd_solver",
    "installed_via_pip": {
        "numpy": "2.5.3",
        "scipy": "1.18.1",
        "matplotlib": "3.11.2",
        "h5py": "3.16.0",
        "pyamg": "5.3.0",
        "scikit-fem": "12.0.2",
        "meshio": "5.3.5",
        "gmsh": "4.15.2 (wheel present but unusable: missing libGLU.so.1)",
    },
    "petsc": {
        "status": "NOT_INSTALLED",
        "attempts": [
            "conda-forge / prefix.dev: blocked by proxy CONNECT 403",
            "pip petsc/petsc4py: no compiler (g++/clang missing system libs)",
            "portable LLVM 18: downloaded but needs libtinfo.so.5",
        ],
        "decision": "PETSc path abandoned within 2h budget; pivot to FEM + scipy SuperLU on reduced DOF",
    },
    "gpu": False,
    "hardware": {"ram_gib": 251, "nproc": 128},
}
Path("/home/daq_user/PlasMEEP/outputs/validation/highres_solver_campaign/phase1/env_manifest.json").write_text(
    json.dumps(out, indent=2) + "\n"
)
print(json.dumps(out, indent=2))
