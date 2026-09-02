"""
Shared MPI launch helpers for PlasMEEP validation harnesses.

Validated configuration (horns_only / full-device):
  OMP_NUM_THREADS=1
  FI_PROVIDER=tcp
  MPICH_CH4_NETMOD=ofi
  mpirun -np 32
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

DEFAULT_MPIRUN = os.environ.get(
    "PLASMEEP_MPIRUN",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/mpirun",
)
DEFAULT_PYTHON = os.environ.get(
    "PLASMEEP_PYTHON",
    "/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python",
)
DEFAULT_RANKS = int(os.environ.get("PLASMEEP_MPI_RANKS", "32"))


def mpi_env(*, extra: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["FI_PROVIDER"] = "tcp"
    env["MPICH_CH4_NETMOD"] = "ofi"
    env.setdefault("UCX_TLS", "tcp,self")
    if extra:
        env.update(extra)
    return env


def mpi_info() -> Dict[str, Any]:
    """Current process MPI metadata for JSON logging."""
    size = 1
    rank = 0
    try:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        size = int(comm.Get_size())
        rank = int(comm.Get_rank())
    except Exception:
        for key, default in (
            ("OMPI_COMM_WORLD_SIZE", "1"),
            ("PMI_SIZE", "1"),
            ("OMPI_COMM_WORLD_RANK", "0"),
            ("PMI_RANK", "0"),
        ):
            if key.endswith("SIZE") and key in os.environ:
                size = int(os.environ[key])
            if key.endswith("RANK") and key in os.environ:
                rank = int(os.environ[key])
    return {
        "ranks": size,
        "rank": rank,
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
        "fi_provider": os.environ.get("FI_PROVIDER", ""),
        "mpich_ch4_netmod": os.environ.get("MPICH_CH4_NETMOD", ""),
        "mpirun": DEFAULT_MPIRUN,
        "python": sys.executable,
    }


def verify_mpi_ranks(
    *,
    ranks: int = DEFAULT_RANKS,
    mpirun: str = DEFAULT_MPIRUN,
    python: str = DEFAULT_PYTHON,
    cwd: Optional[Path | str] = None,
    timeout_s: float = 30.0,
) -> None:
    """Raise if mpirun does not launch the requested number of workers."""
    probe = (
        "from mpi4py import MPI; "
        "import sys; "
        "n=int(MPI.COMM_WORLD.Get_size()); "
        "print(n); "
        "sys.exit(0 if n=="
        f"{ranks}"
        " else 1)"
    )
    cmd = [mpirun, "-np", str(ranks), python, "-c", probe]
    proc = subprocess.run(
        cmd,
        cwd=str(cwd or os.getcwd()),
        env=mpi_env(),
        capture_output=True,
        text=True,
        timeout=timeout_s,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"MPI rank probe failed (expected {ranks}): "
            f"exit={proc.returncode} stdout={proc.stdout!r} stderr={proc.stderr!r}"
        )
    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    if not lines or not any(ln == str(ranks) for ln in lines):
        raise RuntimeError(
            f"MPI rank probe did not report {ranks} workers: stdout={proc.stdout!r}"
        )


def mpirun_cmd(
    script_args: Sequence[str],
    *,
    ranks: int = DEFAULT_RANKS,
    mpirun: str = DEFAULT_MPIRUN,
    python: str = DEFAULT_PYTHON,
) -> List[str]:
    """Wrap a validation script invocation with mpirun."""
    return [mpirun, "-np", str(ranks), python, *script_args]


def run_mpi_script(
    script_args: Sequence[str],
    *,
    ranks: int = DEFAULT_RANKS,
    mpirun: str = DEFAULT_MPIRUN,
    python: str = DEFAULT_PYTHON,
    cwd: Optional[Path | str] = None,
    verify_ranks: bool = True,
    log_path: Optional[Path | str] = None,
) -> subprocess.CompletedProcess:
    """Run a validation script under mpirun with validated env."""
    if verify_ranks:
        verify_mpi_ranks(ranks=ranks, mpirun=mpirun, python=python, cwd=cwd)
    cmd = mpirun_cmd(script_args, ranks=ranks, mpirun=mpirun, python=python)
    env = mpi_env()
    if log_path:
        with open(log_path, "w", encoding="utf-8") as log:
            log.write(f"# {' '.join(cmd)}\n")
            return subprocess.run(
                cmd,
                cwd=str(cwd or os.getcwd()),
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
    return subprocess.run(
        cmd,
        cwd=str(cwd or os.getcwd()),
        env=env,
        check=False,
    )
