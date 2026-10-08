#!/usr/bin/env python3
"""Phase 3 preconditioner screen at low cost (10/15 ppc)."""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "validation" / "highres_solver_campaign" / "phase3"
OUT.mkdir(parents=True, exist_ok=True)
PY = ROOT / ".conda_envs" / "fdfd_solver" / "bin" / "python"
SCRIPT = ROOT / "scripts" / "validation" / "highres_solver_campaign" / "fdfd_shifted_mg.py"
# fallback to plasmeep_mpi if venv missing py path for meep imports — actually prototype
# imports sixport which imports meep. Need plasmeep_mpi python with venv packages.
# Use plasmeep_mpi and PYTHONPATH to venv site-packages for pyamg.
PLAS = Path("/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python")
VENV_SITE = ROOT / ".conda_envs" / "fdfd_solver" / "lib" / "python3.13" / "site-packages"

JOBS = [
    # level1 cheap physics check
    dict(level="level1", ppc=10, method="direct", beta=0.0, timeout=600),
    dict(level="level1", ppc=10, method="shifted_mg", beta=0.5, timeout=600),
    dict(level="level1", ppc=10, method="shifted_ilu", beta=0.5, timeout=600),
    dict(level="level1", ppc=10, method="shifted_pyamg", beta=0.5, timeout=600),
    dict(level="level1", ppc=10, method="shifted_mg", beta=0.3, timeout=600),
    dict(level="level1", ppc=10, method="shifted_mg", beta=0.7, timeout=600),
    # full PMM screen only for methods that look promising — filled dynamically
    dict(level="level3", ppc=10, method="shifted_mg", beta=0.5, timeout=900),
    dict(level="level3", ppc=15, method="shifted_mg", beta=0.5, timeout=1200),
]


def run_job(job: dict) -> dict:
    tag = f"{job['level']}_ppc{job['ppc']}_{job['method']}_b{job['beta']}"
    jout = OUT / f"{tag}.json"
    log = OUT / f"{tag}.log"
    cmd = [
        str(PLAS), str(SCRIPT),
        "--level", job["level"],
        "--ppc", str(job["ppc"]),
        "--method", job["method"],
        "--beta", str(job["beta"]),
        "--timeout", str(job["timeout"]),
        "--maxiter", "200",
        "--json-out", str(jout),
    ]
    env = {"PYTHONPATH": f"{VENV_SITE}:{ROOT / 'outputs' / 'validation' / 'fdfd_feasibility'}:{ROOT / 'scripts' / 'validation'}:{ROOT / 'scripts'}"}
    import os
    full_env = os.environ.copy()
    full_env["PYTHONPATH"] = env["PYTHONPATH"] + ":" + full_env.get("PYTHONPATH", "")
    print(f"RUN {tag}", flush=True)
    t0 = time.time()
    with open(log, "w") as lf:
        p = subprocess.run(cmd, env=full_env, stdout=lf, stderr=subprocess.STDOUT, timeout=job["timeout"] + 120)
    wall = time.time() - t0
    summary = {"tag": tag, "returncode": p.returncode, "wall_s": wall, "json": str(jout)}
    if jout.exists():
        summary["result"] = json.loads(jout.read_text())
    else:
        summary["result"] = None
    print(f"  done rc={p.returncode} wall={wall:.1f}s resid={summary.get('result',{}) and summary['result'].get('true_residual')}", flush=True)
    return summary


def main():
    results = []
    for job in JOBS:
        try:
            results.append(run_job(job))
        except subprocess.TimeoutExpired:
            results.append({"tag": str(job), "error": "timeout"})
        except Exception as e:
            results.append({"tag": str(job), "error": str(e)})
        # early skip: if shifted_mg on level1 fails residual, still try others
        (OUT / "screen_summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps([{k: v for k, v in r.items() if k != "result"} | {"true_residual": (r.get("result") or {}).get("true_residual"), "converged": (r.get("result") or {}).get("converged"), "iters": (r.get("result") or {}).get("iterations"), "total_s": ((r.get("result") or {}).get("timings_s") or {}).get("total")} for r in results], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
