#!/usr/bin/env python3
"""Measure process-parallel A2 scaling at 6/8(+12) workers.

Same 4-frequency Stage-2 probe as BROADBAND_PERFORMANCE_AUDIT:
  4.1, 5.0, 5.65, 5.75 GHz · grade C · Stage-2 q_best · omp=1/worker.

Does not change physics. Does not start Stage-3.
"""
from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
sys_paths = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]
import sys

sys.path[:0] = sys_paths

from a2_efficient_core import eval_bb_A2, rss_gib  # noqa: E402
from multimode_opt_core import orbits_and_masks  # noqa: E402

PROBE_GHZ = [4.1, 5.0, 5.65, 5.75]
OMP_PER_WORKER = 1


def set_parent_threads(n: int = 1):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(n)


def _rss_kib(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except Exception:
        return 0
    return 0


def _tree_rss_kib(root_pid: int) -> int:
    """Sum VmRSS of root + all descendants (KiB)."""
    total = _rss_kib(root_pid)
    try:
        children = []
        # breadth-first via /proc
        queue = [root_pid]
        seen = {root_pid}
        while queue:
            p = queue.pop()
            task = Path(f"/proc/{p}/task")
            if not task.exists():
                continue
            for tid_dir in task.iterdir():
                children_file = tid_dir / "children"
                if not children_file.exists():
                    continue
                for tok in children_file.read_text().split():
                    c = int(tok)
                    if c not in seen:
                        seen.add(c)
                        queue.append(c)
                        total += _rss_kib(c)
    except Exception:
        pass
    return total


def _cpu_jiffies(pid: int) -> int | None:
    try:
        with open(f"/proc/{pid}/stat") as f:
            parts = f.read().split()
        return int(parts[13]) + int(parts[14])
    except Exception:
        return None


def _tree_cpu_jiffies(root_pid: int) -> int:
    total = 0
    queue = [root_pid]
    seen = {root_pid}
    while queue:
        p = queue.pop()
        j = _cpu_jiffies(p)
        if j is not None:
            total += j
        try:
            for tid_dir in Path(f"/proc/{p}/task").iterdir():
                cf = tid_dir / "children"
                if not cf.exists():
                    continue
                for tok in cf.read_text().split():
                    c = int(tok)
                    if c not in seen:
                        seen.add(c)
                        queue.append(c)
        except Exception:
            pass
    return total


class ResourceMonitor:
    """Background peak RSS + average equivalent-core utilization for a process tree."""

    def __init__(self, pid: int, interval_s: float = 0.5):
        self.pid = pid
        self.interval_s = interval_s
        self.peak_rss_kib = 0
        self.samples = 0
        self._stop = threading.Event()
        self._thread = None
        self._clk = os.sysconf(os.sysconf_names["SC_CLK_TCK"])
        self._cpu_integral = 0.0  # cpu-seconds

    def start(self):
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=5)

    def _run(self):
        prev_j = _tree_cpu_jiffies(self.pid)
        prev_t = time.perf_counter()
        while not self._stop.wait(self.interval_s):
            rss = _tree_rss_kib(self.pid)
            if rss > self.peak_rss_kib:
                self.peak_rss_kib = rss
            j = _tree_cpu_jiffies(self.pid)
            t = time.perf_counter()
            if j >= prev_j and t > prev_t:
                self._cpu_integral += (j - prev_j) / float(self._clk)
                self.samples += 1
            prev_j, prev_t = j, t

    def summary(self, wall_s: float) -> dict:
        return {
            "peak_rss_GiB": self.peak_rss_kib / (1024.0 ** 2),
            "cpu_seconds_integrated": self._cpu_integral,
            "equiv_cores_avg": (self._cpu_integral / wall_s) if wall_s > 0 else float("nan"),
            "n_rss_samples": self.samples,
        }


def run_one(grade, q, orbits, freqs_hz, n_workers: int, ser_ref: dict | None):
    set_parent_threads(1)
    mon = ResourceMonitor(os.getpid(), interval_s=0.4)
    mon.start()
    t0 = time.perf_counter()
    rec = eval_bb_A2(
        grade,
        q,
        orbits,
        freqs_hz,
        n_workers=n_workers,
        omp_threads_per_worker=OMP_PER_WORKER,
        profile=False,
    )
    wall = time.perf_counter() - t0
    mon.stop()
    res = mon.summary(wall)
    out = {
        "n_workers": n_workers,
        "omp_threads_per_worker": OMP_PER_WORKER,
        "wall_s": wall,
        "J_BB": rec["J_BB"],
        "grad_norm": rec["grad_norm"],
        "mode": rec.get("mode"),
        "peak_rss_GiB": res["peak_rss_GiB"],
        "equiv_cores_avg": res["equiv_cores_avg"],
        "cpu_seconds_integrated": res["cpu_seconds_integrated"],
        "parent_rss_GiB_end": rss_gib(),
    }
    if ser_ref is not None:
        j_err = abs(rec["J_BB"] - ser_ref["J_BB"]) / max(abs(ser_ref["J_BB"]), 1e-12)
        g_err = float(
            np.linalg.norm(np.asarray(rec["grad_q"]) - np.asarray(ser_ref["grad_q"]))
            / max(np.linalg.norm(ser_ref["grad_q"]), 1e-12)
        )
        out["speedup_vs_serial"] = ser_ref["wall_s"] / wall if wall > 0 else None
        out["rel_err_J_BB"] = j_err
        out["rel_err_grad"] = g_err
        out["max_abs_err_J_per_freq"] = float(
            np.max(
                np.abs(
                    np.asarray(rec["J_per_freq"]) - np.asarray(ser_ref["J_per_freq"])
                )
            )
        )
        out["pass_equivalence"] = bool(
            j_err < 1e-10 and g_err < 1e-9 and out["max_abs_err_J_per_freq"] < 1e-10
        )
    else:
        out["speedup_vs_serial"] = 1.0
        out["rel_err_J_BB"] = 0.0
        out["rel_err_grad"] = 0.0
        out["pass_equivalence"] = True
        out["grad_q"] = np.asarray(rec["grad_q"], float)
        out["J_per_freq"] = rec["J_per_freq"]
    return out, rec


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    set_parent_threads(1)
    opt = json.loads((OUT / "opt_broadband_stage2_A2_C.json").read_text())
    q = np.asarray(opt["q_best"], float)
    orbits, _, _ = orbits_and_masks()
    freqs_hz = [f * 1e9 for f in PROBE_GHZ]
    grade = "C"

    rows = []
    print(f"=== serial baseline n_freq={len(PROBE_GHZ)} {PROBE_GHZ} ===", flush=True)
    ser_out, ser_rec = run_one(grade, q, orbits, freqs_hz, n_workers=1, ser_ref=None)
    ser_ref = {
        "wall_s": ser_out["wall_s"],
        "J_BB": ser_rec["J_BB"],
        "grad_q": ser_rec["grad_q"],
        "J_per_freq": ser_rec["J_per_freq"],
    }
    rows.append(ser_out)
    print(
        f"  wall={ser_out['wall_s']:.1f}s J={ser_out['J_BB']:.6f} "
        f"peakRAM={ser_out['peak_rss_GiB']:.2f}GiB cores≈{ser_out['equiv_cores_avg']:.2f}",
        flush=True,
    )

    # Always re-measure 4 for same-run baseline, then 6 and 8
    worker_list = [4, 6, 8]
    for nw in worker_list:
        print(f"=== n_workers={nw} omp/worker=1 ===", flush=True)
        out, _ = run_one(grade, q, orbits, freqs_hz, n_workers=nw, ser_ref=ser_ref)
        # drop large arrays if any
        out.pop("grad_q", None)
        out.pop("J_per_freq", None)
        rows.append(out)
        print(
            f"  wall={out['wall_s']:.1f}s speedup={out['speedup_vs_serial']:.2f}x "
            f"eq={out['pass_equivalence']} j_err={out['rel_err_J_BB']:.2e} "
            f"g_err={out['rel_err_grad']:.2e} peakRAM={out['peak_rss_GiB']:.2f}GiB "
            f"cores≈{out['equiv_cores_avg']:.2f}",
            flush=True,
        )

    # Optional 12 only if 8 is meaningfully faster than 4 on this probe
    row4 = next(r for r in rows if r["n_workers"] == 4)
    row8 = next(r for r in rows if r["n_workers"] == 8)
    ran_12 = False
    if (
        row8["pass_equivalence"]
        and row8["wall_s"] < 0.92 * row4["wall_s"]
        and row8["peak_rss_GiB"] < 40.0
    ):
        ran_12 = True
        print("=== n_workers=12 omp/worker=1 (optional: 8 beat 4) ===", flush=True)
        out, _ = run_one(grade, q, orbits, freqs_hz, n_workers=12, ser_ref=ser_ref)
        out.pop("grad_q", None)
        out.pop("J_per_freq", None)
        rows.append(out)
        print(
            f"  wall={out['wall_s']:.1f}s speedup={out['speedup_vs_serial']:.2f}x "
            f"eq={out['pass_equivalence']} peakRAM={out['peak_rss_GiB']:.2f}GiB "
            f"cores≈{out['equiv_cores_avg']:.2f}",
            flush=True,
        )
    else:
        print(
            "SKIP n_workers=12: 8-worker run did not show continued useful scaling "
            f"vs 4-worker (8={row8['wall_s']:.1f}s, 4={row4['wall_s']:.1f}s) — "
            "expected on a 4-frequency probe (at most 4 concurrent tasks).",
            flush=True,
        )

    # Fastest measured among equivalent + memory-safe
    candidates = [
        r
        for r in rows
        if r.get("pass_equivalence", False) and r["peak_rss_GiB"] < 40.0
    ]
    best = min(candidates, key=lambda r: r["wall_s"])

    # Stage-3 estimate: historical 121 s/freq × 30 × 83 evals
    s_pf = 281238.1991 / 83 / 28
    n_freq_s3, n_evals = 30, 83
    serial_h = n_evals * n_freq_s3 * s_pf / 3600.0
    # On 4-freq probe, speedup saturates at n_workers>=4; for Stage-3 with 30 freqs,
    # apply measured speedup only for configs that were faster than fewer workers.
    # Conservatively use measured speedup of the recommended config (cannot invent
    # extra parallel efficiency beyond what 4-freq can demonstrate).
    measured_sp = best["speedup_vs_serial"] or 1.0
    # Batch projection for 30 freqs with efficiency from measured 4w/4f
    eff = (row4["speedup_vs_serial"] or 2.93) / 4.0
    proj = {}
    for nw in [1, 4, 6, 8]:
        batches = (n_freq_s3 + nw - 1) // nw
        ideal = n_freq_s3 / batches
        if nw == 1:
            sp = 1.0
            kind = "serial"
        elif nw <= 4:
            # use measured if available
            m = next((r for r in rows if r["n_workers"] == nw), None)
            sp = (m["speedup_vs_serial"] if m else ideal * eff) or 1.0
            kind = "measured_on_4freq_probe"
        else:
            # cannot measure >4 parallel tasks on 4-freq probe; project only
            sp = min(ideal * eff, ideal * 0.85)
            kind = "projected_from_4freq_efficiency"
        proj[str(nw)] = {
            "speedup": sp,
            "wall_hours": serial_h / sp,
            "kind": kind,
            "peak_RAM_GiB_est": 3.2 * nw,
        }

    report = {
        "date": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "probe_GHz": PROBE_GHZ,
        "grade": grade,
        "design": "opt_broadband_stage2_A2_C.json q_best",
        "omp_threads_per_worker": OMP_PER_WORKER,
        "ran_12_workers": ran_12,
        "rows": [{k: v for k, v in r.items() if k not in ("grad_q",)} for r in rows],
        "fastest_measured_safe": {
            "n_workers": best["n_workers"],
            "wall_s": best["wall_s"],
            "speedup_vs_serial": best["speedup_vs_serial"],
            "peak_rss_GiB": best["peak_rss_GiB"],
            "equiv_cores_avg": best["equiv_cores_avg"],
            "pass_equivalence": best["pass_equivalence"],
            "note": (
                "On a 4-frequency probe, n_workers>4 cannot exceed ~4 concurrent "
                "solves; fastest measured is typically n_workers=4."
            ),
        },
        "stage3_estimate": {
            "n_freq": n_freq_s3,
            "n_evals_assumed": n_evals,
            "s_per_freq_historical": s_pf,
            "wall_hours_serial": serial_h,
            "by_n_workers": proj,
            "recommended_launch": {
                "n_workers": int(best["n_workers"]),
                "omp_threads_per_worker": 1,
                "env": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1",
                "basis": "fastest measured numerically-equivalent config on 4-freq probe",
            },
            "status": "NOT STARTED / AWAITING AUTHORIZATION",
        },
    }

    out_json = OUT / "BROADBAND_PARALLEL_SCALING_6_8.json"
    out_json.write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", out_json, flush=True)
    print(
        "RECOMMENDED n_workers=",
        best["n_workers"],
        "speedup=",
        best["speedup_vs_serial"],
        "Stage3_h≈",
        serial_h / (best["speedup_vs_serial"] or 1.0),
        flush=True,
    )


if __name__ == "__main__":
    main()
