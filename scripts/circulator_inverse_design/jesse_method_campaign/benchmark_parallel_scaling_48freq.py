#!/usr/bin/env python3
"""High-worker frequency-parallel scaling on a 48-frequency 4–6 GHz probe.

Measures wall time / speedup / RAM / numerical agreement for
n_workers ∈ {12, 24, 32, 48} with OMP/OPENBLAS/MKL = 1 per worker.

Uses Stage-4 q_best and the validated eval_bb_A2 path. Does NOT change
physics, mesh, objective, derivatives, or tolerances. Does NOT start Stage-5.
"""
from __future__ import annotations

import json
import os
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs" / "circulator_inverse_design" / "jesse_method_campaign"
sys.path[:0] = [
    str(ROOT / "scripts" / "circulator_inverse_design" / "jesse_method_campaign"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "broadband_method"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_optimization"),
    str(ROOT / "scripts" / "circulator_inverse_design" / "multimode_ports"),
    str(ROOT / "scripts" / "circulator_inverse_design"),
    str(ROOT / "scripts" / "validation" / "fem_full_device_validation"),
    str(ROOT / "scripts" / "validation"),
    str(ROOT / "scripts"),
]

from a2_efficient_core import eval_bb_A2, rss_gib  # noqa: E402
from multimode_opt_core import orbits_and_masks  # noqa: E402

# 48 independent frequencies spanning 4–6 GHz so all tested worker counts can
# remain busy (ceil(48/48)=1 batch at the high end).
PROBE_GHZ = [
    4.000, 4.050, 4.100, 4.150, 4.200, 4.250, 4.300, 4.350,
    4.400, 4.450, 4.500, 4.550, 4.600, 4.650, 4.700, 4.725,
    4.750, 4.775, 4.800, 4.850, 4.900, 4.950, 5.000, 5.050,
    5.08125, 5.100, 5.150, 5.200, 5.250, 5.300, 5.350, 5.400,
    5.450, 5.500, 5.550, 5.5875, 5.650, 5.700, 5.750, 5.800,
    5.8125, 5.850, 5.900, 5.950, 6.000,
]
# Pad to exactly 48 with evenly spaced fillers if needed
_extra = [4.025, 4.075, 5.025]
PROBE_GHZ = sorted(set(PROBE_GHZ + _extra))
assert len(PROBE_GHZ) == 48, len(PROBE_GHZ)
assert PROBE_GHZ == sorted(PROBE_GHZ)

WORKER_PLAN = [12, 24, 32, 48]
OMP_PER_WORKER = 1
GRADE = "C"
TAU = 0.05
ALPHA = 0.25


def log(msg: str):
    print(msg, flush=True)


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


def _cpu_jiffies(pid: int) -> int | None:
    try:
        with open(f"/proc/{pid}/stat") as f:
            parts = f.read().split()
        return int(parts[13]) + int(parts[14])
    except Exception:
        return None


def _tree_walk(root_pid: int):
    queue = [root_pid]
    seen = {root_pid}
    while queue:
        p = queue.pop()
        yield p
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


def _tree_rss_kib(root_pid: int) -> int:
    return sum(_rss_kib(p) for p in _tree_walk(root_pid))


def _tree_cpu_jiffies(root_pid: int) -> int:
    total = 0
    for p in _tree_walk(root_pid):
        j = _cpu_jiffies(p)
        if j is not None:
            total += j
    return total


class ResourceMonitor:
    def __init__(self, pid: int, interval_s: float = 0.4):
        self.pid = pid
        self.interval_s = interval_s
        self.peak_rss_kib = 0
        self.samples = 0
        self._stop = threading.Event()
        self._thread = None
        self._clk = os.sysconf(os.sysconf_names["SC_CLK_TCK"])
        self._cpu_integral = 0.0

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


def run_one(grade, q, orbits, freqs_hz, n_workers: int, ref: dict | None):
    set_parent_threads(1)
    mon = ResourceMonitor(os.getpid(), interval_s=0.4)
    mon.start()
    t0 = time.perf_counter()
    rec = eval_bb_A2(
        grade,
        q,
        orbits,
        freqs_hz,
        tau=TAU,
        alpha=ALPHA,
        n_workers=n_workers,
        omp_threads_per_worker=OMP_PER_WORKER,
        profile=True,
    )
    wall = time.perf_counter() - t0
    mon.stop()
    res = mon.summary(wall)
    timings = rec.get("per_freq_timings") or []
    s_pf = (
        float(np.mean([t.get("total", float("nan")) for t in timings]))
        if timings
        else None
    )
    n_freq = len(freqs_hz)
    n_batches = int(np.ceil(n_freq / max(n_workers, 1)))
    serial_est = (s_pf * n_freq) if s_pf is not None else None
    out = {
        "n_workers": n_workers,
        "omp_threads_per_worker": OMP_PER_WORKER,
        "n_freq": n_freq,
        "n_batches": n_batches,
        "wall_s": wall,
        "J_BB": rec["J_BB"],
        "grad_norm": rec["grad_norm"],
        "mode": rec.get("mode"),
        "peak_rss_GiB": res["peak_rss_GiB"],
        "equiv_cores_avg": res["equiv_cores_avg"],
        "cpu_seconds_integrated": res["cpu_seconds_integrated"],
        "parent_rss_GiB_end": rss_gib(),
        "s_per_freq_mean": s_pf,
        "serial_wall_est_s": serial_est,
        "speedup_vs_serial_est": (serial_est / wall) if (serial_est and wall > 0) else None,
        "parallel_efficiency_vs_serial_est": (
            (serial_est / wall) / n_workers
            if (serial_est and wall > 0 and n_workers > 0)
            else None
        ),
    }
    if ref is not None:
        j_err = abs(rec["J_BB"] - ref["J_BB"]) / max(abs(ref["J_BB"]), 1e-12)
        g_err = float(
            np.linalg.norm(np.asarray(rec["grad_q"]) - np.asarray(ref["grad_q"]))
            / max(np.linalg.norm(ref["grad_q"]), 1e-12)
        )
        out["speedup_vs_12workers"] = ref["wall_s"] / wall if wall > 0 else None
        out["rel_err_J_BB"] = j_err
        out["rel_err_grad"] = g_err
        out["max_abs_err_J_per_freq"] = float(
            np.max(np.abs(np.asarray(rec["J_per_freq"]) - np.asarray(ref["J_per_freq"])))
        )
        out["pass_equivalence"] = bool(
            j_err < 1e-10 and g_err < 1e-9 and out["max_abs_err_J_per_freq"] < 1e-10
        )
    else:
        out["speedup_vs_12workers"] = 1.0
        out["rel_err_J_BB"] = 0.0
        out["rel_err_grad"] = 0.0
        out["max_abs_err_J_per_freq"] = 0.0
        out["pass_equivalence"] = True
    return out, rec


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    set_parent_threads(1)

    n_cpu = os.cpu_count() or 0
    log("=" * 72)
    log("48-FREQ HIGH-WORKER PARALLEL SCALING (Stage-5 planning; no optimization)")
    log(f"host_cpus={n_cpu}")
    log(f"probe_GHz ({len(PROBE_GHZ)}): {PROBE_GHZ}")
    log(f"worker_plan: {WORKER_PLAN}")
    log(f"OMP/OPENBLAS/MKL per worker: {OMP_PER_WORKER}")
    log(f"grade={GRADE} tau={TAU} alpha={ALPHA}")
    log("design: opt_broadband_stage4_A2_C.json q_best")
    log("=" * 72)

    opt = json.loads((OUT / "opt_broadband_stage4_A2_C.json").read_text())
    q = np.asarray(opt["q_best"], float)
    orbits, _, _ = orbits_and_masks()
    freqs_hz = [f * 1e9 for f in PROBE_GHZ]

    rows = []
    ref = None
    for nw in WORKER_PLAN:
        log(f"START {nw} workers (batches={int(np.ceil(len(PROBE_GHZ)/nw))})")
        out, rec = run_one(GRADE, q, orbits, freqs_hz, n_workers=nw, ref=ref)
        if nw == 12:
            ref = {
                "wall_s": out["wall_s"],
                "J_BB": rec["J_BB"],
                "grad_q": rec["grad_q"],
                "J_per_freq": rec["J_per_freq"],
            }
        rows.append(out)
        log(
            f"DONE {nw}w: wall={out['wall_s']:.1f}s "
            f"batches={out['n_batches']} "
            f"sp_vs_12={out.get('speedup_vs_12workers')} "
            f"sp_vs_serial_est={out.get('speedup_vs_serial_est')} "
            f"cores≈{out['equiv_cores_avg']:.2f} "
            f"peakRAM={out['peak_rss_GiB']:.2f}GiB "
            f"eq={out['pass_equivalence']} "
            f"j_err={out['rel_err_J_BB']:.2e} g_err={out['rel_err_grad']:.2e}"
        )

    safe = [
        r
        for r in rows
        if r.get("pass_equivalence") and r["peak_rss_GiB"] < 200.0
    ]
    if not safe:
        safe = [r for r in rows if r.get("pass_equivalence")]
    best = min(safe, key=lambda r: r["wall_s"]) if safe else min(rows, key=lambda r: r["wall_s"])

    # Scaling quality vs 12-worker baseline
    row12 = next(r for r in rows if r["n_workers"] == 12)
    scaling_notes = []
    for r in rows:
        if r["n_workers"] == 12:
            continue
        ideal = r["n_workers"] / 12.0
        actual = r.get("speedup_vs_12workers") or float("nan")
        eff = actual / ideal if ideal > 0 else float("nan")
        scaling_notes.append(
            {
                "n_workers": r["n_workers"],
                "ideal_speedup_vs_12": ideal,
                "measured_speedup_vs_12": actual,
                "efficiency_vs_ideal": eff,
                "wall_s": r["wall_s"],
                "peak_rss_GiB": r["peak_rss_GiB"],
                "equiv_cores_avg": r["equiv_cores_avg"],
            }
        )

    report = {
        "date": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "host_cpus": n_cpu,
        "probe_GHz": PROBE_GHZ,
        "n_freq_probe": len(PROBE_GHZ),
        "grade": GRADE,
        "tau": TAU,
        "alpha": ALPHA,
        "design": "opt_broadband_stage4_A2_C.json q_best",
        "omp_threads_per_worker": OMP_PER_WORKER,
        "env": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1",
        "worker_plan": WORKER_PLAN,
        "rows": rows,
        "scaling_vs_12workers": scaling_notes,
        "fastest_measured_safe": {
            "n_workers": best["n_workers"],
            "wall_s": best["wall_s"],
            "n_batches": best["n_batches"],
            "speedup_vs_12workers": best.get("speedup_vs_12workers"),
            "speedup_vs_serial_est": best.get("speedup_vs_serial_est"),
            "peak_rss_GiB": best["peak_rss_GiB"],
            "equiv_cores_avg": best["equiv_cores_avg"],
            "pass_equivalence": best["pass_equivalence"],
            "s_per_freq_mean": best.get("s_per_freq_mean"),
        },
        "baseline_12workers": {
            "wall_s": row12["wall_s"],
            "peak_rss_GiB": row12["peak_rss_GiB"],
            "equiv_cores_avg": row12["equiv_cores_avg"],
            "s_per_freq_mean": row12.get("s_per_freq_mean"),
        },
        "recommendation": {
            "n_workers": int(best["n_workers"]),
            "omp_threads_per_worker": 1,
            "basis": "fastest numerically-equivalent config on 48-freq probe with peak RAM < 200 GiB",
            "note": "Stage-5 NOT STARTED / AWAITING AUTHORIZATION",
        },
        "status": "Stage-5 NOT STARTED / AWAITING AUTHORIZATION",
    }
    out_json = OUT / "BROADBAND_PARALLEL_SCALING_48FREQ.json"
    out_json.write_text(json.dumps(report, indent=2) + "\n")
    log(f"wrote {out_json}")
    log(
        f"RECOMMENDED n_workers={best['n_workers']} "
        f"wall={best['wall_s']:.1f}s "
        f"sp_vs_12={best.get('speedup_vs_12workers')} "
        f"peakRAM={best['peak_rss_GiB']:.1f}GiB"
    )
    log("STAGE-5 STATUS: NOT STARTED / AWAITING AUTHORIZATION")


if __name__ == "__main__":
    main()
