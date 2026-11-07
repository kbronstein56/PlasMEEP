#!/usr/bin/env python3
"""Stage-3-representative frequency-parallel scaling (8–12 freqs).

Uses a subset of the ACTUAL Stage-3 30-frequency set so that 6/8 workers
can run concurrently. Does NOT start Stage-3 optimization.
Does NOT change physics / mesh / objective / solver.
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

# Full Stage-3 set (30) for reference / membership checks
F_STAGE3_GHZ = [
    4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.8, 4.85, 4.9,
    5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.6, 5.625, 5.65, 5.7,
    5.75, 5.8, 5.825, 5.85, 5.9, 5.95, 6.0,
]

# 12-freq representative subset: band coverage + required problem freqs
PROBE_GHZ = [
    4.1, 4.3, 4.5, 4.7, 5.0, 5.2,
    5.4, 5.6, 5.65, 5.75, 5.825, 6.0,
]
assert len(PROBE_GHZ) == 12
assert all(f in F_STAGE3_GHZ for f in PROBE_GHZ)
for req in (4.1, 4.3, 5.65, 5.75, 5.825):
    assert req in PROBE_GHZ

OMP_PER_WORKER = 1
GRADE = "C"


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
        profile=True,
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
        "per_freq_timings": rec.get("per_freq_timings"),
        "s_per_freq_mean": (
            float(np.mean([t.get("total", float("nan")) for t in rec.get("per_freq_timings") or []]))
            if rec.get("per_freq_timings")
            else None
        ),
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
            np.max(np.abs(np.asarray(rec["J_per_freq"]) - np.asarray(ser_ref["J_per_freq"])))
        )
        out["pass_equivalence"] = bool(
            j_err < 1e-10 and g_err < 1e-9 and out["max_abs_err_J_per_freq"] < 1e-10
        )
    else:
        out["speedup_vs_serial"] = 1.0
        out["rel_err_J_BB"] = 0.0
        out["rel_err_grad"] = 0.0
        out["pass_equivalence"] = True
    return out, rec


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    set_parent_threads(1)

    worker_plan = [1, 4, 6, 8]  # 12 appended later only if useful
    n_jobs = len(PROBE_GHZ)

    log("=" * 72)
    log("STAGE-3 SUBSET PARALLEL SCALING BENCHMARK (no Stage-3 optimization)")
    log(f"Exact frequency subset ({n_jobs} GHz): {PROBE_GHZ}")
    log(f"All freqs ⊆ Stage-3 30-set: {all(f in F_STAGE3_GHZ for f in PROBE_GHZ)}")
    log(f"Required problem freqs present: 4.1, 4.3, 5.65, 5.75, 5.825")
    log(f"Worker counts planned: {worker_plan} (+12 only if 8 still scales)")
    log(f"Expected independent jobs per BB eval: {n_jobs}")
    log(f"OMP/OPENBLAS/MKL per worker: {OMP_PER_WORKER}")
    log(f"Grade: {GRADE}; design: Stage-2 q_best")
    log("=" * 72)

    opt = json.loads((OUT / "opt_broadband_stage2_A2_C.json").read_text())
    q = np.asarray(opt["q_best"], float)
    orbits, _, _ = orbits_and_masks()
    freqs_hz = [f * 1e9 for f in PROBE_GHZ]

    rows = []
    ser_ref = None

    for nw in worker_plan:
        label = "serial" if nw == 1 else f"{nw} workers"
        log(f"START {label}")
        out, rec = run_one(
            GRADE, q, orbits, freqs_hz, n_workers=nw, ser_ref=ser_ref
        )
        if nw == 1:
            ser_ref = {
                "wall_s": out["wall_s"],
                "J_BB": rec["J_BB"],
                "grad_q": rec["grad_q"],
                "J_per_freq": rec["J_per_freq"],
            }
        rows.append(out)
        log(
            f"DONE {label}: wall={out['wall_s']:.1f}s "
            f"speedup={out['speedup_vs_serial']:.3f}x "
            f"cores≈{out['equiv_cores_avg']:.2f} "
            f"peakRAM={out['peak_rss_GiB']:.2f}GiB "
            f"eq={out['pass_equivalence']} "
            f"j_err={out['rel_err_J_BB']:.2e} g_err={out['rel_err_grad']:.2e}"
        )

    row4 = next(r for r in rows if r["n_workers"] == 4)
    row8 = next(r for r in rows if r["n_workers"] == 8)
    ran_12 = False
    # Useful scaling: 8 meaningfully beats 4, and we have ≥12 independent freqs
    if (
        n_jobs >= 12
        and row8["pass_equivalence"]
        and row8["wall_s"] < 0.90 * row4["wall_s"]
        and row8["peak_rss_GiB"] < 80.0
    ):
        ran_12 = True
        log("START 12 workers")
        out, _ = run_one(GRADE, q, orbits, freqs_hz, n_workers=12, ser_ref=ser_ref)
        rows.append(out)
        log(
            f"DONE 12 workers: wall={out['wall_s']:.1f}s "
            f"speedup={out['speedup_vs_serial']:.3f}x "
            f"cores≈{out['equiv_cores_avg']:.2f} "
            f"peakRAM={out['peak_rss_GiB']:.2f}GiB "
            f"eq={out['pass_equivalence']}"
        )
    else:
        log(
            "SKIP 12 workers: either n_freq<12 or 8-worker run did not show "
            f"continued useful scaling vs 4 "
            f"(8={row8['wall_s']:.1f}s, 4={row4['wall_s']:.1f}s, "
            f"need 8 < 0.90*4)."
        )

    candidates = [
        r for r in rows if r.get("pass_equivalence") and r["peak_rss_GiB"] < 80.0
    ]
    best = min(candidates, key=lambda r: r["wall_s"])

    # Stage-3 full estimate from historical s/freq and measured speedup
    s_pf = 281238.1991 / 83 / 28
    n_freq_s3, n_evals = 30, 83
    serial_h = n_evals * n_freq_s3 * s_pf / 3600.0
    measured_sp = best["speedup_vs_serial"] or 1.0

    # Also project using probe wall / n_freq * 30 / measured pattern
    by_nw = {}
    for r in rows:
        nw = r["n_workers"]
        sp = r["speedup_vs_serial"] or 1.0
        # For Stage-3 (30 freqs): if probe used fewer freqs than workers,
        # speedup underestimates; if probe freqs >= workers, measured sp applies
        # more directly. Here n_jobs=12 >= 8, so 1/4/6/8 are well-tested.
        # Scale Stage-3 wall by measured speedup of that config.
        by_nw[str(nw)] = {
            "probe_wall_s": r["wall_s"],
            "probe_speedup": sp,
            "probe_peak_rss_GiB": r["peak_rss_GiB"],
            "probe_equiv_cores": r["equiv_cores_avg"],
            "stage3_wall_hours_from_measured_speedup": serial_h / sp,
            "pass_equivalence": r["pass_equivalence"],
            "kind": "measured_on_12freq_stage3_subset",
        }

    report = {
        "date": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "probe_GHz": PROBE_GHZ,
        "n_freq_probe": n_jobs,
        "stage3_full_n_freq": 30,
        "grade": GRADE,
        "design": "opt_broadband_stage2_A2_C.json q_best",
        "omp_threads_per_worker": OMP_PER_WORKER,
        "ran_12_workers": ran_12,
        "rows": rows,
        "fastest_measured_safe": {
            "n_workers": best["n_workers"],
            "wall_s": best["wall_s"],
            "speedup_vs_serial": measured_sp,
            "peak_rss_GiB": best["peak_rss_GiB"],
            "equiv_cores_avg": best["equiv_cores_avg"],
            "pass_equivalence": best["pass_equivalence"],
        },
        "stage3_estimate": {
            "n_freq": n_freq_s3,
            "n_evals_assumed": n_evals,
            "s_per_freq_historical": s_pf,
            "wall_hours_serial": serial_h,
            "by_n_workers": by_nw,
            "recommended_launch": {
                "n_workers": int(best["n_workers"]),
                "omp_threads_per_worker": 1,
                "env": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1",
                "cli": (
                    f"--n-workers {int(best['n_workers'])} --omp-threads-per-worker 1"
                ),
                "basis": (
                    "fastest numerically-equivalent config on 12-freq Stage-3 subset"
                ),
                "expected_wall_hours": serial_h / measured_sp,
                "expected_peak_RAM_GiB": best["peak_rss_GiB"]
                * (min(best["n_workers"], 30) / min(best["n_workers"], n_jobs)),
            },
            "status": "NOT STARTED / AWAITING AUTHORIZATION",
        },
    }
    # Peak RAM for Stage-3: if n_workers busy fully, scale from probe peak
    # when probe already saturated the worker count (n_jobs >= n_workers).
    rec_nw = int(best["n_workers"])
    if n_jobs >= rec_nw:
        report["stage3_estimate"]["recommended_launch"]["expected_peak_RAM_GiB"] = (
            best["peak_rss_GiB"]
        )

    out_json = OUT / "BROADBAND_PARALLEL_SCALING_STAGE3_SUBSET.json"
    out_json.write_text(json.dumps(report, indent=2) + "\n")
    log(f"wrote {out_json}")
    log(
        f"RECOMMENDED n_workers={best['n_workers']} "
        f"speedup={measured_sp:.3f}x "
        f"Stage3≈{serial_h / measured_sp:.1f}h "
        f"peakRAM≈{best['peak_rss_GiB']:.1f}GiB"
    )
    log("STAGE-3 STATUS: NOT STARTED / AWAITING AUTHORIZATION")


if __name__ == "__main__":
    main()
