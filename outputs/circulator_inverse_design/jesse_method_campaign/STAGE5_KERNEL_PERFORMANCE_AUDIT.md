# Stage-5 kernel performance audit

**Status:** COMPLETE  
**Stage-5:** **NOT STARTED / AWAITING AUTHORIZATION**  
**Branch:** `agent/eigenmode-ports`  
**Companion JSON:** `STAGE5_KERNEL_PERFORMANCE_AUDIT.json`  
**Design:** Stage-4 `q_best` (`opt_broadband_stage4_A2_C.json`)  
**Machine env:** `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`

Physics, mesh, A2_hybrid (τ=0.05, α=0.25), ports, B convention, adjoint equations,
bounds, and the Stage-5 **321-frequency / 6.25 MHz** sampling rule were **not** changed.

---

## Decision

### **A — meaningful implementation speedup achieved — use optimized kernel for Stage-5**

| Metric | Reference | Optimized | Result |
|--------|----------:|----------:|--------|
| Single-freq wall (warm) | ~134 s | ~17 s | **~7.8×** |
| 48-freq / 12 workers | **1220.3 s** | **171.9 s** | **7.10×** |
| Numerical J / ∇J / S | — | — | **PASS** (near FP) |
| Revised Stage-5 @74 evals | ~6.3 days | **~24 h** | still multi-hour, far cheaper |

Default production path: `use_fast_kernel=True`  
(`PLASMEEP_A2_FAST_KERNEL=0` restores the reference implementation).

---

## 1. Current (reference) per-frequency kernel profile

Measured on Stage-4 `q_best`, FULL91-C, representative hard frequency **5.08125 GHz**
(and confirmed at 4.5 / 4.73125 / 5.8125). Coarse phases from production `profile=True`:

| Phase | Wall (s) | Share |
|-------|--------:|------:|
| `assemble_build_system` (mesh+ρ+bulk+ports+Robin) | **~85–86** | **~64%** |
| `dA_contraction` (16 orbit assemblies + matvecs) | **~36–37** | **~27%** |
| SuperLU `splu` factor | **~9.0** | **~7%** |
| Adjoint 6 RHS | ~2.1 | ~1.5% |
| Forward 6 RHS | ~0.7 | ~0.5% |
| S + objective | ~0.3–0.5 | — |
| **Total** | **~134** | 100% |

This **confirms** the earlier pre–Stage-3 audit (~63% assembly, ~29% dA, ~6–7% factor)
on the **current** Stage-4/Stage-5 code path.

Finer instrumentation at 5.08125 GHz (`profile_stage5_kernel.py` →
`STAGE5_KERNEL_BASELINE_PROFILE.json`):

| Sub-phase | Wall (s) | Share of total |
|-----------|--------:|---------------:|
| `assemble_port_operator` (24× CSR↔LIL Robin updates) | **57.6** | **42.2%** |
| `dA` sparse assemble (16× full-mesh) | **34.8** | **25.5%** |
| `build_ports` (boundary edge discovery) | **21.2** | **15.5%** |
| SuperLU factor | 8.7 | 6.4% |
| Bulk `assemble_anisotropic` | 5.0 | 3.7% |
| load mesh + masks + ρ + exterior | 3.6 | 2.7% |
| `dA` matvecs | 1.7 | 1.3% |
| forward + adjoint solves + S | ~3.8 | ~2.8% |

So the historical “assembly ~63%” bucket is mostly **port Robin format thrash + port
geometry rebuild**, not the bulk FEM stiffness loop; “dA ~29%” is almost entirely
**re-assembling full-mesh sparse dA** rather than the matvecs.

---

## 2. Implementation speedups (mathematics unchanged)

| Change | What | Why exact |
|--------|------|-----------|
| `MeshGeomCache` | Precompute area, ∇φ (`bx`,`by`), PML stretch, COO index pattern | Same weak form; geometry is mesh-only |
| `DeviceStaticCache` | Cache mesh, material masks, rod masks, exterior, port edge geometry | Frequency-independent for this device topology |
| `assemble_anisotropic_cached` | Use cached geometry; pin empty rows via diagonal COO add (no LIL Python loop) | Identical stiffness/mass; empty rows still get identity |
| `assemble_port_operator_fast` | **One** CSR→LIL→CSR around all Robin updates | Same (−iβ)/(+iα) `c cᵀ` updates as reference |
| `contract_dA_orbit` | Element-local Λ·(dA U) on orbit support only; no CSR `dA` | Same as `assemble_dA_rod` with `mass_scale=0`, `pin_empty=False` |
| Worker initializer | Pre-build device cache when fast kernel enabled | No numeric change; removes cold-start from first task |

**Not changed:** weak form, gyrotropic ρ(s), A2_hybrid, soft-min weights, mesh, solver,
tolerances, ports/B, frequency set, optimizer.

---

## 3. Numerical equivalence (mandatory)

Acceptance: near floating-point for J, full 16-vector ∇q, S, and TE0 fields.

### Single-frequency (4 probes)

| f (GHz) | rel ‖ΔJ‖ | rel ‖Δg‖ | max‖ΔS‖ | max‖ΔU‖ | PASS |
|--------:|---------:|---------:|--------:|--------:|:----:|
| 4.50000 | 3.9e−14 | 4.6e−09 | <1e−12 | <1e−12 | **yes** |
| 4.73125 | 3.2e−15 | 1.1e−11 | <1e−12 | <1e−12 | **yes** |
| 5.08125 | 2.1e−13 | 3.6e−10 | 5.4e−14 | 3.5e−13 | **yes** |
| 5.81250 | 1.6e−14 | 2.3e−12 | <1e−12 | <1e−12 | **yes** |

Isolated bulk-operator check: `max|A_ref − A_fast| = 0` on a full assemble.

### 48-frequency / 12-worker broadband

| | Reference | Fast |
|--|----------:|-----:|
| J_BB | −7.748864969756681 | −7.748864969756545 |
| rel ΔJ_BB | — | **1.8e−14** |
| rel ‖Δ∇q‖ | — | **9.4e−11** |
| Equivalence | — | **PASS** |

Reference path remains selectable via `use_fast_kernel=False` /
`PLASMEEP_A2_FAST_KERNEL=0`.

---

## 4. Single-frequency speedup

Warm fast = production-relevant (cache already built; as in multi-freq workers).

| f (GHz) | ref (s) | fast warm (s) | fast cold (s) | speedup |
|--------:|--------:|--------------:|--------------:|--------:|
| 4.50000 | 134.6 | **17.0** | 43.8 | **7.92×** |
| 4.73125 | 134.5 | **17.2** | 43.7 | **7.83×** |
| 5.08125 | 135.0 | **17.2** | 43.4 | **7.83×** |
| 5.81250 | 133.9 | **17.9** | 44.2 | **7.48×** |
| **mean** | ~134.5 | **~17.3** | ~43.8 | **~7.77×** |

Warm phase split (typical): assemble ~4.5 s, factor ~8.8 s, dA ~0.4 s, solves ~3 s.

Peak RSS (single process): ~3–4 GiB (similar ref vs fast).

---

## 5. Parallel path that matters for Stage-5

Probe: **48** frequencies spanning 4–6 GHz (same class as
`BROADBAND_PARALLEL_SCALING_48FREQ.json`).  
`n_workers=12`, OMP/OpenBLAS/MKL=1.

| Kernel | Wall (s) | J_BB | ‖∇q‖ | Peak note |
|--------|---------:|-----:|-----:|-----------|
| Reference | **1220.3** | −7.74886497 | 1.447784e3 | ~prior ~45 GiB class |
| Fast | **171.9** | −7.74886497 | 1.447784e3 | same worker count |
| **Speedup** | **7.10×** | ΔJ rel 1.8e−14 | Δg rel 9.4e−11 | |

Per-worker wall under load ≈170 s for 4 freqs (~42 s/freq) vs ~17 s isolated warm:
memory-bandwidth contention remains, as in the earlier worker-scaling study. **Keep
`n_workers=12`.**

---

## 6. Revised Stage-5 cost (from measured 48-freq fast wall)

Sampling rule **unchanged**:

```python
freqs = [round(4.0 + i * 0.00625, 6) for i in range(321)]  # n_freq=321
```

Scaling method:  
`t_eval ≈ t_48_fast × ceil(321/12) / ceil(48/12) = 171.9 × 27/4`.

| Quantity | Value |
|----------|------:|
| n_freq | **321** |
| n_workers | **12** |
| batches / BB eval | **27** |
| **wall / BB eval** | **~1160 s ≈ 0.32 h** |
| expected peak RAM | **~45–50 GiB** (12 workers) |
| total @ **50** evals | **~16.1 h** |
| total @ **60** evals | **~19.3 h** |
| total @ **74** evals | **~23.9 h (~1 day)** |

Previous estimate (pre-kernel work): ~2.05 h/eval → ~6.3 days @74.  
**New estimate: ~0.32 h/eval → ~1 day @74.**

---

## 7. Recommendation

**Use the optimized kernel for any future authorized Stage-5 run.**

- Exact frequency rule: full-band 6.25 MHz (`stage5_frequency_set.py`) — unchanged  
- `n_workers=12`, OMP/OpenBLAS/MKL=1  
- Expect ~0.32 h per objective+gradient evaluation, ~24 h for a 74-eval campaign  
- Stage-5 remains **NOT STARTED / AWAITING AUTHORIZATION**

```
STAGE-5: NOT STARTED / AWAITING AUTHORIZATION
```

### Code pointers

- Fast assembly/cache: `scripts/.../jesse_method_campaign/fem_assembly_fast.py`
- Production entry: `a2_efficient_core.eval_and_grad_A2(..., use_fast_kernel=...)`
- Audit driver: `audit_stage5_kernel_performance.py`
- Detailed baseline profiler: `profile_stage5_kernel.py`
