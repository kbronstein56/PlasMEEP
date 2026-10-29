# BROADBAND_PRODUCTION_SPEC_STAGE2_FROZEN

Frozen **2026-10-27** after Objective A2 offline ranking + FD validation.

**Version history**
- 2026-10-27: initial Stage-2 freeze (A2_hybrid, τ=0.05, α=0.25, 28 freqs).
- 2026-10-28: soft-min scale sensitivity accepted — **τ/α unchanged**
  (`A2_SOFTMIN_SCALE_SENSITIVITY.md`). Stage-2 optimization **authorized**
  via `production_broadband_opt_stage2.py` (checkpointed, restartable).
  Do not reopen τ/α unless Stage-2 shows a concrete optimization pathology.

Stage-1 historical freeze remains authoritative for the completed 21-freq run:
`BROADBAND_PRODUCTION_SPEC_FROZEN.md` — **do not overwrite**.

---

## Frozen choices

| Item | Freeze |
|------|--------|
| Symmetry | **C6** (16 vars) |
| Incident mode | **TE0** drive on each of 6 ports |
| Desired-output definition | **Objective A2_hybrid** (mode-preserving): \(J=2\log(T_\mathrm{des}+\varepsilon)-\log(T_\mathrm{rev}+\varepsilon)-1.5\,T_\mathrm{thr}-T_\mathrm{other}-0.25\,R-T_\mathrm{TE1}\), \(\varepsilon=10^{-6}\). No hard IL/isolation dB thresholds. |
| Broadband objective | \(J_{BB}=\alpha\,\mathrm{mean}(J_f)+(1-\alpha)\,\mathrm{softmin}_\tau(J_f)\) |
| τ | **0.05** |
| α | **0.25** |
| Frequency set | **28** points = stage-1 21×100 MHz **plus** `{4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95}` GHz |
| Warm start | `q_best` from `opt_broadband_production_C.json` |
| Optimization mesh | FULL91-**C** |
| Gradient method | Efficient frequency-domain discrete adjoint (`A^\top`, 1 LU + 6+6 RHS) with **state-dependent** A2 channel weights. No JAX. |
| Iteration budget (proposed run) | **30** L-BFGS iters; bounds `[0.05, 1.8]` |
| Expected wall-time | **~55–75 h** (scaled from stage-1 ~44 h × 28/21) |

### Frequency list (fixed during each L-BFGS run)

```
4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.55, 4.6, 4.7, 4.8, 4.85, 4.9,
5.0, 5.1, 5.2, 5.3, 5.35, 5.4, 5.5, 5.6, 5.625, 5.65, 5.7,
5.8, 5.85, 5.9, 5.95, 6.0
```

---

## Evidence pointers

| Gate | Artifact |
|------|----------|
| Alignment audit | `BB_OBJECTIVE_ALIGNMENT_AUDIT.md`, `BB_OBJECTIVE_ALIGNMENT_OFFLINE.json` |
| A vs A2 decision | `BB_OBJECTIVE_A_VS_A2.md` |
| A2 FD (SF + 3-freq soft-min) | `A2_GRADIENT_FD_VALIDATION.json` **PASS** |
| Enrichment | `BB_ENRICHMENT_RECOMMENDATION.md` |
| Implementation | `scripts/.../objective_A2.py` (`A2_FREEZE`) |

---

## Explicit non-goals of this freeze

- Does not reopen τ/α/objective mid-run without concrete pathology evidence.
- Does not auto-start Stage-3 or further enrichment during an active L-BFGS run.
- Does not reopen DtN / gyrotropic absorption validation.
- Does not claim stage-1 spectra meet engineering IL/isolation targets.
- Does not authorize quoting M-mesh single-frequency J as continuum physics.
