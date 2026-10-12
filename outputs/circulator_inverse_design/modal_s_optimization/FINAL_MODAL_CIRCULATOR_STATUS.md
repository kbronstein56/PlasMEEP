# Final modal circulator campaign status

Branch: `agent/eigenmode-ports`  
Pref inverse-design artifacts preserved unchanged.

## Gate board

| Gate | Status |
|------|--------|
| MODAL_S_OPTIMIZED_ONSAGER | **FAIL** |
| MODAL_S_NONUNIFORM_GRADIENT | **PASS** |
| TRUE_MODAL_PEAK_OPT | **NOT_REACHED** (blocked by Onsager FAIL) |
| TRUE_MODAL_FINE_MESH | **NOT_REACHED** |
| TRUE_MODAL_BROADBAND | **NOT_REACHED** |
| FREE91_MODAL | **NOT_ATTEMPTED** |

## Why optimization did not start

Frozen Onsager criterion (`r_F ≤ 5×10⁻³` on fundamental modal S at 3.85 GHz)
**fails** on the Pref sixfold peak with `r_F ≈ 0.040` on C/M/F (mesh-independent).

Root cause: production feed `clear_width=2` supports **propagating TE1** at 3.85 GHz
(`kc₁=π/W < k0`). Fundamental-only S is then an open multimode truncation and
need not obey Onsager; the reaction operator still does (`~10⁻¹³`). Domain PML
(`dpml=2a`) is ≪ `λ₁≈17a`, so TE1 is also poorly absorbed. Remesh-free
mitigations (monitor shifts, two-plane extract, σ scaling, short constrictions)
do not meet the frozen bar. See `MODAL_S_OPTIMIZED_ONSAGER_AUDIT.md`.

Campaign rule: do not optimize through an unvalidated S extractor.

## What did pass

- **Power balance** (Phase 2, FULL91-M): algebraic identity holds; optimized designs
  show ~13–15% PML+higher-mode power vs ~4% on uniform (`MODAL_S_POWER_BALANCE_REPORT.md`).
- **Nonuniform modal-S adjoint** (Phase 3): full quotient-rule FD agreement
  `≤ 5×10⁻⁴` on peak and broadband designs after fine h-sweep
  (`MODAL_S_NONUNIFORM_GRADIENT_VALIDATION.md`).

## Pref peak under modal S (characterization only; Onsager-unclean)

FULL91-M · 3.85 GHz · +0.05 T (from prior eval, unchanged):

| metric | value |
|--------|------:|
| T_des | 0.742 |
| IL | 1.30 dB |
| isolation | 16.0 dB |
| RL | 12.4 dB |
| through | 0.0037 |
| modal Onsager `r_F` | ≈0.040 |

## Recommended path to reopen optimization

Freeze a **new** extractor/geometry package *before* further pass/fail, then meet
the same numerical bar (do not raise 5×10⁻³):

1. Remesh with `clear_width ≤ π/k0 ≈ 1.947` so TE1 is cutoff at 3.85 GHz, **or**
2. Remesh with TE1-capable PML (`dpml` several `λ₁`) + validated 12×12 multimode S,
   and optimize on the fundamental block of that multimode S.

## Commits (this campaign segment)

| SHA | Message |
|-----|---------|
| `90ca515` | Audit optimized-design modal-S Onsager: FAIL frozen 5e-3 gate |
| *(next)* | Phase 2 power balance + Phase 3 nonuniform gradient PASS |

Pushed to `origin/agent/eigenmode-ports` only.
