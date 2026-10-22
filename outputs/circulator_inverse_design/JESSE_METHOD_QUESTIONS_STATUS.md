# Jesse method / science questions — status

Standalone status for the six meeting questions.
Branch `agent/eigenmode-ports`. Continuum SF design: `opt_uniform_C_warm_M.json`.
Do **not** quote M-mesh J=0.723 as physics; use Q–X2 bracket J≈0.626±0.007
(`MESH_CONVERGENCE_MREFINED_DESIGN.md`).

Literature log: `docs/circulator_literature/LITERATURE_LOG.md`.

---

## 1. Broadband sampling (how many points for 4–6 GHz isolation?)

**QUESTION**  
If we want isolation over 4–6 GHz, how many objective frequencies? Too few → sharp spikes?

**WHAT WE KNOW**  
- Band has **no modal-channel switch**: TE0+TE1 throughout; TE2 cut off (~7.49 GHz).
- Soft-min + mean aggregation exists (τ=0.05, α=0.25); weight FD ~5.7e-11.
- Literature (Gedeon 2026 preprint) shows integrated FoMs can hide in-band spikes.

**WHAT HAS BEEN NUMERICALLY VERIFIED**  
- Soft-min weight derivatives. Continuum SF mesh behavior. **Not yet:** dense 4–6 GHz spectrum of M-opt design or grid-sufficiency test.

**WHAT LITERATURE SAYS**  
Hammond (hybrid TD/FD multi-freq TO); Gedeon (autocorrelation smoothness). Neither gives a universal N.

**WHAT IS STILL UNKNOWN**  
Feature widths of *our* isolation/IL/TE1 spectra; N and Δf that pass a frozen grid-sufficiency bar.

**EXACT TEST**  
Staged spectral recon on practical mesh (C/M): 200→100→50/25 MHz local; compare worst iso/IL/RL vs dense audit; verify worst points on F/A. Adaptive enrichment: opt → audit → add worst → restart. Freeze bars **before** judging.

**CURRENT RECOMMENDATION**  
Do **not** claim “11 points.” Use adaptive enrichment; initial guess ~200 MHz then densify. Soft-min alone is insufficient without dense audit.

---

## 2. Time-domain validation

**QUESTION**  
Validate broadband with time domain.

**WHAT WE KNOW**  
Meep/PlasMEEP available; known Meep↔FEM plasma-interface disagreements at B=0.

**VERIFIED**  
FD FEM power balance / Onsager / continuum bracket. TD plan exists (`TIME_DOMAIN_BROADBAND_VALIDATION_PLAN.md`); full-device TD not yet executed for M-opt.

**LITERATURE**  
Hammond 2022: Meep hybrid TO — supports TD as validation/opt path, not automatic ground truth.

**UNKNOWN**  
Magnetized plasma Meep agreement on our FULL91 device over 4–6 GHz.

**EXACT TEST**  
(1) Small Voigt/Faraday Meep vs analytic/FEM. (2) B-sign check. (3) Broadband pulse + DFT/modal extract vs FEM dense spectrum. Report disagreement honestly.

**RECOMMENDATION**  
Meep = independent cross-check. Optimize in FD FEM; validate in TD. Do not tune FEM to Meep.

---

## 3. Differentiable broadband / JAX / autograd

**QUESTION**  
Make broadband differentiable. JAX? Autograd like ceviche?

**WHAT WE KNOW**  
Per-freq J is differentiable via discrete adjoint; BB = soft-min + mean of J(f). Chain rule: ∇J_BB = Σ w_i ∇J(f_i).

**VERIFIED**  
Soft-min weight FD. Efficient single-freq adjoint. **Not yet:** end-to-end multi-frequency FD of ∇J_BB on design coordinates.

**LITERATURE**  
Ceviche uses HIPS autograd with **custom sparse Ax=b VJPs** (= adjoint). Autograd does not differentiate arbitrary SciPy SuperLU without primitives (Hughes workshop FAQ).

**UNKNOWN**  
Whether JAX custom_vjp improves API without wall-clock gain on FULL91.

**EXACT TEST**  
3-frequency end-to-end FD of analytic ∇J_BB vs centered FD. Optional tiny JAX custom_vjp wrapping validated adjoint.

**RECOMMENDATION**  
Keep analytic FD adjoint + soft-min. JAX only as optional wrapper. No FULL91 port.

---

## 4. Gradient method / efficiency (“extra simulations”)

**QUESTION**  
Is there a frequency-domain adjoint? Is another simulation required?

**WHAT WE KNOW / VERIFIED**  
**Production method already is a frequency-domain discrete adjoint** using complex **bilinear** transpose (`A^T`, not casually `A^H`). Efficient path: **1 LU + 6 TE0 fwd + 6 adjoint RHS**. FULL91-C: 1397 s → 1195 s (1.17×). LU ≈ **84%** of runtime; adjoint RHS ~3 s after factor.

**LITERATURE**  
Hughes 2019: forward vs reverse scaling. For N_param=16, forward-mode after shared LU may compete.

**UNKNOWN**  
Measured forward-mode wall time on C6; PARDISO vs SuperLU if available; symbolic-factor reuse across nearby f.

**EXACT TEST**  
Forward-mode bench; solver inventory; document LU dominance for Jesse.

**RECOMMENDATION**  
Explain: adjoint RHS ≠ new factorization. LU dominates; solver improvements matter more than shaving RHS.

---

## 5. Mode preservation

**QUESTION**  
Must we preserve TE0? Can we measure which mode carries power?

**WHAT WE KNOW**  
Current FoM rewards desired TE0 and penalizes TE1. Zhao 2025 / Wang 2021: devices can use TE0↔TE1 conversion.

**VERIFIED**  
TE1 terms in current J. Continuum TE1≈0.0205. **Not yet:** port-power objective comparison opt.

**UNKNOWN**  
Performance sacrificed by TE0 preservation; experimental mode resolution feasibility.

**EXACT TEST**  
Freeze Objective A (mode-preserving) vs B (port-power total); evaluate baselines; short equal-budget opts; spectra of TE0/TE1 desired/reverse. Separate experimental feasibility memo.

**RECOMMENDATION**  
Treat TE0 preservation as a **design choice**. Run controlled A vs B before freezing BB production objective.

---

## 6. C6 vs C3 symmetry

**QUESTION**  
Why sixfold? Wouldn’t threefold make sense? Fan was three ports?

**WHAT WE KNOW**  
C6: 16 vars (center+15×6). Six-port cycle 1→2→…→1 is C6-invariant. Fan/Wang classic device is **3-port → C3**. C6 is a strong constraint, not required by circulator physics.

**VERIFIED**  
C6 orbit geometry PASS. **Not yet:** C3 orbit build + controlled opt.

**LITERATURE**  
Wang & Fan 2005; Li 2022 Y-junction plasma — both C3 three-port.

**UNKNOWN**  
Whether C3 improves J/isolation under same budget; whether 1 (C6) or 2 (C3) representative sources suffice by rotation.

**EXACT TEST**  
Build C3 orbits (~31 vars); map C6 design → C3 warm start; equal-budget C6 vs C3; verify rotational source reduction numerically.

**RECOMMENDATION**  
Answer Jesse: C6 because we imposed 60° symmetry on six ports; Fan is C3 because of three ports. Test C3 empirically before freezing BB symmetry.

---

## Campaign gate

Long 4–6 GHz optimization starts only after `BROADBAND_PRODUCTION_SPEC_FROZEN.md` exists.
