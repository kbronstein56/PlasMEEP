# Solver / Optimization Architecture Audit

**Date:** 2026-09-24  
**Scope:** Architecture decision only — no production changes, no commits.  
**Meep version inspected:** 1.30.0 (`plasmeep_mpi` env)  
**Ceviche inspected:** v0.1.3 wheel source (not installed into production)

---

## Executive recommendation

| Goal | Best near-term architecture |
|---|---|
| **(a) Optimization running soon** | **Meep FDTD + finite differences or coarse Bayesian** on 91 `fp` at **validated 25 ppc**, B fixed — *or* explore Meep adjoint **only if** restricted to B=0 Lorentzian/`MaterialGrid` mapping (see Part 4). Stock Meep adjoint does **not** natively optimize gyrotropic 91-bulb plasmas. |
| **(b) High-resolution optimization** | **Not ready.** FEM is fast but **NOT VALIDATED** vs Meep. Uniform FDFD/Ceviche at 50 ppc is ~8.4M DOFs — loses the local-refinement speed win. |
| **(c) Broadband optimization** | **Meep FDTD(+adjoint if applicable)** wins structurally: one pulse → many DFT frequencies. FEM/Ceviche pay **per frequency** (new `A` + factor). |
| **(d) Final validation** | **Always Meep/PlasMEEP** (trusted reference). |

**Do not treat the ~11 s FEM forward as the architecture winner.** Accuracy gate failed (max ~9 dB B=0 error). Adjoint math for FEM is promising (PoC gradient error ~4×10⁻⁹) but useless for production until ports match Meep.

**Hybrid that makes sense:** keep developing FEM (or anisotropic FDFD) as a *candidate* inner loop **in parallel with** Meep validation; optimize only after B=0 port agreement. Until then, Meep remains both the physics oracle and the only credible optimization engine for gyrotropic plasmas.

---

## Part 1 — What exactly is the FEM?

### 1–2. Built from scratch vs packages

**Mostly custom electromagnetic FEM**, not a full FEM application stack.

| Concern | What is actually used |
|---|---|
| Mesh generation | **Primary meshes:** custom `scipy.spatial.Delaunay` sampler (`fem_delaunay_mesh.py`). **Optional / unused in validated path:** `gmsh` in `fem_mesh_build.py` (import fails here: `OSError` / libGLU) |
| Basis / FE library | **None.** Hand-written **P1 (CG1) triangles**. `skfem` is **installed** but **not imported** by FEM solvers |
| Assembly | **Custom** NumPy vectorized element loops → `scipy.sparse` |
| PML | **Custom** complex coordinate stretch on element centroids |
| Sparse matrices | `scipy.sparse` (COO→CSR/CSC) |
| Factorization / solve | `scipy.sparse.linalg.splu` → **SuperLU** |

### 3–6. gmsh / skfem / SciPy / SuperLU

| Package | Role |
|---|---|
| **gmsh** | Present in code path for OCC meshing; **not used** for current FEM-H artifacts (Delaunay meshes). Runtime import fails in this env. |
| **scikit-fem** | Installed (12.0.2); **unused** by our FEM code. |
| **SciPy** | Core: Delaunay, sparse, SuperLU. |
| **SuperLU** | Yes — via `splu`. Transpose/`H` solves work without refactor (`trans='N'|'T'|'H'`). |

### 7–8. Custom physics code size

Approx. **~1.5k lines** of FEM-specific scripts:

| File | Lines | Role |
|---|---:|---|
| `fem_validated_solver.py` | 659 | Anisotropic ρ, Poynting, numerical-mode source, walls |
| `fem_hz_solver.py` | 394 | Earlier isotropic FEM + \|Hz\|² proxy |
| `fem_delaunay_mesh.py` | 246 | Local refinement mesher |
| `fem_mesh_build.py` | 233 | gmsh path (unused) |
| + diagnostics | ~700 | compare/diagnose/tests |

Plus reuse of PlasMEEP/`sixport_common` for geometry constants and mode JSON.

### 9–10. Feature matrix (implemented vs planned)

| Feature | Status |
|---|---|
| Anisotropic ε / complex off-diag | **Implemented** in validated solver (via ρ=ε⁻¹) |
| Full ρ=ε⁻¹ | **Implemented** |
| PML (stretched) | **Implemented** (heuristic σ) |
| PEC | **ρ→0 metal triangles** (TE-correct). Dirichlet Hz=0 was tried and is **wrong** (PMC) |
| Modal / numerical-mode sources | **Implemented** (nearest-node scatter of Meep mode samples) |
| Guide-normal Poynting ports | **Implemented** |
| Gradients / adjoint | **Not in solver API.** Cheap PoC only (see Part 7) |
| Validated vs Meep | **NOT VALIDATED** (see Part 9) |

### Architecture diagram (FEM)

```
sixport_common geometry + NumericalPortMode JSON
        │
        ▼
Delaunay mesh (.npz) ──► assign ε/ρ (bulbs, quartz, vacuum)
        │                         + ρ→0 horn walls
        ▼
Custom P1 assembly ──► A(ω,p) ∈ ℂ^{n×n}  (scipy.sparse)
        │
        ▼
SuperLU factor(A) ──► solve A x = b_port   (reuse for ports 2..6)
        │
        ▼
E = (i/ω) ρ curl H  ──► ∫ S·n̂ ds  ──► normalize / field-sub P11
```

---

## Part 2 — PlasMEEP architecture (what is reusable)

| Layer | Where | Reusable if solver changes? |
|---|---|---|
| Geometry / horns / 91-bulb array | `sixport_common`, `plasmeep/ports/horn.py`, `PMMI` | **Yes** (coordinates, clear width, bulb centers) |
| Plasma / gyrotropic material model | `plasmeep/lib.py` `Get_Med`, `faraday_benchmark.gyrotropic_drude_eps_eta` | **Yes** (ε(ω,B,fp) formulas) |
| Source construction | `port_formulations`, `numerical_launch`, mode registry | **Partially** — spatial mode profiles yes; Meep `Source` objects no |
| Receivers / flux / modal | `port_formulations`, ports/* | **Partially** — monitor planes/normals yes; Meep DFT monitors no |
| Normalization | `normalize_port`, flux pickle cache | Meep-specific field subtraction; **idea** reusable |
| S-matrix / objective | `simulate_circulator`, `circulator_objective` | **Yes** (operates on power matrix) |
| Optimization | notebook / inverse scripts; **no production adjoint** | Objective definition reusable |
| Meep simulation glue | `Get_Sim`, FDTD run | Meep-only |
| Adjoint | **None in PlasMEEP** | — |

**Observed:** repo has **no** `meep.adjoint.OptimizationProblem` integration for the 91 plasmas.

---

## Part 3 — Can `solve_cw` handle our plasma?

Diagnostics: `solver_architecture/solve_cw_complex_fields.json`  
(with `force_complex_fields=True`; without it, all cases fail `use_real_fields`).

| Material | `solve_cw` result |
|---|---|
| Vacuum / real ε | **OK** |
| Real anisotropic diag / real offdiag | **OK** |
| `D_conductivity` (lossy) | **OK** |
| Complex scalar `epsilon=2+0.1j` in constructor | **TypeError** (API) |
| **Imaginary ε off-diagonal** (`±iη`) | **Hard fail:** `Found non-zero imaginary part of epsilon or mu offdiag` |
| `DrudeSusceptibility` | **Fail:** `dispersive materials are not yet implemented for solve_cw` |
| `GyrotropicDrudeSusceptibility` | **Same dispersive fail** |

### Answers

1. Ordinary nondispersive complex ε: **partial** (loss via conductivity works; complex `epsilon=` constructor does not).
2. Anisotropic real ε: **yes**.
3. Complex off-diagonal: **no** (Meep rejects Im(ε_offdiag)).
4. Fixed-frequency magnetized plasma tensor `[[ε⊥,-iη],[iη,ε⊥]]`: **cannot** be injected as a Medium tensor into this Meep.
5. Direct GyrotropicDrude + `solve_cw`: **no**.
6. Evaluate Drude/gyro at ω and pass fixed tensor: **blocked** by (4).
7. Antisymmetric ±iη via material interface: **not** for static ε; only via **time-domain** gyrotropic susceptibility ODEs.
8. Would require **Meep source changes** (allow complex offdiag and/or CW dispersive) — not a small API tweak.

**Verdict:** Meep **cannot** give a true single-frequency magnetized-plasma CW solve for our tensor without modifying Meep. Time-domain FDTD + GyrotropicDrude remains Meep’s path.

---

## Part 4 — Can Meep’s adjoint optimize our 91 plasmas?

Inspected: `meep.adjoint.OptimizationProblem`, `MaterialGrid` docs in `geom.py` (1.30).

### Native design variables

- `DesignRegion` ↔ **`MaterialGrid` weights** u∈[0,1] interpolating **medium1 ↔ medium2**.
- Supported material types for grids (**documented explicitly**):
  1. frequency-independent isotropic ε (`epsilon_diag` / `epsilon_offdiag` interpolated);
  2. **`LorentzianSusceptibility`** (`sigma` / `sigma_offdiag` interpolated).
- **Not listed:** `DrudeSusceptibility`, **`GyrotropicDrudeSusceptibility`**.

### Checklist

| # | Question | Answer |
|---|---|---|
| 1 | Native DVs | MaterialGrid weights (topology-style), not raw `fp` |
| 2 | MaterialGrid required? | **Yes** for stock adjoint design regions |
| 3 | Interpolate dispersive? | **Lorentzian only** (per docs) |
| 4 | GyrotropicDrude? | **No** (not in supported types) |
| 5 | One param per bulb? | **Possible pattern:** 91 tiny grids/regions — awkward but conceivable for B=0 Lorentzian mapping |
| 6 | Chain ρ→fp→susceptibility→J for 91 | **Not for gyrotropic.** B=0 Lorentzian σ(fp) would need custom mapping outside stock examples → **MAJOR** |
| 7 | Guide-normal / our modal objective | Stock: **`EigenmodeCoefficient`**, FourierFields, LDOS — **not** our custom `num_mode_guide_normal` flux path |
| 8 | Multiple ports | Multiple objective args / multi-objective — **yes in API** |
| 9 | Multiple frequencies one FDTD run | **Yes** (`frequencies=` / fcen,df,nf + DFT) |
| 10 | B≠0 gyrotropic | **Not supported** by MaterialGrid adjoint path |
| 11 | Modifications | New adjoint material model for gyro Drude **or** external FD on Meep FDTD |

### Classification for **exact** 91-plasma gyrotropic problem

**NOT PRACTICAL** with stock Meep adjoint.

| Subcase | Rating |
|---|---|
| B=0, map each bulb vacuum↔Lorentzian via MaterialGrid | **MAJOR CUSTOM WORK** (objective ports + 91 regions + validation) |
| B≠0 gyrotropic 91-fp | **NOT PRACTICAL** without Meep core/adjoint extension |
| Broadband multi-freq DFT adjoint (generic dielectrics) | **SUPPORTED NOW** (but wrong materials for us) |

PlasMEEP: **no existing adjoint derivative path** for `rho`.

---

## Part 5 — `solve_cw` + Meep adjoint together?

**No.**

`OptimizationProblem.forward_run` / `adjoint_run` call `sim.run(... until_after_sources=stop_when_dft_decayed)` — **time-domain DFT**, not `solve_cw`.

`solve_cw` is a separate monochromatic CG-like solver on the Yee system; it is **not** wired as the adjoint engine. Combining them would be new research/engineering, not an installed feature.

---

## Part 6 — Ceviche audit

Source: extracted `ceviche-0.1.3` wheel (not production-installed).

| Question | Answer |
|---|---|
| Hz polarization | **Yes** — `fdfd_hz` |
| Complex ε | **Yes** — `eps_r` array can be complex |
| Anisotropic / ε_xy | **No** in stock `fdfd_hz` — scalar `eps_r` → diagonal 1/ε only |
| Magnetized tensor | **Not stock.** `fdfd_mf_ez` is for **modulation** sidebands, not B∥z gyro plasma |
| Extension invasiveness | Modify `_make_A` for full 2×2 ρ or ε tensor; touch PML averaging — **moderate** (hundreds of lines) but must re-derive AD |
| Autodiff after tensor | **Yes in principle** — assembly already uses `autograd.numpy` + custom `sp_solve` VJPs; new entries must stay in that graph |
| 91 densities as parameters | **Yes naturally** — paint `eps_r` (or tensor fields) from bulb masks; AD through `eps_r` |
| Multi-RHS reuse | Same `A`; Ceviche `sp_solve` can refactor or you cache factorization outside — stock API solves once per `solve()` |
| Nonuniform / FEM refinement | **No** — uniform `dL` grid only |
| 50 ppc unknowns | Domain 30×28 a @ res 100 → **N = 8.4×10⁶** Hz DOFs (vs FEM-H **6.6×10⁵**). Destroys local-refinement advantage; memory/runtime like prior FDFD feasibility conclusions |

**Effort estimate:** 1–3 weeks for anisotropic Hz FDFD + AD + port plumbing; **plus** full Meep validation (same problem FEM has). Not a free adjoint lunch for gyrotropy.

---

## Part 7 — Adjoint for the fast FEM solver

### Discrete problem

\[
A(p,\omega)\,x = b,\qquad
A\in\mathbb{C}^{n\times n},\ p\in\mathbb{R}^{91}.
\]

PoC objective (complex fields, real parameter), \(J=|x_k|^2\):

\[
A^{H}\lambda = e_k\,x_k,
\qquad
\frac{\partial J}{\partial p} = -2\,\mathrm{Re}\big(\lambda^{\dagger} (\partial A/\partial p)\,x\big).
\]

**Measured** on FEM-M vs central difference: **relative error ≈ 3.8×10⁻⁹**  
(`solver_architecture/fem_adjoint_fd_poc.json`).

For a real port-power objective \(J(x,\bar x,p)\), same pattern: one adjoint RHS from \(\partial J/\partial x^*\), then contractions with sparse \(\partial A/\partial p_i\).

### Answers

| # | Answer |
|---|---|
| 1 | Per frequency: **1 factor** + **N_src forward solves** + **1 adjoint solve** (independent of 91) |
| 2 | Reuse LU for adjoint? **Yes** |
| 3 | \(A^H\lambda\) cheap after factor? **Yes** |
| 4 | SuperLU `trans='H'` without refactor? **Yes** |
| 5 | FEM-H bench (`superlu_transpose_bench.json`): factor ~5.6 s; solve N/T/H all **~0.19–0.21 s**, residuals ~1e−13 |
| 6 | Need \(\partial A/\partial p_i\) from \(\partial\rho/\partial f_{p,i}\) on that bulb’s triangles only |
| 7 | Local? **Yes** — each \(p_i\) touches one bulb’s elements |
| 8 | Gradient time estimate (FEM-H, 1 ω): factor 5.6 + 6×0.19 forwards + 0.21 adjoint + 91×(local dA·x + dot) ≪ 1 s → **~7–8 s linear algebra** + assembly/PP. Full candidate with ports ~**15–25 s / frequency** (dominated by assemble/PP today) |

**Caveat:** PoC used \|Hz\|², not validated Poynting S-params. Adjoint for wrong physics is still wrong.

---

## Part 8 — Broadband cost

### What the objective is today

- `circulator_objective(power_matrix)` — **single-frequency** 6×6 power matrix at `fs` (~3.85 GHz).
- Project intent (`PROJECT.md`, validation-status): **broadband** S-parameter behavior — **not yet implemented** as multi-frequency sampling.
- Bandwidth / Nf: **not fixed in code**; agent-kit marks “Broadband objective validated” as TODO.

Assume discrete samples Nf ∈ {1,5,10,20} for planning.

### Meep FDTD (+ adjoint if usable)

- **One** forward pulse can fill **all Nf** DFT monitors.
- **One** adjoint pulse can yield gradients at **all Nf** (stock adjoint).
- Measured: **~12 min** device @ 25 ppc / np32 / rt20 (one source); **~72 min** naive 6 sources.  
  Trusted 50 ppc profile: **~72 min** one source (`trusted_ppc50_profile.json` ~4346 s).
- Adjoint ≈ another full FDTD of similar cost → objective+grad ~**2×** one multi-port strategy (often 1 forward + 1 adjoint with multi-source formulations — depends on implementation).

Broadband is Meep’s structural strength.

### FEM

- **Each frequency:** rebuild ε(ω), assemble, **refactor**.
- Six horns **reuse** factor; adjoint **reuses** factor (`H` solve).
- Using measured FEM-H: ~6 s factor + ~1.2 s (6 solves) + ~0.2 s adjoint ≈ **~8 s LA / freq**, plus material/assemble/PP (~5–10 s) → **~15–25 s / freq**.
- Nf=5 → ~**1.5–2 min**; Nf=10 → ~**3–4 min**; Nf=20 → ~**5–8 min** — **if** validated. Today it is not.

| Nf | FEM obj+grad (proj., if valid) | Meep FDTD 25ppc 6-col + adjoint (order-of-mag) |
|---:|---:|---:|
| 1 | ~0.3–0.5 min | tens of minutes |
| 5 | ~2–3 min | similar to 1 (DFT) |
| 10 | ~3–5 min | similar to 1 |
| 20 | ~6–10 min | similar to 1 |

---

## Part 9 — FEM validation status

From `outputs/validation/fem_meep_validation/FEM_MEEP_VALIDATION.md`:

| Item | Status |
|---|---|
| Poynting flux | Implemented |
| Same port planes / numerical mode | Attempted |
| Normalization / field-sub P11 | Attempted |
| B=0 agreement | **FAIL** — max \|ΔdB\| ≈ **9.3** (P13/P15); P14 ~0.5 dB |
| Mesh convergence | **FAIL** — FEM-M vs H ~2 dB on P14 |
| Full gyrotropic tensor | Code + unit tests; **no** device B≠0 vs Meep |
| B≠0 agreement | **Not run** (gate) |

**Verdict remains NOT VALIDATED.** Speed is irrelevant until this changes.

---

## Part 10 — Decision matrix

| Criterion | A Meep FDTD+adjoint | B Meep solve_cw+adj | C Ceviche FDFD+AD | D Local FEM+custom adj | E Hybrid |
|---|---|---|---|---|---|
| Full gyrotropic plasma | **Yes** (FDTD) | **No** | Needs extension | Code yes / valid no | Meep validates |
| Broadband efficiency | **Excellent** | N/A | Poor (per-ω) | Poor (per-ω) | Meep for BB |
| Single-ω speed @ hi-res | Poor (50ppc ~1h/src) | N/A | Poor (8M DOF) | **Good** (~10–20s) if valid | FEM opt / Meep check |
| Six-port efficiency | 6 FDTD or clever multi | N/A | Multi-RHS OK | **Excellent** (reuse LU) | — |
| Adjoint availability | Stock for MaterialGrid | **No** | Autograd yes | PoC yes; prod no | — |
| 91-fp parameterization | **Hard** (gyro) | No | Natural | Natural | — |
| Local mesh refinement | No (Yee) | No | No | **Yes** | — |
| Hi-res J time | ≫10 min @50ppc | — | Likely ≫10 min | Potentially &lt;5 min | — |
| Hi-res J+∇ time | ~2× FDTD | — | Large | Potentially &lt;10 min | — |
| Dev effort | Low if accept FD/Bayesian | Dead end | Medium–high | High (validate+adjoint ports) | Medium |
| Validation maturity | **Highest** | — | None here | **Failed B=0** | — |
| Risk | Slow opts | Impossible gyro CW | Uniform grid tax | Wrong physics, fast | Coordination |

---

## Part 11 — Hybrid architectures

### Recommended hybrid (if FEM is ever validated)

1. **Inner loop:** FEM (or anisotropic FDFD) forward + discrete adjoint at Nf discrete frequencies, 91 `fp`, fixed B.  
2. **Outer trust:** Meep FDTD broadband + magnetized checks on candidates / final designs.  
3. **Never** replace Meep as the acceptance test.

### Alternative hybrid (optimization sooner)

1. **Optimize in Meep at 25 ppc** with finite differences / Bayesian / evolutionary on 91 `fp` (and maybe B sweep outside).  
2. Use FEM only as a **speed research track** until port errors &lt;~0.25 dB.  
3. Stock Meep adjoint only if you deliberately simplify materials (B=0 Lorentzian topology) — **not** the full scientific target.

### Anti-pattern

Using unvalidated FEM gradients to drive design, then “checking” with Meep once — will optimize numerical artifacts (already seen: wrong PMC walls looked like good P14).

---

## Explicit answers (Part FINAL)

1. **FEM built from:** custom P1 + SciPy Delaunay/sparse/SuperLU; gmsh optional unused; skfem unused.  
2. **Custom FEM physics:** ~1.5k lines (+ diagnostics).  
3. **solve_cw + our gyro tensor:** **No** (no dispersive CW; no Im offdiag ε).  
4. **solve_cw + Meep adjoint:** **No** (adjoint is time-domain DFT).  
5. **Meep adjoint → 91 plasma fp:** **Not natively**; MaterialGrid ≠ fp bulbs.  
6. **Meep adjoint + gyrotropic/dispersive map we need:** **No** for GyrotropicDrude; Lorentzian-only grids.  
7. **Ceviche AD machinery:** **Yes** (autograd + sparse VJPs) for scalar ε.  
8. **Ceviche for off-diag plasma:** needs `_make_A` tensor extension — moderate, must keep AD.  
9. **Uniform Ceviche @50ppc:** **Yes, destroys** FEM’s DOF advantage (~8.4M vs 0.66M).  
10. **FEM adjoint difficulty:** Math straightforward; **code medium**; **validation is the hard part**.  
11. **Reuse LU for adjoint:** **Yes** (`trans='H'`).  
12. **91-bulb gradient time:** LA ~**&lt;1 s** after factor; end-to-end ~**15–25 s/ω** with current PP (proj.).  
13. **Nf=5/10/20:** FEM cost ≈ linear in Nf; Meep DFT ≈ flat in Nf.  
14. **Best for:**  
    - (a) soon: **Meep + non-adjoint opt @25ppc**  
    - (b) hi-res opt: **only after FEM (or other FD) validates**  
    - (c) broadband: **Meep FDTD**  
    - (d) final validation: **Meep/PlasMEEP**  
15. **Implement NEXT (no production leap):**  
    1. **Close FEM↔Meep B=0 port gap** (or formally retire FEM for opt).  
    2. **Define broadband objective** (Nf, band, which S_ij) in writing + Meep DFT harness.  
    3. **Decide B strategy** (fixed vs outer loop) before adjoint investment.  
    4. Only then: FEM port-power adjoint **or** Meep FD/Bayesian pilot on 91 `fp`.

---

## Artifact index

| Path | Content |
|---|---|
| `solve_cw_complex_fields.json` | CW capability matrix |
| `superlu_transpose_bench.json` | FEM-H SuperLU N/T/H timings |
| `fem_adjoint_fd_poc.json` | Adjoint vs FD (rel err ~4e−9) |
| `ceviche_inspect/ceviche/` | Unpacked Ceviche 0.1.3 sources |
| `../fem_meep_validation/FEM_MEEP_VALIDATION.md` | Accuracy NOT VALIDATED |

---

## Bottom line

- **Meep FDTD** is still the only installed path that correctly carries **gyrotropic dispersive** plasma.  
- **Meep adjoint** is excellent for **MaterialGrid / Lorentzian topology**, **not** for our **91 gyrotropic fp** variables.  
- **solve_cw** is a dead end for this plasma model in Meep 1.30.  
- **Ceviche** gives AD but **scalar uniform FDFD** — bad fit for hi-res + gyro without real work, and uniform grids erase FEM’s DOF win.  
- **FEM + custom adjoint** is mathematically attractive and **experimentally adjoint-feasible**, but **physically unvalidated**; do not optimize on it yet.

**Architecture choice until FEM matches Meep:** optimize/validate in **Meep**; treat FEM as a **conditional accelerator** gated on port accuracy.
