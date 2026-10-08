# FDFD GO/NO-GO — Custom single-frequency Hz solver

**Date:** 2026-09-16  
**Verdict: NO-GO** for TRUE 50 ppc / 3–5 minute one-source solves with the tools available on this machine.  
**Formulation note:** Scalar-Hz frequency-domain physics is sound; the blocker is the linear solver stack, not Maxwell reduction.

Artifacts: `outputs/validation/fdfd_feasibility/`  
Prototype (benchmark-only): `fdfd_hz_prototype.py`

---

## 1. Exact frequency-domain physics

### Tensor (Meep / Faraday `e^{-iωt}` convention)

From `faraday_benchmark.gyrotropic_drude_eps_eta` at \(f_s \approx 3.85\,\mathrm{GHz}\) (\(f_s a/c = 0.256844\)):

\[
\varepsilon(f_s,f_p,\gamma,b_z)
=
\begin{pmatrix}
\varepsilon_\perp & -i\eta \\
+i\eta & \varepsilon_\perp
\end{pmatrix}
\]

| bias \(b_z\) (Meep) | \(\varepsilon_\perp\) | \(\eta\) |
|---|---|---|
| 0 | \(-3.31776 + 0.0011215i\) | \(0\) |
| 0.05 | \(-3.48783 + 0.0012575i\) | \(0.87365 - 0.0004717i\) |

Signs match the Faraday benchmark (`eps_+ = eps_perp + eta`, `eps_xy = -i η`).

### Scalar Hz reduction — YES

Polarization: \(H_z, E_x, E_y\) nonzero; \(\partial_z=0\); \(\mu=\mu_0\); B along \(z\).

From \(\nabla\times H = -i\omega\varepsilon_0\varepsilon E\) and \((\nabla\times E)_z = i\omega\mu_0 H_z\):

\[
\nabla\cdot\bigl(\rho\,\nabla H_z\bigr) + k_0^2 H_z = 0,
\qquad
\rho = \varepsilon_{2D}^{-1},
\quad
k_0 = 2\pi f_s\ \text{(Meep units)}.
\]

At \(B=0\): \(\rho = \varepsilon^{-1} I\).  
At \(B\neq 0\): full anisotropic \(\rho\) (9-point stencil). **One scalar unknown per cell.**

Coupled \(E_x,E_y,H_z\) is **not** required for this 2D TE-to-z problem.

PML: complex coordinate stretch \(s_x,s_y = 1 + i\sigma/\omega\) in the divergence form (prototype CFS-like).

---

## 2. Problem size (real six-port domain \(30\times 28\))

| ppc | res | Nx×Ny | scalar unknowns | coupled×3 | nnz (meas./est.) | iterative mem est. | direct LU fill est. |
|---:|---:|---|---:|---:|---:|---:|---:|
| 25 | 50 | 1500×1400 | **2.10M** | 6.3M | 10.2M (meas.) | ~1.6 GiB | ~11 GiB |
| 35 | 70 | 2100×1960 | **4.12M** | 12.3M | 20.0M (meas.) | ~3.2 GiB | ~27 GiB |
| 50 | 100 | 3000×2800 | **8.40M** | 25.2M | ~76M | ~6.5 GiB | ~74 GiB |
| 60 | 120 | 3600×3360 | **12.1M** | 36.3M | ~109M | ~9.4 GiB | ~123 GiB |

**Sparse DIRECT at 50 ppc vs ~252 GiB WSL RAM:**  
Memory *might* fit (~74 GiB fill estimate). **Time does not.** Measured SuperLU at 10 ppc (336k) ≈ 2.5 min factorization → \(N^{1.5}\) projection ≈ **~9.5 hours** at 50 ppc. **Ruled out for 3–5 min. Do not waste time testing SuperLU at 50 ppc.**

---

## 3. Available solvers / hardware

| Item | Status |
|---|---|
| scipy sparse (gmres, bicgstab, splu, spilu, spsolve) | **YES** |
| PETSc / petsc4py | NO |
| hypre | NO |
| PyAMG | NO |
| MUMPS / SuperLU_DIST | NO |
| CuPy | NO |
| NVIDIA GPU (`nvidia-smi`) | **NO GPU** |
| RAM / CPU | ~251 GiB / 128 threads |

### Ranking for ~8M complex Helmholtz/gyrotropic

| Approach | Chance here |
|---|---|
| GMRES / BiCGSTAB / FGMRES alone | Poor (indefinite Helmholtz) |
| Jacobi / ILU0 (scipy) | **Failed** — true residual stalls ~0.05–0.2 |
| Algebraic multigrid (standard) | Weak for high-frequency Helmholtz |
| **Shifted-Laplacian / Helmholtz MG (PETSc/hypre)** | **Best theoretical path** — not installed |
| Sparse direct SuperLU | Works, far too slow at target ppc |
| GPU matrix-free Krylov | N/A — no GPU |

---

## 4. Prototype validation stages

### Level 1 — straight PEC guide, B=0
- Direct SuperLU @ 10 ppc: residual \(10^{-14}\), centered-source energy R/L ≈ **1.055** → operator + PML qualitatively OK.
- Jacobi/ILU GMRES: preconditioned residual drops to ~\(10^{-5}\) while **true** residual remains ~0.05–0.08 after hundreds of iters → **misleading**.

### Level 2 / 3 — horns / full 91-bulb @ 10 ppc (direct)
- Linear solve succeeds (resid ~\(10^{-13}\)).
- Port centers fixed to cell frame (array @ \(n_x/2,n_y/2\)).
- **|Hz|² line proxy ≠ Meep mode-overlap S-parameters.** Cannot claim ≤1 dB agreement.
- Trusted Meep 25 ppc 0.10-fs: P11≈0.168, P1→P2≈0.0137. Proxy ratios at 10 ppc are not comparable on the same scale; Level-3 **physics accuracy vs Meep is NOT validated** (proxy limitation + coarse mesh).
- Per user gate (“do not proceed if multi-dB errors”): **no green light to treat 25/35 iterative timing as physically certified.** Performance conclusions below still stand from solver behavior + direct scaling.

### B ≠ 0 (mandatory smoke)
- Level3 @ 10 ppc, `bias_a=0.05`: solve OK; nnz rises (9-point gyrotropic); port proxies become **asymmetric** vs B=0 (P2/P6 and P3/P5 split). Off-diagonal \(\rho\) is active. Not compared to a trusted Meep B≠0 six-port column (none located at matching ppc).

---

## 5. Performance (measured + projected)

### Measured

| Case | Factor / prec | Solve / sub | Wall 1st RHS | Iters | True resid |
|---|---:|---:|---:|---:|---:|
| L1 direct 10 ppc | ~274 s (spsolve) | — | ~274 s | — | 1e-14 |
| L3 direct 10 ppc | **148 s** (splu) | **0.15 s** | **149 s** | — | 1e-13 |
| L3 2nd RHS (reuse LU) | 0 | **0.16 s** | — | — | 1e-13 |
| L1 BiCGSTAB+shifted 10 ppc | ~0.03 s | 87 s | 87 s | 800 | **0.22 FAIL** |
| L1 GMRES+ILU0 10 ppc | 0.67 s | 324 s | 324 s | 800 | **0.051 FAIL** |

Material/operator @ full geometry (no solve):

| ppc | material | operator | matvec | N |
|---:|---:|---:|---:|---:|
| 25 | 1.9 s | 2.1 s | 49 ms (21/s) | 2.10M |
| 35 | 4.5 s | 3.1 s | 63 ms (16/s) | 4.12M |

### Projection to 50 ppc (do not fudge)

From L3 SuperLU 10 ppc factor \(T_{10}=148\,\mathrm{s}\), assume 2D fill \(T \propto N^{1.5}\):

| ppc | N / N₁₀ | Projected **one-source** factor | vs 3–5 min |
|---:|---:|---:|---|
| 25 | 6.25 | **~71 min** | miss |
| 35 | 12.25 | **~3.3 h** | miss |
| **50** | **25** | **~9.5 h** | **NO-GO** |

If a *convergent* Helmholtz MG needed only ~100 matvecs @ 50 ppc (~4× slower matvec than 25 ≈ 0.2 s): ~20 s matvecs — **would** hit the target. That preconditioner is **not available** and scipy ILU/Jacobi do **not** provide it.

**GO criterion:** projected 50 ppc ≤10 min with path to 3–5 → **FAIL** on available stack.  
**MARGINAL 10–20 min:** not met.  
**NO-GO >20 min or exploding iters:** **met** (hours direct; iterative does not converge).

---

## 6. Multiple RHS / six ports

Same \(A(\rho,\omega)\) for all six sources.

| Work | Reusable? | Evidence |
|---|---|---|
| Sparse operator | YES | rebuild only if ρ/ω changes |
| SuperLU numeric factor | YES | 2nd RHS substitution **0.16 s** vs **148 s** factor @ 10 ppc |
| ILU/Jacobi | YES if used | irrelevant — does not converge |
| Krylov recycling | not available | — |

**Six-source estimate @ 10 ppc:** ≈ 149 + 5×0.16 ≈ **150 s** (factor dominates).  
**@ 50 ppc (direct):** still ≈ one factorization (**~hours**); substitutions negligible.

---

## 7. GPU

No NVIDIA GPU in this WSL environment. No CuPy. **GPU path not testable; would not change the CPU NO-GO by itself without a working Helmholtz precond.**

---

## 8. Answers to the 12 questions

1. **Scalar Hz?** **Yes** (with \(\rho=\varepsilon^{-1}\), including gyrotropic).
2. **50 ppc unknowns?** **8.4M** scalar.
3. **Sparse direct realistic?** Memory maybe; **time no** (~9.5 h projected).
4. **Best iterative method that worked?** **None** among scipy Jacobi/ILU0 + GMRES/BiCGSTAB. Direct SuperLU works but too slow.
5. **Accuracy vs Meep @ 25 ppc?** **Not established.** Level-3 run was 10 ppc with |Hz|² proxy, not mode-overlap S-params. No certified ≤1 dB comparison.
6. **Iterations @ 25 / 35?** Iterative **did not converge** at 10 ppc already; 25/35 Krylov solves not productive to run.
7. **Measured time @ 25 / 35?** Operator+matvec only (table above). Full direct projected **~71 min / ~3.3 h**.
8. **Projected one-source @ 50 ppc?** **~9.5 hours** (direct); iterative unknown/nonconvergent with current preconds.
9. **Projected six-source with reuse?** ≈ **same as one-source** for direct (factor once). Still hours @ 50 ppc.
10. **B ≠ 0 supported?** Operator runs; port asymmetry appears. Not Meep-validated.
11. **GPU improve projection?** **No GPU present.** Even with GPU, need a converging Helmholtz precond first.
12. **TRUE 50 ppc / 3–5 min?** **Clearly infeasible with this environment’s solvers.** Formulation is **plausible only after** installing a real Helmholtz multigrid / shifted-Laplacian stack (CPU). Not “GPU-only.”

---

## Final decision

### **NO-GO**

Do not invest further in Meep FDTD tuning (already ruled out) or in scipy ILU/GMRES for this Helmholtz problem.

### Single next engineering action (if pursuing frequency-domain)

**Install PETSc + hypre (or equivalent) and implement a complex-shifted-Laplacian / Helmholtz multigrid preconditioner for the scalar-Hz operator; then validate proper port mode-overlap S-parameters against trusted Meep 25 ppc (0.10 fs) before any 50 ppc timing claim.**

If that still fails to reach ≤10 min at 50 ppc: next algorithmic change is **nonuniform / adaptive FEM (or domain-decomposition)** to cut DOFs in horns/PML while keeping plasma resolution — not another uniform Cartesian FDFD with weak preconds. Reduced-order models are a later optimization, not the first fix.
