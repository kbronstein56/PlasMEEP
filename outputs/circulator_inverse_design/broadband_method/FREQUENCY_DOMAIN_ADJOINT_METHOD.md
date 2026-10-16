# FREQUENCY_DOMAIN_ADJOINT_METHOD

**Status:** documentation of the *already validated* PlasMEEP multimode derivative.
This is **not** a new physics derivation. The frozen single-frequency objective and
port package (`MULTIMODE_PORT_PACKAGE_FROZEN.md`, `MULTIMODE_SINGLE_FREQ_OBJECTIVE_FROZEN.md`)
are unchanged.

## 1. What we already have

At each frequency \(f\) and design \(s\), the FEM system is

\[
A(f,s)\,x = b.
\]

The direct DtN multimode S-matrix is

\[
S = I + 2i\,Q,\qquad
Q_{mn}=\sqrt{\beta_m\beta_n}\,c_m^\top A^{-1} c_n.
\]

The TE0-incidence objective \(J(s)\) is a smooth function of \(\{|S_{mn}|^2\}\)
(desired / reverse / through / other TE0, reflection, and **all outgoing TE1**).

Differentiating through the linear solve gives the **frequency-domain adjoint**:

\[
\frac{\partial S_{mn}}{\partial s_k}
= -2i\sqrt{\beta_m\beta_n}\,
\bigl(A^{-\top}c_m\bigr)^\top
\Bigl(\frac{\partial A}{\partial s_k}\Bigr)
\bigl(A^{-1}c_n\bigr).
\]

Production code already reuses one SuperLU factor for transpose solves:

```python
lu = splu(A.tocsc())
# forward
U[:, n] = lu.solve(C[:, n])
# transpose-adjoint (same factor — no A.T refactorization)
V[:, m] = lu.solve(C[:, m], trans="T")
```

Validated: forward ≡ adjoint ~1e-13; centered FD on rods/orbits ≤1e-3 relative
(`DIRECT_DTN_ADJOINT_FD_VALIDATION.md`, `C6_TIED_FD_CHECK.md`).

## 2. Method comparison (this problem)

| Method | What it computes | Cost scaling (dominant) | Notes for PlasMEEP |
|--------|------------------|-------------------------|--------------------|
| Finite differences | \(N_\mathrm{param}\) perturbed forward maps | \(O(N_\mathrm{param}\,N_f)\) full assemblies+LU | Used only for spot checks |
| Forward-mode / direct sensitivity | \(A\,(\partial x/\partial s_k)= -(\partial A/\partial s_k)x\) per param | \(O(N_\mathrm{param})\) RHS per frequency after 1 LU | Attractive if \(N_\mathrm{param}\ll N_\mathrm{obj}\) |
| **Frequency-domain adjoint** | one (or few) \(A^\top\lambda=g\) per objective | \(O(N_\mathrm{obj})\) RHS; **independent of \(N_\mathrm{param}\)** after contractions with \(\partial A/\partial s_k\) | **Current production method** |
| JAX / HIPS autograd reverse mode | VJP through the computational graph | Still requires VJP of `solve(A,b)` = adjoint solve | Does **not** remove Maxwell/adjoint work |
| Time-domain / FDTD adjoint | reverse-time adjoint field | Broadband pulse → many frequencies in one pair of runs | Better for *validation* / broadband FoMs; different code path (Meep) |

### Cost sketch for this circulator

Symbols:

- \(N_\mathrm{param}\in\{16,91\}\) (C6 orbits or free rods)
- \(N_s=6\) TE0 source ports (objective averages over them)
- \(N_f\) frequencies in a broadband objective
- \(N_\mathrm{ch}=12\) propagating channels (TE0+TE1), but **incident TE1 is not required for the frozen TE0-drive objective**

After one LU per \((f,s)\):

| Approach | Approx. RHS after LU (per frequency) |
|----------|--------------------------------------|
| FD on 16 orbits | \(2\times16\) full rebuilds (not just RHS) |
| Forward-mode on 16 / 91 | 16 / 91 sensitivity solves (+ field) |
| Baseline adjoint (current `_channel_basis`) | 12 forward + 12 transpose (all channels) + **duplicate** 12 in `direct_S_from_factor` |
| Efficient adjoint (Phase 2 target) | **6 TE0 forward + 6 scalar-objective transpose** |

Reverse-mode AD (JAX/autograd) through `x = A^{-1}b` implements the same
transpose solve; it does not change the asymptotic Maxwell cost.

## 3. Literature map (accurate claims only)

- **Hughes, Williamson, Minkov, Fan, ACS Photonics 2019** (“Forward-Mode Differentiation of Maxwell’s Equations”):
  forward-mode / tangent-linear differentiation of frequency-domain Maxwell;
  competitive when few parameters or many simultaneous objectives. It does **not**
  claim that reverse mode avoids linear solves.

- **Ceviche (HIPS autograd)**: reverse-mode AD for photonic inverse design with
  custom VJPs for sparse solves — mathematically an adjoint method packaged as AD.

For PlasMEEP: we already have an analytic adjoint matching that mathematics.
Migrating to JAX/autograd is a **software** question (maintenance / GPU), not a
new derivative principle.

## 4. Advisor answers (preview)

1. **Is there a frequency-domain adjoint?** Yes — this is the current validated method.
2. **Does AD remove adjoint solves?** No.
3. **Most efficient path for \(N_\mathrm{param}=16\)–\(91\), \(N_s=6\), many \(N_f\)?**
   Frequency-domain adjoint with one LU/frequency, TE0-only forward RHS, and
   scalar-objective transpose RHS (Phase 2), optionally parallelized over \(N_f\)
   (Phase 6).

## 5. What must not change

- Frozen projected measure, DtN S definition, Robin/RHS signs
- Frozen single-frequency weights and \(0.05\le s\le 1.80\) box
- Onsager / power-balance gates for post-opt validation
