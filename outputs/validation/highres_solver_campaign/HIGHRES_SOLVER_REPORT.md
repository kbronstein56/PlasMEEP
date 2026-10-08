# High-Resolution Solver Campaign Report

**Campaign window:** 2026-09-21 (autonomous, ~decision-tree driven)  
**Workspace:** `outputs/validation/highres_solver_campaign/`  
**Scripts:** `scripts/validation/highres_solver_campaign/`  
**Production `plasmeep_mpi`:** not modified. No commits/pushes.

---

## Executive verdict

| Goal | Result |
|---|---|
| Uniform 50 ppc FDFD ≤10 min | **NO-GO** (PETSc unavailable; scipy Helmholtz PCs fail) |
| Locally refined ~0.2 mm FEM ≤10 min | **TIMING TARGET ACHIEVED** |
| Locally refined FEM 3–5 min | **ACHIEVED with margin** (measured **~11 s** one-source) |
| Six-port candidate ≤10 min | **ACHIEVED** (measured estimate **~12 s** with LU reuse) |
| Accurate Meep S-parameter replacement | **NOT YET** — |Hz|² proxy ≠ mode-overlap S; gyrotropic FEM still isotropic-proxy |

**Bottom line:** Cutting unknowns with local refinement + SuperLU on ~0.66M CG1 DOFs is the winning path for *speed*. It is **not** yet proven as a drop-in Meep replacement for inverse design until port modal S-parameters and full anisotropic \(\rho=\varepsilon^{-1}\) are validated.

---

## Phase 0 — Baseline preserved

- Snapshot of `fdfd_feasibility/` prototype + `FDFD_GO_NO_GO.md` under `phase0/`.
- Reproduced Level1 SuperLU @ 10 ppc: true resid ~7e-14, R/L ≈ 1.055 (matches prior baseline).
- Git branch at start: `agent/eigenmode-ports` (untracked validation outputs only).

---

## Phase 1 — Solver environment

| Package | Status |
|---|---|
| Isolated venv `.conda_envs/fdfd_solver` | **YES** (pip) |
| numpy/scipy/matplotlib/h5py/pyamg/scikit-fem/meshio | **YES** |
| PETSc / petsc4py / hypre | **NO** |
| GPU | **NO** |
| System g++ / usable clang | **NO** (`libtinfo.so.5` missing for portable LLVM) |
| gmsh Python wheel | Present but **unusable** (`libGLU.so.1` missing) |
| conda-forge / prefix.dev | **Blocked** (proxy CONNECT 403) |

**Decision (gate):** PETSc install dead-end within budget → skip uniform PETSc-FDFD scaling → **pivot to locally refined FEM**.

Manifest: `phase1/env_manifest.json`.

---

## Phases 2–4 — Uniform FDFD / Helmholtz PC screen

Reused exact scalar-Hz operator from `fdfd_hz_prototype.py`.

| Method | Case | True residual | Notes |
|---|---|---:|---|
| SuperLU direct | L1 10 ppc | ~1e-14 | Works; too slow at high uniform ppc |
| Shifted-Laplacian + geometric MG + GMRES | L1 10 ppc | timeout / failed | Preconditioned residual fell; true residual not driven (timeout left poor solution) |
| Shifted ILU + GMRES | L1 10 ppc | did not finish ≤400 s | Setup alone too heavy |
| Shifted PyAMG + GMRES | L1 10 ppc | did not finish ≤400 s | Setup alone too heavy |

**Uniform-grid GO gate:** projected 50 ppc still hours with direct; iterative path not converging with available tools → **reject uniform Cartesian FDFD** for the inner loop.

---

## Phases 5–6 — Locally refined mesh (no gmsh)

Delaunay generator: `fem_delaunay_mesh.py` (geometry from `PMMI` hex train + `full_horns`).

| Grade | h_crit | Nodes (CG1 DOFs) | vs 8.4M @50 ppc | Median edge | Max air edge |
|---|---:|---:|---:|---:|---:|
| FEM-L | 0.8 mm | 89,062 | **1.1%** | 0.80 mm | 8.5 mm |
| FEM-M | 0.4 mm | 234,307 | **2.8%** | 0.40 mm | 4.9 mm |
| FEM-H | **0.2 mm** | **661,890** | **7.9%** | **0.20 mm** | 3.5 mm |

**Answer to “can local 0.2 mm beat 8.4M?”:** **Yes — ~0.66M unknowns (~12.7× fewer).**

---

## Phases 7–9 — FEM solve, timing, multi-RHS

Solver: `fem_hz_solver.py` — P1 weak form of \(\nabla\cdot((1/\varepsilon)\nabla H_z)+k_0^2 H_z=0\) with PML stretch; PEC nodes stamped along horn polylines; scipy `splu`.

### Timing (MEASURED)

| Case | DOFs | Factor | 1st RHS wall | 2nd RHS | Est. 6-port |
|---|---:|---:|---:|---:|---:|
| FEM-L level1 | 89k | 27 s | 28 s | — | — |
| FEM-L level3 B=0 | 89k | 55 s | 55 s | 0.05 s | ~56 s |
| FEM-M level3 B=0 +PEC | 234k | 28 s | **30 s** | 0.09 s | ~30 s |
| **FEM-H level3 B=0 +PEC** | **662k** | **6.4 s** | **11.3 s** | **0.21 s** | **~12.3 s** |
| FEM-H level3 B=0.05 (ε⊥ proxy) | 662k | 9.4 s | 14.2 s | — | — |

Peak RSS FEM-H: ~1.9 GiB.

### Level1 physics check
FEM-L vacuum domain R/L ≈ **1.002**, true resid ~1e-13 → assembly/PML/direct path sane.

### Port proxies (NOT Meep S-parameters)
FEM-H B=0 |Hz|² proxy (norm to port0):  
P2≈0.013, P3≈0.0004, P4≈0.016, P5≈0.0004, P6≈0.013  

Trusted Meep 25 ppc 0.10-fs: P11≈0.168, P1→P2≈0.0137 (linear power).  
Proxy scale differs; **do not claim ≤1 dB Meep agreement**.

### Multi-RHS
Same \(A\): factorization reusable. Additional RHS ≈ **0.2 s** @ FEM-H. Six-port ≈ one factor + 6×substitution ≈ **12 s**.

---

## Phase 10 — B ≠ 0

Ran FEM-H with `bias_a=0.05` using **isotropic ε_⊥ only** (not full gyrotropic \(\rho=\varepsilon^{-1}\) tensor FEM).  
Solve stable (resid ~1e-12). Port proxies change vs B=0, but this is **not** a validated Faraday gyrotropic operator. Full anisotropic FEM assembly remains required.

---

## Phase 12 — Comparison table

| Method | Resolution | Unknowns | B support | vs Meep S | Setup | 1st RHS | Extra RHS | ~6-port | Mem | Status |
|---|---|---:|---|---|---:|---:|---:|---:|---:|---|
| Meep FDTD | 25 ppc | — | yes | trusted | — | **~6.1 min** | ~same×6 | ~37 min | — | baseline |
| Meep FDTD | 50 ppc | — | yes | trusted | — | **~72 min** | ×6 | hours | — | too slow |
| SciPy FDFD direct | 10 ppc unif. | 0.34M | tensor OK | proxy | ~148 s | ~149 s | 0.16 s | ~150 s | ~1 GiB | slow @hi-res |
| SciPy weak iterative | 10 ppc | 0.34M | — | — | — | fail | — | — | — | stalled |
| PETSc Helmholtz | — | — | — | — | — | — | — | — | — | **unavailable** |
| **FEM-H local 0.2 mm** | local 0.2 mm | **0.66M** | ε⊥ proxy | proxy only | **6.4 s** | **11.3 s** | **0.21 s** | **~12 s** | ~2 GiB | **speed GO; accuracy WIP** |

---

## Answers to required questions

1. **Why SciPy FDFD failed to scale?** Indefinite Helmholtz; Jacobi/ILU/shifted-MG did not reduce **true** residual; SuperLU ~\(N^{1.5}\) → hours at 50 ppc.
2. **Did PETSc + Helmholtz PC fix it?** **No — PETSc could not be installed** in this environment.
3. **Best true residual?** Direct SuperLU/FEM: **~1e-13**. Iterative: not usefully converged.
4. **Fastest accurate uniform-grid solve?** None that is both accurate and fast at ≥25 ppc. Direct @10 ppc ~2.5–4.5 min.
5. **Uniform 50 ppc ≤10 min?** **No.**
6. **Uniform 50 ppc 3–5 min?** **No.**
7. **FEM unknowns for ~0.2 mm interfaces?** **~662k** (7.9% of 8.4M).
8. **FEM accurate vs Meep?** **Not demonstrated** (proxy ≠ S-params; incomplete PEC geometry fidelity; B≠0 not full tensor).
9. **Fastest locally-high-res one-source?** **FEM-H: 11.3 s measured.**
10. **1st vs subsequent RHS?** 11.3 s vs **0.21 s**.
11. **Six-port candidate?** **~12.3 s measured estimate.**
12. **B≠0?** Solves with ε⊥ proxy; **full gyrotropic FEM not done.**
13. **Measured vs projected?** Times above are **measured**. Uniform 50 ppc projections from prior campaign.
14. **Achieved ≤10 min / 3–5 min / ≤10 min six-port?** **Yes / Yes / Yes** on FEM-H **timing**. Accuracy gate still open.
15. **What blocks Meep replacement now?** Modal port S-parameter extraction + anisotropic \(\rho\) FEM + mesh/PEC convergence vs Meep.
16. **Would GPU help?** Not required for the FEM-H timing win; no GPU present.
17. **Advanced Helmholtz DD?** Useful for uniform grids; **local FEM already crossed the time target.**
18. **Single best next engineering step?**  
    **Implement proper numerical-mode port overlap S-parameters on the FEM-H solution and validate P11 & P1→P2…P6 against trusted Meep 25 ppc 0.10-fs within ≤1 dB; then add full gyrotropic \(\rho=\varepsilon^{-1}\) element matrices for B≠0.**

---

## Decision-tree outcome

```
PETSc install fails → FEM
Uniform iterative fails → FEM
FEM-H timing ≤5 min → TARGET ACHIEVED (speed)
FEM Meep S accuracy → NOT YET → next engineering step above
```

---

## Reproduce

```bash
# venv packages already under .conda_envs/fdfd_solver
export PYTHONPATH="$PWD/.conda_envs/fdfd_solver/lib/python3.13/site-packages:$PWD/scripts/validation:$PWD/scripts"

# meshes
/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python -u \
  scripts/validation/highres_solver_campaign/fem_delaunay_mesh.py --grade all

# FEM-H solve
/home/daq_user/miniconda3/envs/plasmeep_mpi/bin/python -u \
  scripts/validation/highres_solver_campaign/fem_hz_solver.py \
  --grade FEM-H --level level3 --method direct --second-rhs \
  --json-out outputs/validation/highres_solver_campaign/phase5_fem/FEM-H_level3_pec.json
```

Artifacts: JSON/logs under `phase0/`, `phase1/`, `phase3/`, `phase5_fem/`.
