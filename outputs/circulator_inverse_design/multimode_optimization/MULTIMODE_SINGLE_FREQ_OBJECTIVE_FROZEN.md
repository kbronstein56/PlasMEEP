# MULTIMODE_SINGLE_FREQ_OBJECTIVE_FROZEN

**Frozen before** the first production multimode inverse-design run.
Do **not** retune weights or change the S definition after seeing results.

Base package: `multimode_ports/MULTIMODE_PORT_PACKAGE_FROZEN.md` (commit `818c982`).
Adjoint validation: `multimode_ports/DIRECT_DTN_ADJOINT_FD_VALIDATION.md`.

## Operating point

| Quantity | Value |
|----------|-------|
| Frequency | 3.85 GHz |
| Bias | \(B_z = +0.05\) T |
| Mesh (production) | FULL91-C (~7 min/eval+grad); M for transfer + short refine (M full opt not practical in budget) |
| Port package | direct DtN 12×12, `measure=projected` |
| Modes in S | TE0 + TE1 propagating (2 per port) |
| Excitation | TE0 incidence on each of 6 horns |
| Design variables | C6 orbits: 1 center + 15 six-rod = **16** `q_m` |
| Bounds | \(0.05 \le s_i = q_m \le 1.80\) (validated in adjoint FD) |

## Scattering definition (unchanged)

\[
S = I + 2i\,Q,\qquad
Q_{mn}=\sqrt{\beta_m\beta_n}\,c_m^\top A^{-1} c_n.
\]

Channel indexing: column/row \(2j\) = TE0 on port \(j\), \(2j+1\) = TE1 on port \(j\).

## Objective (validated in `adjoint_direct_dtn.py`)

For TE0 incidence on source port \(j\) (column \(n=2j\)):

\[
\begin{aligned}
J_j &=
 +1.0\, T_{\mathrm{des,TE0}}
 -2.0\, T_{\mathrm{rev,TE0}}
 -1.5\, T_{\mathrm{thr,TE0}}
 -1.0\, T_{\mathrm{other,TE0}}
 -0.25\, R_{\mathrm{TE0}}
 -1.0\, T_{\mathrm{TE1,total}} \\[4pt]
J &= \mathrm{mean}_{j=0\ldots5}\, J_j
\end{aligned}
\]

with

| Symbol | Definition |
|--------|------------|
| \(T_{\mathrm{des,TE0}}\) | \(\lvert S_{2\,\mathrm{des},\,2j}\rvert^2\), \(\mathrm{des}=(j+1)\bmod 6\) |
| \(T_{\mathrm{rev,TE0}}\) | \(\lvert S_{2\,\mathrm{rev},\,2j}\rvert^2\), \(\mathrm{rev}=(j-1)\bmod 6\) |
| \(T_{\mathrm{thr,TE0}}\) | \(\lvert S_{2\,\mathrm{thr},\,2j}\rvert^2\), \(\mathrm{thr}=(j+3)\bmod 6\) |
| \(T_{\mathrm{other,TE0}}\) | sum of \(\lvert S_{2p,2j}\rvert^2\) over the two remaining ports |
| \(R_{\mathrm{TE0}}\) | \(\lvert S_{2j,2j}\rvert^2\) |
| \(T_{\mathrm{TE1,total}}\) | \(\sum_{p=0}^{5}\lvert S_{2p+1,\,2j}\rvert^2\) (all outgoing TE1) |

Weights (frozen):

```
W_DES=1.0  W_REV=2.0  W_THR=1.5  W_LEAK=1.0  W_REFL=0.25  W_TE1=1.0
```

## What is **not** used

- Pref / soft-source objective
- Inset \(H_z/\partial_n\) Ain/Bout extractor
- Fundamental-only 6×6 S
- Arclength staircase measure

## C6 tying

Verified by `sixfold_orbits.py` / `SIXFOLD_DESIGN_ORBITS.md`:

- \(s_i = q_m\) for \(i\in\) orbit \(m\)
- \(\partial J/\partial q_m = \sum_{i\in m}\partial J/\partial s_i\)

Reuse the validated direct-DtN adjoint; orbit mask = OR of rod masks.

## Engineering metrics (reported, not retuned as gates)

Descriptive reference levels for judging designs (not new pass/fail bars):

- IL \(= -10\log_{10} T_{\mathrm{des}}\) (dB)
- isolation \(= 10\log_{10}(T_{\mathrm{des}}/T_{\mathrm{rev}})\) (dB)
- RL \(= -10\log_{10} R\) (dB)
- Illustrative: IL < 3 dB, isolation > 10 dB, RL > 10 dB

Frozen numerical gates that **do** apply to post-opt validation:

| Gate | Bar |
|------|-----|
| Onsager \(r_F=\|S(+B)-S(-B)^\top\|_F/\|S\|_F\) | ≤ 5×10⁻³ |
| Power balance \(\lvert 1-\sum\|S\|^2-P_{\mathrm{mat}}\rvert\) | ≤ 0.05 per channel |

## Starts (Phase 2 baselines)

A. uniform \(s_i\equiv 1\)  
B. Pref sixfold peak `opt_single_uniform_M_Bp0.0500`  
C. Pref broadband `opt_broadband_C_Bp0.0500`

These Pref designs are **warm starts only**, not final designs under the multimode S.

## Optimizer

- Bounded L-BFGS-B on the 16 orbit variables
- Same objective and bounds for all warm starts
- Keep best-so-far \(q\) separately from the final iterate
