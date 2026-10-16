# BROADBAND_DIFFERENTIABLE_OBJECTIVE

**Frozen before** any production broadband optimization.
Does **not** change the single-frequency frozen objective
(`MULTIMODE_SINGLE_FREQ_OBJECTIVE_FROZEN.md`).

## Target band

\[
4.0\,\mathrm{GHz} \le f \le 6.0\,\mathrm{GHz}.
\]

Mode content in this band (parallel-plate clear width \(W=2a\)):

| Mode | Cutoff | In 4–6 GHz |
|------|--------|------------|
| TE0 | 0 | propagating |
| TE1 | \(\approx 3.747\) GHz | propagating |
| TE2 | \(\approx 7.49\) GHz | evanescent |

→ The same **12-channel TE0+TE1** DtN definition applies at every sample in 4–6 GHz.
No cutoff-crossing treatment is required *inside this band* (unlike 3.75–3.95 GHz).

## Per-frequency score (unchanged physics)

At each \(f\), use the **same** validated TE0-drive multimode score \(J(f;s)\):

\[
\begin{aligned}
J_j(f)
&= +1.0\,T_{\mathrm{des,TE0}}
 -2.0\,T_{\mathrm{rev,TE0}}
 -1.5\,T_{\mathrm{thr,TE0}}
 -1.0\,T_{\mathrm{other,TE0}}
 -0.25\,R_{\mathrm{TE0}}
 -1.0\,T_{\mathrm{TE1,total}},\\
J(f;s) &= \mathrm{mean}_{j=0\ldots5}\,J_j(f).
\end{aligned}
\]

Weights identical to the single-frequency freeze.

## Broadband aggregation (robust, differentiable)

**Do not** optimize a bare arithmetic mean alone (allows deep valleys between samples).

Define a **soft-min** (log-sum-exp) of the negated score so that maximizing soft-min
improves the worst frequency, with a small mean term for average performance:

\[
\begin{aligned}
J_{\mathrm{softmin}}(s)
&= -\tau\log\Bigl(\frac{1}{N_f}\sum_{i=1}^{N_f}
\exp\bigl(-J(f_i;s)/\tau\bigr)\Bigr),\\[4pt]
J_{\mathrm{BB}}(s)
&= (1-\alpha)\,J_{\mathrm{softmin}}(s)
 + \alpha\,\underbrace{\frac{1}{N_f}\sum_i J(f_i;s)}_{J_{\mathrm{mean}}}.
\end{aligned}
\]

### Frozen hyperparameters (before production)

| Symbol | Value | Role |
|--------|------:|------|
| \(\tau\) | \(0.05\) | Soft-min sharpness (same units as \(J\)) |
| \(\alpha\) | \(0.25\) | Mean-performance weight |
| \(1-\alpha\) | \(0.75\) | Worst-frequency emphasis |

Smaller \(\tau\) → closer to \(\min_i J(f_i)\). \(\tau=0.05\) is chosen so that a
\(0.1\) gap in \(J\) between worst and next is amplified without making the
softmax numerically pathological for expected \(J\in[-1,1]\).

**These values are frozen for the sanity campaign.** Retuning after seeing
production optima is out of scope unless a new freeze document is written.

## Analytic gradient (chain rule)

Let \(w_i=\partial J_{\mathrm{BB}}/\partial J(f_i)\). Then

\[
\frac{\partial J_{\mathrm{BB}}}{\partial s}
= \sum_{i=1}^{N_f} w_i\,\frac{\partial J(f_i)}{\partial s},
\]

with

\[
\begin{aligned}
p_i &= \frac{\exp(-J_i/\tau)}{\sum_k\exp(-J_k/\tau)},\\
\frac{\partial J_{\mathrm{softmin}}}{\partial J_i} &= p_i,\\
w_i &= (1-\alpha)\,p_i + \alpha/N_f.
\end{aligned}
\]

Each \(\partial J(f_i)/\partial s\) is the **existing** frequency-domain adjoint
(efficient TE0 form). **No JAX required** for this aggregation.

## Initial vs audit grids

| Grid | Spacing | Role |
|------|---------|------|
| Objective (initial) | 4.0, 4.2, …, 6.0 GHz (**11** pts) | \(J_{\mathrm{BB}}\) during opt |
| Dense audit | 4.0–6.0 every **0.05** GHz (**41** pts) | Find hidden dips |
| Final FD validation | ≤50 MHz; 25 MHz near sharp features | Accept/reject design |

Adaptive enrichment: after a short opt stage, add worst audit frequencies to the
objective set and re-optimize (see `ADAPTIVE_BROADBAND_SAMPLING_REPORT.md`).

## What is not claimed

- Soft-min does not prove continuous-band optimality; density of the audit grid does.
- This freeze does **not** authorize changing single-frequency weights.
