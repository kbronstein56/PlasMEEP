# Circulator objective derivation

Frozen **before** any production optimization. Branch `agent/eigenmode-ports`,
Goal-B closeout `df0cff1`.

## Operating point (initial)

    f = 3.85 GHz
    B = +0.05 T

## Port convention (verified)

Ports are indexed 0…5 counterclockwise, port 0 nearest +x
(`FULL91_GEOMETRY_AUDIT.md`).

Validated FULL91-F `port_power[i,j]` = outward guide-normal Poynting flux at
port `i` for source `j`:

| B | source 0: CCW (port 1) | CW (port 5) | through (port 3) |
|---|------------------------|-------------|------------------|
| +0.05 T | 6.0e-4 | 2.4e-4 | 7.8e-3 |
| −0.05 T | 2.4e-4 | 6.0e-4 | 7.8e-3 |
| 0 | 9.4e-3 | 9.4e-3 | 2.7e-2 |

At **+B** the counterclockwise adjacent port wins over the clockwise
adjacent port. At **−B** the preference reverses. Through still dominates
the unoptimized device; the optimization must raise the desired adjacent
path and suppress through / reverse / leakage.

Classification for source `j` (desired circulation = CCW at +B):

    desired  = (j + 1) mod 6
    reverse  = (j − 1) mod 6
    through  = (j + 3) mod 6
    leak_a   = (j + 2) mod 6
    leak_b   = (j + 4) mod 6
    source   = j

## Observable

Guide-normal outward flux on the production monitor line (vacuum feed),
identical to the validated `guide_normal_flux` / FULL91 `port_power`:

    P_{ij} = ∫ (S · n̂_i) ds
    S_x = ½ Re(E_y H_z*),   S_y = −½ Re(E_x H_z*)

In the vacuum feed this equals

    S·n̂ = ½ Re[ (−i/ω) (∇H_z · n̂) H_z* ]

Monitors sit in vacuum: **no explicit material dependence of P on s**.
All design dependence is through the field `x(s)`.

## Normalization (fixed incident scale)

**Rejected:** dividing by the design-dependent driven-port net flux
`|P_{jj}(s)|`. That can improve a “normalized” figure by reducing accepted
power rather than improving circulation.

**Chosen:** freeze a single positive scale `P_ref` from the **uniform**
production design `s_i ≡ 1` at the same `(f,B)` and the same fixed source
vectors:

    P_ref = (1/6) ∑_j (−P_{jj}(s≡1))

`P_ref` is then **never updated**. Source RHS vectors are fixed (barycentric
injection of the flat feed mode, same as FULL91). Absolute powers with a
frozen scale are equivalent for ranking designs; reporting uses

    t_{ij} = P_{ij} / P_ref

For receiving ports `i≠j`, positive `t_{ij}` is outward transmission.
Source-port outward flux is negative when power enters the cavity; we use

    a_j = (−P_{jj}) / P_ref

as a dimensionless **accepted-power** score (not a true modal |S₁₁|² without
field-subtraction). A mild reward on `a_j` encourages power coupling without
pretending to be a calibrated reflection coefficient.

True S-parameter reflection via straight-guide field subtraction remains
available for **reporting**, not for the optimizer.

## Objective (maximize)

Linear power quantities only (no dB inside the optimizer):

    J = (1/6) ∑_{j=0}^{5} [
          w_des  t_{des,j}
        − w_rev  t_{rev,j}
        − w_thr  t_{thr,j}
        − w_leak (t_{leak_a,j} + t_{leak_b,j})
        + w_acc  a_j
    ]

Frozen weights (chosen from the unoptimized physics: through is the main
competitor; reverse isolation is the circulation signature):

| weight | value | role |
|--------|-------|------|
| `w_des` | 1.0 | raise CCW adjacent |
| `w_rev` | 2.0 | suppress CW adjacent |
| `w_thr` | 1.5 | suppress through |
| `w_leak` | 1.0 | suppress ±120° ports |
| `w_acc` | 0.25 | mild accepted-power reward |

These weights are **not** retuned after seeing optimized results.

The optimizer minimizes `−J` (L-BFGS-B).

## Six-port completeness

Every source `j` appears with equal weight. A design that improves only
port 0 cannot win against a rotationally consistent circulator under
sixfold-tied variables.

## Reporting metrics (dB after the fact)

After optimization, report for each source (and averages / worst-port):

    insertion   = 10 log10(max(t_des, ε))
    isolation   = 10 log10(max(t_des, ε) / max(t_rev, ε))
    through_dB  = 10 log10(max(t_thr, ε))
    leak_dB     = 10 log10(max(t_leak_a + t_leak_b, ε))
