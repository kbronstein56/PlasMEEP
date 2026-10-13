# Direct DtN multimode adjoint vs FD

Grade C, design `opt_single_uniform_M_Bp0.0500`, B=0.05, measure=projected.
J0 = 1.373121e-01. FD rel bar = 0.001.

## Forward ≡ adjoint (machine check)

| rod | g_adj | g_fwd | rel | PASS |
|-----|------:|------:|----:|:----:|
| 0 | -3.9405e-01 | -3.9405e-01 | 7.072e-14 | True |
| 45 | -2.6424e-02 | -2.6424e-02 | 1.488e-13 | True |
| 30 | -3.3058e-02 | -3.3058e-02 | 3.165e-13 | True |

**FORWARD_EQ_ADJOINT: PASS**

## Centered FD

| rod | g_adj | best FD | rel | PASS |
|-----|------:|--------:|----:|:----:|
| 0 | -3.9405e-01 | -3.9405e-01 | 1.402e-06 | True |
| 45 | -2.6424e-02 | -2.6424e-02 | 8.701e-07 | True |
| 30 | -3.3058e-02 | -3.3057e-02 | 4.185e-06 | True |

**DIRECT_DTN_ADJOINT_FD: PASS**

