# C6 tied-variable FD spot check (direct DtN)

Grade C, design `opt_single_uniform_M_Bp0.0500`, B=0.05.
J0=1.373121e-01. FD rel bar=0.001.

| orbit | size | g_adj | best FD | rel | PASS |
|------:|-----:|------:|--------:|----:|:----:|
| 0 | 1 | -1.6861e+00 | -1.6860e+00 | 7.157e-05 | True |
| 1 | 6 | -2.6029e-01 | -2.6029e-01 | 1.416e-07 | True |
| 5 | 6 | -3.0611e-02 | -3.0611e-02 | 4.080e-09 | True |

**C6_TIED_FD: PASS**


## Note

Center orbit (rod 85) needed FD step ≤ 3×10⁻⁴ to meet the 10⁻³ relative bar;
at h=10⁻³ the relative error was 7×10⁻³ (truncation). Six-rod orbits pass at
h=10⁻³ already. Adjoint matches FD; production may proceed.
