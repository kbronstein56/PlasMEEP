# Literature benchmark matrix

A row is a paper reproduction only when the published setup can be rebuilt and a published number can be compared. Figure-only papers are `NOT_REPRODUCIBLE_FROM_PAPER`.

| Paper | Benchmark | Our test | Published value | Our value | Difference | Reproducibility | Status |
|-------|-----------|----------|-----------------|-----------|------------|-----------------|--------|
| Hughes, Williamson, Minkov, Fan, ACS Photonics 6, 3010–3016 (2019), DOI 10.1021/acsphotonics.9b01238 | Forward-mode FDTD sensitivity | conceptual only | no digit table in the available text | — | — | `NOT_REPRODUCIBLE_FROM_PAPER` | conceptual |
| Hughes, Minkov, Williamson, Fan, ACS Photonics 5, 4781–4787 (2018) | Nonlinear adjoint inverse design | not used as a number | device figures | — | — | `NOT_REPRODUCIBLE_FROM_PAPER` | conceptual |
| Minkov et al., ACS Photonics 7, 1729–1741 (2020) | Photonic-crystal autodiff | not used as a number | band figures | — | — | `NOT_REPRODUCIBLE_FROM_PAPER` | conceptual |
| Hamid and Cooray, Adv. Electromagnetics, 30 Dec 2016 | Table I elliptic ferrite echo widths, axial ratio 2, `k0 a = 1` | `hamid_ellipse_benchmark.py`; criteria frozen in `PUBLICATION_BENCHMARK_PASS_CRITERIA.md` | 0.13781 … 0.77898 | finest FEM 0.00062 … 0.0145 (μ_zz = 1 duality as written) | relative error ~0.90–0.998 | reconstructable geometry; Maxwell duality / μ_zz assignment not unique from the paper alone | **FAIL** as paper reproduction. Mie RCS calibration at the same extractor had max relative error `0.021` vs frozen `0.005`, so the ellipse row is also not interpretable under the frozen gate. |
| Platzman/Ozaki; Seshadri magnetized plasma cylinder | analytic magnetized cylinder | already in gyrotropic suite | classic analytic family | FEM vs independent `T_n` worst `1.0e-4` dB | — | analytic family already covered | not a new reproduction |

## Hamid attempt, recorded numbers

Frozen criteria required Mie calibration relative error `≤ 0.005` before trusting the ellipse. Finest circular calibration (`k0 a = 1`, `ε = 4`, h = 0.06):

| φ | analytic σ/λ | FEM | relative |
|--:|-------------:|----:|---------:|
| 0° | 1.26017 | 1.25621 | 0.00315 |
| 45° | 0.72367 | 0.71897 | 0.00649 |
| 90° | 0.11413 | 0.11374 | 0.00343 |
| 180° | 0.01767 | 0.01804 | 0.02105 |

Calibration misses. Ellipse rows are therefore not a PASS even before comparing to Table I. The Table I comparison itself remains off by nearly 100% under the μ = 1 duality insertion written into the criteria. A trial with `μ_zz = 10` changed the widths but did not recover the published table; that trial is diagnostic only and was not used to retune the frozen criteria.

This FAIL does not overturn the independent gyrotropic-cylinder probe agreement already frozen in the gyrotropic campaign.
