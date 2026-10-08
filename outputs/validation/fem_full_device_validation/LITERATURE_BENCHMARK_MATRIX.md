# Literature benchmark matrix

A row is a paper reproduction only when the published setup can be rebuilt and a published number can be compared. Figure-only papers are `NOT_REPRODUCIBLE_FROM_PAPER`. They are conceptual comparisons, not passes.

Our analytic suite already covers Mie, coated T-matrix, Voigt slabs, and Graf clusters. Those are not repeated here as new paper reproductions.

| Paper | Benchmark | Our test | Published value | Our value | Difference | Reproducibility | Status |
|-------|-----------|----------|-----------------|-----------|------------|-----------------|--------|
| Hughes, Williamson, Minkov, Fan, ACS Photonics 6, 3010–3016 (2019), DOI 10.1021/acsphotonics.9b01238, arXiv:1908.10507 | Forward-mode FDTD sensitivity of a dielectric scatterer and a grating coupler | none yet | no table of digits in the paper text available here; results are figures | — | — | `NOT_REPRODUCIBLE_FROM_PAPER` | conceptual only |
| Hughes, Minkov, Williamson, Fan, ACS Photonics 5, 4781–4787 (2018), DOI 10.1021/acsphotonics.8b01522 | Adjoint inverse design of a nonlinear nanophotonic device | not started; blocked until the forward model is validated | device figures, no extractable single-number acceptance table in the abstract record | — | — | `NOT_REPRODUCIBLE_FROM_PAPER` | conceptual only |
| Minkov et al., ACS Photonics 7, 1729–1741 (2020), DOI 10.1021/acsphotonics.0c00327 | Automatic differentiation of photonic-crystal bands | not started | band-structure figures | — | — | `NOT_REPRODUCIBLE_FROM_PAPER` | conceptual only |
| Hamid and Cooray, Advanced Electromagnetics, 30 Dec 2016, “Two-Dimensional Scattering by a Homogeneous Gyrotropic-Type Elliptic Cylinder” | Table I, elliptic ferrite, axial ratio 2, `k0 a = 1`, `μ = (7/9, −5/18, 10)` used by duality as an Hz permittivity | `hamid_ellipse_benchmark.py`, criteria frozen in `PUBLICATION_BENCHMARK_PASS_CRITERIA.md` | six echo widths, 0.13781 to 0.77898 | pending this run | pending | reconstructable from the table and the figure caption | not yet run |
| Monzon and Damaskos, IEEE TAP 34, 1243–1249 (1986); Beker, Umashankar, Taflove, IEEE TAP 37, 1573–1581 (1989) | anisotropic-rod scattering used as a numerical check in later papers | not selected | the open copies found here do not give a complete digit table plus every length in one place | — | — | `NOT_REPRODUCIBLE_FROM_PAPER` from the copies in hand | not used as a pass |
| Platzman and Ozaki, J. Appl. Phys. 31, 1597 (1960); Seshadri, Electron. Lett. 1, 256 (1965) | magnetized plasma cylinder, analytic | already covered by our independent gyrotropic `T_n`, worst FEM error `1.0e-4` dB | the classic papers are not re-digitized here | our `T_n` | — | analytic family already in the gyrotropic suite | not a new reproduction |

## What each Fan paper actually specifies

Hughes 2019 differentiates the FDTD update, not a frequency-domain FEM matrix. The design variables in the two examples are a scatterer permittivity and a grating fill factor. The cost discussion is: one extra forward-mode field per design variable, and that is attractive when the number of outputs exceeds the number of inputs. No mesh, fill factor, or efficiency digit is printed in the text available for this audit.

Hughes 2018 is a time-domain adjoint for a nonlinear Kerr device. The linear frequency-domain adjoint is the same stationarity condition only after the objective is reduced to one real scalar. The paper does not publish a linear FEM Jacobian we can difference.

Minkov 2020 differentiates a plane-wave expansion of a photonic crystal. The state equation is a generalized eigenproblem, not our driven Helmholtz system.

## Equation-level notes kept for later

These are not a claim that our solver matches those papers.

- Frequency-domain forward mode, for `A(p) x = b` and `db = 0`, is `A dx/dp = −(dA/dp) x`. That identity does not depend on FDTD versus FEM. Hughes 2019 reaches the same count of extra solves in the time domain.
- A linear adjoint for one real scalar needs one extra solve with the transpose that the objective’s complex calculus actually produces. `A^T` and `A^H` are not interchangeable. That derivation is not started while the full-device forward gate is open.
