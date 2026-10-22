# Circulator literature log

Persistent log for PlasMEEP / PMM circulator methodological work.
Search date for this seed pass: **2026-10-21**.
Literature informs ideas; it does **not** override validated DtN / gyrotropic / power-balance numerics.

Classification: `CORE` | `DIRECTLY_USEFUL` | `BACKGROUND` | `NOT_APPLICABLE`

---

## LIT-001 — Wang & Fan (2005)

| Field | Value |
|-------|--------|
| Citation | Z. Wang and S. Fan, “Optical circulators in two-dimensional magneto-optical photonic crystals,” *Opt. Lett.* **30**, 1989–1991 (2005). |
| DOI / URL | https://doi.org/10.1364/OL.30.001989 |
| Status | Peer-reviewed journal |
| Device | 3-port MO photonic-crystal cavity circulator |
| Regime | Optical (telecom-scale PC) |
| Ports / modes | 3 waveguides; cavity counter-rotating modes |
| Symmetry | Naturally **C3** (Y-junction) |
| Nonreciprocity | Magneto-optical cavity domain engineering |
| Optimization | Analytic / design guideline (not TO) |
| Solver | Photonic-crystal / modal coupling analysis |
| Why it matters | Canonical answer to Jesse’s C3 vs C6 question: classic circulator literature is three-port → C3. Our six-port cyclic target is C6-compatible but C6 is a *choice*, not physics necessity. |
| Idea to test | Map Fan-style C3 thinking onto six-port as two C3 classes `{1,3,5}` / `{2,4,6}`. |
| Class | **CORE** |
| Code change? | No production change; motivates C3 parameterization study. |

---

## LIT-002 — Mehrpour Bernety et al. (2022) PoP

| Field | Value |
|-------|--------|
| Citation | H. Mehrpour Bernety, L. S. Houriez, J. A. Rodríguez, B. Wang, M. A. Cappelli, “A tunable microwave circulator based on a magnetized plasma as an active gyrotropic element,” *Phys. Plasmas* **29**, 112114 (2022). |
| DOI / URL | https://doi.org/10.1063/5.0123459 |
| Status | Peer-reviewed journal |
| Device | Magnetized plasma microwave circulator |
| Regime | **S/C-band microwave** (same lab / plasma context as ours) |
| Ports | Multiport microwave; plasma GDT as gyrotropic element |
| Nonreciprocity | Anisotropic **ε** of magnetized plasma (not ferrite μ) |
| Validation | Theory + experiment; isolation ≳25 dB reported in band |
| Why it matters | Closest physics sibling to PMM six-port work; validates plasma-ε nonreciprocity experimentally. |
| Idea | Keep experimental isolation/IL targets honest vs continuum FEM, not M-mesh optimistic FoMs. |
| Class | **CORE** |
| Code change? | None. |

---

## LIT-003 — Li et al. (2022) PRApplied Y-junction plasma

| Field | Value |
|-------|--------|
| Citation | F. Li, R. J. Davis, S. M. Kandil, D. F. Sievenpiper, “Electromagnetic Nonreciprocity in a Magnetized Plasma Circulator,” *Phys. Rev. Applied* **18**, 064066 (2022). |
| DOI / URL | https://doi.org/10.1103/PhysRevApplied.18.064066 ; arXiv:2206.12507 |
| Status | Peer-reviewed (+ preprint) |
| Device | E-plane **Y-junction** plasma circulator, WR187 |
| Regime | **4–6 GHz** (exact advisor band) |
| Ports / symmetry | **3-port** → C3 geometry |
| Nonreciprocity | Magnetized plasma column; B direction sets circulation |
| Why it matters | Independent experimental 4–6 GHz plasma circulator; reinforces that three-port devices are C3 and that broadband microwave plasma circulation is experimentally real. |
| Idea | Use as external reality check for isolation bandwidth claims. |
| Class | **CORE** |
| Code change? | None. |

---

## LIT-004 — Hughes, Williamson, Minkov & Fan (2019)

| Field | Value |
|-------|--------|
| Citation | T. W. Hughes, I. A. D. Williamson, M. Minkov, S. Fan, “Forward-Mode Differentiation of Maxwell’s Equations,” *ACS Photonics* **6**, 3010–3016 (2019). |
| DOI / URL | https://doi.org/10.1021/acsphotonics.9b01238 ; arXiv:1908.10507 |
| Status | Peer-reviewed |
| Method | Forward-mode differentiation of Maxwell; compared to adjoint |
| Software | Underpins **ceviche** AD stack |
| Why it matters | Jesse’s JAX/autograd / “extra simulation” question. Forward mode competitive when #parameters ≪ #objectives; adjoint when many parameters / few objectives. |
| Idea | Benchmark forward-mode vs adjoint for our **16**-param C6 problem after shared LU. |
| Class | **DIRECTLY_USEFUL** |
| Code change? | Benchmark only; do not replace production adjoint without win. |

---

## LIT-005 — Ceviche / HIPS autograd (software + Hughes comments)

| Field | Value |
|-------|--------|
| Citation | fancompute/ceviche (GitHub); Hughes workshop FAQ on AD≡reverse-mode adjoint |
| URL | https://github.com/fancompute/ceviche |
| Status | Open-source research software |
| Method | HIPS **autograd** with **custom sparse primitives** for `Ax=b` VJP (= transpose solve) |
| Why it matters | Corrects a common misconception: autograd does **not** magically differentiate SciPy `splu`; ceviche hard-codes sparse-solve adjoints. |
| Idea | Optional JAX `custom_vjp` wrapping *our* validated adjoint (API sugar), not a FULL91 rewrite. |
| Class | **DIRECTLY_USEFUL** |
| Code change? | Tiny custom_vjp prototype only if useful. |

---

## LIT-006 — Hammond et al. (2022) Optics Express

| Field | Value |
|-------|--------|
| Citation | A. M. Hammond, A. Oskooi, M. Chen, Z. Lin, S. G. Johnson, S. E. Ralph, “High-performance hybrid time/frequency-domain topology optimization for large-scale photonics inverse design,” *Opt. Express* **30**, 4467–4491 (2022). |
| DOI / URL | https://doi.org/10.1364/OE.442074 |
| Status | Peer-reviewed (Editor’s Pick) |
| Method | Hybrid **Meep FDTD** + multi-frequency TO; filter-design sources |
| Why it matters | Time-domain broadband optimization / validation path Jesse asked for; Meep as independent check. |
| Idea | Keep Meep for **validation**, FEM FD adjoint for **optimization** unless TD TO clearly wins on our plasma device. |
| Class | **DIRECTLY_USEFUL** |
| Code change? | TD validation plan; no production optimizer swap. |

---

## LIT-007 — Gedeon, Christiansen & Sigmund (2026 preprint)

| Field | Value |
|-------|--------|
| Citation | J. Gedeon, R. E. Christiansen, O. Sigmund, “Spectrally smooth broadband response via autocorrelation-constrained inverse design,” arXiv:2608.13367 (2026). |
| DOI / URL | https://doi.org/10.48550/arXiv.2608.13367 |
| Status | **Preprint** (not yet confirmed journal) |
| Method | Time-domain TO + long-lag **autocorrelation** metric to suppress in-band spikes hidden by integrated FoMs |
| Why it matters | Directly addresses Jesse’s “sharp isolation peaks vs true broadband” worry for aggregate objectives. |
| Idea | Pair soft-min / dense audit with optional spectral-smoothness penalty; adaptive frequency enrichment remains primary FD defense. |
| Class | **DIRECTLY_USEFUL** |
| Code change? | Consider as optional BB regularizer after sampling study; not required for first BB opt. |

---

## LIT-008 — Zhao et al. (2025) ACS Photonics

| Field | Value |
|-------|--------|
| Citation | W. Zhao et al., “Ultrabroadband On-Chip Optical Circulator via Traveling-Wave-Driven Nonreciprocal Mode Conversion,” *ACS Photonics* **12**, 6650–6661 (2025). |
| DOI / URL | https://doi.org/10.1021/acsphotonics.5c01650 |
| Status | Peer-reviewed |
| Device | LNOI EO circulator using **TE0↔TE1** nonreciprocal mode conversion + multiplexers |
| Why it matters | Circulation can *use* mode conversion; TE0 preservation is not universal. Supports Jesse’s mode-preservation question. |
| Idea | Port-power (total desired-port) objective as controlled alternative to TE0-preserving FoM. |
| Class | **CORE** (for mode question) |
| Code change? | Objective-B study only. |

---

## LIT-009 — Wang et al. (2021) PRL modal circulator

| Field | Value |
|-------|--------|
| Citation | J. Wang, J. F. Herrmann, J. D. Witmer, A. H. Safavi-Naeini, S. Fan, “Photonic Modal Circulator Using Temporal Refractive-Index Modulation with Spatial Inversion Symmetry,” *Phys. Rev. Lett.* **126**, 193901 (2021). |
| DOI / URL | https://doi.org/10.1103/PhysRevLett.126.193901 |
| Status | Peer-reviewed |
| Idea | Modes themselves can be circulation channels (modal circulator). |
| Why it matters | Conceptual support that TE1 can carry “information” / routing; measurement of output mode content is scientifically interesting. |
| Class | **BACKGROUND** / **DIRECTLY_USEFUL** for mode narrative |
| Code change? | None; informs experimental mode-resolution discussion. |

---

## Search hygiene

Preferred sources: DOI / publisher / Crossref. Preprints labeled.  
Next pass triggers: new BB objective freeze; C3 results; TD Meep disagreements; new multimode circulator papers 2025–2026.
