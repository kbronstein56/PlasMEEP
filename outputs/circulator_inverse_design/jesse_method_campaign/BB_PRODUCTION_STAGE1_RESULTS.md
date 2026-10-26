# Production broadband optimization — stage-1 results + dense audit

**Do not rerun** `production_broadband_opt.py`.

## Design vector location

**Saved:** `opt_broadband_production_C.json` → field `q_best` (16 C6 orbits).  
Also mirrored in `BB_PRODUCTION_OBJ_GRID_SUMMARY.json`.

## Optimizer outcome

| Quantity | Value |
|----------|------:|
| J_BB initial | **−0.8015** |
| J_BB final | **−0.5290** |
| ΔJ_BB | **+0.273** |
| mean J (obj) | −0.442 |
| worst J (obj) | −0.657 @ **5.7 GHz** |
| nit / stop | 30 / iteration limit |
| elapsed | 158231.5 s (~44 h) |

## Objective-grid metrics (21×100 MHz) — final design

| f | J | T_des | IL | iso | RL | TE1 |
|--:|--:|------:|---:|----:|---:|----:|
| 4.0 | −0.329 | 0.019 | 17.2 | −0.9 | 0.6 | 0.045 |
| 4.1 | −0.428 | 0.000 | 39.9 | −18.8 | 1.6 | 0.220 |
| 4.2 | −0.618 | 0.007 | 21.9 | −2.0 | 3.4 | 0.470 |
| 4.3 | −0.250 | 0.243 | 6.2 | 5.1 | 6.1 | 0.254 |
| 4.4 | −0.338 | 0.136 | 8.8 | 20.8 | 3.7 | 0.341 |
| 4.5 | −0.319 | 0.017 | 17.8 | 8.9 | 0.8 | 0.115 |
| 4.6 | −0.415 | 0.215 | 6.7 | 2.4 | 10.8 | 0.249 |
| 4.7 | −0.480 | 0.219 | 6.6 | 3.9 | 10.6 | 0.455 |
| 4.8 | −0.329 | 0.113 | 9.6 | 1.5 | 5.1 | 0.167 |
| 4.9 | −0.513 | 0.014 | 18.6 | −5.2 | 2.9 | 0.290 |
| 5.0 | −0.136 | 0.364 | 4.6 | 8.8 | 11.8 | 0.372 |
| 5.1 | −0.533 | 0.046 | 14.5 | 3.7 | 7.3 | 0.336 |
| 5.2 | −0.599 | 0.074 | 11.4 | 2.2 | 6.3 | 0.267 |
| 5.3 | −0.555 | 0.044 | 13.6 | −1.0 | 3.4 | 0.318 |
| 5.4 | −0.578 | 0.045 | 13.5 | −1.0 | 5.6 | 0.223 |
| 5.5 | −0.628 | 0.027 | 15.7 | −7.8 | 2.9 | 0.110 |
| 5.6 | −0.409 | 0.284 | 5.5 | 4.2 | 15.6 | 0.263 |
| 5.7 | −0.657 | 0.012 | 19.5 | −8.2 | 3.0 | 0.318 |
| 5.8 | −0.399 | 0.019 | 17.3 | 3.2 | 2.4 | 0.214 |
| 5.9 | −0.281 | 0.016 | 18.0 | −1.4 | 0.4 | 0.017 |
| 6.0 | −0.479 | 0.056 | 12.6 | 10.4 | 4.1 | 0.329 |

Means: T_des≈0.094, IL≈14.3 dB, iso≈1.4 dB, RL≈5.2 dB, TE1≈0.256.

## Dense 50 MHz audit of final design

41 points over 4–6 GHz (`bb_final_audit_50MHz.json`).

| Metric | On 21 obj freqs | On 50 MHz audit | Δ | Bar |
|--------|----------------:|----------------:|--:|-----|
| Worst iso | −18.8 dB @ 4.1 | −18.8 dB @ 4.1 | 0 | PASS |
| Worst IL | 39.9 dB @ 4.1 | 39.9 dB @ 4.1 | 0 | PASS |
| Worst J | −0.657 @ 5.7 | **−0.827 @ 5.65** | 0.17 | **FAIL soft-min** |

## Targeted 25 MHz (around 4.1 and 5.65)

(`bb_final_audit_25MHz_local.json`)

| Metric | Result |
|--------|--------|
| Worst iso / IL | still @ **4.1 GHz** (unchanged) |
| Worst J | **−1.153 @ 5.625 GHz** (Δ ≈ 0.33 vs 50 MHz) |

## Enrichment verdict

**Adaptive enrichment is needed.** Exact priority frequencies to add:

`4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95`

See `BB_ENRICHMENT_RECOMMENDATION.md`. **No automatic re-optimization.**
