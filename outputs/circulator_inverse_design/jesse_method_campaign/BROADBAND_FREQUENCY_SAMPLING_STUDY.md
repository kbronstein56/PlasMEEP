# Broadband frequency sampling study (4–6 GHz)

## Frozen grid-sufficiency bars (before results)

`GRID_SUFFICIENCY_BARS_FROZEN.json`: Δiso < 0.5 dB, ΔIL/RL < 0.2 dB, location ≤100 MHz.

## Results — SF design on FULL91-C

| Grid | Worst iso | @ f | Worst IL | Pass vs finer? |
|------|----------:|-----:|---------:|:--------------:|
| 200 MHz (11) | −24.2 dB | 5.8 | 40.8 dB | **FAIL** vs 100 |
| 100 MHz (21) | −33.1 dB | 5.9 | 43.5 dB | **FAIL** vs 50 |
| 50 MHz (41) | **−52.8 dB** | **5.75** | **58.1 dB** | (finest here) |

200→100: Δiso **8.9 dB**. 100→50: Δiso **19.8 dB**.

Narrow catastrophic feature near **5.75 GHz** is invisible on coarse grids.

## Answer for Jesse

> There is no fixed safe \(N\). For *this* spectrum, even 100 MHz misses a
> −53 dB isolation hole at 5.75 GHz that 50 MHz catches. Production BB uses
> soft-min on a **100 MHz initial objective set**, then **adaptive enrichment**
> from ≤50 MHz (and 25 MHz locally) audits — never changing frequencies mid
> line-search. Final validation is denser still.

See freeze: `BROADBAND_PRODUCTION_SPEC_FROZEN.md`.
