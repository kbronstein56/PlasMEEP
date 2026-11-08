# Stage-3 dense broadband audit

**Design:** `opt_broadband_stage3_A2_C.json` `q_best` (J_BB=−3.659, 30/30)  
**Branch:** `agent/eigenmode-ports`  
**No optimization.** Forward-only A2 metrics, `n_workers=12`, `omp=1`.

## Grids

| Scan | Range | Spacing | n | Artifact |
|------|-------|--------:|--:|----------|
| Full-band dense | 4.00–6.00 GHz | **25 MHz** | 81 | `bb_stage3_dense_25MHz.json` / `.csv` |
| Targeted refine | neighborhoods of worsts + 4.3 / 5.775 | **12.5 MHz** | 160 | `bb_stage3_dense_12p5MHz_local.json` / `.csv` |

25 MHz run seeded 51 frequencies from prior Stage-3 50 MHz + local 25 MHz audits; 30 new nodes solved in parallel.

## Stage-2 vs Stage-3 (intended story)

| Quantity | Stage-2 | Stage-3 |
|----------|--------:|--------:|
| Objective grid | 28 freq, J_BB=**−6.229** | 30 freq (+5.75,+5.825), J_BB=**−3.659** |
| Dense-audit holes motivating Stage-3 | **5.75** A2≈−10.98, iso≈−20.1 dB; **5.825** A2≈−10.71, iso≈−21.3 dB | — |
| Same nodes after Stage-3 | — | **5.75** A2=−2.95, iso=+3.1 dB; **5.825** A2=+1.38, iso=+17.9 dB |
| 4.1 GHz | T_des=0.141, IL=8.5 dB | T_des=0.289, IL=5.4 dB (**improved**) |

**Verdict on Stage-3’s stated goal:** the 5.75 / 5.825 enrichment nodes were **strongly repaired**, and 4.1 GHz was preserved and improved.

## True continuum worsts (this audit)

### On the full-band 25 MHz grid

| Rank | f (GHz) | A2 | IL (dB) | iso (dB) |
|-----:|--------:|---:|--------:|---------:|
| 1 | **4.775** | **−8.587** | **25.90** | **−8.35** |
| 2 | 4.225 | −6.621 | 20.72 | −5.74 |
| 3 | 4.375 | −5.792 | 16.73 | −4.75 |

Previously reported “worst ≈ 4.3 GHz / A2≈−4.57” was an **artifact of the 50 MHz audit** — 25 MHz exposes deeper off-grid corpses, especially **4.775 GHz**.

### After 12.5 MHz local refinement

| Rank | f (GHz) | A2 | IL (dB) | iso (dB) | Comment |
|-----:|--------:|---:|--------:|---------:|---------|
| 1 | **5.8125** | **−12.021** | **31.96** | **−18.67** | Between objective nodes **5.800** and **5.825** |
| 2 | 4.775 | −8.587 | 25.90 | −8.35 | Narrow spike (neighbors much better) |
| 3 | 5.5875 | −8.529 | 23.99 | −10.48 | Mid-band off-grid |
| 4 | 5.7125 | −8.097 | 27.68 | −6.15 | Near old 5.75 neighborhood |

**Critical spectroscopy:** at the Stage-3 objective frequencies themselves:

| f (GHz) | in Stage-3 obj? | A2 | iso (dB) |
|--------:|:---------------:|---:|---------:|
| 5.800 | yes | **+5.53** | **+34.3** |
| **5.8125** | **no** | **−12.02** | **−18.7** |
| 5.825 | yes | +1.38 | +17.9 |

Stage-3 closed the *sampled* holes while leaving a **~12.5 MHz-wide isolation corpse** between enrichment nodes — the same sampling pathology that motivated Stage-2→Stage-3.

5.775 residual (A2≈−4.73) is real but **no longer the dominant** failure once 12.5 MHz sampling is used.

## Are Stage-3’s old holes still repaired?

Yes at the exact enrichment frequencies:

- 5.75 GHz: A2 −10.98→−2.95; iso −20.1→+3.1 dB  
- 5.825 GHz: A2 −10.71→+1.38; iso −21.3→+17.9 dB  

But the continuum between them is **not** repaired.

## Compact data for plotting

- `bb_stage3_dense_25MHz.csv`
- `bb_stage3_dense_12p5MHz_local.csv`
- `BB_STAGE3_DENSE_AUDIT_SUMMARY.json`

## Recommendation

### **C — Stage-3 still contains a substantial hidden broadband failure.**

Evidence:

1. New off-grid isolation/A2 corpse at **5.8125 GHz** (A2≈−12, iso≈−19 dB, IL≈32 dB) comparable to Stage-2’s original 5.75/5.825 failures.  
2. Additional sharp failures at **4.775**, **5.5875**, **5.7125** GHz.  
3. Soft-min on a 30-frequency set again overfit the sampled nodes.

Stage-3 remains a valid *optimization* result on its stated grid, and the 5.75/5.825 / 4.1 goals were met on those nodes. It is **not** yet a validated continuum broadband design.

A narrowly targeted Stage-4 proposal (enrichment only; no launch) is in `BB_STAGE4_PROPOSAL.md`.

**Stage-4 is NOT AUTHORIZED and must not be started.**
