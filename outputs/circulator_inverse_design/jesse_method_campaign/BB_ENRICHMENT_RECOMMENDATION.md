# Adaptive enrichment recommendation (after stage-1 BB opt)

**Do not start another 30-iteration optimization automatically.**

## Enrichment needed?

**YES.**

| Check | Result |
|-------|--------|
| Worst iso on 50 MHz vs 21 obj freqs | **PASS** (both −18.8 dB @ **4.1 GHz**, on-grid) |
| Worst IL on 50 MHz vs 21 obj freqs | **PASS** (both 39.9 dB @ **4.1 GHz**, on-grid) |
| Soft-min / worst-J between-grid | **FAIL** — off-grid hole |

## Critical between-grid failure

| Grid | Worst J_A | Frequency |
|------|----------:|----------:|
| 21×100 MHz objective | −0.657 | 5.7 GHz |
| 50 MHz audit | −0.827 | **5.65 GHz** |
| 25 MHz local | **−1.153** | **5.625 GHz** |

The soft-min objective never saw 5.625 GHz. That point is **worse than every
objective-grid frequency** and must be added before any restart.

## Exact frequencies to add

### Priority set (recommended next objective grid = 21 + these 7)

```
4.55, 4.85, 5.35, 5.625, 5.65, 5.85, 5.95
```

→ new grid size **28** frequencies (still fixed during each L-BFGS run).

### Full between-grid spike list (optional broader enrichment)

```
4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 5.05, 5.15, 5.25, 5.35,
5.625, 5.65, 5.775, 5.85, 5.95
```

## What NOT to change

Frozen production choices remain: C6, Objective A, τ=0.05, α=0.25, mesh C.
Only the **frequency set** should grow via enrichment restart.

## Suggested next step (manual approval)

1. Restart L-BFGS from `q_best` in `opt_broadband_production_C.json`
2. Use priority 28-frequency set (fixed for the whole run)
3. After that run, re-audit at ≤50 MHz; do not auto-chain further opts

Artifacts: `bb_final_audit_50MHz.json`, `bb_final_audit_25MHz_local.json`,
`bb_final_between_grid_report.json`.
