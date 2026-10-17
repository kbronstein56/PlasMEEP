# MULTIMODE_FREQ_CHARACTERIZATION

Design `outputs/circulator_inverse_design/multimode_optimization/opt_uniform_C.json`, grade=C.
TE1 cutoff (analytic parallel-plate): **3.7474 GHz** (W = clear_width × a = 2.0 × 0.02 m).

| f (GHz) | n_prop | β_TE1 (re, im) | T_des | IL | iso | RL | TE1 | note |
|--------:|-------:|----------------|------:|---:|----:|---:|----:|------|
| 3.70 | 1 | (0, 0.2491) | 0.0637 | 12.00 | 0.70 | 5.74 | 0.0000 | TE1 evan. |
| 3.75 | 2 | (0.05846, 0) | 0.0825 | 10.83 | -6.05 | 5.86 | 0.0455 | near cutoff |
| 3.80 | 2 | (0.2641, 0) | 0.0964 | 10.16 | 2.41 | 1.38 | 0.0821 | ok |
| 3.85 | 2 | (0.3701, 0) | 0.8041 | 0.95 | 19.00 | 16.63 | 0.0187 | ok |
| 3.90 | 2 | (0.4528, 0) | 0.3750 | 4.26 | 17.13 | 6.12 | 0.2106 | ok |
| 3.95 | 2 | (0.5234, 0) | 0.0377 | 14.24 | -0.47 | 0.93 | 0.0768 | ok |
| 4.00 | 2 | (0.5864, 0) | 0.0291 | 15.36 | 3.05 | 1.19 | 0.0591 | ok |

## Broadband decision (Phase 10)

**BROADBAND_READY: NEEDS_CUTOFF_DERIVATION**

Single-frequency design lives just above TE1 cutoff. Any band that approaches or crosses ~3.747 GHz changes the number of propagating channels; do not reuse a fixed 12-channel objective without an explicit frozen cutoff treatment.

- Option A: Optimize entirely above TE1 cutoff (safely > 3.7474 GHz), e.g. 3.80–4.00 GHz if TE1 stays propagating and the objective always uses 12 channels.
- Option B: Allow crossing TE1 cutoff with a frequency-dependent propagating-channel objective (count TE1 only when β_TE1 is real). Requires a new frozen derivation before broadband opt.

Broadband optimization itself is deferred to the next campaign.

