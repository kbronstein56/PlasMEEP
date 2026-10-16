# ADAPTIVE_BROADBAND_SAMPLING_REPORT

## Policy (frozen)

**Do not** pick one sparse grid and call it broadband.

1. **Initial objective grid:** \(4.0, 4.2, \ldots, 6.0\) GHz → **11** frequencies.
2. Short optimization stage with frozen \(J_{\mathrm{BB}}\) (soft-min + mean).
3. **Dense audit:** \(4.0\)–\(6.0\) GHz every **0.05** GHz → **41** frequencies (forward \(J(f)\) only).
4. Identify worst \(J(f)\), isolation valleys, IL peaks, TE1 spikes, and criteria violations
   **between** objective samples.
5. **Add** unresolved frequencies to the objective set.
6. Repeat optimize → audit → enrich until:
   - dense-audit worst \(J\) is within a small gap of the objective soft-min, and
   - adding frequencies no longer reveals new hidden failures.
7. Final frequency-domain validation: ≤50 MHz spacing; 25 MHz near sharp features.

## Advisor question 1 — how many frequencies?

**Answer:** Adaptive — not a fixed number asserted in advance.

| Stage | Count | Role |
|-------|------:|------|
| Initial objective | 11 | Optimization |
| Dense audit | 41+ | Discovery of hidden dips |
| Enriched objective | 11 + (added worst pts) | Re-opt |
| Final FD validation | ~41–80 | Accept/reject |

The number that was **ultimately necessary** is recorded below after the sanity
campaign (fill-in section).

## Results (filled as campaign proceeds)

### Iteration 0 — initial grid only

- Objective frequencies: 11 (4.0–6.0 / 0.2 GHz)
- Status: awaiting Phase 7 sanity optimization

### Iteration 1 — after first dense audit

- TBD: worst audit \(f\), Δ vs nearest objective sample, frequencies added

### Convergence

- TBD: final objective count, final audit worst \(J\), stop reason
