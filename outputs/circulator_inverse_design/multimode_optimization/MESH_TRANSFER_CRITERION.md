# Mesh-transfer criterion (frozen BEFORE M→F judgment)

Applied to the **M-optimized** design evaluated on F (not the C optimum
transferred through M/F).

## Pass (call mesh-transfer OK)

All of the following on TE0-drive means at 3.85 GHz, B=+0.05 T:

| metric | |Δ(F − M)| bar |
|--------|----------------|
| J | ≤ 0.05 |
| T_des | ≤ 0.04 |
| IL_dB | ≤ 0.40 dB |
| isolation_dB | ≤ 1.5 dB |
| RL_dB | ≤ 1.5 dB |
| TE1 | ≤ 0.015 |
| through | ≤ 0.010 |
| other TE0 | ≤ 0.030 |

## Fail

Any bar exceeded → **MESH_CONVERGENCE remains FAIL**. Do not declare
mesh-converged. Next step is local interface refinement / finer reference,
not more optimizer iterations alone.

## Note

C→M→F of the C optimum already failed continuum (ΔJ ~ 0.22 per step). This
stricter M→F check after M refinement asks whether **M-optima** are stable on F.
