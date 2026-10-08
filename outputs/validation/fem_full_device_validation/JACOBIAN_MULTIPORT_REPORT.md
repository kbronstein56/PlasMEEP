# Reduced multiport Jacobian validation

One plasma disk between a fixed source line and a receiving monitor line.
Differentiated observables:

- complex port amplitude `a = Σ w_i H_z(x_i)`
- port phase `arg(a)`
- guide-normal Poynting power through the monitor line
- scalar objective `J = |a|²` with `dJ/ds = 2 Re(a* da/ds)`

Vacuum at the monitor, so `da/ds` is linear in `dx/ds` (no local `dρ` term).
Power uses the complex product rule on `(E, H)` along the line.

Script: `jacobian_multiport.py`  
Data: `jacobian_multiport.json`  
Frozen relative tolerance ≤ 1e-4 with truncation visibility.

## Results

| B (T) | field | amp | power | \|a\|² | phase | pass |
|-------|-------|-----|-------|--------|-------|------|
| 0     | 2.9e-8 | 2.9e-8 | 9.3e-8 | 1.1e-7 | 2.4e-8 | yes |
| +0.05 | 5.0e-8 | 3.4e-8 | 5.7e-8 | 6.0e-8 | 1.3e-8 | yes |

Real and imaginary parts of `da/ds` agree separately at the same order.

## Verdict

**MULTIPORT_JACOBIAN: PASS**
