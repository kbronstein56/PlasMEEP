# Exact Meep `num_mode_guide_normal` source / receiver specification

Derived from `scripts/validation/port_formulations.py`,
`plasmeep/ports/numerical_launch.py`, `plasmeep/ports/numerical_mode.py`,
and `sixport_common.normalize_port` / device column assembly.

Time convention: Meep DFT / Faraday validation use \(e^{-i\omega t}\).

---

## Source

1. Load cached `NumericalPortMode` via `get_numerical_mode(port, res, fs_a, …)`.
2. Mode fields: complex samples `hz_complex[k]` on offsets `offsets_a[k]` along the
   **port tangent** \(t = (-n_y, n_x)\) through `horn[source_center]`.
3. L2-normalize: \(\phi = H_z / \|H_z\|_2\).
4. Discrete source weights:
   \[
   a_k = \sqrt{P_0}\,\overline{\phi_k},\qquad \sum_k |a_k|^2 = P_0
   \]
   (default \(P_0=1\) before incident normalization).
5. Meep places point `Hz` sources at
   \[
   x_k = x_{\mathrm{src}} + s_k\, t
   \]
   with `GaussianSource(frequency=fs_a, fwidth=source_df)`, `source_df = 0.10 * fs_a`.

Absolute amplitude is fixed by the **straight-feed reference** (below), not by \(P_0\).

---

## Receiver (guide-normal flux)

`make_guide_normal_flux_spec(center, outward_n̂, n_points=31, span_factor=0.96)`:

- Center = `monitor_center_for_port` (not source center).
- Outward unit normal \(n̂\) = `effective_port_dir(port)`.
- Tangent \(t = (-n_y,n_x)\).
- Span \(L = 0.96 \times\) `clear_width`.
- Sample offsets \(s_j = \mathrm{linspace}(-L/2,+L/2,31)\).
- Trapezoid weights \(w_j\) with half-weight endpoints; `weight_profile="uniform"`.
- Meep FluxRegions are axis-aligned point monitors; weight of \(X\)-flux is \(n_x w_j\),
  of \(Y\)-flux is \(n_y w_j\). Combined:

\[
P_{\mathrm{raw}} = \sum_j w_j\, (n_x S_x + n_y S_y)\big|_{x_{\mathrm{mon}}+s_j t}
\]

i.e. discrete \(\displaystyle\int S\cdot n̂\,ds\) across the feed cross-section.

Meep’s DFT flux equals the continuum time-average Poynting for the 2D TE (\(H_z,E_x,E_y\)) case:

\[
S_x = \tfrac12\mathrm{Re}(E_y H_z^*),\qquad
S_y = -\tfrac12\mathrm{Re}(E_x H_z^*).
\]

(Confirmed also in `extract_dft_sdotn_power`.)

Positive \(P_{\mathrm{raw}}\) = power in the \(+n̂\) (outward) direction.
`make_guide_normal_flux_spec` returns `sign=1.0`; `extract_flux_powers` multiplies by that sign.

---

## Incident normalization + reflection subtraction

**Reference run** (`normalize_port`): straight PEC feed through the source line,
same formulation sources, flux monitor at the **same** measure center / outward normal.
Records:

- \(P_{\mathrm{inc}} = |P_{\mathrm{raw,ref}}|\)
- Meep `get_flux_data(monitor)` for later subtraction

**Device run** (source port \(i\)):

1. Place flux monitors on all output ports with outward normals.
2. At the **source** port monitor only: `load_minus_flux_data(mon_i, incident_flux_data_i)`
   (subtracts the reference incident field contribution so diagonal ≈ reflection).
3. After ringdown: \(P_{\mathrm{raw},j}\) from each monitor.
4. Normalized column:

\[
T_{ji} = \frac{P_{\mathrm{raw},j}}{P_{\mathrm{inc},i}}
\]

---

## What FEM must reproduce

| Quantity | FEM implementation |
|---|---|
| Source shape | Same \(\phi_k\), \(s_k\), \(x_{\mathrm{src}}\), \(t\) → soft RHS on nodes near line |
| Frequency | Single \(\omega = 2\pi f_s\) (dispersive \(\varepsilon(\omega)\)) |
| Receiver | Reconstruct \(E\) from \(H_z\), integrate \(S\cdot n̂\) on same line |
| \(P_{\mathrm{inc}}\) | Same straight-feed geometry + same flux definition |
| Reflection | Subtract reference incident flux contribution at source port (FEM: \(P_{\mathrm{raw},src}-P_{\mathrm{inc}}\) if reference is pure outgoing incident with matching phase/amplitude — see Phase 5 notes) |

### FEM incident subtraction note

Meep subtracts **complex DFT flux data** (field-level), not just a scalar.
For single-frequency FEM, an equivalent approach is:

1. Solve straight-feed reference → \(P_{\mathrm{inc}}\) and complex \(H_z^{\mathrm{inc}}\) on the source monitor line.
2. On the device, compute \(P_{\mathrm{raw}}\) with the same quadrature.
3. For the source port only, form the **scattered** power via
   \(P_{\mathrm{scat}} = P_{\mathrm{device,src}} - P_{\mathrm{inc}}\)
   **only if** amplitudes/phases of incident match. Safer FEM approach matching Meep intent:
   - Scale device source so reference outgoing power = \(P_{\mathrm{inc}}\) of Meep, **or**
   - Compute FEM \(P_{\mathrm{inc}}\) independently and report \(T_{ji}=P_{\mathrm{raw},j}/P_{\mathrm{inc}}^{\mathrm{FEM}}\),
     with source-port reflection \(R = (P_{\mathrm{raw,src}}-P_{\mathrm{inc}})/P_{\mathrm{inc}}\) when
     the reference and device share the same source weights and the monitor sees the same incident wave.

We implement the scalar power form \(R=(P_{\mathrm{src}}-P_{\mathrm{inc}})/P_{\mathrm{inc}}\) after verifying
on a straight guide that \(P_{\mathrm{src}}\approx P_{\mathrm{inc}}\) (matched launch → \(R\approx 0\)).

---

## Coordinates

All centers from `horn_for_port` / `monitor_center_for_port` are in **Meep cell coordinates**
when `set_geometry_context` is active with the device array placement
(array at `(nx_ports/2, ny_ports/2)` already baked into horn geometry in sixport).
