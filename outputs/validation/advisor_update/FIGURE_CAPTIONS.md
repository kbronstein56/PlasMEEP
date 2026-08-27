# Advisor-update figure captions

## fig1_faraday_benchmark.png
Homogeneous Faraday rotation benchmark using the PlasMEEP `GyrotropicDrudeSusceptibility`
path (fs=5 GHz, fp=2 GHz, |B|=0.05 T). Left: measured κ at B=0 / +B / −B versus independent
cold-plasma theory. Right: |κ| relative error at res=32 (~0.09%) and res=64 (~0.02%).
PASS criteria use magnitude and B-reversal; absolute Stokes/theory sign convention is opposite.

## fig2_reciprocity_vs_resolution.png
B=0 six-port reciprocity residual for the worst orientation class (P1↔P2, axis↔diagonal)
versus Meep resolution. Baseline notebook-faithful Hz-line sources peak at **1.64 dB** at
res=64; TE1 launch reduces that to **1.03 dB** but remains above the 0.2 dB gate.
Diagonal↔diagonal control pair P2↔P3 stays ~0 at res=32 for both formulations.

## fig3_incident_power_uniformity.png
Ratio of launched incident power for axis-aligned ports to diagonal (60°) ports.
Baseline retains a ~10–14% orientation bias at all resolutions. TE1 cos-profile launch
brings the ratio within a few percent of unity (inside the shaded ±3% band).

## Fig. 4 — Measurement formulation screen
`fig4_measurement_formulation_screen.png`

Cheap (res=32, rt=40) P1↔P2 residuals for orientation-aware receivers under TE1 (and baseline) launch.
None of the guide-normal / DFT S·n variants reach the 0.2 dB gate or beat baseline+axis at this setting;
dense sampling matches sparse. Supports stopping measurement escalation in favor of PEC/grid / reference-guide work.

