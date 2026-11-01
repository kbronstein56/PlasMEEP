# Presentation visual references (conventions only — no copyrighted figure copies)

This note captures **visual conventions** useful for PlasMEEP talks. It does **not** reproduce copyrighted paper figures.

## Device / circulator schematics

- Show **six ports** on a hex layout with clear **P1…P6** labels.
- Draw **B** as out-of-page / into-page with an explicit **sign** (B_PLUS vs B_MINUS).
- Circulation arrows must match the code convention for that B sign (here: B_PLUS ⇒ j→j+1).
- Prefer a **SCHEMATIC** badge whenever rods/horns are idealized.

## Field visualizations (when real data exist)

- Prefer **Re[Hz e^{−iωt}]** cyclic movies or static **complex** maps (mag+phase), not fake |H| flicker.
- Use a diverging colormap (e.g. RdBu) centered at 0 for instantaneous Re[Hz].
- Always annotate: frequency, B sign, source port, design tag, solver.
- Side-by-side **uniform vs optimized** on a **shared color scale**.

## Spectra / optimization

- Frequency on the horizontal axis (GHz); keep 4–6 GHz band limits visible.
- Mark **objective-grid nodes** vs **audit densification** distinctly (enrichment story).
- Plot **worst-frequency migration** and **4.1 GHz recovery** as first-class callouts for Stage-2.
- FoM curves: higher J_BB better; label soft-min / A2 explicitly.

## Validation communication

- Ladder / traffic-light statuses only from **saved audits** (PASS / FAIL / OPEN).
- Separate material Faraday checks from full-device Meep↔FEM claims.
- Do not imply B≠0 full-device Meep↔FEM agreement — it is **OPEN / never run**.

## Mode content

- TE0 vs TE1: show transverse structure across the **40 mm** guide width.
- In multimode S talks, distinguish TE0-preserving objectives from total port-power objectives.

## Typography / slide hygiene

- Large labels readable at 16:9 projection distance.
- One message per scene; badge data provenance every time.
- Avoid decorative clutter that looks like data.
