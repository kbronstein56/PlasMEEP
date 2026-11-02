# Presentation media QC

## Environment

- Worktree: `/home/daq_user/PlasMEEP-viz` (`agent/circulator-viz`)
- Production worktree `/home/daq_user/PlasMEEP` is **read-only** for this agent
- Conda env: `plasmeep_mpi` (matplotlib + ffmpeg)
- FEM field movies consume **precomputed** production exports only (no solves in viz)

## Scientific QC checklist

| Check | Result |
|-------|--------|
| Schematic vs data vs analytic vs simulation labeled on frames | Required by `stamp_kind()` |
| Port numbering P1…P6 matches `sixport_common.port_dirs` | Used in 01 schematic + field meta |
| B sign: B_PLUS ⇒ desired = (j+1)%6 | Documented in 01 / field meta |
| Stage-2 videos use real `accepted_iter` / `start_eval` + `q` from history | 04, 09 |
| Final J_BB quoted from `opt_broadband_stage2_A2_C.json` | ≈ −6.229 |
| 4.1 GHz recovery numbers from history `at_4p1` | T_des ~3e−4→~0.14; IL ~40→~8.5 dB |
| Worst-frequency migration from history | 4.1 → 5.65 GHz (final) |
| Enrichment story uses Stage-1 spectral audits + documented 5.625 hole | 06 |
| Validation ladder statuses copied from overview script / summary | 08 — no invented PASS/FAIL |
| TE0/TE1 movie uses analytic 40 mm profiles; animates Re[Hz e^{−iωt}] | 07 |
| Device EM fields from production Direct-DtN exports | **02 / 03 real** — not fabricated |

## Field movies (02 / 03) — provenance

Exports written by production `export_presentation_fields.py`:

```
…/jesse_method_campaign/presentation_fields/
  uniform_P1_4p100GHz_Bplus/
  stage2_best_P1_4p100GHz_Bplus/
  stage2_best_P1_5p650GHz_Bplus/
  stage2_best_P1_5p750GHz_Bplus/
  INDEX.json, README.md
```

Animation:

```text
frame(t) = Re[ (Hz_real + 1j*Hz_imag) * exp(-1j*2π*f*t) ]
```

| Clip | Shows |
|------|--------|
| 02 | Uniform vs Stage-2 best at **4.1 GHz**, source P1, B+ |
| 03a | Stage-2 best circulation field at **5.65 GHz** |
| 03b | Same q_best: **5.65 GHz** vs audit hole **5.75 GHz** |

Legacy stub posters remain under `02_field_movies_stubs/` / `03_field_movies_stubs/` for history only.

## Media QC (format)

- 1920×1080, 30 fps, H.264 `yuv420p`, no audio (`-an`)
- Poster PNG + GIF preview for built IDs
- Frame dumps under `*/frames/` are ephemeral (gitignore)

## Known limitations (honest)

- 01 is geometry/circulation **schematic**, not a field solve.
- 05 Stage-2 comparison uses **objective-grid history** metrics (28 freqs), not a new dense spectral audit of the Stage-2 final design (Stage-2 50/25 MHz audits exist separately in production).
- 07 is **analytic guide modes**, not multimode DtN eigenfields from the FEM port solver.
- 08 OPEN/FAIL rows remain OPEN/FAIL — plasma B=0 Meep↔FEM and B≠0 Meep↔FEM are unresolved.
- Field grids are **linear interpolations** of nodal FEM Hz onto a regular Cartesian lattice for imshow; nodal stride dumps preserve provenance.
