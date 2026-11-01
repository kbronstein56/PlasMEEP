# Presentation media QC

## Environment

- Worktree: `/home/daq_user/PlasMEEP-viz` (`agent/circulator-viz`)
- Production worktree `/home/daq_user/PlasMEEP` is **read-only** for this agent
- Conda env: `plasmeep_mpi` (matplotlib + ffmpeg)
- No FEM / Meep solves executed for media builds

## Scientific QC checklist

| Check | Result |
|-------|--------|
| Schematic vs data vs analytic vs simulation labeled on frames | Required by `stamp_kind()` |
| Port numbering P1…P6 matches `sixport_common.port_dirs` | Used in 01 schematic |
| B sign: B_PLUS ⇒ desired = (j+1)%6 | Documented in 01 / stubs |
| Stage-2 videos use real `accepted_iter` / `start_eval` + `q` from history | 04, 09 |
| Final J_BB quoted from `opt_broadband_stage2_A2_C.json` | ≈ −6.229 |
| 4.1 GHz recovery numbers from history `at_4p1` | T_des ~3e−4→~0.14; IL ~40→~8.5 dB |
| Worst-frequency migration from history | 4.1 → 5.65 GHz (final) |
| Enrichment story uses Stage-1 spectral audits + documented 5.625 hole | 06 |
| Validation ladder statuses copied from overview script / summary | 08 — no invented PASS/FAIL |
| TE0/TE1 movie uses analytic 40 mm profiles; animates Re[Hz e^{−iωt}] | 07 |
| No fabricated device EM fields | 02/03 stubs only |

## Awaiting Stage-2 / production field exports (blocks 02 & 03)

Needed complex **Hz** dumps (real+imag or complex), shared XY grid, meta:

1. **Uniform start**, source **P1**, at least one frequency (prefer **4.1 GHz** and final worst ~**5.65 GHz**).
2. **Stage-2 best** (`q_best` in `opt_broadband_stage2_A2_C.json`), same source/frequencies.
3. Optional: **B_MINUS** companion for circulation-sense / Onsager illustration.

Suggested layout (example):

```
outputs/presentation_media/_field_exports/
  uniform_P1_4p1GHz/{Hz_real.npy,Hz_imag.npy,x.npy,y.npy,meta.json}
  stage2best_P1_4p1GHz/...
  stage2best_P1_5p65GHz/...
```

`meta.json` must include: `f_GHz`, `B_tesla`, `source_port`, `design_tag`, `solver`, `mesh_grade`, `units`.

Animation formula once exports exist:

```text
frame(t) = Re[ (Hz_real + 1j*Hz_imag) * exp(-1j*2π*f*t) ]
```

Do **not** ship |Hz| oscillating movies.

## Media QC (format)

- 1920×1080, 30 fps, H.264 `yuv420p`, no audio (`-an`)
- Poster PNG + GIF preview for built IDs
- Frame dumps under `*/frames/` are ephemeral (gitignore candidates)

## Known limitations (honest)

- 01 is geometry/circulation **schematic**, not a field solve.
- 05 Stage-2 comparison uses **objective-grid history** metrics (28 freqs), not a new dense spectral audit of the Stage-2 final design (unless/until a Stage-2 dense audit JSON is added).
- 07 is **analytic guide modes**, not multimode DtN eigenfields from the FEM port solver.
- 08 OPEN/FAIL rows remain OPEN/FAIL — plasma B=0 Meep↔FEM and B≠0 Meep↔FEM are unresolved.
