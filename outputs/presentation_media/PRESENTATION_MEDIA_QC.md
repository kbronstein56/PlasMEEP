# Presentation media QC

## Environment

- Worktree: `/home/daq_user/PlasMEEP-viz` (`agent/circulator-viz`)
- Production worktree `/home/daq_user/PlasMEEP` is **read-only** for viz renders
- Conda env: `plasmeep_mpi` (matplotlib + ffmpeg)
- FEM field movies consume **precomputed** production exports only

## Stage-3 playback QC (2026-11-08)

| Check | Result |
|-------|--------|
| Field movies ~12–15 s @ 30 fps | **10/11/12 = 13.0 s**; sweep/evolution = 14.0 s |
| ~2 visual phase periods + start/end holds | `FIELD_N_PERIODS=2`, `FRAMES_PER_PERIOD=180`, `HOLD=15` |
| GIF loops indefinitely (`ffmpeg -loop 0`) | yes for 10–13, 18 |
| Poster PNG present for each important clip | yes |
| Shared color scale on S2 vs S3 pairs | yes (`shared_scale=True`) |
| Animate `Re[Hz e^{-iωt}]` not `|Hz|` | yes |
| Titles / colorbars readable on posters | inspected 10 poster + fig 13 |
| Ports / B sign from export meta | P1, B_PLUS |
| Static figures tell story without video | 13–18 PNGs in `stage3_figures/` |

## Scientific QC checklist

| Check | Result |
|-------|--------|
| Schematic vs data vs analytic vs simulation labeled | `stamp_kind()` |
| Port numbering P1…P6 | device summary + field meta |
| B_PLUS ⇒ desired=(j+1)%6 | documented |
| Stage-3 hole repair metrics match dense audit | 5.75 / 5.825 bars |
| New hidden failures marked on spectral figures | 5.8125, 4.775 |
| Dense-audit recommendation | **C** (proposal only Stage-4) |

## Field exports consumed

```
…/presentation_fields/
  stage2_best_P1_{4p100,5p650,5p750,5p825}GHz_Bplus/
  stage3_best_P1_{4p100,4p300,4p775,5p750,5p812,5p825}GHz_Bplus/
```

Note: `5.8125` GHz export stem rounds to `5p812` via `:.3f`; meta.json retains exact `f_GHz=5.8125`.

## Known follow-ups

- Optionally rebuild legacy 02/03 with the new slow schedule (script updated; large re-encode).
- Optional field movie at 5.8125 / 4.775 for the Stage-4 motivation slide.
