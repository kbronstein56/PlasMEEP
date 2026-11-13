# Presentation media catalog (PlasMEEP-viz)

**Worktree:** `/home/daq_user/PlasMEEP-viz` · branch `agent/circulator-viz`  
**Default format:** 16:9 · 1920×1080 · 30 fps · H.264 MP4 · no audio · PNG poster · **looping GIF** (`-loop 0`)  
**Field movies:** ~13 s (2 visual RF periods + 0.5 s holds) · animate `Re[Hz e^{-iωt}]`  
**Source policy:** production spectra / histories / design vectors + **production FEM field exports** — **no FEM/Meep solves inside this worktree**.

## Legacy 01–09

| ID | File | Kind | Status | Script |
|----|------|------|--------|--------|
| 01 | `01_device_circulation_overview/` | **SCHEMATIC** | built | `render_01_device_circulation_overview.py` |
| 02 | `02_uniform_vs_optimized_field_movie/` | **SIMULATION** | rebuilt slower | `render_02_03_real_field_movies.py` |
| 03a/b | `03_circulation_field_movie/` | **SIMULATION** | rebuilt slower | `render_02_03_real_field_movies.py` |
| 04 | `04_inverse_design_evolution/` | **DATA** | built | `render_stage2_progress.py` |
| 05 | `05_broadband_frequency_sweep/` | **DATA** | built | `render_05_broadband_frequency_sweep.py` |
| 06 | `06_adaptive_frequency_enrichment/` | **DATA** | built | `render_06_adaptive_frequency_enrichment.py` |
| 07 | `07_TE0_TE1_modes/` | **ANALYTIC** | built | `render_07_TE0_TE1_modes.py` |
| 08 | `08_validation_ladder/` | **DATA** | built | `render_08_validation_ladder.py` |
| 09 | `09_stage2_progress_latest/` | **DATA** | built | `render_stage2_progress.py` |

## Stage-3 package (primary)

| ID | Artifact | Kind | Demonstrates |
|----|----------|------|--------------|
| 10 | `10_stage3_hole_repair_5p75/` MP4+GIF+poster | **SIMULATION** | Stage-2 vs Stage-3 real Hz @ **5.75 GHz** |
| 11 | `11_stage3_hole_repair_5p825/` MP4+GIF+poster | **SIMULATION** | Stage-2 vs Stage-3 real Hz @ **5.825 GHz** |
| 12 | `12_stage3_4p1_improvement/` MP4+GIF+poster | **SIMULATION** | Stage-2 vs Stage-3 real Hz @ **4.1 GHz** |
| 13 | `13_stage2_vs_stage3_broadband_sweep/` + `stage3_figures/13_*.png` | **DATA** | Spectral S2→S3 comparison; marks repaired + new holes |
| 14 | `stage3_figures/14_stage3_hole_repair_5p75.png` | **DATA** | Quantitative 5.75 repair bars |
| 15 | `stage3_figures/15_stage3_hole_repair_5p825.png` | **DATA** | Quantitative 5.825 repair bars |
| 16 | `stage3_figures/16_stage3_4p1_improvement.png` | **DATA** | 4.1 GHz improvement |
| 17 | `stage3_figures/17_stage3_device_summary.png` | **DATA** | Lab-meeting / poster summary |
| 18 | `18_stage3_optimization_evolution/` + PNG | **DATA** | Stage-2→Stage-3 checkpoint evolution |

Scripts: `render_stage3_field_movies.py`, `render_stage3_static_figures.py`, `render_stage3_spectral_sweep.py`, `render_stage3_evolution.py`.

## Stage-4 package

| ID | Artifact | Kind | Demonstrates |
|----|----------|------|--------------|
| 19 | `19_stage3_vs_stage4_broadband_sweep/` + `stage4_figures/19_*.png` | **DATA** | Spectral S3→S4; marks improved notches + new holes |
| 20 | `20_stage4_notch_5p8125/` MP4+GIF+poster + PNG | **SIMULATION** / **DATA** | 5.8125 before/after (improved, not fully repaired) |
| 21 | `21_stage4_notch_4p775/` MP4+GIF+poster + PNG | **SIMULATION** / **DATA** | 4.775 before/after (node repaired) |
| 22 | `22_stage4_new_hole_5p08125/` MP4+GIF+poster + PNG | **SIMULATION** / **DATA** | True dense-audit worst @ **5.08125 GHz** |
| 23 | `stage4_figures/23_stage4_device_summary.png` | **DATA** | Final Stage-4 device / spectrum summary |
| 24 | `24_stage4_new_hole_4p725/` + `stage4_figures/24_*.png` | **SIMULATION** / **DATA** | Secondary new hole @ 4.725 GHz |

Scripts: `render_stage4_field_movies.py`, `render_stage4_static_figures.py`, `render_stage4_spectral_sweep.py`.

**Scientific note:** Stage-4 improved 5.8125/4.775 on-grid but dense audit found new off-grid corpses at 5.08125 and 4.725. Do not treat Stage-4 as continuum-complete.

## Production field exports used

`/home/daq_user/PlasMEEP/outputs/circulator_inverse_design/jesse_method_campaign/presentation_fields/`

- Stage-2: 4.1, 5.65, 5.75, 5.825 GHz  
- Stage-3: 4.1, 4.3, 4.725, 4.775, 5.08125, 5.75, 5.8125, 5.825 GHz  
- Stage-4: 4.725, 4.775, 5.08125, 5.8125 GHz  

## Playback / motion policy (Stage-3 update)

- Target duration **12–15 s** @ 30 fps for substantive field movies.  
- Default: **2** visual periods × **180** frames/period + **15**-frame holds → **13.0 s**.  
- GIF previews: 6 fps, ≤72 frames, **`-loop 0`** (infinite).  
- Always ship **MP4 + looping GIF + PNG poster** for important results.

## Provenance badges

Every frame: **SCHEMATIC / DATA / ANALYTIC / SIMULATION**.
