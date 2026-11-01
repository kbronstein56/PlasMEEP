# Presentation media catalog (PlasMEEP-viz)

**Worktree:** `/home/daq_user/PlasMEEP-viz` · branch `agent/circulator-viz`  
**Default format:** 16:9 · 1920×1080 · 30 fps · H.264 MP4 · no audio · PNG poster · GIF preview where practical  
**Source policy:** existing spectra / histories / design vectors / presentation figures only — **no FEM/Meep solves** in this pipeline.

| ID | File | Kind | Status | Script |
|----|------|------|--------|--------|
| 01 | `01_device_circulation_overview/01_device_circulation_overview.mp4` | **SCHEMATIC** | built | `render_01_device_circulation_overview.py` |
| 02 | `02_field_movies_stubs/` | SIMULATION (awaiting exports) | **stub only** | `render_02_03_field_movie_stubs.py` |
| 03 | `03_field_movies_stubs/` | SIMULATION (awaiting exports) | **stub only** | `render_02_03_field_movie_stubs.py` |
| 04 | `04_inverse_design_evolution/04_inverse_design_evolution.mp4` | **DATA** (Stage-2 accepted checkpoints) | built | `render_stage2_progress.py` |
| 05 | `05_broadband_frequency_sweep/05_broadband_frequency_sweep_stage1.mp4` (+ stage2 comparison alias) | **DATA** | built | `render_05_broadband_frequency_sweep.py` |
| 06 | `06_adaptive_frequency_enrichment/06_adaptive_frequency_enrichment.mp4` | **DATA** (Stage-1 audits) | built | `render_06_adaptive_frequency_enrichment.py` |
| 07 | `07_TE0_TE1_modes/07_TE0_TE1_modes.mp4` | **ANALYTIC** (40 mm guide) | built | `render_07_TE0_TE1_modes.py` |
| 08 | `08_validation_ladder/08_validation_ladder.mp4` | **DATA** (documented ladder) | built | `render_08_validation_ladder.py` |
| 09 | `09_stage2_progress_latest/09_stage2_progress_latest.mp4` (+ `.png`) | **DATA** (final Stage-2) | built | `render_stage2_progress.py` |

## Read-only production inputs

Under `/home/daq_user/PlasMEEP/outputs/circulator_inverse_design/jesse_method_campaign/`:

- `bb_stage2_history.json`, `opt_broadband_stage2_A2_C.json`, `bb_stage2_checkpoint.json`
- `opt_broadband_production_C.json`, `bb_production_history.json`
- `spectral_*.json`, `bb_final_audit_*.json`, `spectral_grid_comparison.json`
- enrichment docs: `BB_ENRICHMENT_RECOMMENDATION.md`, `BROADBAND_FREQUENCY_SAMPLING_STUDY.md`

Presentation figures: `/home/daq_user/PlasMEEP/outputs/validation/fem_meep_validation/presentation_figures/`

## Conventions (must stay consistent)

- Ports **P1…P6** at 0°, 60°, … (code `port_dirs`).
- For source `j` and **B_PLUS**: desired `(j+1)%6`, reverse `(j−1)%6`.
- Animate **Re[Hz exp(−iωt)]** for complex fields; never oscillate **|H|**.
- Badge every frame: SCHEMATIC / DATA / ANALYTIC / SIMULATION.

## Rebuild

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate plasmeep_mpi
cd /home/daq_user/PlasMEEP-viz
python scripts/presentation_media/run_all_presentation_media.py
```

Large MP4s and `frames/` dumps are **not** meant for git; commit scripts, docs, posters, and small GIFs.
