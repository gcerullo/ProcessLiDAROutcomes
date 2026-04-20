# SBE LiDAR sensitivity (Nature revision)

**Numbered scripts `01`–`06`** are the intended **run order**; files named **`sbe_*.R`** are **library-style** helpers (sourced by other scripts, not meant as stand-alone drivers).

R scripts for **Sabah Biodiversity Experiment (SBE)**–related **LiDAR sensitivity**: paired CHMs (2013 vs 2020) on four-hectare units, **implied carbon change** from mean height, **megatree share** change across height thresholds, and **composite figures** (maps + bars). Run everything from the **ProcessLiDAROutcomes project root** (`setwd` / RStudio project folder).

## Configuration

- `../config.R` — NR2 paths (`Outputs/NR2`, `Inputs`, figures). Sourced by each driver script.

## Pipeline (typical order)

1. **`01_prep_polygons_and_rasters.R`** — Prepares plot vectors / rasters when building from raw assets (see script header).
2. **`02_determine_lidar_change_2013_2020.R`** — Extracts CHM statistics per plot; writes **`Outputs/NR2/Sabah_logging_recovery_data_LiDAR.csv`** (downstream scripts read this).
3. **`03_lidar_SBE_carbon_validation.R`** — Implied Delta C (Mg C ha-1 yr-1) from mean height change (5.2 Mg C ha-1 m-1, annualised over 7 yr); stratum means and 95% t-CIs, **`lidar_SBE_class_summary.csv`**, figures.
4. **`04_lidar_SBE_megatree_sensitivity.R`** — Share of pixels above 45 / 50 / 55 m; annualised Delta p per year; stratum means and 95% t-CIs, **`lidar_SBE_megatree_sensitivity_summary.csv`**, **`lidar_SBE_megatree_plot_level_long.csv`**, figures.
5. **`05_sbe_lidar_maps_bars_composite.R`** — Three-site CHM maps + LiDAR carbon and megatree (50 m) bars + class legend; optional **my model** overlay. Requires **patchwork**. **`Outputs/NR2/figures/SBE_LiDAR_maps_and_barplots_composite.png`** (+ PDF).
6. **`06_mapping_repeat_LIDAR.R`** — Map-only three-site figure + CHM legend (optional; same helpers as 05).

## Shared helpers (sourced, not run alone)

| File | Role |
|------|------|
| `sbe_lidar_plot_class_labels.R` | Plot classes, palettes, bar/map class lists, megatree epoch length (`sbe_lidar_megatree_delta_years()`). |
| `sbe_lidar_site_maps.R` | CHM panels, overlap filter, bottom class legend, CHM colour-bar column. |
| `sbe_lidar_megatree_from_plot_level.R` | Extra megatree / pooled carbon rows from **`lidar_SBE_plot_level.csv`** when summaries omit strata. |
| `sbe_carbon_recovery_model_for_composite.R` | Reads quick-compare or annual carbon summary; builds model overlay points for the composite. |

## Inputs / outputs (short)

- **CHMs & plots:** under `RawData/SBE_data/` (see `02_determine_lidar_change_2013_2020.R`).
- **NR2 tables & figures:** under `Outputs/NR2/` and `Outputs/NR2/figures/`.
- **Optional model overlay:** `Inputs/carbon_recovery__quick_compare_model_vs_philipson.csv` (fallback: `Inputs/carbon_recovery__annual_acd_inc__summary.csv`).

For the wider Nature-revision tree, see the repository **`README.md`** (`Scripts/Nature_Revision`).
