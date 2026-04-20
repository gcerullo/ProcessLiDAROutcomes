# Nature Revision (NR2) — R workflow

This folder contains the **Nature revision** analysis scripts for megatree (tall-canopy) models, logging–recovery scenarios, and **Sabah Biodiversity Experiment (SBE)** LiDAR sensitivity work. Paths are written for the **`ProcessLiDAROutcomes`** project: set the R working directory to the **repository root** before sourcing anything (same assumption as RStudio “Session → Set Working Directory → To Project Source”).

The **`old/`** subfolder holds superseded one-off scripts; they are **not** part of the active pipeline and are ignored below.

---

## Configuration

**`config.R`** defines shared locations and helpers:

| Symbol | Role |
|--------|------|
| `nr2_config` | `raw_data_dir`, `input_dir`, `output_dir` (`Outputs/NR2`), `figures_dir`, `models_dir` |
| `nr2_input_path(...)` | Files under `Inputs/` |
| `nr2_output_path(...)` | Files under `Outputs/NR2/` |
| `nr2_raw_path(...)` | Files under `RawData/TreeHeightDanum/` |

It creates `Outputs/NR2` (and figures/models subfolders) if missing.

---

## One-shot driver

**`run_all.R`** sources, in order:

1. **`config.R`**
2. **`01_full_NAture_megatrees_model.R`** — habitat-specific zero-inflated beta-binomial megatree models at 45, 50, and 55 m; posterior summaries, PPCs, and saved **`.rds`** models under `Outputs/NR2/models/`.
3. **`02_full_Nature_process_scenarios.R`** — propagates megatree posterior draws through fixed scenarios (`Inputs/FixedScenarioParmams.R`, `MasterAllScenarios.rds`, carbon tables); writes combined scenario performance tables to **`Outputs/NR2/`**.
4. **`SBE_sensitivity_analysis/01_prep_polygons_and_rasters.R`** — builds combined SBE plot polygons and optional aligned RIL CHM GeoTIFFs under **`RawData/SBE_data/`** (requires plot vectors and, for RIL, rasters under `chms/ril_source/`).
5. **`SBE_sensitivity_analysis/02_determine_lidar_change_2013_2020.R`** — extracts 2013/2020 CHM statistics per plot/hex for SBE, Danum, and RIL; writes **`Outputs/NR2/Sabah_logging_recovery_data_LiDAR.csv`** (needs CHM stacks under `RawData/SBE_data/chms/` and outputs from step 4 where applicable).
6. **`SBE_sensitivity_analysis/03_lidar_SBE_carbon_validation.R`** — implied carbon sensitivity from mean canopy height (2013 vs 2020); tables and figures under **`Outputs/NR2/`** (including e.g. **`sbe_uncertainty_carbon.csv`**).
7. **`SBE_sensitivity_analysis/04_lidar_SBE_megatree_sensitivity.R`** — megatree *pixel share* sensitivity at 45/50/55 m (\(p = n(\mathrm{CH}>h)/N_\mathrm{total}\)); CSVs and figures under **`Outputs/NR2/`**.

From the project root:

```r
source(file.path("Scripts", "Nature_Revision", "run_all.R"))
```

**Note:** Items 4–7 in the list above need **Sabah LiDAR inputs** under `RawData/SBE_data/`. If those are absent, the SBE block will fail while the megatree + scenario block (items 2–3) may still succeed if their prerequisites exist.

---

## Main scripts (by file)

### `01_full_NAture_megatrees_model.R`

Fits **brms** zero-inflated beta-binomial models per habitat and height threshold. Reads hectare-scale CHM summary **`.rds`** objects from **`Outputs/NR2/`** (e.g. `primary1haCHMvals.rds`, `onceLogged_vals1haCHMvals.rds`, `restored_vals1haCHMvals.rds`, `twiceLogged_vals1haCHMvals.rds`). These are produced elsewhere in the wider project, not by `run_all.R` alone. Writes model fits, summaries, habitat-level draws, and PPC outputs under **`Outputs/NR2/`** and **`Outputs/NR2/models/`**.

### `02_full_Nature_process_scenarios.R`

Joins scenario timelines to carbon parameters and megatree posterior structure; exports long-format scenario × threshold performance (**`.rds`** and **`.csv`** in **`Outputs/NR2/`**).

### `03_full_Nature_predictive_scenarios.R`

**Not** sourced by `run_all.R`. Optional **predictive** scenario path: averages over multiple representative hectare-like units per posterior draw before aggregation (see header comments inside the file). Uses saved models from **`Outputs/NR2/models/`** and **`Inputs/`** scenario helpers. Run manually when you need that middle-ground predictive workflow.

### `SBE_sensitivity_analysis/`

| Script | Purpose |
|--------|--------|
| `01_prep_polygons_and_rasters.R` | Polygon prep and RIL CHM alignment; **Toby Jackson** authorship in header. |
| `02_determine_lidar_change_2013_2020.R` | Raster extraction and plot-level summaries; **Toby Jackson** authorship in header. |
| `03_lidar_SBE_carbon_validation.R` | Carbon-from-height sensitivity and NR2-facing tables/figures. RIL is kept in by default, coded as `once_logged_ril`, with a `model_habitat` column (`once_logged`, etc.) for comparison to once-logged model predictions. |
| `04_lidar_SBE_megatree_sensitivity.R` | Megatree share and Δ\(p\) sensitivity vs 2013/2020; same RIL / `model_habitat` behaviour as the carbon script. |
| `05_sbe_lidar_maps_bars_composite.R` | Map row (three sites + CHM scale), then grey carbon + megatree (50 m) bar rows and bottom class legend. Not in `run_all.R`. |
| `06_mapping_repeat_LIDAR.R` | **Three** map panels (Primary Danum, SBE, RIL) + CHM legend: 2020 CHM greyscale, outlines use **`sbe_palette_by_display()`**; polygons with **no overlapping CHM pixels are dropped** (e.g. empty RIL tiles). Requires **`patchwork`** and (recommended) **`stars`**. Not in `run_all.R`. |
| `sbe_lidar_plot_class_labels.R` | Shared `plot_class` / display names / megatree habitat strings — sourced by the carbon, megatree, map, and composite scripts (not run on its own). |
| `sbe_lidar_site_maps.R` | Helpers to load SBE + Danum + RIL polygons with aligned labels and to build one LiDAR + outlines `ggplot` panel per site. |

SBE scripts assume **`source(..., "config.R")`**-compatible paths from the **project root** (or call `config.R` first, as each runnable script does).

---

## Dependencies (high level)

- **NR2 core:** `tidyverse`, `brms`, `bayesplot`, `tidyr`, `ggplot2`, `data.table`, `dplyr`, `ggpubr`, `stringr`, `cowplot`, `purrr`, etc. (see each script’s `library()` calls).
- **SBE prep / extraction:** `terra`, `sf`, `dplyr`, `tidyr`, `magrittr`, and additional packages in `01_prep_polygons_and_rasters.R` (`readxl`, `MetBrewer`, `patchwork`, …).
- **SBE maps / composite figure:** `stars` (recommended for `geom_stars` alignment), `patchwork` (three-site map row and optional full composite with barplots).

---

## Outputs

Primary artefacts live under **`Outputs/NR2/`** (and **`Outputs/NR2/figures/`** where scripts write plots). Sabah extraction targets **`Sabah_logging_recovery_data_LiDAR.csv`** in that same NR2 output tree so the sensitivity scripts share one canonical table.

---

## `old/`

Legacy NR2 scripts kept for reference only; **`run_all.R` does not source them.** Use the numbered `01`–`03` and `SBE_sensitivity_analysis/` scripts above for replication.
