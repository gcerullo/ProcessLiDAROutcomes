# Prepare SBE plot polygons and (optionally) RIL CHM GeoTIFFs for `Sabah_recovery_for_Gianluca.R`
# Written by Toby Jackson, 17th April 2026 (tobydjackson@gmail.com)
#
# Original paths pointed at a local Dropbox tree; this version uses the project layout:
#   RawData/SBE_data/plots/   — vectors (.fgb, .shp)
#   RawData/SBE_data/chms/   — outputs: ril_chms/*.tiff ; place source RIL rasters under chms/ril_source/
#
# Run from ProcessLiDAROutcomes project root, then run `02_determine_lidar_change_2013_2020.R` in this folder.

library(tidyverse)
library(magrittr)
library(readxl)
library(ggplot2)
library(MetBrewer)
library(stats)
library(terra)
library(sf)
library(patchwork)

sbe <- function(...) file.path("RawData", "SBE_data", ...)

# --- helper (must be defined before use) ------------------------------------
filter_polygons_by_coverage <- function(A, B, threshold = 0.95,
                                        id_col = NULL,
                                        make_valid = TRUE,
                                        return_coverage = TRUE) {
  if (make_valid) {
    A <- st_make_valid(A)
    B <- st_make_valid(B)
  }
  if (st_crs(A) != st_crs(B)) {
    B <- st_transform(B, st_crs(A))
  }
  if (is.null(id_col)) {
    A <- A %>% mutate(.tmp_id = row_number())
    id_col <- ".tmp_id"
  }
  A <- A %>%
    mutate(.area_A = as.numeric(st_area(geometry)))
  B_union <- st_union(B)
  coverage <- st_intersection(A %>% select(all_of(id_col), .area_A), B_union) %>%
    mutate(.area_overlap = as.numeric(st_area(geometry))) %>%
    st_drop_geometry() %>%
    group_by(across(all_of(id_col)), .area_A) %>%
    summarise(.area_overlap = sum(.area_overlap), .groups = "drop") %>%
    mutate(prop_covered = .area_overlap / .area_A)
  A_out <- A %>%
    left_join(coverage %>% select(all_of(id_col), prop_covered), by = id_col) %>%
    mutate(prop_covered = coalesce(prop_covered, 0))
  A_out <- A_out %>%
    filter(prop_covered >= threshold)
  if (!return_coverage) {
    A_out <- A_out %>% select(-prop_covered)
  }
  if (id_col == ".tmp_id") {
    A_out <- A_out %>% select(-.tmp_id)
  }
  A_out
}

# --- mirror prior Dropbox layout under RawData/SBE_data/plots ---------------
# Place these files (or update paths) after copying from your SBE data share:
#   plots/Plot_locations/SBE_whole_area_4ha_grid.fgb
#   plots/Plot_locations/SBE_plots_corrected_with_GPS_2025/SBE_plots_corrected_with_GPS_2025.shp
#   plots/logging/USFR_compart_sf_TJ/USFR_compart_sf_TJ.shp   (or keep your subfolder names)

plot_loc <- function(...) sbe("plots", "Plot_locations", ...)

sbe_4ha <- terra::vect(plot_loc("SBE_whole_area_4ha_grid.fgb")) %>%
  terra::project("epsg:32650") %>%
  sf::st_as_sf()

sbe_logging <- terra::vect(
  sbe("plots", "logging", "USFR_compart_sf_TJ", "USFR_compart_sf_TJ.shp")
) %>%
  terra::project("epsg:32650") %>%
  sf::st_as_sf()

sbe_once_logged <- sbe_logging %>% dplyr::filter(is.na(logg_2) == 1) %>% st_union()
sbe_twice_logged <- sbe_logging %>% dplyr::filter(is.na(logg_2) == 0) %>% st_union()

plots_sbe <- terra::vect(
  plot_loc("SBE_plots_corrected_with_GPS_2025", "SBE_plots_corrected_with_GPS_2025.shp")
) %>%
  terra::project("epsg:32650") %>%
  sf::st_as_sf()

sbe_4ha_not_experiment <- sbe_4ha[lengths(st_intersects(sbe_4ha, plots_sbe)) == 0, ]
sbe_4ha_twice_logged <- filter_polygons_by_coverage(sbe_4ha_not_experiment, sbe_twice_logged, threshold = 0.9)
sbe_4ha_once_logged <- filter_polygons_by_coverage(sbe_4ha_not_experiment, sbe_once_logged, threshold = 0.9)

ggplot() +
  geom_sf(data = sbe_4ha_once_logged, fill = "blue", alpha = 0.5) +
  geom_sf(data = sbe_4ha_twice_logged, fill = "red", alpha = 0.5) +
  geom_sf(data = plots_sbe, fill = "grey", alpha = 0.5)

sbe_4ha_once_logged$class <- "once logged - 1960"
sbe_4ha_once_logged$restoration <- "no"
sbe_4ha_twice_logged$class <- "twice logged - 2012"
sbe_4ha_twice_logged$restoration <- "no"
plots_sbe$class <- "once logged - 1960"

plots_sbe$restoration <- NA
plots_sbe$restoration[plots_sbe$treatment == "0_spp"] <- "no"
plots_sbe$restoration[plots_sbe$treatment == "1_spp"] <- "Enrichment 1 species - 2002"
plots_sbe$restoration[plots_sbe$treatment == "4_spp"] <- "Enrichment 4 species - 2002"
plots_sbe$restoration[plots_sbe$treatment == "16_spp"] <- "Enrichment 16 species - 2002"
plots_sbe$restoration[plots_sbe$treatment == "16_cut" & plots_sbe$block == "South"] <- "Liana cutting - 2011"
plots_sbe$restoration[plots_sbe$treatment == "16_cut" & plots_sbe$block == "North"] <- "Liana cutting - 2014"

sbe_4ha_combined <- bind_rows(sbe_4ha_once_logged, sbe_4ha_twice_logged, plots_sbe)
sbe_4ha_combined %<>% select(class, restoration, geometry)
sbe_4ha_combined$id <- seq_len(nrow(sbe_4ha_combined))

out_fgb <- sbe("plots", "SBE_4ha_combined.fgb")
dir.create(dirname(out_fgb), recursive = TRUE, showWarnings = FALSE)
sf::st_write(sbe_4ha_combined, out_fgb, delete_dsn = TRUE, quiet = TRUE)
message("Wrote ", normalizePath(out_fgb, winslash = "/", mustWork = FALSE))

# --- RIL CHMs: read source GeoTIFFs, write aligned stack to chms/ril_chms -----
# Place source files under RawData/SBE_data/chms/ril_source/ (same basenames as below), or edit paths.
# Expected layout (copy from your ALS folder):  ril_source/2013/*.tif  ril_source/2020/*.tif
ril_src <- sbe("chms", "ril_source")
ril_2013_chm <- terra::rast(
  file.path(ril_src, "2013", "chm_lspikefree_multi3.1_slope1.75_offset2.1.tif")
)
ril_2020_chm <- terra::rast(
  file.path(ril_src, "2020", "chm_lspikefree_multi3.1_slope1.75_offset2.1.tif")
)
ril_2020_dtm <- terra::rast(file.path(ril_src, "2020", "dtm.tif"))

ril_2013_cropped <- crop(ril_2013_chm, ril_2020_chm)
ril_2013_resampled <- resample(ril_2013_cropped, ril_2020_chm, method = "bilinear")
ril_chm_diff <- ril_2020_chm - ril_2013_resampled

names(ril_2013_resampled) <- "chm_2013"
names(ril_2020_chm) <- "chm_2020"
names(ril_2020_dtm) <- "dtm"
names(ril_chm_diff) <- "chm_diff"

ril_out <- sbe("chms", "ril_chms")
dir.create(ril_out, recursive = TRUE, showWarnings = FALSE)

terra::writeRaster(ril_2013_resampled, filename = file.path(ril_out, "chm_2013.tiff"), overwrite = TRUE)
terra::writeRaster(ril_2020_chm, filename = file.path(ril_out, "chm_2020.tiff"), overwrite = TRUE)
terra::writeRaster(ril_2020_dtm, filename = file.path(ril_out, "dtm.tiff"), overwrite = TRUE)
terra::writeRaster(ril_chm_diff, filename = file.path(ril_out, "chm_diff.tiff"), overwrite = TRUE)
message("Wrote RIL CHM stack under ", normalizePath(ril_out, winslash = "/", mustWork = FALSE))
