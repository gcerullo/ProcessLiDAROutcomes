
# This script extracts canopy height values from rasters for 4 ha plots in Sabah
# The aim is to assess recovery after logging and restoration
# Written by Toby Jackson, 17th April 2026 (tobydjackson@gmail.com)
# (Renamed from `Sabah_recovery_for_Gianluca.R`; run `01_prep_polygons_and_rasters.R` first when building from raw inputs.)
#
# Paths: run from the ProcessLiDAROutcomes project root (RStudio / setwd).
# Inputs:  RawData/SBE_data/chms/{sbe_chms,dan_chms,ril_chms}/*.tiff
#          RawData/SBE_data/plots/{vector assets — see paths below}
# Output:  Outputs/NR2/Sabah_logging_recovery_data_LiDAR.csv  (NR2 sensitivity scripts read this)

source(file.path("Scripts", "Nature_Revision", "config.R"))

library(terra)
library(sf)
library(dplyr)
library(magrittr)
library(tidyr)

sbe <- function(...) file.path("RawData", "SBE_data", ...)

# SBE ########

# Load in rasters as stack
sbe_rasters <- terra::rast(list.files(sbe("chms", "sbe_chms"), pattern = "\\.tiff?$", full.names = TRUE))

# Load in 4 ha plots and add ID column
sbe_4ha <- terra::vect(sbe("plots", "SBE_4ha_combined.fgb"))
sbe_4ha$ID <- seq_len(nrow(sbe_4ha))

# Extract raster pixels by each plot
dt_sbe <- terra::extract(sbe_rasters, sbe_4ha)

# Re-join to plot data by ID column
dt_sbe <- left_join(st_drop_geometry(st_as_sf(sbe_4ha)), dt_sbe, by = "ID")

# Add name and reduce columns
dt_sbe$site <- "Sabah Biodiversity Experiment"
dt_sbe %<>% select(site, class, restoration, ID, dtm, chm_2013, chm_2020, chm_diff) %>% drop_na()

# Summarize
dts_sbe <- dt_sbe %>%
  group_by(site, class, restoration, ID) %>%
  summarize(
    elevation = mean(dtm, na.rm = TRUE),
    total_pixels = n(),
    canopy_height_2013_mean = mean(chm_2013, na.rm = TRUE),
    pixels_over_5m_in_2013 = sum(chm_2013 > 5),
    pixels_over_45m_in_2013 = sum(chm_2013 > 45),
    pixels_over_50m_in_2013 = sum(chm_2013 > 50),
    pixels_over_55m_in_2013 = sum(chm_2013 > 55),
    canopy_height_2020_mean = mean(chm_2020, na.rm = TRUE),
    pixels_over_5m_in_2020 = sum(chm_2020 > 5),
    pixels_over_45m_in_2020 = sum(chm_2020 > 45),
    pixels_over_50m_in_2020 = sum(chm_2020 > 50),
    pixels_over_55m_in_2020 = sum(chm_2020 > 55),
    .groups = "drop"
  )


# Danum ###########
dan_rasters <- terra::rast(list.files(sbe("chms", "dan_chms"), pattern = "\\.tiff?$", full.names = TRUE))
dan_4ha <- terra::vect(sbe("plots", "Dan_4ha_grid_clipped.shp"))
dan_4ha$ID <- seq_len(nrow(dan_4ha))
dt_dan <- terra::extract(dan_rasters, dan_4ha)
dt_dan <- left_join(st_drop_geometry(st_as_sf(dan_4ha)), dt_dan, by = "ID")
dt_dan %<>% mutate(site = "Danum", class = "Primary forest", restoration = "none")
dt_dan %<>% select(site, class, restoration, ID, dtm, chm_2013, chm_2020, chm_diff) %>% drop_na()

dts_dan <- dt_dan %>%
  group_by(site, class, restoration, ID) %>%
  summarize(
    elevation = mean(dtm, na.rm = TRUE),
    total_pixels = n(),
    canopy_height_2013_mean = mean(chm_2013, na.rm = TRUE),
    pixels_over_5m_in_2013 = sum(chm_2013 > 5),
    pixels_over_45m_in_2013 = sum(chm_2013 > 45),
    pixels_over_50m_in_2013 = sum(chm_2013 > 50),
    pixels_over_55m_in_2013 = sum(chm_2013 > 55),
    canopy_height_2020_mean = mean(chm_2020, na.rm = TRUE),
    pixels_over_5m_in_2020 = sum(chm_2020 > 5),
    pixels_over_45m_in_2020 = sum(chm_2020 > 45),
    pixels_over_50m_in_2020 = sum(chm_2020 > 50),
    pixels_over_55m_in_2020 = sum(chm_2020 > 55),
    .groups = "drop"
  )


# Reduced Impact Logging ######
ril_rasters <- terra::rast(list.files(sbe("chms", "ril_chms"), pattern = "\\.tiff?$", full.names = TRUE))
ril_4ha <- sf::st_read(sbe("plots", "RIL_whole_area_4ha_grid.fgb"))
ril_4ha$ID <- seq_len(nrow(ril_4ha))
dt_ril <- terra::extract(ril_rasters, ril_4ha)
dt_ril <- left_join(st_drop_geometry(st_as_sf(ril_4ha)), dt_ril, by = "ID")
dt_ril %<>% mutate(site = "Reduced impact logging", class = "NA", restoration = "none")
dt_ril %<>% select(site, class, restoration, ID, dtm, chm_2013, chm_2020, chm_diff) %>% drop_na()

dts_ril <- dt_ril %>%
  group_by(site, class, restoration, ID) %>%
  summarize(
    elevation = mean(dtm, na.rm = TRUE),
    total_pixels = n(),
    canopy_height_2013_mean = mean(chm_2013, na.rm = TRUE),
    pixels_over_5m_in_2013 = sum(chm_2013 > 5),
    pixels_over_45m_in_2013 = sum(chm_2013 > 45),
    pixels_over_50m_in_2013 = sum(chm_2013 > 50),
    pixels_over_55m_in_2013 = sum(chm_2013 > 55),
    canopy_height_2020_mean = mean(chm_2020, na.rm = TRUE),
    pixels_over_5m_in_2020 = sum(chm_2020 > 5),
    pixels_over_45m_in_2020 = sum(chm_2020 > 45),
    pixels_over_50m_in_2020 = sum(chm_2020 > 50),
    pixels_over_55m_in_2020 = sum(chm_2020 > 55),
    .groups = "drop"
  )

dts <- bind_rows(dts_sbe, dts_dan, dts_ril)

out_csv <- nr2_output_path("Sabah_logging_recovery_data_LiDAR.csv")
readr::write_csv(dts, file = out_csv)
message("Wrote ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))
