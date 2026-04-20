# =============================================================================
# Shared LiDAR / SBE plot-class labels (carbon, megatree, maps)
# -----------------------------------------------------------------------------
# Source from ProcessLiDAROutcomes project root. Defines plot_class codes and
# display strings aligned across 03_lidar_SBE_carbon_validation.R,
# 04_lidar_SBE_megatree_sensitivity.R, and 06_mapping_repeat_LIDAR.R.
# =============================================================================

#' @param site_w,class_w,restoration_w Trimmed character vectors (same length)
sbe_plot_class_from_w <- function(site_w, class_w, restoration_w) {
  dplyr::case_when(
    site_w == "Reduced impact logging" & restoration_w == "none" &
      (class_w == "" | stringr::str_detect(class_w, stringr::regex("^na$", ignore_case = TRUE))) ~ "once_logged_ril",
    site_w == "Danum" &
      stringr::str_detect(class_w, stringr::regex("^primary\\s*forest$", ignore_case = TRUE)) &
      restoration_w == "none" ~ "primary",
    site_w == "Sabah Biodiversity Experiment" &
      stringr::str_detect(class_w, stringr::regex("twice\\s*logged", ignore_case = TRUE)) &
      restoration_w == "no" ~ "twice_logged",
    site_w == "Sabah Biodiversity Experiment" &
      stringr::str_detect(class_w, stringr::regex("once\\s*logged", ignore_case = TRUE)) &
      restoration_w == "no" ~ "sbe_once_logged_control",
    site_w == "Sabah Biodiversity Experiment" &
      stringr::str_detect(class_w, stringr::regex("once\\s*logged", ignore_case = TRUE)) &
      stringr::str_detect(restoration_w, "^Enrichment") ~ "enrichment",
    site_w == "Sabah Biodiversity Experiment" &
      stringr::str_detect(class_w, stringr::regex("once\\s*logged", ignore_case = TRUE)) &
      stringr::str_detect(restoration_w, "^Liana") ~ "liana_cutting",
    TRUE ~ "other"
  )
}

#' Carbon / map legend labels (Title Case; matches class_summary display)
sbe_plot_class_display <- function(plot_class) {
  dplyr::case_when(
    plot_class == "once_logged_ril" ~ "Once-logged (RIL)",
    plot_class == "sbe_once_logged_control" ~ "Once-logged (SBE control)",
    plot_class == "sbe_once_logged_control_plus_enrichment" ~ "Once-logged (control + enrichment pooled)",
    plot_class == "once_logged_enrichment_ril_pooled" ~ "Once-logged (RIL and enrichment) pooled",
    plot_class == "enrichment" ~ "Once-logged (enrichment)",
    plot_class == "liana_cutting" ~ "Restored (liana cut)",
    plot_class == "twice_logged" ~ "Twice logged (SBE)",
    plot_class == "primary" ~ "Primary (Danum)",
    plot_class == "other" ~ "Other",
    TRUE ~ as.character(plot_class)
  )
}

#' Megatree / barplot habitat strings (must match 04_lidar_SBE_megatree_sensitivity)
sbe_megatree_habitat_from_plot_class <- function(plot_class) {
  dplyr::case_when(
    plot_class == "liana_cutting" ~ "restored (liana cut)",
    plot_class == "sbe_once_logged_control" ~ "once-logged (control)",
    plot_class == "enrichment" ~ "once-logged (enrichment)",
    plot_class == "once_logged_ril" ~ "once-logged (RIL)",
    plot_class == "primary" ~ "primary (Danum)",
    plot_class == "twice_logged" ~ "twice logged (SBE)",
    TRUE ~ "other / unclassed"
  )
}

sbe_megatree_habitat_levels <- function() {
  c(
    "primary (Danum)", "once-logged (enrichment)", "once-logged (control)",
    "once-logged (RIL)", "restored (liana cut)", "twice logged (SBE)",
    "other / unclassed"
  )
}

#' Map megatree habitat label to carbon-style display (for composite legends)
sbe_habitat_to_display <- function(habitat) {
  dplyr::case_when(
    habitat == "primary (Danum)" ~ "Primary (Danum)",
    habitat == "once-logged (control)" ~ "Once-logged (SBE control)",
    habitat == "once-logged (enrichment)" ~ "Once-logged (enrichment)",
    habitat == "once-logged (RIL)" ~ "Once-logged (RIL)",
    habitat == "restored (liana cut)" ~ "Restored (liana cut)",
    habitat == "twice logged (SBE)" ~ "Twice logged (SBE)",
    habitat == "once-logged (enrichment + RIL)" ~ "Once-logged (RIL and enrichment) pooled",
    habitat == "once-logged (enrichment + RIL pooled)" ~ "Once-logged (RIL and enrichment) pooled",
    habitat == "other / unclassed" ~ "Other",
    TRUE ~ as.character(habitat)
  )
}

#' Strata shown on maps (polygons only; no pooled footprint; no SBE control).
sbe_figure_map_classes <- function() {
  c(
    "Primary (Danum)",
    "Once-logged (enrichment)",
    "Restored (liana cut)",
    "Twice logged (SBE)",
    "Once-logged (RIL)"
  )
}

# Back-compat name.
sbe_figure_six_classes <- function() {
  sbe_figure_map_classes()
}

#' Barplot + bottom legend strata (map classes plus pooled enrichment+RIL; no SBE control).
sbe_figure_bar_classes <- function() {
  c(
    sbe_figure_map_classes(),
    sbe_plot_class_display("once_logged_enrichment_ril_pooled")
  )
}

# Back-compat (same as bar-class list).
sbe_figure_seven_classes <- function() {
  sbe_figure_bar_classes()
}

#' Map `plot_class` to display used on maps (same as `sbe_plot_class_display`;
#' pooled codes are dropped by `map_display_limits` if not listed).
sbe_map_display_from_plot_class <- function(plot_class) {
  sbe_plot_class_display(plot_class)
}

sbe_palette_map <- function() {
  sbe_palette_by_display()[sbe_figure_map_classes()]
}

sbe_palette_bar <- function() {
  sbe_palette_by_display()[sbe_figure_bar_classes()]
}

# Back-compat.
sbe_palette_six <- function() {
  sbe_palette_map()
}

# Back-compat.
sbe_palette_seven <- function() {
  sbe_palette_bar()
}

#' Ordered display levels for colour scales (maps + composite bars)
sbe_display_level_order <- function() {
  c(
    "Primary (Danum)",
    "Twice logged (SBE)",
    "Once-logged (SBE control)",
    "Once-logged (enrichment)",
    "Restored (liana cut)",
    "Once-logged (RIL)",
    "Once-logged (control + enrichment pooled)",
    "Once-logged (RIL and enrichment) pooled",
    "Once-logged (enrichment + RIL)",
    "Other"
  )
}

#' Years between LiDAR epochs (2013 → 2020); used to express megatree Δp per year.
sbe_lidar_megatree_delta_years <- function() {
  7L
}

#' Okabe–Ito palette named by `sbe_plot_class_display` levels (subset allowed)
sbe_palette_by_display <- function() {
  ord <- sbe_display_level_order()
  pal <- c(
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
    "#D55E00", "#000000", "#F0E442", "#999999", "#CCCCCC"
  )
  stats::setNames(rep(pal, length.out = length(ord)), ord)
}
