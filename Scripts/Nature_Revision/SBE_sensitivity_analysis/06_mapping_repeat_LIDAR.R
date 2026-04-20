# =============================================================================
# Maps: four-hectare units vs LiDAR CHM (Primary Danum, SBE, RIL)
# -----------------------------------------------------------------------------
# One multi-panel figure with the same map-class display names as
# 03_lidar_SBE_carbon_validation.R (see sbe_lidar_plot_class_labels.R).
#
# Run from project root. Requires package `patchwork`.
# Polygons with no overlapping CHM pixels are dropped per site.
#
# Writes: Outputs/NR2/figures/SBE_LiDAR_three_site_maps.png (+ .pdf)
# =============================================================================

source(file.path("Scripts", "Nature_Revision", "config.R"))
source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_lidar_site_maps.R"))

library(ggplot2)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install package `patchwork` to build the three-site map figure.")
}
suppressPackageStartupMessages(library(patchwork))

pol <- sbe_load_site_polygons_for_maps()
pal_map <- sbe_palette_map()
map_classes <- sbe_figure_map_classes()

chm_dan <- sbe_resolve_chm_path("dan_chms")
chm_sbe <- sbe_resolve_chm_path("sbe_chms")
chm_ril <- sbe_resolve_chm_path("ril_chms")

map_target_cells <- 650L
map_title_size <- 9.2

p_dan <- if (!is.null(pol$danum) && nrow(pol$danum) > 0L) {
  sbe_ggplot_lidar_panel(
    pol$danum,
    chm_dan,
    "Primary (Danum)",
    target_cells = map_target_cells,
    pal_named = pal_map,
    map_display_limits = map_classes,
    plot_title_size = map_title_size,
    square_extent = TRUE,
    axis_title_y = "Northing (m)",
    axis_title_x = NULL,
    show_axis_text_x = TRUE,
    show_axis_text_y = TRUE
  )
} else {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Danum 4 ha grid not available", size = 3) +
    ggplot2::theme_void()
}

p_sbe <- sbe_ggplot_lidar_panel(
  pol$sbe,
  chm_sbe,
  "Sabah Biodiversity\nExperiment (SBE)",
  target_cells = map_target_cells,
  pal_named = pal_map,
  map_display_limits = map_classes,
  plot_title_size = map_title_size,
  square_extent = TRUE,
  axis_title_x = "Easting (m)",
  axis_title_y = NULL,
  show_axis_text_x = TRUE,
  show_axis_text_y = TRUE
)

p_ril <- if (!is.null(pol$ril) && nrow(pol$ril) > 0L) {
  sbe_ggplot_lidar_panel(
    pol$ril,
    chm_ril,
    "Reduced-impact\nlogging (RIL)",
    target_cells = map_target_cells,
    pal_named = pal_map,
    map_display_limits = map_classes,
    plot_title_size = map_title_size,
    square_extent = TRUE,
    axis_title_x = NULL,
    axis_title_y = NULL,
    show_axis_text_x = TRUE,
    show_axis_text_y = TRUE
  )
} else {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "RIL 4 ha grid not available", size = 3) +
    ggplot2::theme_void()
}

p_chm_leg <- sbe_chm_legend_column_ggplot()
p_class_leg <- sbe_class_legend_bottom_ggplot(pal_map, map_classes)

p_all <- (p_dan | p_sbe | p_ril | p_chm_leg) / p_class_leg +
  patchwork::plot_layout(heights = c(1, 0.22), widths = c(1, 1, 1, 0.14))

fig_dir <- nr2_config$figures_dir
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fig_png <- file.path(fig_dir, "SBE_LiDAR_three_site_maps.png")
fig_pdf <- file.path(fig_dir, "SBE_LiDAR_three_site_maps.pdf")

h_in <- 5.05
w_in <- 12.4

ggsave(fig_png, p_all, width = w_in, height = h_in, dpi = 200, bg = "white")
tryCatch(
  ggsave(fig_pdf, p_all, width = w_in, height = h_in, units = "in", device = grDevices::cairo_pdf),
  error = function(e) message("PDF not written: ", conditionMessage(e))
)

message("\nWrote:\n  ", normalizePath(fig_png, winslash = "/", mustWork = FALSE))
message("  ", normalizePath(fig_pdf, winslash = "/", mustWork = FALSE))
