# =============================================================================
# Composite: row1 = three-site LiDAR maps (+ CHM scale); row2 = carbon + megatree
# bars; row3 = land-class legend (horizontal, full width, bottom).
# -----------------------------------------------------------------------------
# Run from project root *after*:
#   lidar_SBE_class_summary.csv, lidar_SBE_megatree_sensitivity_summary.csv,
#   lidar_SBE_plot_level.csv (for RIL + pooled megatree rows if missing from
#   megatree summary), CHMs + polygons.
#
# Bar y-axis uses `sbe_figure_bar_classes()` (map strata + pooled enrichment+RIL).
# Maps use `sbe_figure_map_classes()` (SBE control polygons omitted).
# Optional carbon overlay: `Inputs/carbon_recovery__quick_compare_model_vs_philipson.csv`
# (my_model columns), else `Inputs/carbon_recovery__annual_acd_inc__summary.csv`.
# Requires package `patchwork`.
#
# Writes: Outputs/NR2/figures/SBE_LiDAR_maps_and_barplots_composite.png (+ .pdf)
# =============================================================================

source(file.path("Scripts", "Nature_Revision", "config.R"))
source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_lidar_site_maps.R"))
source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_lidar_megatree_from_plot_level.R"))
source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_carbon_recovery_model_for_composite.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install package `patchwork` to build the composite figure.")
}
suppressPackageStartupMessages(library(patchwork))

map_classes <- sbe_figure_map_classes()
bar_classes <- sbe_figure_bar_classes()
pal_map <- sbe_palette_map()
pal_bar <- sbe_palette_bar()

allowed_carbon_pc <- c(
  "primary",
  "enrichment",
  "liana_cutting",
  "twice_logged",
  "once_logged_ril",
  "once_logged_enrichment_ril_pooled"
)

pol <- sbe_load_site_polygons_for_maps()

chm_dan <- sbe_resolve_chm_path("dan_chms")
chm_sbe <- sbe_resolve_chm_path("sbe_chms")
chm_ril <- sbe_resolve_chm_path("ril_chms")

map_target_cells <- 650L
map_title_size <- 9.2

# Per-site square extent (tight zoom); do not share one large window across sites
# or Danum/SBE shrink inside fixed patchwork panel sizes.
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
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Danum grid missing") +
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
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "RIL grid missing") +
    ggplot2::theme_void()
}

p_chm_leg <- sbe_chm_legend_column_ggplot()
p_class_leg <- sbe_class_legend_bottom_ggplot(pal_bar, bar_classes)

maps_row <- (p_dan | p_sbe | p_ril | p_chm_leg) +
  patchwork::plot_layout(ncol = 4L, widths = c(1, 1, 1, 0.14)) &
  ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 2, 4))

# --- Carbon ------------------------------------------------------------------
cs_path <- nr2_output_path("lidar_SBE_class_summary.csv")
if (!file.exists(cs_path)) {
  stop("Missing ", cs_path, " — run 03_lidar_SBE_carbon_validation.R first.")
}
cs_raw <- readr::read_csv(cs_path, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(.data$mean_implied_delta_carbon_MgC_ha_yr))

if ("plot_class" %in% names(cs_raw)) {
  cs <- cs_raw %>%
    dplyr::mutate(
      plot_class_display = sbe_plot_class_display(as.character(.data$plot_class))
    ) %>%
    dplyr::filter(.data$plot_class %in% .env$allowed_carbon_pc)
} else {
  cs <- cs_raw %>%
    dplyr::mutate(plot_class_display = as.character(.data$plot_class_display)) %>%
    dplyr::filter(.data$plot_class_display %in% .env$bar_classes)
}

if (!any(cs$plot_class == "once_logged_enrichment_ril_pooled", na.rm = TRUE)) {
  pr <- sbe_carbon_pooled_enrichment_ril_row(cs)
  if (!is.null(pr)) {
    cs <- dplyr::bind_rows(cs, pr)
  }
}

cs <- cs %>%
  dplyr::filter(.data$plot_class_display %in% .env$bar_classes) %>%
  dplyr::mutate(
    ymin = .data$lwr_95_MgC_ha_yr,
    ymax = .data$upr_95_MgC_ha_yr,
    cls = factor(.data$plot_class_display, levels = rev(.env$bar_classes))
  )

model_pts <- sbe_carbon_model_points_for_composite(bar_classes)
if (nrow(model_pts) > 0L) {
  model_pts <- sbe_carbon_model_apply_x_offset(model_pts)
}
model_col <- sbe_carbon_model_overlay_colour()

# --- Megatree (50 m) + plot-level fill-ins -----------------------------------
mg_path <- nr2_output_path("lidar_SBE_megatree_sensitivity_summary.csv")
if (!file.exists(mg_path)) {
  stop("Missing ", mg_path, " — run 04_lidar_SBE_megatree_sensitivity.R first.")
}
mg_raw <- readr::read_csv(mg_path, show_col_types = FALSE)

mg50_base <- mg_raw %>%
  dplyr::filter(.data$height_m == 50L, !is.na(.data$mean_delta_prop)) %>%
  dplyr::mutate(
    disp = if ("habitat_display" %in% names(mg_raw)) {
      as.character(.data$habitat_display)
    } else {
      sbe_habitat_to_display(as.character(.data$habitat))
    }
  )

plot_level_path <- nr2_output_path("lidar_SBE_plot_level.csv")
mg_extra <- sbe_megatree_extra_rows_from_plot_level(plot_level_path)
if (nrow(mg_extra) > 0L) {
  have_h <- unique(as.character(mg50_base$habitat))
  add <- mg_extra %>%
    dplyr::filter(!.data$habitat %in% have_h) %>%
    dplyr::mutate(disp = sbe_habitat_to_display(as.character(.data$habitat)))
  if (nrow(add) > 0L) {
    mg50_base <- dplyr::bind_rows(mg50_base, add)
  }
}

mg7 <- mg50_base %>%
  dplyr::filter(.data$disp %in% .env$bar_classes) %>%
  dplyr::distinct(.data$disp, .keep_all = TRUE) %>%
  dplyr::mutate(cls = factor(.data$disp, levels = rev(.env$bar_classes)))

n_for_label <- function(short_name) {
  nc <- cs$n_plots[match(short_name, cs$plot_class_display)]
  nm <- mg7$n_plots[match(short_name, mg7$disp)]
  nc <- if (length(nc) == 0L || is.na(nc)) 0L else as.integer(nc[[1L]])
  nm <- if (length(nm) == 0L || is.na(nm)) 0L else as.integer(nm[[1L]])
  as.integer(max(nc, nm, na.rm = TRUE))
}

y_labels <- stats::setNames(
  vapply(bar_classes, function(s) paste0(s, "\n(n = ", n_for_label(s), ", 4 ha plots)"), character(1L)),
  bar_classes
)

theme_bar <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "plain", size = 10.2, hjust = 0),
    legend.position = "none",
    plot.margin = ggplot2::margin(2, 4, 2, 4)
  )

theme_carbon <- if (nrow(model_pts) > 0L) {
  theme_bar +
    ggplot2::theme(
      legend.position = c(0.99, 0.02),
      legend.justification = c(1, 0),
      legend.background = ggplot2::element_rect(fill = "white", colour = "grey70", linewidth = 0.25),
      legend.margin = ggplot2::margin(3, 5, 3, 5),
      legend.text = ggplot2::element_text(size = 7.5),
      legend.key.size = grid::unit(0.75, "lines")
    )
} else {
  theme_bar
}

p_carbon <- ggplot2::ggplot(cs, ggplot2::aes(x = .data$mean_implied_delta_carbon_MgC_ha_yr, y = .data$cls)) +
  ggplot2::geom_col(width = 0.72, fill = "#8a8a8a", colour = "grey35", linewidth = 0.25, alpha = 0.95) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      x = .data$mean_implied_delta_carbon_MgC_ha_yr,
      y = .data$cls,
      xmin = .data$ymin,
      xmax = .data$ymax
    ),
    orientation = "y",
    width = 0.2,
    linewidth = 0.32,
    colour = "grey25",
    inherit.aes = FALSE
  )

if (nrow(model_pts) > 0L) {
  pos_y <- ggplot2::position_nudge(y = sbe_carbon_model_y_nudge())
  p_carbon <- p_carbon +
    ggplot2::geom_errorbar(
      data = model_pts,
      ggplot2::aes(
        x = .data$mean_plot,
        y = .data$cls,
        xmin = .data$ymin_plot,
        xmax = .data$ymax_plot,
        colour = "My model"
      ),
      orientation = "y",
      inherit.aes = FALSE,
      width = 0.22,
      linewidth = 0.55,
      position = pos_y
    ) +
    ggplot2::geom_point(
      data = model_pts,
      ggplot2::aes(x = .data$mean_plot, y = .data$cls, colour = "My model"),
      inherit.aes = FALSE,
      size = 2.65,
      shape = 16,
      position = pos_y
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c("My model" = model_col),
      drop = FALSE
    )
}

p_carbon <- p_carbon +
  ggplot2::scale_y_discrete(
    drop = FALSE,
    labels = function(brk) vapply(as.character(brk), function(b) y_labels[[b]], character(1L))
  ) +
  ggplot2::labs(
    title = "LiDAR-derived implied carbon change",
    y = "Four-hectare unit",
    x = expression("Implied" ~ Delta * C ~ "(" * Mg ~ C ~ ha^{-1} ~ yr^{-1} * ")")
  ) +
  theme_carbon +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7.5))

p_mega <- ggplot2::ggplot(mg7, ggplot2::aes(x = .data$mean_delta_prop, y = .data$cls)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.25) +
  ggplot2::geom_col(width = 0.72, fill = "#8a8a8a", colour = "grey35", linewidth = 0.25, alpha = 0.95) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      x = .data$mean_delta_prop,
      y = .data$cls,
      xmin = .data$lwr_95,
      xmax = .data$upr_95
    ),
    orientation = "y",
    width = 0.2,
    linewidth = 0.32,
    colour = "grey25",
    inherit.aes = FALSE
  ) +
  ggplot2::scale_y_discrete(drop = FALSE, labels = NULL) +
  ggplot2::labs(
    title = "LiDAR megatree share change (p) per year, 50 m threshold",
    y = NULL,
    x = "Mean \u0394p per year ((p2020 \u2212 p2013) / 7)"
  ) +
  theme_bar +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )

if (requireNamespace("scales", quietly = TRUE)) {
  p_mega <- p_mega +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.0001, trim = TRUE))
} else {
  fmt_dp <- function(z) format(z, scientific = FALSE, trim = TRUE, digits = 5)
  p_mega <- p_mega + ggplot2::scale_x_continuous(labels = fmt_dp)
}

bars_row <- (p_carbon | p_mega) +
  patchwork::plot_layout(ncol = 2L, widths = c(1.12, 1)) &
  ggplot2::theme(plot.margin = ggplot2::margin(2, 4, 1, 4))

p_all <- maps_row / bars_row / p_class_leg +
  patchwork::plot_layout(heights = c(1, 0.88, 0.11))

fig_dir <- nr2_config$figures_dir
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
fig_png <- file.path(fig_dir, "SBE_LiDAR_maps_and_barplots_composite.png")
fig_pdf <- file.path(fig_dir, "SBE_LiDAR_maps_and_barplots_composite.pdf")

ggsave(fig_png, p_all, width = 12.8, height = 9.35, dpi = 600, bg = "white")
tryCatch(
  ggsave(fig_pdf, p_all, width = 12.8, height = 9.35, units = "in", device = grDevices::cairo_pdf),
  error = function(e) message("PDF not written: ", conditionMessage(e))
)

message("\nWrote:\n  ", normalizePath(fig_png, winslash = "/", mustWork = FALSE))
message("  ", normalizePath(fig_pdf, winslash = "/", mustWork = FALSE))
