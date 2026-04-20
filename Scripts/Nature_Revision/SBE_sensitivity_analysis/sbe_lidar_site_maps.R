# =============================================================================
# Reusable LiDAR + 4 ha polygon map panels (SBE, Danum primary, RIL)
# -----------------------------------------------------------------------------
# Source after setting wd to project root. Auto-sources
# sbe_lidar_plot_class_labels.R if needed.
# =============================================================================

if (!exists("sbe_plot_class_from_w", inherits = TRUE)) {
  source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_lidar_plot_class_labels.R"))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

sbe_raw_path <- function(...) {
  file.path("RawData", "SBE_data", ...)
}

#' First CHM GeoTIFF in chms/<subdir>, preferring filenames matching chm_2020.
sbe_resolve_chm_path <- function(subdir) {
  d <- sbe_raw_path("chms", subdir)
  if (!dir.exists(d)) {
    return(NA_character_)
  }
  fl <- sort(list.files(d, pattern = "\\.tiff?$", full.names = TRUE))
  if (length(fl) == 0L) {
    return(NA_character_)
  }
  hit <- grepl("chm_2020", basename(fl), ignore.case = TRUE)
  if (any(hit)) {
    return(fl[[which(hit)[[1]]]])
  }
  fl[[1L]]
}

#' Load sf polygons with `plot_class`, `plot_class_display`, `panel` columns.
sbe_load_site_polygons_for_maps <- function() {
  sbe_fgb <- sbe_raw_path("plots", "SBE_4ha_combined.fgb")
  if (!file.exists(sbe_fgb)) {
    stop("Missing ", sbe_fgb, " — run 01_prep_polygons_and_rasters.R")
  }
  out <- list()

  out$sbe <- sf::st_read(sbe_fgb, quiet = TRUE) %>%
    sf::st_make_valid() %>%
    dplyr::mutate(
      site_w = "Sabah Biodiversity Experiment",
      class_w = stringr::str_squish(as.character(class)),
      restoration_w = stringr::str_squish(as.character(restoration)),
      panel = "Sabah Biodiversity Experiment (SBE)"
    ) %>%
    dplyr::mutate(
      plot_class = sbe_plot_class_from_w(site_w, class_w, restoration_w),
      plot_class_display = sbe_plot_class_display(plot_class)
    ) %>%
    dplyr::select(panel, plot_class, plot_class_display, geometry)

  dan_shp <- sbe_raw_path("plots", "Dan_4ha_grid_clipped.shp")
  if (file.exists(dan_shp)) {
    out$danum <- sf::st_read(dan_shp, quiet = TRUE) %>%
      sf::st_make_valid() %>%
      dplyr::mutate(
        site_w = "Danum",
        class_w = "Primary forest",
        restoration_w = "none",
        panel = "Primary (Danum)"
      ) %>%
      dplyr::mutate(
        plot_class = sbe_plot_class_from_w(site_w, class_w, restoration_w),
        plot_class_display = sbe_plot_class_display(plot_class)
      ) %>%
      dplyr::select(panel, plot_class, plot_class_display, geometry)
  } else {
    out$danum <- NULL
    message("Danum shapefile not found: ", dan_shp)
  }

  ril_fgb <- sbe_raw_path("plots", "RIL_whole_area_4ha_grid.fgb")
  if (file.exists(ril_fgb)) {
    out$ril <- sf::st_read(ril_fgb, quiet = TRUE) %>%
      sf::st_make_valid() %>%
      dplyr::mutate(
        site_w = "Reduced impact logging",
        class_w = "NA",
        restoration_w = "none",
        panel = "Reduced impact logging (RIL)"
      ) %>%
      dplyr::mutate(
        plot_class = sbe_plot_class_from_w(site_w, class_w, restoration_w),
        plot_class_display = sbe_plot_class_display(plot_class)
      ) %>%
      dplyr::select(panel, plot_class, plot_class_display, geometry)
  } else {
    out$ril <- NULL
    message("RIL grid not found: ", ril_fgb)
  }

  out
}

#' Keep only polygons that overlap at least one non-NA CHM pixel (same CRS as `r_chm`).
sbe_filter_plots_chm_overlap <- function(plots_sf, r_chm) {
  if (is.null(plots_sf) || nrow(plots_sf) == 0L) {
    return(plots_sf)
  }
  r1 <- if (terra::nlyr(r_chm) > 1L) r_chm[[1L]] else r_chm
  v <- terra::vect(sf::st_transform(plots_sf, terra::crs(r1)))
  ex <- tryCatch(
    as.data.frame(terra::extract(r1, v, cells = TRUE)),
    error = function(e) NULL
  )
  if (is.null(ex) || nrow(ex) == 0L) {
    message("Could not raster-extract for CHM overlap filter; plotting all polygons.")
    return(plots_sf)
  }
  lyr <- names(r1)[1L]
  if (!lyr %in% names(ex)) {
    cand <- setdiff(names(ex), c("ID", "cell", "x", "y"))
    lyr <- if (length(cand) > 0L) cand[[1L]] else return(plots_sf)
  }
  idcol <- if ("ID" %in% names(ex)) {
    "ID"
  } else if ("id" %in% names(ex)) {
    "id"
  } else {
    message("Unexpected extract() columns; skip CHM overlap filter.")
    return(plots_sf)
  }
  keep_ids <- ex %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(idcol))) %>%
    dplyr::summarise(ok = any(!is.na(.data[[lyr]])), .groups = "drop") %>%
    dplyr::filter(.data$ok) %>%
    dplyr::pull(dplyr::all_of(idcol))
  n0 <- nrow(plots_sf)
  out <- plots_sf[sort(unique(as.integer(keep_ids))), ]
  nd <- n0 - nrow(out)
  if (nd > 0L) {
    message("Dropped ", nd, " of ", n0, " polygons with no overlapping CHM pixels (", lyr, ").")
  }
  out
}

#' Largest axis span (m) for a square map extent around figure-class plots that overlap the CHM.
#' Use `max()` across sites and pass as `global_square_m` so all panels share the same scale.
sbe_map_square_side_m <- function(
    plots_sf,
    chm_path,
    map_display_limits = sbe_figure_map_classes()
) {
  if (is.null(plots_sf) || nrow(plots_sf) == 0L) {
    return(NA_real_)
  }
  if (is.na(chm_path) || !nzchar(chm_path) || !file.exists(chm_path)) {
    return(NA_real_)
  }
  r <- terra::rast(chm_path)
  ps <- sf::st_transform(plots_sf, sf::st_crs(r))
  ps <- sbe_filter_plots_chm_overlap(ps, r)
  if (nrow(ps) == 0L) {
    return(NA_real_)
  }
  ps <- ps %>%
    dplyr::mutate(
      .map_disp = sbe_map_display_from_plot_class(as.character(.data$plot_class))
    ) %>%
    dplyr::filter(.data$.map_disp %in% .env$map_display_limits)
  if (nrow(ps) == 0L) {
    return(NA_real_)
  }
  bb <- sf::st_bbox(ps)
  w <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  max(w, h)
}

#' Horizontal legend for land classes (fill swatches).
sbe_class_legend_bottom_ggplot <- function(
    pal_named = sbe_palette_bar(),
    classes = sbe_figure_bar_classes()
) {
  d <- data.frame(cls = factor(classes, levels = classes), yy = 1)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$cls, y = .data$yy, fill = .data$cls)) +
    ggplot2::geom_point(shape = 22, size = 3.6, stroke = 0, alpha = 0) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = pal_named[classes],
      limits = classes,
      breaks = classes,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        nrow = 2L,
        direction = "horizontal",
        label.position = "right",
        label.theme = ggplot2::element_text(size = 7),
        keywidth = ggplot2::unit(3, "mm"),
        keyheight = ggplot2::unit(3, "mm"),
        override.aes = list(alpha = 1, colour = NA)
      )
    ) +
    ggplot2::coord_cartesian(ylim = c(0.998, 1.002), expand = FALSE, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.just = "center",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(0, 6, 0, 6)
    )
}

#' Single CHM colourbar (no fill raster beside it — avoids a “double grey bar”).
sbe_chm_legend_column_ggplot <- function() {
  y <- seq(0, 60, length.out = 200L)
  d <- data.frame(y = y, chm = y)
  ggplot2::ggplot(d, ggplot2::aes(x = 1, y = .data$y, color = .data$chm)) +
    ggplot2::geom_point(alpha = 0, stroke = 0, size = 0.01) +
    ggplot2::scale_color_gradientn(
      colours = grDevices::gray.colors(8, start = 0.97, end = 0.22),
      limits = c(0, 60),
      name = "CHM 2020 (m)",
      na.value = NA_character_,
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(4.5, "mm"),
        barheight = ggplot2::unit(5.8, "cm"),
        title.position = "top",
        frame.colour = "grey45",
        ticks.colour = "grey35"
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(1, 1), expand = FALSE, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.margin = ggplot2::margin(2, 2, 2, 4),
      plot.margin = ggplot2::margin(4, 0, 4, 0)
    )
}

#' One ggplot panel: greyscale CHM, polygon outlines by class colour,
#' dashed valid-data footprint. `plots_sf` must have `plot_class` and
#' `geometry` (and usually `plot_class_display` from loaders).
#' When `global_square_m` is set (same value for all panels), each map uses a
#' square extent of that side in metres, centred on the site, so metres per
#' cm match across Primary, SBE, and RIL.
sbe_ggplot_lidar_panel <- function(
    plots_sf,
    chm_path,
    panel_title,
    pad_frac = 0.12,
    target_cells = 1800L,
    pal_named = sbe_palette_map(),
    map_display_limits = sbe_figure_map_classes(),
    show_chm_fill_legend = FALSE,
    axis_title_x = NULL,
    axis_title_y = NULL,
    show_axis_text_x = TRUE,
    show_axis_text_y = TRUE,
    plot_title_size = 10.5,
    square_extent = FALSE,
    global_square_m = NULL
) {
  if (is.na(chm_path) || !nzchar(chm_path) || !file.exists(chm_path)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = paste0("Missing CHM:\n", panel_title), size = 3.2
        ) +
        ggplot2::theme_void()
    )
  }

  r <- terra::rast(chm_path)
  names(r) <- "chm_m"
  plots_sf <- sf::st_transform(plots_sf, sf::st_crs(r))
  plots_sf <- sbe_filter_plots_chm_overlap(plots_sf, r)
  if (nrow(plots_sf) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = paste0("No plots overlap CHM:\n", panel_title), size = 3.2
        ) +
        ggplot2::theme_void()
    )
  }

  if (!"plot_class" %in% names(plots_sf)) {
    stop("sbe_ggplot_lidar_panel: plots_sf must have column plot_class")
  }

  plots_sf <- plots_sf %>%
    dplyr::mutate(
      .map_disp = sbe_map_display_from_plot_class(as.character(.data$plot_class))
    ) %>%
    dplyr::filter(.data$.map_disp %in% .env$map_display_limits) %>%
    dplyr::mutate(
      plot_class_display = factor(.data$.map_disp, levels = map_display_limits)
    ) %>%
    dplyr::select(-".map_disp")

  if (nrow(plots_sf) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = paste0("No figure-class plots in panel:\n", panel_title), size = 3.2
        ) +
        ggplot2::theme_void()
    )
  }

  bb <- sf::st_bbox(plots_sf)
  gsq <- if (!is.null(global_square_m)) global_square_m[[1L]] else NA_real_
  if (!is.null(global_square_m) && length(gsq) == 1L && is.finite(gsq) && gsq > 0) {
    xc <- mean(as.numeric(bb[c("xmin", "xmax")]))
    yc <- mean(as.numeric(bb[c("ymin", "ymax")]))
    hlf <- as.numeric(gsq) / 2
    bb[["xmin"]] <- xc - hlf
    bb[["xmax"]] <- xc + hlf
    bb[["ymin"]] <- yc - hlf
    bb[["ymax"]] <- yc + hlf
  } else if (isTRUE(square_extent)) {
    w <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
    h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
    side <- max(w, h)
    xc <- mean(as.numeric(bb[c("xmin", "xmax")]))
    yc <- mean(as.numeric(bb[c("ymin", "ymax")]))
    bb[["xmin"]] <- xc - side / 2
    bb[["xmax"]] <- xc + side / 2
    bb[["ymin"]] <- yc - side / 2
    bb[["ymax"]] <- yc + side / 2
  }
  pad <- max(diff(bb[c("xmin", "xmax")]), diff(bb[c("ymin", "ymax")])) * pad_frac
  ext_crop <- terra::ext(
    bb[["xmin"]] - pad, bb[["xmax"]] + pad,
    bb[["ymin"]] - pad, bb[["ymax"]] + pad
  )
  r_crop <- tryCatch(terra::crop(r, ext_crop, snap = "out"), error = function(e) r)

  nc <- ncol(r_crop)
  nr <- nrow(r_crop)
  fact <- max(1L, as.integer(ceiling(max(nc, nr) / target_cells)))
  r_d <- if (fact > 1L) {
    terra::aggregate(r_crop, fact = fact, fun = "mean", na.rm = TRUE)
  } else {
    r_crop
  }

  # Cap raster size so non-stars `geom_raster` + `coord_sf` does not expand to
  # millions of `geom_rect` draws (very slow). Few halving steps are cheap.
  max_cells <- 12000L
  while (terra::ncell(r_d) > max_cells) {
    r_d <- terra::aggregate(r_d, fact = 2L, fun = "mean", na.rm = TRUE)
  }

  foot <- terra::ifel(is.na(r_d[[1L]]), NA, 1)
  foot_poly <- terra::as.polygons(foot, dissolve = TRUE) %>%
    sf::st_as_sf() %>%
    sf::st_transform(sf::st_crs(plots_sf))

  xlim_sf <- c(bb[["xmin"]] - pad, bb[["xmax"]] + pad)
  ylim_sf <- c(bb[["ymin"]] - pad, bb[["ymax"]] + pad)

  ord <- map_display_limits
  pal_vals <- unname(pal_named[ord])
  names(pal_vals) <- ord

  chm_guide <- if (isTRUE(show_chm_fill_legend)) {
    ggplot2::guide_colorbar(
      barwidth = ggplot2::unit(3, "mm"),
      barheight = ggplot2::unit(22, "mm"),
      title.position = "top"
    )
  } else {
    "none"
  }

  ax_x <- if (is.null(axis_title_x)) {
    ggplot2::element_blank()
  } else {
    ggplot2::element_text(size = 8)
  }
  ax_y <- if (is.null(axis_title_y)) {
    ggplot2::element_blank()
  } else {
    ggplot2::element_text(size = 8)
  }
  txt_x <- if (isTRUE(show_axis_text_x)) {
    ggplot2::element_text(size = 7, angle = 28, hjust = 1, vjust = 1, margin = ggplot2::margin(t = 3))
  } else {
    ggplot2::element_blank()
  }
  txt_y <- if (isTRUE(show_axis_text_y)) {
    ggplot2::element_text(size = 7)
  } else {
    ggplot2::element_blank()
  }

  common_theme <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = plot_title_size, hjust = 0),
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.2),
      axis.title.x = ax_x,
      axis.title.y = ax_y,
      axis.text.x = txt_x,
      axis.text.y = txt_y,
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  if (requireNamespace("stars", quietly = TRUE)) {
    st_r <- stars::st_as_stars(r_d)
    p <- ggplot2::ggplot() +
      stars::geom_stars(data = st_r, ggplot2::aes(fill = chm_m)) +
      ggplot2::scale_fill_gradientn(
        colours = grDevices::gray.colors(8, start = 0.97, end = 0.25),
        name = "CHM 2020\n(m)",
        na.value = NA_character_,
        guide = chm_guide
      ) +
      ggplot2::geom_sf(
        data = plots_sf,
        ggplot2::aes(colour = plot_class_display),
        fill = NA,
        linewidth = 0.45,
        lineend = "round",
        linejoin = "mitre"
      ) +
      ggplot2::geom_sf(
        data = foot_poly,
        fill = NA,
        colour = "#009E73",
        linewidth = 0.28,
        linetype = "dashed",
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = pal_vals,
        limits = ord,
        name = NULL,
        drop = FALSE,
        na.value = "grey50",
        guide = "none"
      ) +
      ggplot2::coord_sf(
        expand = FALSE,
        crs = sf::st_crs(plots_sf),
        default_crs = sf::st_crs(plots_sf),
        xlim = xlim_sf,
        ylim = ylim_sf
      ) +
      ggplot2::labs(title = panel_title, x = axis_title_x, y = axis_title_y) +
      common_theme
  } else {
    chm_df <- as.data.frame(r_d, xy = TRUE, na.rm = TRUE)
    names(chm_df) <- c("x", "y", "chm_m")
    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = chm_df,
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$chm_m),
        interpolate = FALSE
      ) +
      ggplot2::scale_fill_gradientn(
        colours = grDevices::gray.colors(8, start = 0.97, end = 0.25),
        name = "CHM 2020\n(m)",
        na.value = NA_character_,
        guide = chm_guide
      ) +
      ggplot2::geom_sf(
        data = plots_sf,
        ggplot2::aes(colour = plot_class_display),
        fill = NA,
        linewidth = 0.45
      ) +
      ggplot2::geom_sf(
        data = foot_poly,
        fill = NA,
        colour = "#009E73",
        linewidth = 0.28,
        linetype = "dashed",
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = pal_vals,
        limits = ord,
        name = NULL,
        drop = FALSE,
        guide = "none"
      ) +
      ggplot2::coord_sf(
        crs = sf::st_crs(plots_sf),
        expand = FALSE,
        default_crs = sf::st_crs(plots_sf),
        xlim = xlim_sf,
        ylim = ylim_sf
      ) +
      ggplot2::labs(title = panel_title, x = axis_title_x, y = axis_title_y) +
      common_theme
  }
  p
}
