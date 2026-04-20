# =============================================================================
# Megatree Δp (50 m) + carbon pooled row from plot-level / class summaries
# -----------------------------------------------------------------------------
# Reads Outputs/NR2/lidar_SBE_plot_level.csv for habitats missing from
# lidar_SBE_megatree_sensitivity_summary.csv (RIL, pooled enrichment+RIL).
#
# Source from project root (uses dplyr; config optional if paths passed in).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

if (!exists("sbe_lidar_megatree_delta_years", inherits = TRUE)) {
  source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_lidar_plot_class_labels.R"))
}

#' Per-plot megatree share change at 50 m, **per year**:
#'   (p_2020 − p_2013) / `sbe_lidar_megatree_delta_years()` (same years as LiDAR epochs).
sbe_plot_level_add_delta_prop_50m <- function(plot_df) {
  yrs <- as.numeric(sbe_lidar_megatree_delta_years())
  if (!is.finite(yrs) || yrs <= 0) {
    yrs <- 7
  }
  plot_df %>%
    dplyr::mutate(
      p2013 = .data$pixels_over_50m_in_2013 / .data$total_pixels,
      p2020 = .data$pixels_over_50m_in_2020 / .data$total_pixels,
      delta_prop = (.data$p2020 - .data$p2013) / yrs
    )
}

#' One-row summary: mean delta_prop (per yr), normal approx 95% CI on the mean.
sbe_megatree_one_row <- function(df, habitat_label, height_m = 50L) {
  x <- df$delta_prop
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(
      tibble::tibble(
        habitat = habitat_label,
        height_m = as.integer(height_m),
        n_plots = 0L,
        mean_delta_prop = NA_real_,
        sd_delta = NA_real_,
        lwr_95 = NA_real_,
        upr_95 = NA_real_
      )
    )
  }
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    s <- 0
  }
  se <- s / sqrt(n)
  tibble::tibble(
    habitat = habitat_label,
    height_m = as.integer(height_m),
    n_plots = as.integer(n),
    mean_delta_prop = m,
    sd_delta = s,
    lwr_95 = m - 1.96 * se,
    upr_95 = m + 1.96 * se
  )
}

#' Extra megatree summary rows from plot-level CSV (RIL + pooled enrichment+RIL).
#'
#' @param plot_level_path Path to lidar_SBE_plot_level.csv
#' @return Tibble with columns habitat, height_m, n_plots, mean_delta_prop,
#'   sd_delta, lwr_95, upr_95
sbe_megatree_extra_rows_from_plot_level <- function(plot_level_path) {
  if (!file.exists(plot_level_path)) {
    return(tibble::tibble())
  }
  pl <- readr::read_csv(plot_level_path, show_col_types = FALSE)
  req <- c(
    "plot_class", "total_pixels", "pixels_over_50m_in_2013", "pixels_over_50m_in_2020"
  )
  if (!all(req %in% names(pl))) {
    warning("lidar_SBE_plot_level.csv missing required columns; skipping megatree extras.")
    return(tibble::tibble())
  }

  pl <- sbe_plot_level_add_delta_prop_50m(pl)

  ril <- pl %>% dplyr::filter(.data$plot_class == "once_logged_ril")
  pool <- pl %>% dplyr::filter(.data$plot_class %in% c("enrichment", "once_logged_ril"))

  dplyr::bind_rows(
    sbe_megatree_one_row(ril, "once-logged (RIL)"),
    sbe_megatree_one_row(pool, "once-logged (enrichment + RIL pooled)")
  )
}

#' Pooled carbon row (enrichment + RIL) from two class-summary rows.
#'
#' Uses a weighted mean of stratum means and an independence approximation for
#' the standard error. Returns NULL if either stratum is missing.
sbe_carbon_pooled_enrichment_ril_row <- function(cs_df) {
  a <- cs_df %>%
    dplyr::filter(.data$plot_class %in% c("enrichment", "once_logged_ril")) %>%
    dplyr::filter(!is.na(.data$mean_implied_delta_carbon_MgC_ha_yr), !is.na(.data$n_plots)) %>%
    dplyr::mutate(
      sd_MgC_ha_yr = ifelse(is.na(.data$sd_MgC_ha_yr), 0, .data$sd_MgC_ha_yr)
    )
  if (nrow(a) < 2L) {
    return(NULL)
  }
  N <- sum(a$n_plots)
  mp <- sum(a$n_plots * a$mean_implied_delta_carbon_MgC_ha_yr) / N
  se_p <- sqrt(sum((a$n_plots / N)^2 * (a$sd_MgC_ha_yr^2 / a$n_plots)))
  tibble::tibble(
    plot_class = "once_logged_enrichment_ril_pooled",
    n_plots = as.integer(N),
    mean_implied_delta_carbon_MgC_ha_yr = mp,
    sd_MgC_ha_yr = NA_real_,
    lwr_95_MgC_ha_yr = mp - 1.96 * se_p,
    upr_95_MgC_ha_yr = mp + 1.96 * se_p,
    plot_class_display = sbe_plot_class_display("once_logged_enrichment_ril_pooled")
  )
}
