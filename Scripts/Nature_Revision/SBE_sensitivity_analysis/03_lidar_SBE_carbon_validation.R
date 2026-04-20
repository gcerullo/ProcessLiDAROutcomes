# =============================================================================
# SENSITIVITY: implied forest carbon from LiDAR mean canopy height (2013 vs 2020)
# ------------------------------------------------------------------------------
# Uses: Outputs/NR2/Sabah_logging_recovery_data_LiDAR.csv (from
#   `Scripts/Nature_Revision/SBE_sensitivity_analysis/02_determine_lidar_change_2013_2020.R`).
# Assumption: 5.2 Mg C ha^-1 of stock change for each 1 m change in mean canopy
#   height, annualised by the 7-year window 
# Outputs: Outputs/NR2 (tables + figures) — see messages at end of the script.
# Related: megatree height-threshold *sensitivity* (n(>h) / N_total) is
#   `04_lidar_SBE_megatree_sensitivity.R` in this folder.
# =============================================================================

source(file.path("Scripts", "Nature_Revision", "config.R"))
source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_lidar_plot_class_labels.R"))

library(tidyverse)

# Built by `02_determine_lidar_change_2013_2020.R` in this folder → Outputs/NR2 (same basename as legacy RawData copy)
raw_path <- nr2_output_path("Sabah_logging_recovery_data_LiDAR.csv")

carbon_per_m_height_MgC_ha <- 5.2
years_elapsed <- as.numeric(2020L - 2013L)
# Keep RIL in summaries; classify as once-logged stratum for comparison to once-logged model predictions.
exclude_ril <- FALSE

required_cols <- c(
  "site", "class", "restoration",
  "canopy_height_2013_mean", "canopy_height_2020_mean"
)

d <- read_csv(raw_path, show_col_types = FALSE)
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0L) {
  stop("Input CSV missing columns: ", paste(missing_cols, collapse = ", "))
}

d <- d %>%
  mutate(
    site_w = str_trim(coalesce(as.character(site), "")),
    class_w = str_trim(coalesce(as.character(class), "")),
    restoration_w = str_trim(coalesce(as.character(restoration), "")),
    mean_heigt_diff = `canopy_height_2020_mean` - `canopy_height_2013_mean`,
    implied_delta_carbon_MgC_ha = mean_heigt_diff * carbon_per_m_height_MgC_ha,
    implied_delta_carbon_MgC_ha_yr = implied_delta_carbon_MgC_ha / years_elapsed,
    plot_class = sbe_plot_class_from_w(site_w, class_w, restoration_w),
    # Align with habitat labels used in megatree / Phillips-style once-logged predictions
    model_habitat = case_when(
      plot_class == "once_logged_ril" ~ "once_logged",
      plot_class == "sbe_once_logged_control" ~ "once_logged",
      plot_class == "enrichment" ~ "once_logged",
      plot_class == "twice_logged" ~ "twice_logged",
      plot_class == "primary" ~ "primary",
      plot_class == "liana_cutting" ~ "restored",
      TRUE ~ NA_character_
    )
  )

viol_primary <- d %>%
  filter(
    restoration_w == "none",
    str_detect(class_w, regex("primary\\s*forest", ignore_case = TRUE)),
    site_w != "Danum"
  )
if (nrow(viol_primary) > 0) {
  stop(
    "Primary forest + restoration 'none' must only occur at Danum. Offending rows:\n",
    paste(utils::capture.output(print(viol_primary %>% select(site, class, restoration, ID))), collapse = "\n")
  )
}

d <- d %>%
  select(-site_w, -class_w, -restoration_w)

if (exclude_ril) {
  n_ril_excluded <- sum(d$plot_class == "once_logged_ril", na.rm = TRUE)
  d <- d %>%
    filter(plot_class != "once_logged_ril")
  message("Excluding RIL / once_logged_ril (n = ", n_ril_excluded, " plots) from summaries, figures, and exported plot table.")
}

pool_control_enrichment <- c("sbe_once_logged_control", "enrichment")

# 95% two-sided CIs for the per-class *mean* (one-sample t, df = n-1, within-class plot SD)
add_mean_ci_95 <- function(summ) {
  summ %>%
    mutate(
      ci_95_marg = if_else(
        n_plots < 2L,
        NA_real_,
        stats::qt(0.975, df = n_plots - 1L) *
          (sd_MgC_ha_yr / sqrt(as.double(n_plots)))
      ),
      lwr_95_MgC_ha_yr = if_else(
        n_plots < 2L,
        NA_real_,
        mean_implied_delta_carbon_MgC_ha_yr - ci_95_marg
      ),
      upr_95_MgC_ha_yr = if_else(
        n_plots < 2L,
        NA_real_,
        mean_implied_delta_carbon_MgC_ha_yr + ci_95_marg
      )
    ) %>%
    select(-ci_95_marg)
}

class_summary <- d %>%
  group_by(plot_class) %>%
  summarise(
    n_plots = n(),
    mean_implied_delta_carbon_MgC_ha_yr = mean(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE),
    sd_MgC_ha_yr = if_else(
      n_plots < 2L,
      NA_real_,
      sd(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  add_mean_ci_95()

combo_summary <- d %>%
  filter(plot_class %in% pool_control_enrichment) %>%
  summarise(
    plot_class = "sbe_once_logged_control_plus_enrichment",
    n_plots = n(),
    mean_implied_delta_carbon_MgC_ha_yr = mean(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE),
    sd_MgC_ha_yr = if_else(
      n() < 2L,
      NA_real_,
      sd(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE)
    )
  ) %>%
  add_mean_ci_95()

class_summary <- bind_rows(class_summary, combo_summary) %>%
  mutate(
    plot_class_display = sbe_plot_class_display(plot_class),
    plot_class_lab = forcats::fct_reorder(plot_class_display, mean_implied_delta_carbon_MgC_ha_yr)
  )

message("Per-plot data: nrow = ", nrow(d))
print(
  d %>%
    select(site, class, restoration, plot_class, model_habitat, ID, mean_heigt_diff, implied_delta_carbon_MgC_ha_yr) %>%
    head(8L)
)

message("\nMean and 95% CI of implied Mg C ha^-1 yr^-1 (t on mean) by plot_class:")
print(class_summary %>% select(-plot_class_lab), n = 100)

plot_dt <- class_summary %>%
  mutate(
    plot_class_lab_n = forcats::fct_reorder(
      paste0(plot_class_display, " (n = ", n_plots, ", 4 ha plots)"),
      mean_implied_delta_carbon_MgC_ha_yr
    ),
    ymin = lwr_95_MgC_ha_yr,
    ymax = upr_95_MgC_ha_yr
  )

p <- ggplot(plot_dt, aes(x = plot_class_lab_n, y = mean_implied_delta_carbon_MgC_ha_yr)) +
  geom_col(width = 0.72, fill = "#8a8a8a", colour = "grey35", linewidth = 0.25, alpha = 0.95) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.22, linewidth = 0.35, colour = "grey25") +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Implied" ~ Delta * C ~ "(" * Mg ~ C ~ ha^{-1} ~ yr^{-1} * ")"),
    title = "LiDAR-derived implied carbon change",
    subtitle = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", size = 12),
    plot.subtitle = element_blank(),
    axis.text.y = element_text(size = 8.2)
  )

fig_file <- file.path(nr2_config$figures_dir, "lidar_SBE_implied_MgC_ha_yr_by_plot_class.png")
ggsave(fig_file, p, width = 9.2, height = 6.5, dpi = 300, bg = "white")

# Quick comparison: mean with 95% CI on a common Mg C ha^-1 yr^-1 scale (point + range, zero line)
p_compare <- ggplot(plot_dt, aes(x = forcats::fct_rev(plot_class_lab_n), y = mean_implied_delta_carbon_MgC_ha_yr)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
  geom_pointrange(
    aes(ymin = ymin, ymax = ymax),
    colour = "grey35",
    linewidth = 0.55,
    size = 2.2
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Implied" ~ Delta * C ~ "(" * Mg ~ C ~ ha^{-1} ~ yr^{-1} * ")"),
    title = "LiDAR-derived implied carbon change",
    subtitle = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", size = 12),
    plot.subtitle = element_blank(),
    axis.text.y = element_text(size = 8.2)
  )

fig_compare <- file.path(nr2_config$figures_dir, "lidar_SBE_MgC_ha_yr_comparison.png")
ggsave(fig_compare, p_compare, width = 9.2, height = 6.5, dpi = 300, bg = "white")

# Alternate view: drop SBE once-logged control and the control+enrichment pooled row;
# combine enrichment + RIL into one "once-logged (enrichment + RIL pooled)" class.
pool_enrich_ril <- c("enrichment", "once_logged_ril")
alt_per_class <- d %>%
  filter(!plot_class %in% c("sbe_once_logged_control", pool_enrich_ril)) %>%
  group_by(plot_class) %>%
  summarise(
    n_plots = n(),
    mean_implied_delta_carbon_MgC_ha_yr = mean(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE),
    sd_MgC_ha_yr = if_else(
      n_plots < 2L,
      NA_real_,
      sd(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  add_mean_ci_95()

alt_pooled_enrich_ril <- d %>%
  filter(plot_class %in% pool_enrich_ril) %>%
  summarise(
    plot_class = "once_logged_enrichment_ril_pooled",
    n_plots = n(),
    mean_implied_delta_carbon_MgC_ha_yr = mean(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE),
    sd_MgC_ha_yr = if_else(
      n() < 2L,
      NA_real_,
      sd(implied_delta_carbon_MgC_ha_yr, na.rm = TRUE)
    )
  ) %>%
  add_mean_ci_95()

class_summary_alt <- bind_rows(alt_per_class, alt_pooled_enrich_ril) %>%
  mutate(
    plot_class_display = sbe_plot_class_display(plot_class),
    plot_class_lab = forcats::fct_reorder(plot_class_display, mean_implied_delta_carbon_MgC_ha_yr)
  )

plot_dt_alt <- class_summary_alt %>%
  mutate(
    plot_class_lab_n = forcats::fct_reorder(
      paste0(plot_class_display, " (n = ", n_plots, ", 4 ha plots)"),
      mean_implied_delta_carbon_MgC_ha_yr
    ),
    ymin = lwr_95_MgC_ha_yr,
    ymax = upr_95_MgC_ha_yr
  )

p_alt <- ggplot(plot_dt_alt, aes(x = plot_class_lab_n, y = mean_implied_delta_carbon_MgC_ha_yr)) +
  geom_col(width = 0.72, fill = "#8a8a8a", colour = "grey35", linewidth = 0.25, alpha = 0.95) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.22, linewidth = 0.35, colour = "grey25") +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Implied" ~ Delta * C ~ "(" * Mg ~ C ~ ha^{-1} ~ yr^{-1} * ")"),
    title = "LiDAR-derived implied carbon change",
    subtitle = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", size = 12),
    plot.subtitle = element_blank(),
    axis.text.y = element_text(size = 8.2)
  )

fig_alt <- file.path(
  nr2_config$figures_dir,
  "lidar_SBE_implied_MgC_ha_yr_enrichment_RIL_pooled_no_control.png"
)
ggsave(fig_alt, p_alt, width = 9.2, height = 6.5, dpi = 300, bg = "white")

write_csv(d, nr2_output_path("lidar_SBE_plot_level.csv"))
write_csv(
  class_summary %>% select(-plot_class_lab),
  nr2_output_path("lidar_SBE_class_summary.csv")
)
write_csv(
  class_summary %>% select(
    plot_class, plot_class_display, n_plots, mean_implied_delta_carbon_MgC_ha_yr, sd_MgC_ha_yr,
    lwr_95_MgC_ha_yr, upr_95_MgC_ha_yr
  ),
  nr2_output_path("sbe_uncertainty_carbon.csv")
)
write_csv(
  class_summary_alt %>% select(-plot_class_lab),
  nr2_output_path("lidar_SBE_class_summary_enrichment_RIL_pooled_no_control.csv")
)

message("\nWrote:\n  ", normalizePath(nr2_output_path("lidar_SBE_plot_level.csv"), mustWork = FALSE))
message("  ", normalizePath(nr2_output_path("lidar_SBE_class_summary.csv"), mustWork = FALSE))
message("  ", normalizePath(nr2_output_path("sbe_uncertainty_carbon.csv"), mustWork = FALSE))
message("  ", normalizePath(nr2_output_path("lidar_SBE_class_summary_enrichment_RIL_pooled_no_control.csv"), mustWork = FALSE))
message("  ", normalizePath(fig_file, mustWork = FALSE))
message("  ", normalizePath(fig_compare, mustWork = FALSE))
message("  ", normalizePath(fig_alt, mustWork = FALSE))
