# LiDAR height change (2013 vs 2020) and implied carbon rate using 5.2 Mg C ha^-1 per m canopy height.
# SBE / Sabah logging–recovery plot table.

source(file.path("Scripts", "Nature_Revision", "config.R"))

library(tidyverse)

raw_path <- file.path("RawData", "Sabah_logging_recovery_data_LiDAR.csv")

carbon_per_m_height_MgC_ha <- 5.2
years_elapsed <- as.numeric(2020L - 2013L)
exclude_ril <- TRUE

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
    plot_class = case_when(
      site_w == "Reduced impact logging" & class_w == "" & restoration_w == "none" ~ "RIL",
      site_w == "Danum" &
        str_detect(class_w, regex("^primary\\s*forest$", ignore_case = TRUE)) &
        restoration_w == "none" ~ "primary",
      site_w == "Sabah Biodiversity Experiment" &
        str_detect(class_w, regex("twice\\s*logged", ignore_case = TRUE)) &
        restoration_w == "no" ~ "twice_logged",
      site_w == "Sabah Biodiversity Experiment" &
        str_detect(class_w, regex("once\\s*logged", ignore_case = TRUE)) &
        restoration_w == "no" ~ "sbe_once_logged_control",
      site_w == "Sabah Biodiversity Experiment" &
        str_detect(class_w, regex("once\\s*logged", ignore_case = TRUE)) &
        str_detect(restoration_w, "^Enrichment") ~ "enrichment",
      site_w == "Sabah Biodiversity Experiment" &
        str_detect(class_w, regex("once\\s*logged", ignore_case = TRUE)) &
        str_detect(restoration_w, "^Liana") ~ "liana_cutting",
      TRUE ~ "other"
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
  n_ril_excluded <- sum(d$plot_class == "RIL", na.rm = TRUE)
  d <- d %>%
    filter(plot_class != "RIL")
  message("Excluding RIL (n = ", n_ril_excluded, " plots) from summaries, figures, and exported plot table.")
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
    plot_class_lab = forcats::fct_reorder(plot_class, mean_implied_delta_carbon_MgC_ha_yr)
  )

message("Per-plot data: nrow = ", nrow(d))
print(
  d %>%
    select(site, class, restoration, plot_class, ID, mean_heigt_diff, implied_delta_carbon_MgC_ha_yr) %>%
    head(8L)
)

message("\nMean and 95% CI of implied Mg C ha^-1 yr^-1 (t on mean) by plot_class:")
print(class_summary %>% select(-plot_class_lab), n = 100)

plot_dt <- class_summary %>%
  mutate(
    ymin = lwr_95_MgC_ha_yr,
    ymax = upr_95_MgC_ha_yr
  )

p <- ggplot(plot_dt, aes(x = plot_class_lab, y = mean_implied_delta_carbon_MgC_ha_yr)) +
  geom_col(width = 0.72, fill = "#2E6F40", alpha = 0.85) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.22, linewidth = 0.35) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Implied" ~ Delta * C ~ "(" * Mg ~ C ~ ha^{-1} ~ yr^{-1} * ")"),
    title = "LiDAR-derived implied carbon change rate by plot class",
    subtitle = if (exclude_ril) {
      "RIL excluded; pooled control+enrichment; error bars: 95% CI of class mean; 5.2 Mg C ha\u207b\u00b9 per m; 2013–2020 / 7 yr"
    } else {
      "Pooled control+enrichment; error bars: 95% CI of class mean; 5.2 Mg C ha\u207b\u00b9 per m; 2013–2020 / 7 yr"
    }
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

fig_file <- file.path(nr2_config$figures_dir, "lidar_SBE_implied_MgC_ha_yr_by_plot_class.png")
ggsave(fig_file, p, width = 8, height = 6.2, dpi = 300, bg = "white")

# Quick comparison: mean with 95% CI on a common Mg C ha^-1 yr^-1 scale (point + range, zero line)
p_compare <- ggplot(plot_dt, aes(x = forcats::fct_rev(plot_class_lab), y = mean_implied_delta_carbon_MgC_ha_yr)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
  geom_pointrange(
    aes(ymin = ymin, ymax = ymax),
    colour = "#134d3b",
    linewidth = 0.55,
    size = 2.2
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Implied" ~ Delta * C ~ "(" * Mg ~ C ~ ha^{-1} ~ yr^{-1} * ")"),
    title = "Compare implied carbon rates across plot classes",
    subtitle = if (exclude_ril) {
      "RIL excluded; 95% CI of class mean; pooled control+enrichment in addition to separate classes"
    } else {
      "95% CI of class mean; pooled control+enrichment in addition to separate classes"
    }
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9.5)
  )

fig_compare <- file.path(nr2_config$figures_dir, "lidar_SBE_MgC_ha_yr_comparison.png")
ggsave(fig_compare, p_compare, width = 8, height = 6.2, dpi = 300, bg = "white")

write_csv(d, nr2_output_path("lidar_SBE_plot_level.csv"))
write_csv(
  class_summary %>% select(-plot_class_lab),
  nr2_output_path("lidar_SBE_class_summary.csv")
)
write_csv(
  class_summary %>% select(
    plot_class, n_plots, mean_implied_delta_carbon_MgC_ha_yr, sd_MgC_ha_yr,
    lwr_95_MgC_ha_yr, upr_95_MgC_ha_yr
  ),
  nr2_output_path("sbe_uncertainty_carbon.csv")
)

message("\nWrote:\n  ", normalizePath(nr2_output_path("lidar_SBE_plot_level.csv"), mustWork = FALSE))
message("  ", normalizePath(nr2_output_path("lidar_SBE_class_summary.csv"), mustWork = FALSE))
message("  ", normalizePath(nr2_output_path("sbe_uncertainty_carbon.csv"), mustWork = FALSE))
message("  ", normalizePath(fig_file, mustWork = FALSE))
message("  ", normalizePath(fig_compare, mustWork = FALSE))
