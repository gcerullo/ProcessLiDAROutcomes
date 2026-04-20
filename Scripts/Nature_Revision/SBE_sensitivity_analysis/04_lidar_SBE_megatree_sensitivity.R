# =============================================================================
# SENSITIVITY: megatree (tall canopy) area share from LiDAR (2013 vs 2020)
# ------------------------------------------------------------------------------
# Proportion of megatree pixels at height h:  p = n(CH > h) / N_total, with
#   N_total = `total_pixels` (fixed plot/hex sample size in the table). The same
#   N is used for 2013 and 2020 so the year-on-year *change* in p is not driven
#   by a moving canopy mask in the denominator.
# For each of 45, 50, 55 m: delta = (p(2020) - p(2013)) / years  (change in share per year).
#
# plot_class rules match `03_lidar_SBE_carbon_validation.R`. For *plotting* we
#   label SBE liana as "restored (liana cut)"; once-logged enrichment and
#   control stay separate as before.
# exclude_ril: same behaviour as the carbon validation script.
# =============================================================================

source(file.path("Scripts", "Nature_Revision", "config.R"))
source(file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis", "sbe_lidar_plot_class_labels.R"))

library(tidyverse)

# Built by `02_determine_lidar_change_2013_2020.R` in this folder → Outputs/NR2
raw_path <- nr2_output_path("Sabah_logging_recovery_data_LiDAR.csv")
height_thresholds_m <- c(45L, 50L, 55L)
exclude_ril <- FALSE
megatree_delta_years <- as.numeric(sbe_lidar_megatree_delta_years())

required_base <- c("site", "class", "restoration", "ID", "total_pixels")
hcols <- c(
  "pixels_over_45m_in_2013", "pixels_over_50m_in_2013", "pixels_over_55m_in_2013",
  "pixels_over_45m_in_2020", "pixels_over_50m_in_2020", "pixels_over_55m_in_2020"
)
required_all <- c(required_base, hcols)

d0 <- read_csv(raw_path, show_col_types = FALSE)
if (length(setdiff(required_all, names(d0))) > 0L) {
  miss <- setdiff(required_all, names(d0))
  stop("Input CSV missing columns: ", paste(miss, collapse = ", "))
}

d <- d0 %>%
  mutate(
    site_w = str_trim(coalesce(as.character(site), "")),
    class_w = str_trim(coalesce(as.character(class), "")),
    restoration_w = str_trim(coalesce(as.character(restoration), "")),
    plot_class = sbe_plot_class_from_w(site_w, class_w, restoration_w)
  ) %>%
  mutate(
    habitat = factor(
      sbe_megatree_habitat_from_plot_class(plot_class),
      levels = sbe_megatree_habitat_levels()
    ),
    model_habitat = case_when(
      plot_class == "once_logged_ril" ~ "once_logged",
      plot_class == "sbe_once_logged_control" ~ "once_logged",
      plot_class == "enrichment" ~ "once_logged",
      plot_class == "twice_logged" ~ "twice_logged",
      plot_class == "primary" ~ "primary",
      plot_class == "liana_cutting" ~ "restored",
      TRUE ~ NA_character_
    ),
    habitat_display = factor(
      sbe_habitat_to_display(as.character(habitat)),
      levels = sbe_display_level_order()
    )
  )

viol <- d %>%
  filter(
    restoration_w == "none",
    str_detect(class_w, regex("primary\\s*forest", ignore_case = TRUE)),
    site_w != "Danum"
  )
if (nrow(viol) > 0) {
  stop("Primary + none outside Danum (see `03_lidar_SBE_carbon_validation.R`).")
}
d <- d %>% select(-site_w, -class_w, -restoration_w)

if (exclude_ril) {
  d <- d %>% filter(plot_class != "once_logged_ril")
}

# Proportion: megatree pixels as a share of the full sample (total valid pixels in plot)
to_prop <- function(nh, ntot) {
  ntot <- as.numeric(ntot)
  ok <- ntot > 0 & is.finite(ntot)
  p <- as.numeric(nh) / ntot
  p[!ok] <- NA_real_
  p
}

# Long: one row per (plot) × (height in m in {45, 50, 55})
d_long <- map_dfr(
  height_thresholds_m,
  function(h) {
    col13 <- paste0("pixels_over_", h, "m_in_2013")
    col20 <- paste0("pixels_over_", h, "m_in_2020")
    d %>%
      mutate(
        height_m = h,
        prop_megatree_2013 = to_prop(.data[[col13]], .data$total_pixels),
        prop_megatree_2020 = to_prop(.data[[col20]], .data$total_pixels),
        delta_prop_megatree = (prop_megatree_2020 - prop_megatree_2013) / megatree_delta_years
      ) %>%
      select(
        site, class, restoration, ID, total_pixels, plot_class, habitat, habitat_display,
        model_habitat, height_m,
        prop_megatree_2013, prop_megatree_2020, delta_prop_megatree
      )
  }
)

# 95% one-sample t-CI of mean(Δp per year) by habitat and height
add_mean_ci_95_delta <- function(summ) {
  summ %>%
    mutate(
      ci_95 = if_else(
        n_plots < 2L, NA_real_,
        stats::qt(0.975, df = n_plots - 1L) *
          (sd_delta / sqrt(as.double(n_plots)))
      ),
      lwr_95 = if_else(
        n_plots < 2L, NA_real_, mean_delta_prop - ci_95
      ),
      upr_95 = if_else(
        n_plots < 2L, NA_real_, mean_delta_prop + ci_95
      )
    ) %>%
    select(-ci_95)
}

megatree_sensitivity_summary <- d_long %>%
  group_by(habitat, habitat_display, height_m) %>%
  summarise(
    n_plots = n(),
    mean_delta_prop = mean(delta_prop_megatree, na.rm = TRUE),
    sd_delta = if_else(
      n_plots < 2L, NA_real_, sd(delta_prop_megatree, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  add_mean_ci_95_delta() %>%
  filter(n_plots > 0L) %>%
  mutate(
    model_habitat = case_when(
      habitat == "once-logged (RIL)" ~ "once_logged",
      habitat == "once-logged (control)" ~ "once_logged",
      habitat == "once-logged (enrichment)" ~ "once_logged",
      habitat == "twice logged (SBE)" ~ "twice_logged",
      habitat == "primary (Danum)" ~ "primary",
      habitat == "restored (liana cut)" ~ "restored",
      TRUE ~ NA_character_
    )
  )

h_lab <- tibble::tibble(
  height_m = c(45L, 50L, 55L),
  pnl = c(
    "45 m (p = n(>h) / N_total)", "50 m (p = n(>h) / N_total)", "55 m (p = n(>h) / N_total)"
  )
)

plot_sum <- megatree_sensitivity_summary %>%
  left_join(h_lab, by = "height_m") %>%
  mutate(
    pnl = factor(
      pnl, levels = h_lab$pnl
    ),
    habitat = forcats::fct_drop(habitat),
    habitat_display = forcats::fct_drop(habitat_display)
  )

pal_mega <- sbe_palette_by_display()[levels(plot_sum$habitat_display)]

p <- ggplot(plot_sum, aes(x = habitat_display, y = mean_delta_prop, fill = habitat_display)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_col(color = "grey20", width = 0.78, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = lwr_95, ymax = upr_95), width = 0.2, colour = "grey20", linewidth = 0.3
  ) +
  coord_flip() +
  facet_wrap(vars(pnl), ncol = 1L, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8.2),
    legend.position = "none"
  ) +
  scale_x_discrete(limits = rev) +
  scale_fill_manual(values = pal_mega, na.value = "grey70") +
  labs(
    x = NULL,
    y = "Change in p per year ((p2020 \u2212 p2013) / 7 yr); p = n(>h) / N_total",
    title = "Sensitivity: annualized change in LiDAR megatree share (N = total pixels), 2013\u20132020",
    subtitle = if (exclude_ril) {
      "RIL excluded. Liana \u2192 restored. p = n(>h)/N_total (same N both years from table). Mean(\u0394p per yr) per habitat; 95% t-CI. h \u2208 {45, 50, 55} m."
    } else {
      "RIL included as once-logged (RIL); model_habitat column for overlay vs once-logged models. Liana \u2192 restored. p = n(>h)/N_total. Mean(\u0394p/7 yr); 95% t-CI. Facets: 45, 50, 55 m."
    }
  )

fig_mega <- file.path(
  nr2_config$figures_dir, "lidar_SBE_megatree_sensitivity_deltap_by_habitat.png"
)
ggsave(fig_mega, p, width = 9, height = 8.2, dpi = 300, bg = "white")

# Shared y-scale: compare 45, 50, 55 m directly
ps_m <- dplyr::filter(plot_sum, !is.na(.data$mean_delta_prop))
pal_mega2 <- sbe_palette_by_display()[as.character(unique(ps_m$habitat_display))]

p2 <- ggplot(
  ps_m,
  aes(x = habitat_display, y = mean_delta_prop, fill = habitat_display)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.2) +
  geom_col(color = "grey20", width = 0.78, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = lwr_95, ymax = upr_95), width = 0.2, colour = "grey20", linewidth = 0.25
  ) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  facet_wrap(~height_m, labeller = label_both) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(face = "bold")) +
  scale_fill_manual(values = pal_mega2, na.value = "grey70") +
  labs(
    y = "Change in p per year ((p2020 \u2212 p2013) / 7 yr; p = n(>h) / N_total)",
    x = NULL,
    title = "Megatree sensitivity: one y-scale, one panel per height h (m)"
  ) +
  labs(subtitle = if (exclude_ril) "RIL excluded" else "RIL included as once-logged (RIL)")
fig2 <- file.path(
  nr2_config$figures_dir, "lidar_SBE_megatree_sensitivity_deltap_shared_scale.png"
)
ggsave(fig2, p2, width = 8.2, height = 5, dpi = 300, bg = "white")

write_csv(d_long, nr2_output_path("lidar_SBE_megatree_plot_level_long.csv"))
write_csv(
  megatree_sensitivity_summary, nr2_output_path("lidar_SBE_megatree_sensitivity_summary.csv")
)

message("Megatree sensitivity: long table rows = ", nrow(d_long))
print(megatree_sensitivity_summary, n = 30)

message("\nWrote:\n  ", normalizePath(nr2_output_path("lidar_SBE_megatree_plot_level_long.csv"), mustWork = FALSE))
message("  ", normalizePath(nr2_output_path("lidar_SBE_megatree_sensitivity_summary.csv"), mustWork = FALSE))
message("  ", normalizePath(fig_mega, mustWork = FALSE))
message("  ", normalizePath(fig2, mustWork = FALSE))
