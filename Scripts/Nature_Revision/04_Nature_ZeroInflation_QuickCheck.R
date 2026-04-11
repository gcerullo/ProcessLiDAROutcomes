# Quick iterative model checks for zero-heavy megatree data.
# This script uses a subset of the data to compare zero-friendly model families
# and produce habitat-level PPC outputs quickly.

source(file.path("Scripts", "Nature_Revision", "config.R"))

library(tidyverse)
library(brms)
library(bayesplot)
library(loo)

# -----------------------------
# Quick-run settings (modify when model testing)
# -----------------------------
set.seed(123)
height_thresholds <- c(50)        # keep narrow for fast iteration
subset_id_fraction <- 1        # sample fraction of hex IDs
min_cells_filter <- 6000

fit_chains <- 4
fit_iter <- 2000
fit_warmup <- 1000
fit_cores <- 4
ppc_draws <- 150

common_ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

# -----------------------------
# Build modelling dataset
# -----------------------------
primary_vals <- readRDS(nr2_output_path("primary1haCHMvals.rds")) %>%
  rename(CH = Primary) %>%
  mutate(habitat = "primary")

once_vals <- readRDS(nr2_output_path("onceLogged_vals1haCHMvals.rds")) %>%
  rename(CH = Band_1) %>%
  mutate(habitat = "once-logged")

restored_vals <- readRDS(nr2_output_path("restored_vals1haCHMvals.rds")) %>%
  rename(CH = Restored) %>%
  mutate(habitat = "restored")

twice_vals <- readRDS(nr2_output_path("twiceLogged_vals1haCHMvals.rds")) %>%
  rename(CH = Band_1) %>%
  mutate(habitat = "twice-logged")

all_vals <- bind_rows(primary_vals, once_vals, restored_vals, twice_vals) %>%
  group_by(ID) %>%
  mutate(num_cells = n()) %>%
  ungroup()

CH_thresh <- function(x, max_h) {
  x %>%
    group_by(ID, habitat, num_cells) %>%
    summarise(bigTrees = sum(CH > max_h), .groups = "drop") %>%
    mutate(height_filt = max_h, propBigTrees = bigTrees / num_cells)
}

df_all <- bind_rows(lapply(height_thresholds, function(h) CH_thresh(all_vals, h))) %>%
  filter(num_cells > min_cells_filter) %>%
  mutate(
    habitat = factor(habitat, levels = c("primary", "once-logged", "restored", "twice-logged")),
    prop_raw = pmin(pmax(propBigTrees, 0), 1),
    prop_adj = pmin(pmax(propBigTrees, 1e-6), 1 - 1e-6)
  )

# Keep the full quick-check dataset for this single-model run.
# If you want to speed things up again later, swap this back to an ID subset.
df <- df_all

# Quick exploratory check of the adjusted proportion distribution
prop_adj_hist <- ggplot(df, aes(x = prop_adj)) +
  geom_histogram(bins = 40, color = "white", fill = "grey40") +
  facet_grid(habitat ~ height_filt, scales = "free_y") +
  labs(
    title = "Distribution of adjusted megatree proportions",
    x = "Adjusted megatree proportion (prop_adj)",
    y = "Count"
  ) +
  theme_minimal(base_size = 11)

# ggsave(
#   filename = file.path(nr2_config$figures_dir, "quick_prop_adj_hist_by_height_habitat.png"),
#   plot = prop_adj_hist,
#   width = 10,
#   height = 8,
#   dpi = 300
# )
# 
# ggsave(
#   filename = file.path(nr2_config$figures_dir, "quick_prop_adj_hist_by_height_habitat.pdf"),
#   plot = prop_adj_hist,
#   width = 10,
#   height = 8
# )

zero_check <- df %>%
  group_by(habitat) %>%
  summarise(
    n = n(),
    zero_prop = mean(prop_raw == 0),
    one_prop = mean(prop_raw == 1),
    mean_prop = mean(prop_raw),
    sd_prop = sd(prop_raw),
    .groups = "drop"
  )

# -----------------------------
# zero inflated model
# -----------------------------

priors_zibb <- c(
  prior(normal(0, 2), class = "b"),
  prior(normal(0, 3), class = "Intercept"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "b", dpar = "phi"),
  prior(normal(0, 2), class = "b", dpar = "zi"),
  prior(logistic(0, 1), class = "Intercept", dpar = "zi")
)

# Zero-inflated beta-binomial on the count response at the 50 m threshold.
# Habitat-specific mean, dispersion, and zero inflation are all allowed to
# differ, which should better capture the very different spread and zero mass
# across habitats at a single threshold.

model_zibb <- brm(
  formula = bf(
    bigTrees | trials(num_cells) ~ habitat,
    phi ~ habitat,
    zi ~ habitat
  ),
  data = df,
  family = zero_inflated_beta_binomial(),
  prior = priors_zibb,
  chains = fit_chains,
  iter = fit_iter,
  warmup = fit_warmup,
  cores = fit_cores,
  seed = 2002,
  backend = "cmdstanr",
  control = common_ctrl
)
saveRDS(model_zibb, file.path(nr2_config$models_dir, "nr2_zibb_re.rds"))


# -----------------------------
# PPC and fit summaries for the single ZIBB model
# -----------------------------
make_ppc_pdf <- function(model_obj, model_name, dat) {
  pdf_path <- file.path(nr2_config$figures_dir, paste0("quick_ppc_by_habitat_", model_name, ".pdf"))
  yrep_counts <- posterior_predict(model_obj, newdata = dat, ndraws = ppc_draws, re_formula = NA)
  yrep_prop <- sweep(yrep_counts, 2, dat$num_cells, "/")

  pdf(pdf_path, width = 8, height = 6)
  habitats <- levels(dat$habitat)
  for (h in habitats) {
    idx <- which(dat$habitat == h)
    yrep_h <- yrep_prop[, idx, drop = FALSE]
    keep <- which(rowSums(is.na(yrep_h)) == 0)
    if (length(keep) < 5) next
    yrep_ok <- yrep_h[keep, , drop = FALSE]
    yrep_show <- yrep_ok[seq_len(min(50, nrow(yrep_ok))), , drop = FALSE]
    y_obs <- dat$prop_raw[idx]

    p1 <- bayesplot::ppc_dens_overlay(y = y_obs, yrep = yrep_show)
    print(p1 + ggplot2::ggtitle(paste(model_name, "- Density overlay -", h)))
    p2 <- bayesplot::ppc_hist(y = y_obs, yrep = yrep_show)
    print(p2 + ggplot2::ggtitle(paste(model_name, "- Histogram -", h)))
    p3 <- bayesplot::ppc_scatter_avg(y = y_obs, yrep = yrep_show)
    print(p3 + ggplot2::ggtitle(paste(model_name, "- Scatter observed vs predicted -", h)))
    p4 <- bayesplot::ppc_stat(y = y_obs, yrep = yrep_show, stat = "mean")
    print(p4 + ggplot2::ggtitle(paste(model_name, "- Predicted means -", h)))
    p5 <- bayesplot::ppc_stat(y = y_obs, yrep = yrep_show, stat = "sd")
    print(p5 + ggplot2::ggtitle(paste(model_name, "- Predicted SDs -", h)))
    p6 <- bayesplot::ppc_stat(y = y_obs, yrep = yrep_show, stat = function(x) mean(x == 0))
    print(p6 + ggplot2::ggtitle(paste(model_name, "- Predicted zero proportions -", h)))
  }
  dev.off()
}

make_fit_summary <- function(model_obj, model_name, dat) {
  yrep_counts <- posterior_predict(model_obj, newdata = dat, re_formula = NA)
  yrep_prop <- sweep(yrep_counts, 2, dat$num_cells, "/")
  idx_by_habitat <- split(seq_len(nrow(dat)), as.character(dat$habitat))

  purrr::map_dfr(names(idx_by_habitat), function(hab_name) {
    idx <- idx_by_habitat[[hab_name]]
    y_obs <- dat$prop_raw[idx]
    yrep_h <- yrep_prop[, idx, drop = FALSE]
    pred_means <- rowMeans(yrep_h)
    pred_sds <- apply(yrep_h, 1, sd)
    pred_zero_props <- apply(yrep_h, 1, function(x) mean(x == 0))

    tibble::tibble(
      model = model_name,
      habitat = hab_name,
      obs_mean = mean(y_obs),
      obs_sd = sd(y_obs),
      obs_zero_prop = mean(y_obs == 0),
      pred_mean_median = median(pred_means),
      pred_mean_lwr = quantile(pred_means, 0.025),
      pred_mean_upr = quantile(pred_means, 0.975),
      pred_sd_median = median(pred_sds),
      pred_sd_lwr = quantile(pred_sds, 0.025),
      pred_sd_upr = quantile(pred_sds, 0.975),
      pred_zero_prop_median = median(pred_zero_props),
      pred_zero_prop_lwr = quantile(pred_zero_props, 0.025),
      pred_zero_prop_upr = quantile(pred_zero_props, 0.975)
    )
  })
}

make_ppc_pdf(model_zibb, "zibb_re_50m", df)

fit_summary_zibb <- make_fit_summary(model_zibb, "zibb_re", df)

saveRDS(fit_summary_zibb, nr2_output_path("model_fit_summary_by_habitat.rds"))

message("Quick zero-inflated beta-binomial model check complete.")
message("Check Outputs/NR2/figures/quick_prop_adj_hist_by_height_habitat.png for the faceted prop_adj histogram.")
message("Check Outputs/NR2/quick_model_fit_summary_by_habitat.csv and quick_ppc_by_habitat_zibb_re.pdf")
