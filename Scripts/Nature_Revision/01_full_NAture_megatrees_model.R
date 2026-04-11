# Fit full habitat-specific zero-inflated beta-binomial models for multiple
# megatree thresholds and export threshold-specific posterior predictive checks.

source(file.path("Scripts", "Nature_Revision", "config.R"))

library(tidyverse)
library(brms)
library(bayesplot)

# -----------------------------
# Full-run settings
# -----------------------------
set.seed(123)
height_thresholds <- c(45, 50, 55)
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

build_threshold_data <- function(max_h) {
  CH_thresh(all_vals, max_h) %>%
    filter(num_cells > min_cells_filter) %>%
    mutate(
      habitat = factor(habitat, levels = c("primary", "once-logged", "restored", "twice-logged")),
      prop_raw = pmin(pmax(propBigTrees, 0), 1),
      prop_adj = pmin(pmax(propBigTrees, 1e-6), 1 - 1e-6)
    )
}

priors_zibb <- c(
  prior(normal(0, 2), class = "b"),
  prior(normal(0, 3), class = "Intercept"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "b", dpar = "phi"),
  prior(normal(0, 2), class = "b", dpar = "zi"),
  prior(logistic(0, 1), class = "Intercept", dpar = "zi")
)

# -----------------------------
# PPC and fit summary helpers
# -----------------------------
make_ppc_pdf <- function(model_obj, model_name, dat) {
  pdf_path <- file.path(nr2_config$figures_dir, paste0("ppc_by_habitat_", model_name, ".pdf"))
  yrep_counts <- posterior_predict(model_obj, newdata = dat, ndraws = ppc_draws)
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
  yrep_counts <- posterior_predict(model_obj, newdata = dat)
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
      height_filt = unique(dat$height_filt),
      habitat = hab_name,
      obs_mean = mean(y_obs),
      obs_sd = sd(y_obs),
      obs_zero_prop = mean(y_obs == 0),
      pred_mean_median = median(pred_means),
      pred_mean_lwr95 = quantile(pred_means, 0.025),
      pred_mean_upr95 = quantile(pred_means, 0.975),
      pred_mean_lwr80 = quantile(pred_means, 0.10),
      pred_mean_upr80 = quantile(pred_means, 0.90),
      pred_sd_median = median(pred_sds),
      pred_sd_lwr95 = quantile(pred_sds, 0.025),
      pred_sd_upr95 = quantile(pred_sds, 0.975),
      pred_sd_lwr80 = quantile(pred_sds, 0.10),
      pred_sd_upr80 = quantile(pred_sds, 0.90),
      pred_zero_prop_median = median(pred_zero_props),
      pred_zero_prop_lwr95 = quantile(pred_zero_props, 0.025),
      pred_zero_prop_upr95 = quantile(pred_zero_props, 0.975),
      pred_zero_prop_lwr80 = quantile(pred_zero_props, 0.10),
      pred_zero_prop_upr80 = quantile(pred_zero_props, 0.90)
    )
  })
}

make_prediction_outputs <- function(model_obj, model_name, dat, max_h) {
  representative_trials <- as.integer(round(median(dat$num_cells)))
  newdata <- tibble::tibble(
    habitat = factor(
      c("primary", "once-logged", "restored", "twice-logged"),
      levels = levels(dat$habitat)
    ),
    num_cells = representative_trials,
    height_filt = max_h
  )

  epreds_counts <- posterior_epred(model_obj, newdata = newdata)
  epreds_prop <- epreds_counts / representative_trials

  pred_summary <- as.data.frame(epreds_prop) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(cols = -draw, names_to = "obs_id", values_to = "estimate") %>%
    mutate(obs_id = as.integer(gsub("V", "", obs_id))) %>%
    group_by(obs_id) %>%
    summarise(
      mean = mean(estimate),
      lower95 = quantile(estimate, 0.025),
      upper95 = quantile(estimate, 0.975),
      lower80 = quantile(estimate, 0.10),
      upper80 = quantile(estimate, 0.90),
      .groups = "drop"
    ) %>%
    left_join(newdata %>% mutate(obs_id = row_number()), by = "obs_id") %>%
    select(habitat, height_filt, mean, lower95, upper95, lower80, upper80)

  thin_interval <- max(1, floor(nrow(epreds_prop) / 500))
  thinned_draws <- seq(1, nrow(epreds_prop), by = thin_interval)
  epreds_thinned <- epreds_prop[thinned_draws, , drop = FALSE]

  epreds_long <- as.data.frame(epreds_thinned) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(cols = -draw, names_to = "obs_id", values_to = "estimate") %>%
    mutate(obs_id = as.integer(gsub("V", "", obs_id))) %>%
    left_join(newdata %>% mutate(obs_id = row_number()), by = "obs_id") %>%
    select(draw, habitat, height_filt, estimate)

  summary_stem <- paste0("nature_megatrees_prediction_summary_", max_h, "m")
  draws_stem <- paste0("nature_megatrees_habitat_draws_", max_h, "m")

  write.csv(
    pred_summary,
    nr2_output_path(paste0(summary_stem, ".csv")),
    row.names = FALSE
  )
  saveRDS(
    pred_summary,
    nr2_output_path(paste0(summary_stem, ".rds"))
  )
  saveRDS(
    epreds_long,
    nr2_output_path(paste0(draws_stem, ".rds"))
  )

  list(
    pred_summary = pred_summary,
    epreds_long = epreds_long
  )
}

# -----------------------------
# Fit one model per threshold
# -----------------------------
fit_one_threshold <- function(max_h) {
  message("Fitting habitat-specific ZIBB model for threshold ", max_h, " m")
  dat <- build_threshold_data(max_h)
  model_name <- paste0("nature_megatrees_zibb_", max_h, "m")

  model_obj <- brm(
    formula = bf(
      bigTrees | trials(num_cells) ~ habitat,
      phi ~ habitat,
      zi ~ habitat
    ),
    data = dat,
    family = zero_inflated_beta_binomial(),
    prior = priors_zibb,
    chains = fit_chains,
    iter = fit_iter,
    warmup = fit_warmup,
    cores = fit_cores,
    seed = 2000 + max_h,
    backend = "cmdstanr",
    control = common_ctrl
  )

  saveRDS(
    model_obj,
    file.path(nr2_config$models_dir, paste0(model_name, ".rds"))
  )

  make_ppc_pdf(model_obj, model_name, dat)
  fit_summary <- make_fit_summary(model_obj, model_name, dat)
  prediction_outputs <- make_prediction_outputs(model_obj, model_name, dat, max_h)

  list(
    fit_summary = fit_summary,
    pred_summary = prediction_outputs$pred_summary,
    epreds_long = prediction_outputs$epreds_long
  )
}

threshold_results <- purrr::map(height_thresholds, fit_one_threshold)
fit_summaries_all <- purrr::map_dfr(threshold_results, "fit_summary")
pred_summaries_all <- purrr::map_dfr(threshold_results, "pred_summary")
epreds_long_all <- purrr::map_dfr(threshold_results, "epreds_long")

write.csv(
  fit_summaries_all,
  nr2_output_path("full_nature_megatrees_model_fit_summary_by_habitat.csv"),
  row.names = FALSE
)
saveRDS(
  fit_summaries_all,
  nr2_output_path("full_nature_megatrees_model_fit_summary_by_habitat.rds")
)

write.csv(
  pred_summaries_all,
  nr2_output_path("full_nature_megatrees_prediction_summary_by_habitat.csv"),
  row.names = FALSE
)
saveRDS(
  pred_summaries_all,
  nr2_output_path("full_nature_megatrees_prediction_summary_by_habitat.rds")
)
saveRDS(
  epreds_long_all,
  nr2_output_path("full_nature_megatrees_habitat_draws_all_thresholds.rds")
)

prediction_plot <- pred_summaries_all %>%
  mutate(height_filt = factor(height_filt, levels = height_thresholds)) %>%
  ggplot(aes(x = habitat, y = mean)) +
  geom_col(fill = "grey60", width = 0.7) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2, linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower80, ymax = upper80), width = 0.35, linewidth = 0.9) +
  facet_wrap(~ height_filt, scales = "free_y") +
  labs(
    y = "Proportion of canopy > height threshold",
    x = "Habitat type"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  filename = file.path(nr2_config$figures_dir, "nature_megatrees_predictions_by_threshold.pdf"),
  plot = prediction_plot,
  width = 11,
  height = 6,
  units = "in"
)

message("Full Nature megatree models complete.")
message("Check Outputs/NR2/models for threshold-specific model files.")
message("Check Outputs/NR2/figures for threshold-specific PPC PDFs.")
message("Check Outputs/NR2/figures/nature_megatrees_predictions_by_threshold.pdf for habitat predictions across thresholds.")
