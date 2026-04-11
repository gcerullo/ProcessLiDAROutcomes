# Calculate megatrees per habitat and fit baseline Bayesian model.

source(file.path("Scripts", "Nature_Revision", "config.R"))

library(terra)
library(tidyverse)
library(data.table)
library(ggpubr)
library(brms)
library(tidybayes)

max_height <- c(45, 50, 55, 60, 65, 70)

hex <- vect(nr2_raw_path("Sabah100HexShapefile.shp"))
primary <- rast(nr2_raw_path("Primary.tif"))
once_logged <- rast(nr2_raw_path("Once_Logged.tif"))
restore <- rast(nr2_raw_path("Restored.tif"))
twice_logged <- rast(nr2_raw_path("Twice_Logged.tif"))

# Optional one-time extraction from rasters to 1ha hex cells:
# primary_vals <- terra::extract(primary, hex) %>% na.omit()
# once_vals <- terra::extract(once_logged, hex) %>% na.omit()
# restored_vals <- terra::extract(restore, hex) %>% na.omit()
# twice_vals <- terra::extract(twice_logged, hex) %>% na.omit()
# saveRDS(primary_vals, nr2_output_path("primary1haCHMvals.rds"))
# saveRDS(once_vals, nr2_output_path("onceLogged_vals1haCHMvals.rds"))
# saveRDS(restored_vals, nr2_output_path("restored_vals1haCHMvals.rds"))
# saveRDS(twice_vals, nr2_output_path("twiceLogged_vals1haCHMvals.rds"))

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

df <- bind_rows(lapply(max_height, function(h) CH_thresh(all_vals, h))) %>%
  filter(num_cells > 6000) %>%
  mutate(
    habitat = factor(habitat, levels = c("primary", "once-logged", "restored", "twice-logged")),
    prop_adj = pmin(pmax(propBigTrees, 1e-6), 1 - 1e-6)
  )

bayes_model <- brm(
  formula = prop_adj ~ habitat * height_filt,
  data = df,
  family = Beta(),
  prior = c(
    prior(normal(0, 2.5), class = "b"),
    prior(normal(0, 5), class = "Intercept")
  ),
  chains = 4,
  iter = 4000,
  cores = 4,
  seed = 123,
  backend = "cmdstanr"
)

saveRDS(bayes_model, file.path(nr2_config$models_dir, "bayesian_megatree_model_22_10.rds"))

newdata <- expand.grid(
  habitat = levels(df$habitat),
  height_filt = unique(df$height_filt)
)

epreds <- posterior_epred(bayes_model, newdata = newdata, re_formula = NA)
thin_interval <- max(1, floor(nrow(epreds) / 500))
thinned_draws <- seq(1, nrow(epreds), by = thin_interval)
epreds_thinned <- epreds[thinned_draws, , drop = FALSE]

epreds_long <- as.data.frame(epreds_thinned) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(cols = -draw, names_to = "obs_id", values_to = "estimate") %>%
  mutate(obs_id = as.integer(gsub("V", "", obs_id))) %>%
  left_join(newdata %>% mutate(obs_id = row_number()), by = "obs_id") %>%
  select(draw, habitat, height_filt, estimate)

saveRDS(epreds_long, nr2_output_path("megatrees_draws_per_hab.rds"))

# Habitat-level fit summary and PPC by habitat
ppreds <- posterior_predict(bayes_model, newdata = df, re_formula = NA)
habitat_index <- split(seq_len(nrow(df)), as.character(df$habitat))

habitat_fit_summary <- purrr::map_dfr(names(habitat_index), function(hab_name) {
  idx <- habitat_index[[hab_name]]
  observed_vals <- df$prop_adj[idx]
  pred_subset <- ppreds[, idx, drop = FALSE]

  pred_means <- rowMeans(pred_subset)
  pred_sds <- apply(pred_subset, 1, sd)

  tibble(
    habitat = hab_name,
    obs_mean = mean(observed_vals),
    obs_sd = sd(observed_vals),
    pred_mean_median = median(pred_means),
    pred_mean_lwr = quantile(pred_means, 0.025),
    pred_mean_upr = quantile(pred_means, 0.975),
    pred_sd_median = median(pred_sds),
    pred_sd_lwr = quantile(pred_sds, 0.025),
    pred_sd_upr = quantile(pred_sds, 0.975)
  )
})

saveRDS(habitat_fit_summary, nr2_output_path("model_fit_by_habitat_summary.rds"))
write.csv(habitat_fit_summary, nr2_output_path("model_fit_by_habitat_summary.csv"), row.names = FALSE)

habitats <- levels(df$habitat)
pdf(file.path(nr2_config$figures_dir, "posterior_pp_checks_by_habitat.pdf"), width = 8, height = 6)
for (h in habitats) {
  habitat_data <- df %>% filter(habitat == h)
  yrep <- posterior_predict(bayes_model, newdata = habitat_data, ndraws = 200, re_formula = NA)
  keep <- which(rowSums(is.na(yrep)) == 0)
  if (length(keep) < 5) next
  yrep_ok <- yrep[keep, , drop = FALSE]
  yrep_show <- yrep_ok[seq_len(min(50, nrow(yrep_ok))), , drop = FALSE]

  p1 <- bayesplot::ppc_dens_overlay(y = habitat_data$prop_adj, yrep = yrep_show)
  print(p1 + ggtitle(paste("Density overlay -", h)))

  p2 <- bayesplot::ppc_hist(y = habitat_data$prop_adj, yrep = yrep_show)
  print(p2 + ggtitle(paste("Histogram -", h)))

  p3 <- bayesplot::ppc_scatter_avg(y = habitat_data$prop_adj, yrep = yrep_show)
  print(p3 + ggtitle(paste("Scatter observed vs predicted -", h)))

  p4 <- bayesplot::ppc_stat(y = habitat_data$prop_adj, yrep = yrep_show, stat = "mean")
  print(p4 + ggtitle(paste("Predicted means -", h)))

  p5 <- bayesplot::ppc_stat(y = habitat_data$prop_adj, yrep = yrep_show, stat = "sd")
  print(p5 + ggtitle(paste("Predicted SDs -", h)))
}
dev.off()

message("NR2 script 01 complete.")
