# Fit and compare improved NR2 Bayesian models, then run habitat PPCs.

source(file.path("Scripts", "Nature_Revision", "config.R"))

library(tidyverse)
library(brms)
library(bayesplot)
library(loo)
library(ggpubr)

max_height <- c(45, 50, 55, 60, 65, 70)

x <- df %>% filter(height_filt == 60)
hist(x$propBigTrees)
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
   #prevent exact 0s and 1s
     prop_adj = pmin(pmax(propBigTrees, 1e-6), 1 - 1e-6),
    height_sc = as.numeric(scale(height_filt))
  )

priors_beta <- c(
  prior(normal(0, 2), class = "b"),
  prior(normal(0, 3), class = "Intercept")
)

priors_count <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 2), class = "Intercept")
)

common_ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

# Toggle behaviour:
# TRUE  -> fit models from scratch and save to Outputs/NR2/models
# FALSE -> skip fitting and rely on previously saved model files
fit_models <- FALSE

# Candidate 1: simple baseline Beta model on proportions.
# We keep this as the reference model to test whether added complexity is justified.
# Candidate 2: Beta model with random intercepts by hexagon (ID) and habitat/height-specific dispersion.
# This was selected because repeated measures within ID violate independence and PPC suggested
# variance mismatch by habitat that phi can absorb.
# Candidate 3: beta-binomial model on counts (bigTrees out of num_cells) with random ID effects.
# This is included because the response is naturally count-with-denominator data and contains
# many near-zero outcomes, where count likelihoods often calibrate tails better than Beta on proportions.
if (isTRUE(fit_models)) {
  model_beta_base <- brm(
    formula = prop_adj ~ habitat * height_sc,
    data = df,
    family = Beta(),
    prior = priors_beta,
    chains = 4,
    iter = 3000,
    cores = 4,
    seed = 123,
    backend = "cmdstanr",
    control = common_ctrl
  )
  saveRDS(model_beta_base, file.path(nr2_config$models_dir, "nr2_model_beta_base.rds"))

  model_beta_phi_re <- brm(
    formula = bf(
      prop_adj ~ habitat * height_sc + (1 | ID),
      phi ~ habitat + height_sc
    ),
    data = df,
    family = Beta(),
    prior = priors_beta,
    chains = 4,
    iter = 3000,
    cores = 4,
    seed = 124,
    backend = "cmdstanr",
    control = common_ctrl
  )
  saveRDS(model_beta_phi_re, file.path(nr2_config$models_dir, "nr2_model_beta_phi_re.rds"))

  model_beta_binomial_re <- brm(
    formula = bigTrees | trials(num_cells) ~ habitat * height_sc + (1 | ID),
    data = df,
    family = beta_binomial(),
    prior = priors_count,
    chains = 4,
    iter = 3000,
    cores = 4,
    seed = 125,
    backend = "cmdstanr",
    control = common_ctrl
  )
  saveRDS(model_beta_binomial_re, file.path(nr2_config$models_dir, "nr2_model_beta_binomial_re.rds"))
}

# Lower section: read in model outputs when needed (or after fitting, to standardise workflow).
model_paths <- list(
  beta_base = file.path(nr2_config$models_dir, "nr2_model_beta_base.rds"),
  beta_phi_re = file.path(nr2_config$models_dir, "nr2_model_beta_phi_re.rds")#,
 # beta_binomial_re = file.path(nr2_config$models_dir, "nr2_model_beta_binomial_re.rds")
)

missing_models <- names(model_paths)[!file.exists(unlist(model_paths))]
if (length(missing_models) > 0) {
  stop("Missing model files in Outputs/NR2/models: ", paste(missing_models, collapse = ", "))
}

models <- purrr::imap(model_paths, function(path, model_name) {
  message("Loading model: ", model_name)
  readRDS(path)
})

# Use LOO to rank candidates by expected out-of-sample predictive performance.
loos <- lapply(models, loo)
loo_comp <- loo_compare(loos)
loo_comp_df <- as.data.frame(loo_comp) %>%
  rownames_to_column("model") %>%
  arrange(elpd_diff)

write.csv(loo_comp_df, nr2_output_path("model_comparison_loo.csv"), row.names = FALSE)
saveRDS(loo_comp_df, nr2_output_path("model_comparison_loo.rds"))

# Save the top-ranked model so downstream diagnostics always use the same selection rule.
best_model_name <- rownames(loo_comp)[1]
best_model <- models[[best_model_name]]
saveRDS(best_model, file.path(nr2_config$models_dir, "nr2_best_model.rds"))

to_proportion_yrep <- function(yrep_raw, dat, family_name) {
  if (grepl("binomial", family_name, ignore.case = TRUE)) {
    denom <- matrix(dat$num_cells, nrow = nrow(yrep_raw), ncol = ncol(yrep_raw), byrow = TRUE)
    return(yrep_raw / denom)
  }
  yrep_raw
}

habitat_index <- split(seq_len(nrow(df)), as.character(df$habitat))
make_model_fit_summary <- function(model_obj, model_name, dat, idx_by_habitat) {
  family_name <- family(model_obj)$family
  ppreds_raw <- posterior_predict(model_obj, newdata = dat, re_formula = NA)
  ppreds_prop <- to_proportion_yrep(ppreds_raw, dat, family_name)

  purrr::map_dfr(names(idx_by_habitat), function(hab_name) {
    idx <- idx_by_habitat[[hab_name]]
    observed_vals <- dat$propBigTrees[idx]
    pred_subset <- ppreds_prop[, idx, drop = FALSE]
    pred_means <- rowMeans(pred_subset)
    pred_sds <- apply(pred_subset, 1, sd)

    tibble::tibble(
      model = model_name,
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
}

make_ppc_pdf <- function(model_obj, model_name, dat, idx_by_habitat) {
  family_name <- family(model_obj)$family
  ppreds_raw <- posterior_predict(model_obj, newdata = dat, ndraws = 200, re_formula = NA)
  ppreds_prop <- to_proportion_yrep(ppreds_raw, dat, family_name)

  pdf_path <- file.path(
    nr2_config$figures_dir,
    paste0("improved_posterior_pp_checks_by_habitat_", model_name, ".pdf")
  )
  pdf(pdf_path, width = 8, height = 6)

  for (h in names(idx_by_habitat)) {
    idx <- idx_by_habitat[[h]]
    yrep <- ppreds_prop[, idx, drop = FALSE]
    keep <- which(rowSums(is.na(yrep)) == 0)
    if (length(keep) < 5) next

    yrep_ok <- yrep[keep, , drop = FALSE]
    yrep_show <- yrep_ok[seq_len(min(50, nrow(yrep_ok))), , drop = FALSE]
    y_obs <- dat$propBigTrees[idx]

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
  }

  dev.off()
  message("Saved PPC PDF: ", pdf_path)
}

all_model_fit_summaries <- purrr::imap_dfr(models, function(model_obj, model_name) {
  make_model_fit_summary(model_obj, model_name, df, habitat_index)
})

write.csv(
  all_model_fit_summaries,
  nr2_output_path("improved_model_fit_by_habitat_summary_all_models.csv"),
  row.names = FALSE
)
saveRDS(
  all_model_fit_summaries,
  nr2_output_path("improved_model_fit_by_habitat_summary_all_models.rds")
)

purrr::iwalk(models, function(model_obj, model_name) {
  make_ppc_pdf(model_obj, model_name, df, habitat_index)
})

message("NR2 script 03 complete. Best model by LOO: ", best_model_name)
