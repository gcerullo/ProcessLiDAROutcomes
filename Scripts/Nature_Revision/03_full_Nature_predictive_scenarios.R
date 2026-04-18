source(file.path("Scripts", "Nature_Revision", "config.R"))

library(tidyverse)
library(brms)
library(data.table)
library(purrr)

source(nr2_input_path("FixedScenarioParmams.R"))

set.seed(123)
height_thresholds <- c(45, 50, 55)

# Middle-ground predictive setting:
# instead of propagating one fully predictive hectare per habitat, average over a
# user-defined number of representative hectare-like units within each posterior draw.
# This keeps posterior predictive variation, but targets the habitat-average signal
# that is more appropriate for scenario-scale aggregation.
effective_hectares_per_habitat <- 25L

log_progress <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message("[", timestamp, "] ", paste0(..., collapse = ""))
}

existing_model_specs <- tibble::tribble(
  ~height_filt, ~model_label, ~model_file,
  45, "nature_megatrees_zibb_45m", file.path(nr2_config$models_dir, "nature_megatrees_zibb_45m.rds"),
  50, "nature_megatrees_zibb_50m", file.path(nr2_config$models_dir, "nature_megatrees_zibb_50m.rds"),
  55, "nature_megatrees_zibb_55m", file.path(nr2_config$models_dir, "nature_megatrees_zibb_55m.rds")
)

prepare_saved_model_data <- function(dat, height_filt) {
  dat <- tibble::as_tibble(dat)

  if (!"propBigTrees" %in% names(dat) && all(c("bigTrees", "num_cells") %in% names(dat))) {
    dat <- dat %>% mutate(propBigTrees = bigTrees / num_cells)
  }

  if (!"height_filt" %in% names(dat)) {
    dat <- dat %>% mutate(height_filt = height_filt)
  }

  dat %>%
    mutate(
      habitat = factor(habitat, levels = c("primary", "once-logged", "restored", "twice-logged"))
    )
}

build_prediction_newdata <- function(dat, height_filt, effective_hectares) {
  representative_trials <- as.integer(round(median(dat$num_cells)))

  habitat_rows <- tibble::tibble(
    habitat = factor(
      c("primary", "once-logged", "restored", "twice-logged"),
      levels = levels(dat$habitat)
    ),
    num_cells = representative_trials,
    height_filt = height_filt
  )

  # Replicate each habitat row so posterior_predict() simulates several
  # representative hectare-like observations per habitat and posterior draw.
  newdata <- tidyr::crossing(
    habitat_draw_rep = seq_len(effective_hectares),
    habitat_rows
  )

  if ("height_sc" %in% names(dat)) {
    newdata <- newdata %>% mutate(height_sc = unique(dat$height_sc)[1])
  }

  if ("ID" %in% names(dat)) {
    newdata <- newdata %>% mutate(ID = dat$ID[1])
  }

  newdata
}

average_predictive_draws <- function(draw_matrix, newdata) {
  # Convert predictive counts/proportions into one habitat-average value per
  # posterior draw by averaging across the replicated hectare-like units.
  as.data.frame(draw_matrix) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(cols = -draw, names_to = "obs_id", values_to = "estimate") %>%
    mutate(obs_id = as.integer(gsub("V", "", obs_id))) %>%
    left_join(newdata %>% mutate(obs_id = row_number()), by = "obs_id") %>%
    group_by(draw, habitat, height_filt, num_cells) %>%
    summarise(
      effective_hectares_per_habitat = dplyr::n(),
      estimate = mean(estimate),
      .groups = "drop"
    )
}

summarize_averaged_draws <- function(averaged_draws) {
  averaged_draws %>%
    group_by(habitat, height_filt, num_cells, effective_hectares_per_habitat) %>%
    summarise(
      mean = mean(estimate),
      lower95 = quantile(estimate, 0.025),
      upper95 = quantile(estimate, 0.975),
      lower80 = quantile(estimate, 0.10),
      upper80 = quantile(estimate, 0.90),
      .groups = "drop"
    ) %>%
    mutate(uncertainty_stream = "posterior_predictive_habitat_average", .before = 1) %>%
    select(
      uncertainty_stream, habitat, height_filt, num_cells,
      effective_hectares_per_habitat, mean, lower95, upper95, lower80, upper80
    )
}

thin_averaged_draws <- function(averaged_draws) {
  draw_ids <- sort(unique(averaged_draws$draw))
  thin_interval <- max(1, floor(length(draw_ids) / 500))
  keep_draws <- draw_ids[seq(1, length(draw_ids), by = thin_interval)]

  averaged_draws %>%
    filter(draw %in% keep_draws) %>%
    mutate(
      draw = match(draw, keep_draws),
      uncertainty_stream = "posterior_predictive_habitat_average",
      .after = draw
    ) %>%
    select(
      draw, uncertainty_stream, habitat, height_filt, num_cells,
      effective_hectares_per_habitat, estimate
    )
}

make_predictive_outputs <- function(height_filt, model_label, model_file) {
  if (!file.exists(model_file)) {
    stop("Expected model file not found: ", model_file)
  }

  log_progress("Generating predictive draws for saved model ", model_label, " (threshold ", height_filt, " m)")
  model_obj <- readRDS(model_file)

  if (is.null(model_obj$data)) {
    stop("Saved model does not contain data: ", model_label)
  }

  dat <- prepare_saved_model_data(model_obj$data, height_filt = height_filt)
  newdata <- build_prediction_newdata(
    dat,
    height_filt = height_filt,
    effective_hectares = effective_hectares_per_habitat
  )
  log_progress(
    "Running posterior_predict() for ", model_label, " on ", nrow(newdata),
    " replicated habitat rows (", effective_hectares_per_habitat, " representative hectares per habitat)"
  )

  predictive_counts <- posterior_predict(model_obj, newdata = newdata, re_formula = NA)
  predictive_prop <- sweep(predictive_counts, 2, newdata$num_cells, "/")
  log_progress(
    "Finished posterior_predict() for ", model_label,
    "; averaging predictive draws within habitat before scenario propagation"
  )

  averaged_draws <- average_predictive_draws(predictive_prop, newdata)
  pred_summary <- summarize_averaged_draws(averaged_draws)
  predictive_draws <- thin_averaged_draws(averaged_draws)

  summary_stem <- paste0(model_label, "_predictive_prediction_summary")
  draws_stem <- paste0(model_label, "_predictive_habitat_draws")

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
    predictive_draws,
    nr2_output_path(paste0(draws_stem, ".rds"))
  )
  log_progress("Saved predictive habitat outputs for ", model_label)

  list(
    pred_summary = pred_summary,
    predictive_draws = predictive_draws
  )
}

scenarios <- readRDS(nr2_input_path("MasterAllScenarios.rds"))
hab_carbon <- read.csv(nr2_input_path("allHabCarbon_60yr_withDelays.csv"))
hab_carbon <- as.data.table(hab_carbon)
log_progress("Loaded ", length(scenarios), " scenarios and habitat carbon lookup table")
log_progress(
  "Using posterior predictive habitat averages with ",
  effective_hectares_per_habitat,
  " representative hectares per habitat and posterior draw"
)

age_fun <- function(x) {
  x <- as.data.table(x)
  result <- x[hab_carbon, on = .(habitat = habitat, original_habitat = original_habitat),
              allow.cartesian = TRUE, nomatch = 0]
  result <- result[harvest_delay == 1]
  result
}

scenarios <- lapply(scenarios, age_fun)

cols_to_remove <- c(
  "full_carbon", "full_carbon_lwr", "full_carbon_upr",
  "ACD", "lwr_ACD", "upr_ACD", "X"
)

scenarios <- scenarios %>%
  map(~ .x %>% select(-any_of(cols_to_remove)))
log_progress("Prepared scenario inputs for predictive propagation")

function_scenario_60yr_predictive <- function(x, megatrees) {
  result <- x %>%
    left_join(megatrees, by = "functional_habitat", relationship = "many-to-many") %>%
    as.data.table()

  result[, megatree_hab_year := prop_megatrees * num_parcels]
  result <- result[, .(megatree_hab_year = sum(megatree_hab_year)),
                   by = .(scenarioName, scenarioStart, index, draw, production_target, true_year, original_habitat, habitat)]

  result <- result[, .(landscape_megatree_yr = sum(megatree_hab_year)),
                   by = .(scenarioName, scenarioStart, index, draw, production_target, true_year)]

  result[, .(megatree_60yr = sum(landscape_megatree_yr) / 1000),
         by = .(scenarioName, scenarioStart, index, draw, production_target)]
}

summarize_megatrees_predictive <- function(dt, max_h) {
  dt[, .(
    height_filt = max_h,
    effective_hectares_per_habitat = unique(effective_hectares_per_habitat)[1],
    landscape_prop = median(megatree_60yr, na.rm = TRUE),
    landscape_prop_lwr95 = quantile(megatree_60yr, 0.025, na.rm = TRUE),
    landscape_prop_upr95 = quantile(megatree_60yr, 0.975, na.rm = TRUE),
    landscape_prop_lwr80 = quantile(megatree_60yr, 0.10, na.rm = TRUE),
    landscape_prop_upr80 = quantile(megatree_60yr, 0.90, na.rm = TRUE),
    primary_SL_prop = median(primary_SL_prop, na.rm = TRUE),
    primary_SL_prop_lwr95 = quantile(primary_SL_prop, 0.025, na.rm = TRUE),
    primary_SL_prop_upr95 = quantile(primary_SL_prop, 0.975, na.rm = TRUE),
    primary_SL_prop_lwr80 = quantile(primary_SL_prop, 0.10, na.rm = TRUE),
    primary_SL_prop_upr80 = quantile(primary_SL_prop, 0.90, na.rm = TRUE)
  ), by = .(scenarioName, scenarioStart, index, production_target)]
}

process_one_threshold_predictive <- function(max_h, megatrees) {
  log_progress("Propagating predictive draws through scenarios for threshold ", max_h, " m")

  starting_landscape_megatrees <- megatrees %>%
    filter(functional_habitat == "primary") %>%
    mutate(primary_SL_prop = prop_megatrees * 61) %>%
    select(primary_SL_prop, draw)

  total_scenarios <- length(scenarios)
  outcomes_summary <- purrr::map2_dfr(scenarios, seq_along(scenarios), function(x, scenario_idx) {
    scenario_name <- if ("scenarioName" %in% names(x)) unique(x$scenarioName)[1] else paste0("scenario_", scenario_idx)
    log_progress(
      "Threshold ", max_h, " m: starting scenario ", scenario_idx, " of ", total_scenarios,
      " (", scenario_name, ")"
    )

    scenario_summary <- function_scenario_60yr_predictive(x, megatrees) %>%
      left_join(starting_landscape_megatrees, by = "draw") %>%
      as.data.table() %>%
      summarize_megatrees_predictive(max_h = max_h) %>%
      tibble::as_tibble()

    log_progress(
      "Threshold ", max_h, " m: finished scenario ", scenario_idx, " of ", total_scenarios,
      " (", scenario_name, ")"
    )

    scenario_summary
  })

  output <- outcomes_summary %>%
    mutate(
      outcome = "megatrees",
      uncertainty_stream = "posterior_predictive_habitat_average",
      .before = 1
    ) %>%
    select(
      uncertainty_stream, outcome, height_filt, effective_hectares_per_habitat,
      index, production_target, scenarioName, scenarioStart,
      landscape_prop, landscape_prop_lwr95, landscape_prop_upr95,
      landscape_prop_lwr80, landscape_prop_upr80,
      primary_SL_prop, primary_SL_prop_lwr95, primary_SL_prop_upr95,
      primary_SL_prop_lwr80, primary_SL_prop_upr80
    )

  saveRDS(
    output,
    nr2_output_path(paste0("full_nature_scenario_megatree_performance_predictive_", max_h, "m.rds"))
  )
  write.csv(
    output,
    nr2_output_path(paste0("full_nature_scenario_megatree_performance_predictive_", max_h, "m.csv")),
    row.names = FALSE
  )
  log_progress("Saved predictive scenario outputs for threshold ", max_h, " m")

  output
}

log_progress("Starting predictive habitat-output rebuild")
predictive_model_results <- purrr::pmap(
  existing_model_specs,
  make_predictive_outputs
)
log_progress("Finished predictive habitat-output rebuild")

predictive_pred_summaries <- purrr::map_dfr(predictive_model_results, "pred_summary")
predictive_draws_all <- purrr::map_dfr(predictive_model_results, "predictive_draws")

write.csv(
  predictive_pred_summaries,
  nr2_output_path("full_nature_megatrees_predictive_prediction_summary_by_habitat.csv"),
  row.names = FALSE
)
saveRDS(
  predictive_pred_summaries,
  nr2_output_path("full_nature_megatrees_predictive_prediction_summary_by_habitat.rds")
)
saveRDS(
  predictive_draws_all,
  nr2_output_path("full_nature_megatrees_predictive_habitat_draws_all_thresholds.rds")
)
log_progress("Saved combined predictive habitat outputs across all thresholds")

log_progress("Starting predictive scenario propagation across all thresholds")
scenario_predictive_outputs <- purrr::map2_dfr(
  existing_model_specs$height_filt,
  predictive_model_results,
  ~ process_one_threshold_predictive(
    max_h = .x,
    megatrees = .y$predictive_draws %>%
      rename(
        prop_megatrees = estimate,
        functional_habitat = habitat
      ) %>%
      filter(!is.na(draw))
  )
)
log_progress("Finished predictive scenario propagation across all thresholds")

saveRDS(
  scenario_predictive_outputs,
  nr2_output_path("full_nature_scenario_megatree_performance_predictive_all_thresholds.rds")
)
write.csv(
  scenario_predictive_outputs,
  nr2_output_path("full_nature_scenario_megatree_performance_predictive_all_thresholds.csv"),
  row.names = FALSE
)

log_progress("Standalone predictive-draw workflow complete")
log_progress("Check Outputs/NR2 for files with the '_predictive_' naming convention")
