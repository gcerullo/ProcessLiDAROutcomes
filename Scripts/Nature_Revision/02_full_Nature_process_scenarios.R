# Propagate threshold-specific posterior megatree draws through scenarios.

source(file.path("Scripts", "Nature_Revision", "config.R"))

library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr)
library(cowplot)
library(purrr)

source(nr2_input_path("FixedScenarioParmams.R"))

height_thresholds <- c(45, 50, 55)

scenarios <- readRDS(nr2_input_path("MasterAllScenarios.rds"))
hab_carbon <- read.csv(nr2_input_path("allHabCarbon_60yr_withDelays.csv"))
hab_carbon <- as.data.table(hab_carbon)

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

function_scenario_60yr_uncertainty <- function(x, megatrees) {
  result <- x %>% left_join(megatrees, by = "functional_habitat", relationship = "many-to-many")

  result[, megatree_hab_year := prop_megatrees * num_parcels]
  result <- result[, .(megatree_hab_year = sum(megatree_hab_year)),
                   by = .(scenarioName, scenarioStart, index, draw, production_target, true_year, original_habitat, habitat)]

  result <- result[, .(landscape_megatree_yr = sum(megatree_hab_year)),
                   by = .(scenarioName, scenarioStart, index, draw, production_target, true_year)]

  result <- result[, .(megatree_60yr = sum(landscape_megatree_yr) / 1000),
                   by = .(scenarioName, scenarioStart, index, draw, production_target)]

  result
}

summarize_megatrees <- function(dt, max_h) {
  dt[, .(
    height_filt = max_h,
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

process_one_threshold <- function(max_h) {
  message("Processing scenarios for threshold ", max_h, " m")

  draws_file <- paste0("nature_megatrees_zibb_", max_h, "m_habitat_draws.rds")

  megatrees <- readRDS(
    nr2_output_path(draws_file)
  ) %>%
    rename(
      prop_megatrees = estimate,
      functional_habitat = habitat
    ) %>%
    filter(!is.na(draw))

  outcomes <- map(scenarios, function(x) function_scenario_60yr_uncertainty(x, megatrees))

  starting_landscape_megatrees <- megatrees %>%
    filter(functional_habitat == "primary") %>%
    mutate(primary_SL_prop = prop_megatrees * 61) %>%
    select(primary_SL_prop, draw)

  add_SL_fun <- function(x) {
    x %>% left_join(starting_landscape_megatrees, by = "draw")
  }

  outcomes <- map(outcomes, add_SL_fun)

  outcomes_summary <- map(outcomes, ~ summarize_megatrees(.x, max_h)) %>%
    rbindlist()

  output <- outcomes_summary %>%
    select(
      height_filt, index, production_target, scenarioName, scenarioStart,
      landscape_prop, landscape_prop_lwr95, landscape_prop_upr95,
      landscape_prop_lwr80, landscape_prop_upr80,
      primary_SL_prop, primary_SL_prop_lwr95, primary_SL_prop_upr95,
      primary_SL_prop_lwr80, primary_SL_prop_upr80
    ) %>%
    cbind(outcome = "megatrees")

  saveRDS(
    output,
    nr2_output_path(paste0("full_nature_scenario_megatree_performance_", max_h, "m.rds"))
  )
  write.csv(
    output,
    nr2_output_path(paste0("full_nature_scenario_megatree_performance_", max_h, "m.csv")),
    row.names = FALSE
  )

  output
}

all_threshold_outputs <- map_dfr(height_thresholds, process_one_threshold)

saveRDS(
  all_threshold_outputs,
  nr2_output_path("full_nature_scenario_megatree_performance_all_thresholds.rds")
)
write.csv(
  all_threshold_outputs,
  nr2_output_path("full_nature_scenario_megatree_performance_all_thresholds.csv"),
  row.names = FALSE
)

message("Full Nature scenario processing complete.")
message("Check Outputs/NR2 for threshold-specific and combined scenario outputs.")
