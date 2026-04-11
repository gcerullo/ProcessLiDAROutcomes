# Propagate posterior megatree draws through scenarios.

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

megaTreeThreshold <- 50

scenarios <- readRDS(nr2_input_path("MasterAllScenarios.rds"))
hab_carbon <- read.csv(nr2_input_path("allHabCarbon_60yr_withDelays.csv"))
megatrees <- readRDS(nr2_output_path("megatrees_draws_per_hab.rds"))

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

megatrees <- megatrees %>%
  filter(height_filt == megaTreeThreshold) %>%
  rename(
    prop_megatrees = estimate,
    functional_habitat = habitat
  ) %>%
  filter(!is.na(draw))

function_scenario_60yr_uncertainty <- function(x) {
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

outcomes <- map(scenarios, function_scenario_60yr_uncertainty)

starting_landscape_megatrees <- megatrees %>%
  filter(functional_habitat == "primary") %>%
  mutate(primary_SL_prop = prop_megatrees * 61) %>%
  select(primary_SL_prop, draw)

add_SL_fun <- function(x) {
  x %>% left_join(starting_landscape_megatrees, by = "draw")
}

outcomes <- map(outcomes, add_SL_fun)

summarize_megatrees <- function(dt) {
  dt[, .(
    landscape_prop = median(megatree_60yr, na.rm = TRUE),
    landscape_prop_lwr = quantile(megatree_60yr, 0.025, na.rm = TRUE),
    landscape_prop_upr = quantile(megatree_60yr, 0.975, na.rm = TRUE),
    primary_SL_prop = median(primary_SL_prop, na.rm = TRUE),
    primary_SL_prop_lwr = quantile(primary_SL_prop, 0.025, na.rm = TRUE),
    primary_SL_prop_upr = quantile(primary_SL_prop, 0.975, na.rm = TRUE)
  ), by = .(scenarioName, scenarioStart, index, production_target)]
}

outcomes_summary <- map(outcomes, summarize_megatrees) %>%
  rbindlist()

output <- outcomes_summary %>%
  select(
    index, production_target, scenarioName, scenarioStart,
    landscape_prop, landscape_prop_lwr, landscape_prop_upr,
    primary_SL_prop, primary_SL_prop_lwr, primary_SL_prop_upr
  ) %>%
  cbind(outcome = "megatrees")

saveRDS(output, nr2_output_path("MasterMegatreePerformance_with_uncertainty.rds"))

message("NR2 script 02 complete.")
