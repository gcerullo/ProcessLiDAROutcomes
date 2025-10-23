# propagate posterior draws thru scenarios 

library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr) 
library(cowplot)
library(purrr)

#Read in Inputs ####
#read in the scenario parametres containing conversion factors for converting from point to parcel/entire landscape  
source('Inputs/FixedScenarioParmams.R')

#Define params ####

#define megatree threshold (45 to 70 in increments of 5)
megaTreeThreshold <- 50

#Read in Data ####
#----------------read in scenarios -------------------------------
#yield matched scenarios where 1/30th of plantation conversion happens annually
scenarios <- readRDS("Inputs/MasterAllScenarios.rds")
scenario_composition <- rbindlist(scenarios, use.names=TRUE)

#read in hab carbon through time - this dataframe allows us to add a temporal gradient
hab_carbon <-read.csv("Inputs/allHabCarbon_60yr_withDelays.csv")

#read in posterior draws of height by raster (Output from CalculateMegatreesPerHab.R)
megatrees <- readRDS("Outputs/megatrees_draws_per_hab.rds") 

names(scenario_composition)
names(hab_carbon)

#convert hab_harbon to data-table 
hab_carbon <- as.data.table(hab_carbon)

x <- scenarios[[1]]

age_fun <- function(x) {
  x <- as.data.table(x)
  
  # Perform the left join
  result <- x[hab_carbon, on = .(habitat = habitat, original_habitat = original_habitat), 
              allow.cartesian = TRUE, nomatch = 0]
  
  # Filter for harvest_delay == 1 for simplicity and efficiency of runtime and because we assume that 
  #megatrees recover equivalently to deaths, do delay dynamic is unnecessary
  result <- result[harvest_delay == 1]
  
  return(result)
}

#code starts here ####
# ------Add temporal carbon data to scenarios -----------------
scenarios <- lapply(scenarios,age_fun)


#retain only necessary columns to speed up downstream processing
names(scenarios[[1]])

cols_to_remove <- c(
  "full_carbon", "full_carbon_lwr", "full_carbon_upr",
  "ACD", "lwr_ACD", "upr_ACD","X"
)

scenarios <- scenarios %>%
  map(~ .x %>% select(-any_of(cols_to_remove)))

#filter only pre-definded threshold
megatrees <- megatrees %>%  filter(height_filt == megaTreeThreshold) %>% 
  rename(prop_megatrees = estimate, 
         functional_habitat = habitat) %>% 
  filter(!is.na(draw))

#remind ourselves of scenario structure
x <- scenarios[[1]]
y <- x %>% filter(index == "11858") #remember what a single scenario is
#so joining add 500 posterior draw for each functional habitat type and year
m <- y %>% left_join(megatrees) 
megatrees %>% group_by( functional_habitat) %>%count()

#APPLY FUNCTION FOR ESTIMATING LANDSCAPE MEGATREES #####

#function to calculate megatree coverage
function_scenario_60yr_uncertainty <-  function(x){
  
  #add 500 posterior draws for each row of scenario 
  result <- x %>%  left_join(megatrees, by = "functional_habitat", relationship = "many-to-many")
  
    # Step 2: Group and summarize mhatrees for specific year and habitat transition
  #[for each true year and habitat transition, calculate megatree years combined across the staggered
  # harvesting schedule - note this is only really needed if we consider multiple harvest delays, but included here incase applying multiple harvest delays) 
  
  #calculate the area of each megatree coverage per year
  result[, megatree_hab_year := prop_megatrees * num_parcels]
  result <- result[, .(megatree_hab_year = sum(megatree_hab_year)), 
                        by = .(scenarioName, scenarioStart,index, draw, production_target, true_year, original_habitat, habitat)]
  
  # Step 3: Calculate landscape megatrees per year
  #[Across habitat type transitions (e.g for ALL hab_parcel transitions) in a scenario, calculate occupancy for a given year]
  result <- result[, .(landscape_megatree_yr = sum(megatree_hab_year)), 
                   by = .(scenarioName, scenarioStart,index, draw, production_target, true_year)]
  
  # Step 4: Calculate scenario's landscape megatrees across years
  result <- result[, .(megatree_60yr = sum(landscape_megatree_yr)/1000), 
                   by = .(scenarioName, scenarioStart,index, draw, production_target)]
  
  
  return(result)
}


#test_scenario_run <- scenarios[[1]]
#outcomes <- function_scenario_60yr_uncertainty(test_scenario_run)
outcomes <- map(scenarios, function_scenario_60yr_uncertainty)

#---------  Calculate  starting landscape megatree proportion -------

# #calculate the proportion of megatrees in an old-growth  starting landscape 
starting_landscape_megatrees <- megatrees %>% filter(functional_habitat == "primary") %>%
  mutate(primary_SL_prop = prop_megatrees*61) %>% select(primary_SL_prop, draw)


#join the primary SL prop coverage to to each scenario  

#add the starting landscape megatrees to each scenarios
add_SL_fun <- function(x){
  x %>% left_join(starting_landscape_megatrees, by = "draw")
}

outcomes <- map(outcomes, add_SL_fun)


#now we can calculate using the posterior draws the uncertainty in our calculation

summarize_megatrees <- function(dt) {
  dt[, .(
    landscape_prop = median(megatree_60yr, na.rm = TRUE),
    landscape_prop_lwr    = quantile(megatree_60yr, 0.025, na.rm = TRUE),
    landscape_prop_upr    = quantile(megatree_60yr, 0.975, na.rm = TRUE),
    
    primary_SL_prop    = median(primary_SL_prop, na.rm = TRUE),
    primary_SL_prop_lwr       = quantile(primary_SL_prop, 0.025, na.rm = TRUE),
    primary_SL_prop_upr       = quantile(primary_SL_prop, 0.975, na.rm = TRUE)
  ), by = .(scenarioName, scenarioStart, index, production_target)]
}

outcomes_summary <- map(outcomes, summarize_megatrees)
outcomes_summary <- rbindlist(outcomes_summary)

#-----EXPORT OUTCOME PERFORMANCE for consolidated figure of all outcomes -----
getwd()
names(outcomes_summary)
output <- outcomes_summary %>% select(index, production_target, scenarioName,scenarioStart,
                                 landscape_prop,landscape_prop_lwr, landscape_prop_upr,
                                 primary_SL_prop,primary_SL_prop_lwr,primary_SL_prop_upr) %>% cbind(outcome = "megatrees")
saveRDS(output, "Outputs/nature_output/MasterMegatreePerformance_with_uncertainty.rds")

