#10.06.24
#Assess the LiDAR outcomes of different scenarios 

library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr) 
library(cowplot)

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

#read in raw values showing megatrees per ha in diff habitat types (Output from CalculateMegatreesPerHab.R)
megatrees <- read.csv("Outputs/rawMegatreesPerHa.csv") 

names(scenario_composition)
names(hab_carbon)

#convert hab_harbon to data-table 
hab_carbon <- as.data.table(hab_carbon)

age_fun <- function(x){
  x <- as.data.table(x)
  # Perform the left join operation in data table 
  result <- x[hab_carbon, on = .(habitat = habitat, original_habitat = original_habitat), allow.cartesian = TRUE, nomatch = 0]
}

#code starts here ####
# ------Add temporal carbon data to scenairos -----------------
scenarios <- lapply(scenarios,age_fun)

#filer only hexagons for which most 1m cells were captured
megatrees <- megatrees %>% filter(num_cells > 8650)
megatrees <- megatrees %>%  filter(height_filt == megaTreeThreshold) 

test_megatrees <- megatrees %>% group_by(habitat) %>% 
  summarise(mean = mean(propBigTrees), 
            sd = sd(propBigTrees)) %>% 
  mutate(functional_habitat = habitat)

#fit a binomial regression that estimates the proportion of big tree cells (successes) per number of 
#sampled LiDAR cells
y= cbind(megatrees$bigTrees, megatrees$num_cells-megatrees$bigTrees)

#estimate the proportion of cells that are MEGATREES for of a given habitat type 
#megatreeModel <-  glm(y ~ habitat-1, family = "binomial", data = megatrees)
#summary(megatreeModel)
# Fit model with quasibinomial to account for overdispersion
megatreeModel <- glm(y ~ habitat-1, family = "quasibinomial", data = megatrees)
summary(megatreeModel)

#extract the values for each habitat type  
model_predictions<-data.frame(habitat = unique(megatrees$habitat))

# #for each habitat type, predict the proportion of each habitat type that has tall trees, and the standard error
# megatrees_prediction <- as.data.frame(predict(megatreeModel, newdata = model_predictions, type = 'response', se.fit = TRUE)) %>%  
#   cbind((habitat = unique(megatrees$habitat)))  %>%  
#   select(-residual.scale) %>% 
#   rename(prop_megatrees = fit, 
#          se_megatrees = se.fit,
#          habitat = 3) %>% 
#   #multiply SE by 1.96 to get the 5 and 95th percentile
#   mutate(upr_megatrees=  prop_megatrees + se_megatrees * 1.96,
#          lwr_megatrees= prop_megatrees- se_megatrees * 1.96) %>%
#   select(prop_megatrees, habitat, upr_megatrees,lwr_megatrees) %>%
#    rename(functional_habitat = habitat)
# 
# #------- calculate the amount of megatrees per scenario landscape ---------
# #get number of staggered harvests to define harvest window (this must match harvests)
# J <- scenarios[[12]] 
# harvest_window <- J$harvest_delay %>% unique %>% length()
# 
# ##temporal megatrees function 
# scenario_megatree_fun <- function(x){
#   x %>%  left_join(megatrees_prediction, by = "functional_habitat") %>%  
#     #1. assuming 1/30th of of each habitat type is applied to each harvesting delay schedule
#     #, calculate the total ACD for a given habitat type in a given year
#     #NB- if there is no habitat transition, then don't need to divide by harvest window 
#     mutate(
#       prop_megatrees_stag =  prop_megatrees * num_parcels / harvest_window,
#       lwr_megatrees_stag = lwr_megatrees * num_parcels / harvest_window,
#       upr_megatrees_stag = upr_megatrees * num_parcels / harvest_window)  %>% 
#     
#     
#   #2. for each true year and habitat transition, calculate megatrees combined across the staggered
#   #harvesting schedule (i.e. the megatrees in a given habitat transition for a given year) 
#   group_by(index,production_target, original_habitat, habitat, true_year) %>%  
#     mutate(hab_megatree_year = sum(prop_megatrees_stag,na.rm = TRUE), 
#            hab_megatree_year_lwr= sum(lwr_megatrees_stag,na.rm = TRUE), 
#            hab_megatree_year_upr = sum(upr_megatrees_stag,na.rm = TRUE)) %>%  ungroup %>%  
#     
#     #select a single harvest delay worth of data, as we now have calculated megatrees across harvesting schedules
#     filter(harvest_delay == 15) %>% select(-harvest_delay) %>% 
#     
#     
#     #4. Across habitat type transitions (e.g for all hab_parcel transitions) in a scenario, calculate megatree for a given year
#     group_by(index, production_target, true_year) %>%  
#     mutate(scen_megatree_year = sum(hab_megatree_year), 
#            scen_megatree_year_lwr= sum(hab_megatree_year_lwr), 
#            scen_megatree_year_upr = sum(hab_megatree_year_upr)) %>%  ungroup() %>% 
#     
#     
#     #now we make sure we only have one row for each scenario and year, showing scen_megatrees_year
#     select(index, production_target,scenarioName,scenarioStart, true_year, 
#            scen_megatree_year,scen_megatree_year_lwr, scen_megatree_year_upr) %>%  
#     group_by(true_year,index,production_target) %>%  slice(1) %>% 
#     ungroup() %>%  
#     
#     #calculate the number of parcels across the entire scenario landscape covered in megatrees
#     #[/1000 then to give the proortion of the landscape covered by big trees]
#     
#     #deforested,albizia and eucalyptus have NA megatree values, so we use na.rm =TRUE to ignore these when computing landscape sums 
#     #this sum = megatree years (e.g. a value of 60 = each year being covering entire canopy)
#     group_by(index, production_target) %>% 
#     mutate(landscape_prop = sum(scen_megatree_year)/1000, 
#            landscape_prop_lwr = sum(scen_megatree_year_lwr)/1000, 
#            landscape_prop_upr = sum(scen_megatree_year_upr)/1000) %>%  
#     ungroup() %>% 
#     
#     
#     #now we make sure we only have one row for each scenario, showing scen_megatree_year
#     select(index, production_target,scenarioName,scenarioStart, 
#            landscape_prop,landscape_prop_lwr, landscape_prop_upr) %>%  
#     group_by(production_target,index) %>%  slice(1) %>% 
#     ungroup()
# }
# 
#for each habitat type, predict the proportion of each habitat type that has tall trees, and the standard error
megatrees_prediction <- as.data.frame(predict(megatreeModel, newdata = model_predictions, type = 'response', se.fit = TRUE)) %>%  
  cbind((habitat = unique(megatrees$habitat)))  %>%  
  select(-residual.scale) %>% 
  rename(prop_megatrees = fit, 
         se_megatrees = se.fit,
         habitat = 3) %>% 
  # Store variance instead of SE for proper error propagation
  mutate(var_megatrees = se_megatrees^2,
         upr_megatrees = prop_megatrees + se_megatrees * 1.96,
         lwr_megatrees = prop_megatrees - se_megatrees * 1.96) %>%
  select(prop_megatrees, habitat, se_megatrees, var_megatrees, upr_megatrees, lwr_megatrees) %>%
  rename(functional_habitat = habitat)

#------- calculate the amount of megatrees per scenario landscape ---------
#get number of staggered harvests to define harvest window (this must match harvests)
J <- scenarios[[12]] 
harvest_window <- J$harvest_delay %>% unique %>% length()

##temporal megatrees function 
scenario_megatree_fun <- function(x){
  x %>%  left_join(megatrees_prediction, by = "functional_habitat") %>%  
    #1. assuming 1/30th of of each habitat type is applied to each harvesting delay schedule
    #, calculate the total ACD for a given habitat type in a given year
    #NB- if there is no habitat transition, then don't need to divide by harvest window 
    mutate(
      # Scale mean by num_parcels/harvest_window
      prop_megatrees_stag = prop_megatrees * num_parcels / harvest_window,
      # Scale variance correctly: Var(aX) = a^2 * Var(X)
      var_megatrees_stag = var_megatrees * (num_parcels / harvest_window)^2) %>% 
    
    #2. for each true year and habitat transition, calculate megatrees combined across the staggered
    #harvesting schedule (i.e. the megatrees in a given habitat transition for a given year) 
    group_by(index, production_target, original_habitat, habitat, true_year) %>%  
    summarize(
      # Sum the means
      hab_megatree_year = sum(prop_megatrees_stag, na.rm = TRUE),
      # Sum the variances (proper error propagation for sum of independent variables)
      hab_megatree_year_var = sum(var_megatrees_stag, na.rm = TRUE),
      # Other variables to keep
      scenarioName = first(scenarioName),
      scenarioStart = first(scenarioStart),
      .groups = "drop") %>%
    # Calculate SE and confidence intervals after summing
    mutate(
      hab_megatree_year_se = sqrt(hab_megatree_year_var),
      hab_megatree_year_lwr = hab_megatree_year - 1.96 * hab_megatree_year_se,
      hab_megatree_year_upr = hab_megatree_year + 1.96 * hab_megatree_year_se) %>%
    
    #4. Across habitat type transitions (e.g for all hab_parcel transitions) in a scenario, calculate megatree for a given year
    group_by(index, production_target, true_year) %>%  
    summarize(
      scen_megatree_year = sum(hab_megatree_year, na.rm = TRUE),
      # Sum variances for proper error propagation
      scen_megatree_year_var = sum(hab_megatree_year_var, na.rm = TRUE),
      scenarioName = first(scenarioName),
      scenarioStart = first(scenarioStart),
      .groups = "drop") %>%
    # Calculate SE and confidence intervals after summing
    mutate(
      scen_megatree_year_se = sqrt(scen_megatree_year_var),
      scen_megatree_year_lwr = scen_megatree_year - 1.96 * scen_megatree_year_se,
      scen_megatree_year_upr = scen_megatree_year + 1.96 * scen_megatree_year_se) %>%
    
    #calculate the number of parcels across the entire scenario landscape covered in megatrees
    #[/1000 then to give the proportion of the landscape covered by big trees]
    group_by(index, production_target) %>% 
    mutate(
      landscape_prop = sum(scen_megatree_year)/1000,
      # Calculate total variance across all years for landscape proportion
      landscape_var = sum(scen_megatree_year_var)/1000^2) %>%
    # Calculate SE and confidence intervals after summing and scaling
    mutate(
      landscape_se = sqrt(landscape_var),
      landscape_prop_lwr = landscape_prop - 1.96 * landscape_se,
      landscape_prop_upr = landscape_prop + 1.96 * landscape_se) %>%
    ungroup() %>% 
    
    #now we make sure we only have one row for each scenario, showing scen_megatree_year
    select(index, production_target, scenarioName, scenarioStart, 
           landscape_prop, landscape_prop_lwr, landscape_prop_upr) %>%  
    group_by(production_target, index) %>%  
    slice(1) %>% 
    ungroup()
}

#############################
#calculated proportional megatrees through time 
outcomes<- lapply(scenarios, scenario_megatree_fun)

#---------  Calculate  starting landscape megatree proportion -------

# #calculate the proportion of megatrees in an old-growth  starting landscape 
# starting_landscape_megatrees <- megatrees_prediction %>% filter(functional_habitat == "primary") %>% 
#   mutate(primary_SL_prop = prop_megatrees *61, 
#          primary_SL_prop_lwr = lwr_megatrees*61,
#          primary_SL_prop_upr = upr_megatrees*61) %>% select(primary_SL_prop,primary_SL_prop_lwr,primary_SL_prop_upr)
# 

#calculate the proportion of megatrees in an old-growth starting landscape 
starting_landscape_megatrees <- megatrees_prediction %>% 
  filter(functional_habitat == "primary") %>% 
  # Calculate variance for proper error propagation
  mutate(
    # Store variance for error propagation
    var_megatrees = se_megatrees^2,
    # Scale mean by 61
    primary_SL_prop = prop_megatrees * 61,
    # Scale variance correctly: Var(aX) = a^2 * Var(X)
    primary_SL_var = var_megatrees * 61^2,
    # Calculate SE after scaling
    primary_SL_se = sqrt(primary_SL_var),
    # Calculate confidence intervals using scaled SE
    primary_SL_prop_lwr = primary_SL_prop - 1.96 * primary_SL_se,
    primary_SL_prop_upr = primary_SL_prop + 1.96 * primary_SL_se
  ) %>% 
  select(primary_SL_prop, primary_SL_prop_lwr, primary_SL_prop_upr)

#join the primary SL prop coverage to to each scenario  

#add the starting landscape megatrees to each scenarios
add_SL_fun <- function(x){
  x %>% cbind(starting_landscape_megatrees)
}

outcomes <- lapply(outcomes, add_SL_fun)
outcomes_df <- rbindlist(outcomes,use.names=TRUE)



#-----EXPORT OUTCOME PERFORMANCE for consolidated figure of all outcomes -----
getwd()
names(outcomes_df)
output <- outcomes_df %>% select(index, production_target, scenarioName,scenarioStart,
                                 landscape_prop,landscape_prop_lwr, landscape_prop_upr,
                                 primary_SL_prop,primary_SL_prop_lwr,primary_SL_prop_upr) %>% cbind(outcome = "megatrees")
saveRDS(output, "Outputs/MasterMegatreePerformance_with_uncertainty.rds")

