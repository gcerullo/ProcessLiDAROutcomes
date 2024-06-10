#10.06.24
#Assess the LiDA outcomes of different scenarios 

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


#Read in Data ####
#----------------read in scenarios -------------------------------
#NEEDS TO BE UPDATED WITH CORRECT YIELDS !!!!!!!!!
#yield matched scenarios where 1/30th of plantation conversion happens annually
scenarios <- readRDS("Inputs/allScenariosStaggered.rds")
scenario_composition <- rbindlist(scenarios, use.names=TRUE)

#read in hab carbon through time - this dataframe allows us to add a temporal gradient
hab_carbon <-read.csv("Inputs/allHabCarbon_60yr_withDelays.csv")

# ------Add temporal carbon data to scenairos -----------------
names(scenario_composition)
names(hab_carbon)

#convert hab_harbon to data-table 
hab_carbon <- as.data.table(hab_carbon)

age_fun <- function(x){
  x <- as.data.table(x)
  # Perform the left join operation in data table 
  result <- x[hab_carbon, on = .(habitat = habitat, original_habitat = original_habitat), allow.cartesian = TRUE, nomatch = 0]
}

scenarios <- lapply(scenarios,age_fun)




