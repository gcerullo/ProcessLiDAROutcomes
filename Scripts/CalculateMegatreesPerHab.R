#10.06.24
#Calculate the number of megatrees per habitat type from LiDAR data covering the concession 

library(terra)
library(tidyverse)
library(data.table)
library(ggpubr)

#define inputs ####

#define proportion of cells surpassing a height threshold per 1ha hexagon
max_height <- c(45, 50, 55, 60, 65, 70)

#read in a hexagon grid of 1ha cells for Sabah  
hex<- vect("RawData/TreeHeightDanum/Sabah100HexShapefile.shp")


# read in data ####
#read in rasters of CTH for different habitats 
primary <- rast("RawData/TreeHeightDanum/Primary.tif")
once_logged <- rast("RawData/TreeHeightDanum/Once_Logged.tif")
restore <- rast("RawData/TreeHeightDanum/Restored.tif")
twice_logged <- rast("RawData/TreeHeightDanum/Twice_Logged.tif")

plot(primary)
plot(once_logged)
plot(restore)
plot(twice_logged)

#correct for primary forest cloud cover ####

#there is a cloud in the primary forest data (dark green smudge on the western side of the raster) 
#that could bias tree height estimates (looks like super tall tree cluster) and that must therefore be removed

# Create an Extent object rectangled around the cloud
cloud_ext <- ext(c(584475, 584670, 548050, 548350))

# View the cloud
cloud <- crop(primary, cloud_ext)
plot(cloud)

#set the values within the cloud extent to NA 
primary[cloud_ext] <- NA
plot(primary)


#assign each 1m raster cell to the ID of a hexagon grid ####

#ONLY NEED TO RUN ONCE: 
# #extract values for each 1ha hexagon [do once]
# primary_vals <- terra::extract(primary,hex) %>%  na.omit()
# once_vals <- terra::extract(once_logged,hex) %>%  na.omit()
# restored_vals <- terra::extract(restore,hex) %>%  na.omit()
#twice_vals <- terra::extract(twice_logged,hex) %>%  na.omit()
# 
# 
# #save outputs 
# saveRDS(primary_vals, "Outputs/primary1haCHMvals.rds")
# saveRDS(once_vals, "Outputs/onceLogged_vals1haCHMvals.rds")
# saveRDS(restored_vals, "Outputs/restored_vals1haCHMvals.rds")
# saveRDS(twice_vals, "Outputs/twiceLogged_vals1haCHMvals.rds")

#CAN reSTART FROM HERE: 
primary_vals <- readRDS("Outputs/primary1haCHMvals.rds") %>% rename(CH = Primary) %>% cbind(habitat = "primary")
once_vals <-  readRDS("Outputs/onceLogged_vals1haCHMvals.rds")  %>% rename(CH = Band_1) %>% cbind(habitat = "once-logged")
restored_vals <- readRDS("Outputs/restored_vals1haCHMvals.rds")  %>% rename(CH = Restored) %>%  cbind(habitat = "restored")
twice_vals <- readRDS("Outputs/twiceLogged_vals1haCHMvals.rds") %>% rename(CH = Band_1) %>% cbind(habitat= "twice-logged")

#calculate number of 1m cells per hex grid
num_cells <- function(x){
  x %>% group_by(ID) %>% mutate(num_cells = n()) %>% ungroup()
}

primary_vals <- num_cells(primary_vals)
once_vals <- num_cells(once_vals)
restored_vals <- num_cells(restored_vals)
twice_vals <- num_cells(twice_vals)

# Define the function CH_thresh, which uses different max height thresholds to determine the proportion of the canopy 
#surpassing a given height threshold in each 1ha hexagon
CH_thresh <- function(x, max_height) {
  canopy_prop <- x %>%
    group_by(ID) %>%
    filter(CH > max_height) %>%
    mutate(bigTrees = n()) %>%
    mutate(propBigTrees = bigTrees / num_cells) %>%
    slice(1) %>%
    ungroup()
  canopy_prop$height_filt <- max_height
  return(canopy_prop)
}


# Create a list to store the results
CH_list <- list()

# Apply CH_thresh for each value in max_height and each dataset
datasets <- list(primary_vals, once_vals, restored_vals, twice_vals)
dataset_names <- c("CH_P", "CH_1", "CH_R", "CH_2")

for (i in seq_along(datasets)) {
  CH_list[[dataset_names[i]]] <- lapply(max_height, function(x) {
    CH_thresh(datasets[[i]], x)
  })
}

# Combine the results into a single dataframe for plotting
df <- do.call(rbind, CH_list)
df <- rbindlist(df) %>%   mutate(habitat = factor(habitat, levels = c("primary", "once-logged", "restored", "twice-logged")))


# Calculate mean and standard deviation for canopy height
result_df <- df %>%
  group_by(habitat, height_filt) %>% 
  summarize(
    meanProp = mean(propBigTrees),
    se = sd(propBigTrees)/sqrt(n()),
    sd = sd(propBigTrees), 
    minProp_se = meanProp - se, 
    maxProp_se = meanProp + se,
    minProp_sd = meanProp - sd, 
    minProp_sd = meanProp - sd, 
    
    
  ) %>%  
  #reorder 
  mutate(habitat = factor(habitat, levels = c("primary", "once-logged", "restored", "twice-logged")))



#plot proportion of canopy over height threshold####

df %>%  group_by(habitat, height_filt) %>% 
  
  ggplot() +
  #plt raw data
  geom_point(data = df, aes(x = habitat, y = propBigTrees), alpha = 0.1, colour = "grey90", fill = "transparent",
             position = position_jitter(width = 0.1, height = 0.1),  size =1) + 
  #plt std error
  geom_errorbar(data = result_df, aes(x = habitat, y = meanProp, ymin = minProp_se, ymax = maxProp_se, colour = habitat), width = 0.2) +
  #plt mean
  geom_point(data = result_df, aes(habitat, meanProp, colour= habitat), size =2)+
  #  ylim(0,1)+
  facet_wrap(~ height_filt, labeller = label_parsed, ncol = 2, scales = "free_y") +
  scale_fill_viridis_d(alpha = 0.9)+
  labs(x = "Habitat", y = "Proportion of canopy/ha over height threshold") +
  theme_pubr(base_size = 18)+
  theme(axis.text.x = element_blank(), 
        legend.position = "bottom")



#save outputs ####
write.csv(df, "Outputs/rawMegatreesPerHa.csv")
write.csv(result_df, "Outputs/ProportionalCanopyHeightByHabitat.csv")











