#10.06.24
#Calculate the number of megatrees per habitat type from LiDAR data covering the concession 

library(terra)
library(tidyverse)
library(data.table)
library(ggpubr)
library(terra)
library(brms)
library(tidybayes)


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

#### Load extracted CHM values per 1ha hex ####
primary_vals <- readRDS("Outputs/primary1haCHMvals.rds") %>%
  rename(CH = Primary) %>%
  mutate(habitat = "primary")

once_vals <- readRDS("Outputs/onceLogged_vals1haCHMvals.rds") %>%
  rename(CH = Band_1) %>%
  mutate(habitat = "once-logged")

restored_vals <- readRDS("Outputs/restored_vals1haCHMvals.rds") %>%
  rename(CH = Restored) %>%
  mutate(habitat = "restored")

twice_vals <- readRDS("Outputs/twiceLogged_vals1haCHMvals.rds") %>%
  rename(CH = Band_1) %>%
  mutate(habitat = "twice-logged")

# Combine all habitat datasets
all_vals <- bind_rows(primary_vals, once_vals, restored_vals, twice_vals)

# Count number of 1m² cells per hex ID
all_vals <- all_vals %>%
  group_by(ID) %>%
  mutate(num_cells = n()) %>%
  ungroup()

#### Function to calculate big tree counts per height threshold ####
CH_thresh <- function(x, max_height) {
  x %>%
    group_by(ID, habitat, num_cells) %>%
    summarise(
      bigTrees = sum(CH > max_height),
      .groups = "drop"
    ) %>%
    mutate(
      height_filt = max_height,
      propBigTrees = bigTrees / num_cells
    )
}

# Apply for all height thresholds
df_list <- lapply(max_height, function(h) CH_thresh(all_vals, h))
df <- bind_rows(df_list)

# Ensure correct factor order
df$habitat <- factor(df$habitat, levels = c("primary", "once-logged", "restored", "twice-logged"))

#filer only hexagons for which most 1m cells were captured
df <- df %>% filter(num_cells > 6000)

# #prepare to fit a beta regression (which can't handle exact 0-1 values, so transform as recommended from Smithson & Verkuilen, 2006
range(df$propBigTrees)
df <- df %>%
  mutate(prop_adj = propBigTrees+0.000001)

# # now fit the model
# bayes_model_binom <- brm(
#   formula = bigTrees | trials(num_cells) ~ habitat * height_filt,
#   data = df,
#   family = binomial(),
#   prior = c(
#     prior(normal(0, 1), class = "b"),
#     prior(normal(logit_mean, 1), class = "Intercept")
#   ),
#   chains = 4, iter = 4000, cores = 4, seed = 123
# )

#### Fit Bayesian beta regression model ####
bayes_model <- brm(
  formula = prop_adj ~ habitat * height_filt, # +(1|ID),
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

# Save model
saveRDS(bayes_model, "Outputs/bayesian_megatree_model_22_10.rds")
saveRDS(bayes_model_re, "Outputs/bayesian_megatree_model_with_site_random_effect.rds")

#read in models 
#bayes_model <- readRDS("Outputs/bayesian_megatree_model.rds")
bayes_model <- readRDS("Outputs/bayesian_megatree_model_22_10.rds")
#bayes_model <- readRDS("Outputs/bayesian_megatree_model_with_site_random_effect.rds")

#### Posterior Predictions ####
# Create new data for prediction
newdata <- expand.grid(
  habitat = levels(df$habitat),
  height_filt = unique(df$height_filt)
)

# model checks ####

# 1. Overlayed density of predicted vs observed proportions
pp_check(bayes_model, type = "dens_overlay") +
  ggtitle("Posterior Predictive Check: Density Overlay")

# 2. Histogram of test statistics (e.g. mean proportion of megatrees)
pp_check(bayes_model, type = "stat", stat = "mean") +
  ggtitle("PPC: Distribution of Mean Proportions") # looks like prediction are a bit underestimated, but nothing too serious(looking at the values)

# 3. PPC scatterplot: predicted vs observed (binned)
pp_check(bayes_model, type = "scatter_avg") +
  ggtitle("PPC: Average Predicted vs Observed")

# Add fitted draws (posterior predictive mean) for each row
df_fitted <- df %>%
  add_fitted_draws(bayes_model, re_formula = NA) # ignore random effects if you want population-level predictions

df_summary <- df_fitted %>%
  group_by(habitat) %>%
  summarise(
    obs_mean = mean(prop_adj),
    pred_mean = mean(.value),
    pred_lower = quantile(.value, 0.025),
    pred_upper = quantile(.value, 0.975)
  )

df_summary

ggplot(df_summary, aes(x = habitat)) +
  geom_point(aes(y = obs_mean), color = "blue", size = 3) +
  geom_point(aes(y = pred_mean), color = "red", size = 3) +
  geom_errorbar(aes(ymin = pred_lower, ymax = pred_upper), width = 0.2, color = "red") +
  ylab("Mean proportion of big trees") +
  ggtitle("Observed vs predicted by habitat") +
  theme_minimal()


# Get posterior expected values ####
epreds <- posterior_epred(bayes_model, newdata = newdata, re_formula = NA)

# Summarize posterior predictions
pred_summary <- as.data.frame(epreds) %>%
  pivot_longer(cols = everything(), names_to = "row", values_to = "estimate") %>%
  group_by(row) %>%
  summarise(
    mean = mean(estimate),
    lower = quantile(estimate, 0.025),
    upper = quantile(estimate, 0.975)
  ) %>%
  bind_cols(newdata)


#Extract posterior draws for habitat and height filt

# Get full posterior predictions

# Thinning: take every 16th row (8000 / 16 = 500)
thin_interval <- floor(nrow(epreds) / 500)
thinned_draws <- seq(1, nrow(epreds), by = thin_interval)

epreds_thinned <- epreds[thinned_draws, ]  # 500 x 24


epreds_long <- as.data.frame(epreds_thinned) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols = -draw,
    names_to = "obs_id",
    values_to = "estimate"
  ) %>%
  mutate(obs_id = as.integer(gsub("V", "", obs_id))) %>%
  left_join(
    newdata %>% mutate(obs_id = row_number()),
    by = "obs_id"
  ) %>%
  select(draw, habitat, height_filt, estimate)

epreds_long %>%
  count(habitat, height_filt)

#### Plot Bayesian Results ####
#decide if I want to filter a give canopy height only
plot<-pred_summary %>% 
  filter(height_filt %in% c(45, 50)) %>% 
  ggplot( aes(x = habitat, y = mean)) +
  geom_col(fill = "grey60", width = 0.7) +
  
  geom_col(position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, position = position_dodge(0.7)) +
  facet_wrap(~ height_filt, scales = "free_y") +
  labs(
    y = "Proportion of canopy > height threshold",
    x = "Habitat type"
  ) +
  #scale_fill_viridis_d(option = "D", end = 0.85) +
  theme_pubr(base_size = 16) +
  theme(legend.position = "bottom")

plot
#### Save Outputs ####

saveRDS(epreds_long, "Outputs/megatrees_draws_per_hab.rds")

width <-11
height <- 6

ggsave(plot, 
       filename = "Figures/megatrees_by_height.pdf",
       width =  width, #in pixels 
       height = height,
       units = "in")

