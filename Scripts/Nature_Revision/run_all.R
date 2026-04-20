# Run the full Nature Revision (NR2) workflow.

source(file.path("Scripts", "Nature_Revision", "config.R"))

message("Running NR2 script 01: fit megatree models and export posterior outputs")
source(file.path("Scripts", "Nature_Revision", "01_full_NAture_megatrees_model.R"))

message("Running NR2 script 02: propagate uncertainty through scenarios")
source(file.path("Scripts", "Nature_Revision", "02_full_Nature_process_scenarios.R"))

sbe_dir <- file.path("Scripts", "Nature_Revision", "SBE_sensitivity_analysis")
message("Running Sabah (SBE) LiDAR pipeline: prep polygons/rasters (01)")
source(file.path(sbe_dir, "01_prep_polygons_and_rasters.R"))
message("Running Sabah (SBE) LiDAR pipeline: extract plot-level CHM metrics (02)")
source(file.path(sbe_dir, "02_determine_lidar_change_2013_2020.R"))
message("Running Sabah (SBE) LiDAR sensitivity: implied carbon from mean height (03)")
source(file.path(sbe_dir, "03_lidar_SBE_carbon_validation.R"))
message("Running Sabah (SBE) LiDAR sensitivity: megatree pixel share (45/50/55 m) (04)")
source(file.path(sbe_dir, "04_lidar_SBE_megatree_sensitivity.R"))

message("NR2 workflow complete (megatree models, scenarios, and SBE LiDAR sensitivity).")
message("Optional next step: source Scripts/Nature_Revision/03_full_Nature_predictive_scenarios.R for predictive scenario aggregation.")
message("Archived legacy scripts are in Scripts/Nature_Revision/old (not used by run_all.R).")
