# Run the full Nature Revision (NR2) workflow.

source(file.path("Scripts", "Nature_Revision", "config.R"))

message("Running NR2 script 01: fit megatree models and export posterior outputs")
source(file.path("Scripts", "Nature_Revision", "01_full_NAture_megatrees_model.R"))

message("Running NR2 script 02: propagate uncertainty through scenarios")
source(file.path("Scripts", "Nature_Revision", "02_full_Nature_process_scenarios.R"))

message("NR2 core workflow complete.")
message("Archived legacy scripts are now in Scripts/Nature_Revision/old.")
