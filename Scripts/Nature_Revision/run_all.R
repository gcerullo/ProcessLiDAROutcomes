# Run the full Nature Revision (NR2) workflow.

source(file.path("Scripts", "Nature_Revision", "config.R"))

message("Running NR2 script 01: estimate habitat-level megatree proportions")
source(file.path("Scripts", "Nature_Revision", "01_Nature_calclulate_megatrees_per_hab.R"))

message("Running NR2 script 02: propagate uncertainty through scenarios")
source(file.path("Scripts", "Nature_Revision", "02_Nature_Process_Scenario_Megatree_Outcomes.R"))

message("NR2 core workflow complete.")
message("Optional: run script 03 for model comparison and improved PPC diagnostics.")
