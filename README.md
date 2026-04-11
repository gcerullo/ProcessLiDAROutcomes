# ProcessLiDAROutcomes

Code for estimating proportional megatree coverage from LiDAR and propagating uncertainty through scenario outcomes.

## Legacy scripts

### `Scripts/01_CalculateMegatreesPerHab.R`
Calculates habitat-specific canopy proportions above multiple height thresholds (45-70 m).

### `Scripts/02_ProcessScenarioMegatreeOutcomes.R`
Propagates habitat-level megatree estimates through scenario trajectories.

## Nature Revision workflow (NR2)

All NR2 scripts are in `Scripts/Nature_Revision`.

### `Scripts/Nature_Revision/config.R`
- Defines shared paths for NR2.
- Creates:
  - `Outputs/NR2`
  - `Outputs/NR2/figures`
  - `Outputs/NR2/models`

### `Scripts/Nature_Revision/01_Nature_calclulate_megatrees_per_hab.R`
- Reads LiDAR canopy rasters and 1 ha hex grid.
- Computes megatree proportions by habitat and height threshold.
- Fits the baseline Bayesian beta model.
- Produces habitat-level posterior predictive diagnostics.

Inputs:
- `RawData/TreeHeightDanum/Sabah100HexShapefile.shp`
- `RawData/TreeHeightDanum/Primary.tif`
- `RawData/TreeHeightDanum/Once_Logged.tif`
- `RawData/TreeHeightDanum/Restored.tif`
- `RawData/TreeHeightDanum/Twice_Logged.tif`
- Cached extracted NR2 inputs (if skipping extraction):
  - `Outputs/NR2/primary1haCHMvals.rds`
  - `Outputs/NR2/onceLogged_vals1haCHMvals.rds`
  - `Outputs/NR2/restored_vals1haCHMvals.rds`
  - `Outputs/NR2/twiceLogged_vals1haCHMvals.rds`

Outputs:
- `Outputs/NR2/models/bayesian_megatree_model_22_10.rds`
- `Outputs/NR2/megatrees_draws_per_hab.rds`
- `Outputs/NR2/model_fit_by_habitat_summary.csv`
- `Outputs/NR2/model_fit_by_habitat_summary.rds`
- `Outputs/NR2/figures/posterior_pp_checks_by_habitat.pdf`

### `Scripts/Nature_Revision/02_Nature_Process_Scenario_Megatree_Outcomes.R`
- Joins posterior megatree draws onto scenario trajectories.
- Summarises 60-year megatree performance and uncertainty.

Inputs:
- `Inputs/FixedScenarioParmams.R`
- `Inputs/MasterAllScenarios.rds`
- `Inputs/allHabCarbon_60yr_withDelays.csv`
- `Outputs/NR2/megatrees_draws_per_hab.rds`

Outputs:
- `Outputs/NR2/MasterMegatreePerformance_with_uncertainty.rds`

### `Scripts/Nature_Revision/03_Nature_Improved_Megatree_Model.R`
- Fits improved model candidates:
  - `beta_base`
  - `beta_phi_re`
  - `beta_binomial_re`
- Compares models using LOO and saves the best model.
- Produces habitat-level PPC diagnostics for the selected model.

Outputs:
- `Outputs/NR2/models/nr2_model_beta_base.rds`
- `Outputs/NR2/models/nr2_model_beta_phi_re.rds`
- `Outputs/NR2/models/nr2_model_beta_binomial_re.rds`
- `Outputs/NR2/models/nr2_best_model.rds`
- `Outputs/NR2/model_comparison_loo.csv`
- `Outputs/NR2/model_comparison_loo.rds`
- `Outputs/NR2/improved_model_fit_by_habitat_summary.csv`
- `Outputs/NR2/improved_model_fit_by_habitat_summary.rds`
- `Outputs/NR2/figures/improved_posterior_pp_checks_by_habitat.pdf`

### `Scripts/Nature_Revision/run_all.R`
Runs the core NR2 pipeline:
1. `01_Nature_calclulate_megatrees_per_hab.R`
2. `02_Nature_Process_Scenario_Megatree_Outcomes.R`

Run from project root:
- `source("Scripts/Nature_Revision/run_all.R")`
