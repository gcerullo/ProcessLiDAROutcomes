# Nature Revision (NR2) configuration

nr2_config <- list(
  raw_data_dir = file.path("RawData", "TreeHeightDanum"),
  input_dir = "Inputs",
  output_dir = file.path("Outputs", "NR2"),
  figures_dir = file.path("Outputs", "NR2", "figures"),
  models_dir = file.path("Outputs", "NR2", "models")
)

dir.create(nr2_config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nr2_config$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nr2_config$models_dir, recursive = TRUE, showWarnings = FALSE)

nr2_input_path <- function(...) {
  file.path(nr2_config$input_dir, ...)
}

nr2_raw_path <- function(...) {
  file.path(nr2_config$raw_data_dir, ...)
}

nr2_output_path <- function(...) {
  file.path(nr2_config$output_dir, ...)
}
