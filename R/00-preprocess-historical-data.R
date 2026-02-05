# =============================================================================
# Pre-process Historical Data with MPA Spatial Joins
# =============================================================================
# This script performs spatial joins between historical survey data and MPA
# boundaries, which can be slow with old GDAL libraries. Only needed if
# historical data have not been processed and saved yet.

source(here::here("R", "00-setup.R"))

cleaned_data_dir <- here::here("data-generated", "cleaned-species-data")
output_dir <- here::here("data-generated", "historical-data-processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# TODO: will need to update this if we add new species and surveys
# Species to process
species_list <- c(
  "yelloweye rockfish",
  "north pacific spiny dogfish",
  "lingcod",
  "quillback rockfish",
  "pacific halibut"
)

# Surveys to process
surveys <- c("HBLL-OUT-N", "HBLL-OUT-S", "HBLL-INS-N")

# Load MPA spatial data (once)
simple_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa.rds"))

# Create task grid
task_grid <- expand.grid(
  species = species_list,
  survey = surveys,
  stringsAsFactors = FALSE
)

# Process each combination
results <- purrr::pmap(task_grid, function(species, survey) {
  message("Species: ", species)
  message("Survey: ", survey)

  input_file <- file.path(cleaned_data_dir, paste0(sp_to_hyphens(species), "-", survey, ".rds"))
  if (!file.exists(input_file)) {
    message("Input file not found: ", input_file)
    return(NULL)
  }

  output_file <- file.path(output_dir, paste0(sp_to_hyphens(species), "-", survey, "-processed.rds"))
  if (file.exists(output_file)) {
    message("Already processed (delete cache to reprocess)")
    return(output_file)
  }

  hdat0 <- readRDS(input_file)

  hdat_sf <- XY_to_sf(hdat0, crs_to = st_crs(simple_mpa))

  hdat_joined <- st_join(hdat_sf, simple_mpa, join = st_within)

  hdat_processed <- hdat_joined |>
    mutate(
      restricted = ifelse(is.na(uid), 0, 1),
      last_sampled_year = max(year),
      year_covariate = 0,
      historical = TRUE
    ) |>
    st_drop_geometry()
  # Save processed data
  message("  Saving to file: ", output_file)
  saveRDS(hdat_processed, output_file)

  return(output_file)
})
