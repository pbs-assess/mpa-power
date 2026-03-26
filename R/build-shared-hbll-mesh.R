# =============================================================================
# Build one shared HBLL mesh from historical data for a single species
# =============================================================================

source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# Configuration
species_name <- "yelloweye rockfish"
survey_keep <- c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")
mesh_cutoff <- 10

# Paths
mesh_dir <- here::here("data-generated", "mesh-cache")
dir.create(mesh_dir, recursive = TRUE, showWarnings = FALSE)

sp <- sp_to_hyphens(species_name)
mesh_file <- file.path(mesh_dir, "HBLL-combined-mesh.rds")

# Load inputs
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
  sf::st_drop_geometry() |>
  dplyr::select(X, Y, uid, restricted)

sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

# Prepare historical HBLL data and reduce to unique coordinates for mesh building
mesh_xy <- sp_dat0 |>
  dplyr::filter(survey_abbrev %in% survey_keep) |>
  prep_hbll_data(bait_counts = bait_counts, restricted_df = historical_locations) |>
  dplyr::distinct(X, Y)

mesh <- make_mesh(mesh_xy, xy_cols = c("X", "Y"), cutoff = mesh_cutoff)

saveRDS(mesh, mesh_file)
