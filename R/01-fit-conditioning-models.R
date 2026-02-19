# =============================================================================
# Fit conditioning models
# =============================================================================

# Load setup and functions
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# devtools::load_all("~/R_DFO/sdmTMB") # need betabinomial branch

library(tidyr)
library(patchwork)
library(digest)

# =============================================================================
# Configuration
# =============================================================================
USE_PARALLEL <- FALSE
N_WORKERS <- NULL

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 40 #NULL
}

if (Sys.info()['user'] == "jilliandunic") {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 8
}

# Setup directories
fit_dir <- here::here("data-generated", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
cleaned_data_dir <- here::here("data-generated", "cleaned-species-data")
dir.create(cleaned_data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa.rds"))

# Fitting parameters
check_cache <- TRUE
silent <- TRUE

# -----------------------------------------------------------------------------
# Run fits
# -----------------------------------------------------------------------------

# Single species (for testing)
# test <- fit_species("yelloweye rockfish", save_cleaned_data = FALSE)

# Species list
sp_list <- c(
  "yelloweye rockfish",
  "north pacific spiny dogfish",
  "lingcod",
  "quillback rockfish",
  "pacific halibut"
)

# # Setup parallel processing
map_fn <- setup_parallel(USE_PARALLEL, N_WORKERS)

# Fit all species
if (USE_PARALLEL) {
  all_fits <- map_fn(sp_list, fit_species,
                     check_cache = check_cache,
                     silent = silent,
                     save_cleaned_data = TRUE,
                     .options = furrr::furrr_options(seed = TRUE))
} else {
  all_fits <- map_fn(sp_list, fit_species,
                     check_cache = check_cache,
                     silent = silent,
                     save_cleaned_data = TRUE)
}

# Reset to sequential
future::plan(future::sequential)
