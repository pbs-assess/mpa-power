# =============================================================================
# Shared run configuration for stages 01–04
# Edit here; all scripts source this file.
# =============================================================================

# --- Run tag (determines ms_dir in 00-setup.R) --------------------------------
# run_tag <- "no-depth"
run_tag <- "ms"       # resdoc outputs
# run_tag <- "0-phi=1000"  # minimize observation error

# Set to NULL for canonical run_tags ("ms", "no-depth") that build everything
# from scratch. Set to an existing run_tag to symlink its Stage 1 outputs
# (mesh, conditioning fits, core geography/depth files) into a new diagnostic
# run_tag instead of regenerating them -- see R/link-run-tag.R.
BASE_RUN_TAG <- "ms"

# --- Parallelism --------------------------------------------------------------
USE_PARALLEL <- TRUE#TRUE
N_WORKERS    <- 8L

if (Sys.info()['user'] %in% c("dunic", "anderson")) N_WORKERS <- 20L
if (Sys.info()['user'] == "jilliandunic")           N_WORKERS <- 8L


# --- Stage 01 (01-fit-conditioning-models.R) ----------------------------------
# Fitting parameters -- these don't work right now
check_cache <- TRUE
silent <- FALSE

# CONDITIONING_FORMULA     <- catch_prop ~ 0 + fyear + restricted
# CONDITIONING_FORMULA_TAG <- "fyear-restricted"
CONDITIONING_FORMULA     <- catch_prop ~ 0 + fyear + restricted + log_depth + I(log_depth^2)
CONDITIONING_FORMULA_TAG <- "fyear-restricted-depth"

ALL_SPECIES <- c(
  "yelloweye rockfish",
  "north pacific spiny dogfish",
  "lingcod",
  "quillback rockfish",
  "pacific halibut",
  "canary rockfish",
  "silvergray rockfish"
)

FIT_SP_LIST <- ALL_SPECIES
# FIT_SP_LIST <- c("lingcod") # Set to ALL_SPECIES to run all species

# --- Stage 02 (02-generate-simulated-data.R) ----------------------------------
# SIM_SP_LIST <- ALL_SPECIES
SIM_SP_LIST <- c("yelloweye rockfish", "lingcod") # Set to ALL_SPECIES to run all species
SIM_NREPS      <- 220L
SIM_REPLICATES <- 1:1
# SIM_TOTAL_INCREASES <- NULL   # NULL = all (1.05, 1.10, 1.25, 1.5); subset e.g. c(1.25) to narrow
SIM_TOTAL_INCREASES <- c(1.25)

SIM_FORMULA_SCENARIOS <- tribble(
  ~formula_scenario, ~formula,
  # "standard", list(~ 1 + log_depth + I(log_depth^2) + restricted * year_covariate)
  "no_depth", list(~ 1 + restricted * year_covariate)
)


# --- Shared filters (used in both 03 and 04) ---------------------------------
FILTER_SPECIES       <- "lingcod"
FILTER_SURVEY        <- NULL
FILTER_MPA_TREND     <- 1.009    # 25% recovery; use 1.0164 for 50%
FILTER_AR1_SCENARIO  <- "fitted_AR1"
FILTER_TIME_SCENARIO <- "thirty_years"
FILTER_REPLICATES    <- 1:1

# --- Stage 03 (03-sample-simulated.R) ----------------------------------------
RUN_NON_BOOTSTRAP_PLANS <- TRUE

# --- Stage 04 (04-fit-simulation.R) ------------------------------------------
# FORMULA <- catch_prop ~ 0 + fyear + restricted + year_covariate +
#   restricted:future_year_covariate +
#   log_depth + I(log_depth^2)
FORMULA <- catch_prop ~ 0 + fyear + restricted + year_covariate +
  restricted:future_year_covariate
TREND_PARAM <- "restricted:future_year_covariate"

FILTER_PLAN             <- "status quo"
EVALUATION_YEARS        <- c(2030, 2034, 2038, 2042, 2046)
FILTER_EVALUATION_YEARS <- NULL   # NULL = all; subset e.g. c(2038, 2046) to narrow
RUN_DEFENSIVE_CHECKS <- FALSE #TRUE # if TRUE

# =============================================================================
# Testing/Debugging options
# =============================================================================
SAVE_TEST_FITS <- TRUE  # Save test fit objects for inspection (if RUN_DEFENSIVE_CHECKS)
TEST_FITS_DIR <- here::here("data-generated", "test-fits")
