# =============================================================================
# Shared run configuration for stages 03 (sample) and 04 (fit)
# Edit here; both scripts source this file.
# =============================================================================

# --- Parallelism --------------------------------------------------------------
USE_PARALLEL <- TRUE
N_WORKERS    <- 8L

if (Sys.info()['user'] %in% c("dunic", "anderson")) N_WORKERS <- 50L
if (Sys.info()['user'] == "jilliandunic")           N_WORKERS <- 8L
# Add server username here, e.g.:
# if (Sys.info()['user'] == "server_user")          N_WORKERS <- 40L

# --- Shared filters (used in both 03 and 04) ---------------------------------
FILTER_SPECIES       <- "yelloweye rockfish"
FILTER_SURVEY        <- NULL
FILTER_MPA_TREND     <- 1.009    # 25% recovery; use 1.0164 for 50%
FILTER_AR1_SCENARIO  <- "fitted_AR1"
FILTER_TIME_SCENARIO <- "thirty_years"
FILTER_REPLICATES    <- 1:100

# --- Stage 03 (03-sample-simulated.R) ----------------------------------------
RUN_NON_BOOTSTRAP_PLANS <- TRUE

# --- Stage 04 (04-fit-simulation.R) ------------------------------------------
FILTER_PLAN             <- "status quo"
EVALUATION_YEARS        <- c(2030, 2034, 2038, 2042, 2046)
FILTER_EVALUATION_YEARS <- NULL   # NULL = all; subset e.g. c(2038, 2046) to narrow
RUN_DEFENSIVE_CHECKS <- TRUE
