# Species data cache
if (Sys.info()['user'] == "jilliandunic") synopsis_cache <- "~/R_DFO/gfsynopsis-2024-data/report/data-cache-2025-03"
if (Sys.info()['user'] == "seananderson") synopsis_cache <- "../gfsynopsis-2024/report/data-cache-2025-03"

# so that there is a place to put some of the data dependencies
dir.create(here::here("data-generated", "spatial"), recursive = TRUE, showWarnings = FALSE)

hbll_ssids <- c(22, 36, 39, 40)
syn_ssids <- c(1, 3, 4, 16)

survey_lu <- tibble::tibble(
  survey_abbrev = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG",
                    "HBLL OUT N", "HBLL OUT S",
                    "HBLL INS N", "HBLL INS S"),
  survey_series_id = c(1, 3, 4, 16, 22, 36, 39, 40)
)

if (!file.exists(here::here("data-generated", "spatial", "comm-ll-draft-activity-status.rds"))) {
  source(here::here("R", "01-prepare-spatial-data.R"))
} else {
  message("REMINDER to JD - update/resend comm-ll-draft-activity-status.rds once final shapefile available")
}
# Something like this could be helpful to add later
# make load data file - sean included this which was smart
# dir.create("data-generated", showWarnings = FALSE)
# dir.create("figs", showWarnings = FALSE)

# if (Sys.info()[["user"]] != "seananderson") {
#   stop("This file does not need to be run; all outputs have been cached and commited in Git.")
# }

