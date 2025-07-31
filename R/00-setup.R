# Species data cache
if (Sys.info()['user'] == "jilliandunic") synopsis_cache <- "~/R_DFO/gfsynopsis-2024-data/report/data-cache-2025-03"

hbll_ssids <- c(22, 36, 39, 40)
syn_ssids <- c(1, 3, 4, 16)

survey_lu <- tibble::tibble(
  survey_abbrev = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG",
                    "HBLL OUT N", "HBLL OUT S",
                    "HBLL INS N", "HBLL INS S"),
  survey_series_id = c(1, 3, 4, 16, 22, 36, 39, 40)
)

# Something like this could be helpful to add later
# make load data file - sean included this which was smart
# dir.create("data-generated", showWarnings = FALSE)
# dir.create("figs", showWarnings = FALSE)

# if (Sys.info()[["user"]] != "seananderson") {
#   stop("This file does not need to be run; all outputs have been cached and commited in Git.")
# }
