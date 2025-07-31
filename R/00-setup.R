# Species data cache
if (Sys.info()['user'] == "jilliandunic") synopsis_cache <- "~/R_DFO/gfsynopsis-2024-data/report/data-cache-2025-03"

hbll_ssids <- c(22, 36, 39, 40)
syn_ssids <- c(1, 3, 4, 16)

survey_lu <- tibble(
  survey_abbrev = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG",
                    "HBLL OUT N", "HBLL OUT S",
                    "HBLL INS N", "HBLL INS S"),
  survey_series_id = c(1, 3, 4, 16, 22, 36, 39, 40)
)
