# Prepare HBLL survey allocations
# Based on sampling-notes.Rmd

# Notes about the sampling for the HBLL:
# - they are randomly sampled within depth AND spatial strata
# - OUT is spatially stratified by PFMA
# - INS is spatially stratified by N/S (within INS N, and INS S)


# Load settings and utilities
source('R/00-utils.R')

# Full grouping table
gt <- readRDS(here::here("data-raw", "grouping-table.rds")) |>
  select(-strata_depth_label) |>
  mutate(max_depth_m = as.numeric(max_depth_m)) |>
  mutate(strata_depth = case_when(
    is.na(min_depth_m) & is.na(max_depth_m) ~ NA,
    depth_operator == "> MIN_DEPTH_M and <= MAX_DEPTH_M" ~ paste0(min_depth_m + 1, "-", max_depth_m),
    depth_operator == "> MIN_DEPTH_M and < MAX_DEPTH_M" ~ paste0(min_depth_m + 1, "-", max_depth_m - 1),
    depth_operator == ">= MIN_DEPTH_M and < MAX_DEPTH_M" ~ paste0(min_depth_m, "-", max_depth_m - 1),
    depth_operator == ">= MIN_DEPTH_M and <= MAX_DEPTH_M" ~ paste0(min_depth_m, "-", max_depth_m),
    depth_operator == "> MIN_DEPTH_M" ~ paste0(min_depth_m + 1, "-", max_depth_m),
    depth_operator == ">= MIN_DEPTH_M" ~ paste0(min_depth_m, "-", max_depth_m),
    is.na(depth_operator) ~ paste0(min_depth_m, "-", max_depth_m),
    TRUE ~ paste0(min_depth_m, "-", max_depth_m)  # fallback
  ))
saveRDS(gt, file.path("data-generated", "grouping-table-clean.rds"))

# Load strata lookup data
strata0 <- readRDS(here::here("data-raw", "strata-lookup.rds")) |>
  rename(strata_depth = strata_depth_label) # TODO - rename this in the draft gfdata function

# Prepare HBLL strata -----
hbll_strata <- strata0 |> filter(survey_series_id %in% c(22, 36, 39, 40)) |>
  mutate(strata_depth = ifelse(
      grepl("<=", depth_operator),
      paste0(min_depth_m, "-", max_depth_m),
      paste0(min_depth_m, "-", max_depth_m - 1)
    ))

survey_lu <- tibble::tribble(
  ~survey_series_id, ~survey_abbrev,
  22, "HBLL OUT N",
  36, "HBLL OUT S",
  39, "HBLL INS N",
  40, "HBLL INS S"
)

# Create HBLL Inside survey sampling allocation
# Source: Williams and Haggarty 2022
hbll_ins_allocation <- data.frame(
  survey_abbrev = c(rep("HBLL INS N", 4), rep("HBLL INS S", 16)),
  pfma = as.character(c(12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 28, 28, 29, 29)),
  strata_depth = c("40-70 m", "71-100 m", "40-70 m", "71-100 m",
                   "40-70 m", "71-100 m", "40-70 m", "71-100 m", "40-70 m", "71-100 m",
                   "40-70 m", "71-100 m", "40-70 m", "71-100 m", "40-70 m", "71-100 m",
                   "40-70 m", "71-100 m", "40-70 m", "71-100 m"),
  relative_allocation = c(36.92, 16.92, 29.23, 16.92,
                          7.14, 7.14, 11.43, 11.43, 10.00, 8.57, 8.57, 2.86,
                          7.14, 4.29, 1.43, 1.43, 5.71, 4.29, 4.29, 4.29),
  n_blocks_70 = c(26, 12, 20, 12, 5, 5, 8, 8, 7, 6, 6, 2, 5, 3, 1, 1, 4, 3, 3, 3) # match format of survey briefings; target allocation is 70 blocks
)

# Source: HBLL2025_allocation.xlsx
# These allocation files are found on the network drive: GFSurveys/longline_hook/HBLL/inside/YYYY/allocation/
hbll_ins_allocation <- tibble::tribble(
    ~survey_series_id, ~grouping_code, ~pfma, ~strata_depth, ~relative_allocation, ~n_blocks_70, ~survey_abbrev,
    39, 279, "12", "40-70",  36.92, 26, "HBLL INS N",
    39, 280, "12", "71-100", 16.92, 12, "HBLL INS N",
    39, 281, "13", "40-70",  29.23, 20, "HBLL INS N",
    39, 282, "13", "71-100", 16.92, 12, "HBLL INS N",
    40, 283, "14", "40-70",   7.14,  5, "HBLL INS S",
    40, 284, "14", "71-100",  7.14,  5, "HBLL INS S",
    40, 285, "15", "40-70",  11.43,  8, "HBLL INS S",
    40, 286, "15", "71-100", 11.43,  8, "HBLL INS S",
    40, 287, "16", "40-70",  10.00,  7, "HBLL INS S",
    40, 288, "16", "71-100",  8.57,  6, "HBLL INS S",
    40, 289, "17", "40-70",   8.57,  6, "HBLL INS S",
    40, 290, "17", "71-100",  2.86,  2, "HBLL INS S",
    40, 291, "18", "40-70",   7.14,  5, "HBLL INS S",
    40, 292, "18", "71-100",  4.29,  3, "HBLL INS S",
    40, 293, "19", "40-70",   1.43,  1, "HBLL INS S",
    40, 294, "19", "71-100",  1.43,  1, "HBLL INS S",
    40, 295, "28", "40-70",   5.71,  4, "HBLL INS S",
    40, 296, "28", "71-100",  4.29,  3, "HBLL INS S",
    40, 297, "29", "40-70",   4.29,  3, "HBLL INS S",
    40, 298, "29", "71-100",  4.29,  3, "HBLL INS S"
  ) |>
  select(survey_series_id, grouping_code, n_blocks_70) |>
  left_join(hbll_strata, by = c("survey_series_id", "grouping_code"))

hbll_ins_allocation

# Source: HBLL OUT S 2024 briefing document and 2020 RFP: "Appendix_F_SOW_RFP_to_vessels_final_2020_south.pdf"
# These briefing documents are found on the network drive: GFSurveys/longline_hook/HBLL/outside/YYYY/documents/
# depth_stratum = grouping_depth_id (strata table)
# pfma = grouping_spatial_id (strata table)
hbll_out_s_2024 <- tibble::tribble(
    ~pfma, ~depth_stratum, ~depth_m_min, ~depth_m_max, ~n_blocks_174, ~n_blocks_198,
    "3C",  1,              20,           70,            14,            16,
    "3C",  2,              71,           150,           11,            12,
    "3C",  3,              151,          260,           2,             2,
    "3D",  1,              20,           70,            21,            24,
    "3D",  2,              71,           150,           25,            29,
    "3D",  3,              151,          260,           8,             9,
    "4B",  1,              20,           70,            1,             1,
    "4B",  2,              71,           150,           1,             1,
    "4B",  3,              151,          260,           0,             0,
    "5A",  1,              20,           70,            16,            19,
    "5A",  2,              71,           150,           27,            31,
    "5A",  3,              151,          260,           12,            14,
    "5B",  1,              20,           70,            8,             9,
    "5B",  2,              71,           150,           20,            22,
    "5B",  3,              151,          260,           8,             9
  ) |>
  mutate(survey_abbrev = rep("HBLL OUT S", 15), survey_series_id = 36)

hbll_out_n_2025 <- tibble::tribble(
    ~pfma, ~depth_stratum, ~depth_m_min, ~depth_m_max, ~n_active_blocks, ~n_blocks_198, ~n_blocks_174, ~n_blocks_150,
    "5B",  1,              20,           70,            24,               2,             2,             2,
    "5B",  2,              70,           150,           54,               3,             3,             2,
    "5B",  3,              150,          260,           85,               5,             4,             4,
    "5C",  1,              20,           70,            310,              20,            18,            15,
    "5C",  2,              70,           150,           647,              43,            37,            32,
    "5C",  3,              150,          260,           425,              27,            24,            20,
    "5D",  1,              20,           70,            270,              18,            16,            14,
    "5D",  2,              70,           150,           383,              24,            21,            18,
    "5D",  3,              150,          260,           146,              9,             8,             7,
    "5E",  1,              20,           70,            203,              13,            11,            10,
    "5E",  2,              70,           150,           276,              18,            16,            14,
    "5E",  3,              150,          260,           258,              16,            14,            12
  ) |>
  mutate(survey_abbrev = rep("HBLL OUT N", 12), survey_series_id = 22)

hbll_out_allocations <- bind_rows(hbll_out_s_2024, hbll_out_n_2025) |>
  select(survey_series_id, grouping_spatial_id = pfma, grouping_depth_id = depth_stratum, n_blocks_174, n_blocks_198, n_blocks_150) |>
  mutate(grouping_depth_id = as.character(grouping_depth_id)) |>
  left_join(hbll_strata, by = c("survey_series_id", "grouping_spatial_id", "grouping_depth_id"))


hbll_allocations <- bind_rows(hbll_ins_allocation, hbll_out_allocations) |>
  mutate(allocation = ifelse(survey_series_id %in% c(39, 40), n_blocks_70, n_blocks_198)) |>
  left_join(survey_lu, by = "survey_series_id")

# Save the allocations
saveRDS(hbll_allocations, hbll_allocations_file)
