# Prepare HBLL survey allocations
# Based on sampling-notes.Rmd

library(dplyr)

# Load settings and utilities
source('R/00-setup.R')
source('R/00-utils.R')

# Load strata lookup data
strata0 <- readRDS(here::here("data-raw", "strata-lookup.rds")) |>
  rename(strata_depth = strata_depth_label) # TODO - rename this in the draft gfdata function

# Prepare HBLL Inside strata
hbll_ins_strata <- strata0 |> filter(survey_series_id %in% c(39, 40)) |>
  mutate(strata_depth = ifelse(
      grepl("<=", depth_operator),
      paste0(min_depth_m, "-", max_depth_m),
      paste0(min_depth_m, "-", max_depth_m - 1)
    ))

# Create HBLL Inside survey sampling allocation
hbll_ins_allocation <- data.frame(
  Survey = c(rep("North", 4), rep("South", 16)),
  PFMA = as.character(c(12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 28, 28, 29, 29)),
  Depth_Strata = c("40-70 m", "71-100 m", "40-70 m", "71-100 m",
                   "40-70 m", "71-100 m", "40-70 m", "71-100 m", "40-70 m", "71-100 m",
                   "40-70 m", "71-100 m", "40-70 m", "71-100 m", "40-70 m", "71-100 m",
                   "40-70 m", "71-100 m", "40-70 m", "71-100 m"),
  Relative_Allocation = c(36.92, 16.92, 29.23, 16.92,
                          7.14, 7.14, 11.43, 11.43, 10.00, 8.57, 8.57, 2.86,
                          7.14, 4.29, 1.43, 1.43, 5.71, 4.29, 4.29, 4.29),
  Allocation_2021 = c(26, 12, 20, 12, 5, 5, 8, 8, 7, 6, 6, 2, 5, 3, 1, 1, 4, 3, 3, 3)
)

# Prepare HBLL Outside strata
hbll_out_strata <- strata0 |> filter(survey_series_id %in% c(22, 36)) |>
  mutate(strata_depth = ifelse(
      grepl("<=", depth_operator),
      paste0(min_depth_m, "-", max_depth_m),
      paste0(min_depth_m, "-", max_depth_m - 1)
    ))

# Create HBLL Outside survey sampling allocation from Doherty 2019 Table 2.7
hbll_out_allocation <- data.frame(
  Survey = c(rep("North", 12), rep("South", 12)),
  Area = c("5E", "5E", "5E", "5D", "5D", "5D", "5C", "5C", "5C", "5B", "5B", "5B",
           "5B", "5B", "5B", "5A4B", "5A4B", "5A4B", "3D", "3D", "3D", "3C", "3C", "3C"),
  Strata = c(1:12, 1:12),
  Depth_Strata = c("20-70", "71-150", "151-260", "20-70", "71-150", "151-260",
                   "20-70", "71-150", "151-260", "20-70", "71-150", "151-260",
                   "20-70", "71-150", "151-260", "20-70", "71-150", "151-260",
                   "20-70", "71-150", "151-260", "20-70", "71-150", "151-260"),
  Area_km2 = c(808, 1124, 1032, 1104, 1520, 584, 1276, 2604, 1716, 96, 216, 340,
                520, 1296, 496, 1132, 1832, 824, 1376, 1664, 528, 924, 704, 100),
  Historical_Allocation = c(6, 9, 9, 8, 13, 5, 9, 21, 13, 1, 2, 3,
                           4, 11, 5, 9, 16, 9, 11, 14, 6, 7, 6, 1),
  Area_Optimal = c(7, 9, 8, 9, 12, 5, 10, 21, 14, 1, 2, 3,
                    5, 11, 4, 10, 16, 7, 12, 15, 5, 8, 6, 1),
  Optimal_Allocation = c(7, 19, 27, 2, 8, 1, 3, 15, 6, 2, 4, 7,
                         2, 13, 7, 4, 19, 19, 6, 21, 7, 2, 1, 0)
)

# Clean and combine allocations
hbll_ins_clean_allocation <- hbll_ins_allocation |>
  mutate(Depth_Strata = gsub(" m", "", Depth_Strata)) |>
  left_join(hbll_ins_strata, by = c("PFMA" = "grouping_spatial_id", "Depth_Strata" = "strata_depth")) |>
  select(survey = Survey, pfma = PFMA, strata_depth = Depth_Strata,
         relative_allocation = Relative_Allocation, allocation_2021 = Allocation_2021,
         survey_series_id, grouping_code, min_depth_m, max_depth_m)

hbll_out_clean_allocation <- hbll_out_allocation |>
  mutate(Depth_Strata = gsub("260", "259", Depth_Strata)) |> # this matches GFBio's depth strata
  left_join(
    hbll_out_strata |>
      # Keep 5A and 4B separate, but create a grouping column for allocation lookup
      mutate(allocation_group = if_else(grouping_spatial_id %in% c("5A", "4B"), "5A4B", grouping_spatial_id)) |>
      mutate(strata_depth = ifelse(
        grepl("<=", depth_operator),
        paste0(min_depth_m, "-", max_depth_m),
        paste0(min_depth_m, "-", max_depth_m - 1)
      )) |>
      left_join(survey_lu, by = "survey_series_id") |>
      mutate(Survey = case_when(
        survey_abbrev == "HBLL OUT N" ~ "North",
        survey_abbrev == "HBLL OUT S" ~ "South",
        TRUE ~ survey_abbrev
      )),
    by = c("Survey" = "Survey", "Area" = "allocation_group", "Depth_Strata" = "strata_depth")
  ) |>
  select(survey = Survey, pfma = grouping_spatial_id, strata_depth = Depth_Strata,
         area_optimal_allocation = Area_Optimal, historical_allocation = Historical_Allocation,
         survey_series_id, grouping_code, min_depth_m, max_depth_m)

# Combine all HBLL allocations
hbll_allocations <- bind_rows(
  hbll_ins_clean_allocation |> mutate(allocation = allocation_2021),
  hbll_out_clean_allocation |> mutate(allocation = historical_allocation)
) |>
  select(survey_series_id, pfma, grouping_code, strata_depth, allocation) |>
  left_join(survey_lu, by = "survey_series_id")

# Save the allocations
saveRDS(hbll_allocations, file.path("data-generated", "hbll-allocations.rds"))
