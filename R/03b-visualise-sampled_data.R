sample_dir <- here::here("data-generated", "sampled-data")
ye_samps <- readRDS(file.path(sample_dir, "yelloweye-rockfish-all-sampled.rds"))

# stash historical data here for now
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets
sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
  prep_hbll_data(bait_counts = bait_counts) |>
  mutate(
    obs_id = factor(row_number()),
    catch_prop = catch_count / hook_count,
    log_depth = log(depth_m)
  )
# Prepare historical data for comparison and future modelling
historical <- sp_dat |>
  XY_to_sf(crs_to = 32609) |>
  st_join(simple_mpa |> st_transform(crs = 32609), join = st_within) |>
  mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
  st_join(hbll_grid_poly |> select(block_id, grouping_code), join = st_within) |>
  st_drop_geometry() |>
  select(ssid, survey_abbrev, year, fishing_event_id, latitude, longitude, X, Y,
    block_id, fe_grouping_code = grouping_code.x, grouping_code = grouping_code.y,
    depth_m, offset, hook_count,
    catch_count, restricted) |>
  mutate(after = 0) |>
  left_join(hbll_allocations,
    by = c("survey_abbrev", "grouping_code", "ssid" = "survey_series_id")) |>
  mutate(observed = catch_count, replicate = 0, d = "historical") |>
  group_by(survey_abbrev) |>
  mutate(year_counter = year - min(year) + 1) |>
  ungroup()