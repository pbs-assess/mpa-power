library(dplyr)
library(ggplot2)


# Survey data depths -----------------------------------------------------------

sp_dat_depths0 <- readRDS(here::here("data-generated", "spatial", "hbll-survey-dem-depths.rds"))

sp_dat_depths <-sp_dat_depths0 |>
  group_by(survey_abbrev, fishing_event_id) |>
  mutate(dem_depth = mean(c(depth_start, depth_end), na.rm = TRUE))

filter(sp_dat_depths, is.na(dem_depth))

# Compare ship-based depth measurements with DEM extracted depths
sp_dat_depths |>
  group_by(survey_abbrev, fishing_event_id) |>
  mutate(dem_depth = mean(c(depth_start, depth_end), na.rm = TRUE)) |>
  ggplot() +
  aes(x = dem_depth, y = depth_m) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_rect(data = strata_rects, aes(xmin = depth_min, xmax = depth_max, ymin = -Inf, ymax = Inf, fill = strata_depth_label), alpha = 0.3, colour = NA, inherit.aes = FALSE) +
  # geom_rect(data = strata_rects, aes(ymin = depth_min, ymax = depth_max, xmin = -Inf, xmax = Inf, fill = strata_depth_label), alpha = 0.3, colour = NA, inherit.aes = FALSE) +
  theme_minimal() +
  facet_wrap(~ survey_abbrev, ncol = 2, scales = "free", dir = "v")

sp_dat_depths |>
  group_by(survey_abbrev, fishing_event_id) |>
  mutate(dem_depth = mean(c(depth_start, depth_end), na.rm = TRUE)) |>
  ggplot() +
  aes(x = dem_depth, y = depth_m) +
  geom_rect(data = survey_shading,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "grey90", inherit.aes = FALSE) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_vline(data = strata_lines, aes(xintercept = boundary),
    linetype = "dashed", colour = "grey30") +
  facet_wrap(~ survey_abbrev, ncol = 2, scales = "free", dir = "v") +
  labs(x = "DEM depth — mean of start/end (m)", y = "Ship depth (m)")

depth_summaries_sf <- hbll_grid_poly_transformed |>
  left_join(depth_summaries, by = c("block_id", "survey_series_id")) |>
  arrange(survey_abbrev, min_depth_m) |>
  mutate(strata_depth_label = factor(strata_depth_label, levels = unique(strata_depth_label)))


# DEM data visualisations ------------------------------------------------------
depth_summaries <- readRDS(here::here("data-generated", "spatial", "hbll-depth-summaries.rds"))
depth_summaries_sf <- hbll_grid_poly_transformed |>
  left_join(depth_summaries, by = c("block_id", "survey_series_id")) |>
  arrange(survey_abbrev, min_depth_m) |>
  mutate(strata_depth_label = factor(strata_depth_label, levels = unique(strata_depth_label)))

# HBLL grid polygons ----
hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S", "HBLL INS N", "HBLL INS S"))
# Transform HBLL grid to match DEM CRS
hbll_grid_poly_transformed <- hbll_grid_poly |>
  st_transform(crs = terra::crs(dem)) |>
  left_join(strata, by = "grouping_code")

strata_rects <- distinct(hbll_grid_poly_transformed, survey_abbrev, survey_series_id, strata_depth_label) |>
  tidyr::separate(strata_depth_label, into = c("depth_min", "depth_max"),
    sep = " - ", convert = TRUE, remove = FALSE)

strata_ranges <- depth_summaries_sf |>
  st_drop_geometry() |>
  distinct(survey_abbrev, strata_depth_label, min_depth_m, max_depth_m) |>
  filter(!is.na(min_depth_m), !is.na(max_depth_m)) |>
  rename(dem_stratum = strata_depth_label, strata_min = min_depth_m, strata_max = max_depth_m) |>
  arrange(survey_abbrev, strata_min) |>
  mutate(dem_stratum = factor(dem_stratum, levels = unique(dem_stratum)))


# Comparison of the DEM and OG grid depths
depth_summaries_sf |>
  st_drop_geometry() |>
  select(survey_abbrev, strata_depth_label, depth_m, depth_dem_mean) |>
  tidyr::pivot_longer(c(depth_m, depth_dem_mean),
    names_to = "source", values_to = "depth") |>
  mutate(source = recode(source,
    depth_m = "Original grid", depth_dem_mean = "DEM mean")) |>
  ggplot(aes(x = strata_depth_label, y = depth, fill = source)) +
  geom_boxplot(outlier.alpha = 0.2) +
  # geom_violin(alpha = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(data = strata_ranges, aes(yintercept = strata_min), linetype = "dashed") +
  geom_hline(data = strata_ranges, aes(yintercept = strata_max), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey80") +
  scale_fill_viridis_d(option = "turbo", begin = 0.15, end = 0.7) +
  facet_wrap(~ survey_abbrev, scales = "free") +
  labs(x = "Depth stratum", y = "Depth (m)", fill = NULL)


# Histogram of depth summaries - with target depth strata highlighted
survey_shading <- strata_ranges |>
  group_by(survey_abbrev) |>
  summarise(lo = min(strata_min), hi = max(strata_max), .groups = "drop") |>
  tidyr::pivot_longer(c(lo, hi)) |>
  mutate(xmin = if_else(name == "lo", -Inf, value),
         xmax = if_else(name == "lo", value, Inf))

strata_lines <- strata_ranges |>
  tidyr::pivot_longer(c(strata_min, strata_max), values_to = "boundary") |>
  distinct(survey_abbrev, boundary)

depth_summaries_sf |>
    st_drop_geometry() |>
    select(survey_abbrev, depth_dem_mean, depth_m) |>
    tidyr::pivot_longer(c(depth_dem_mean, depth_m),
      names_to = "source", values_to = "depth") |>
    mutate(source = recode(source,
      depth_dem_mean = "DEM mean", depth_m = "Original grid")) |>
    ggplot(aes(x = depth)) +
    geom_rect(data = survey_shading,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "grey90", inherit.aes = FALSE) +
    geom_histogram(bins = 40, fill = "grey50") +
    geom_vline(data = strata_lines, aes(xintercept = boundary),
      linetype = "dashed", colour = "grey30") +
    facet_grid(source ~ survey_abbrev, scales = "free_x") +
    labs(x = "Depth (m)", y = "Count")

# Compare difference between DEM and original grid depths
depth_summaries_sf |>
  st_drop_geometry() |>
  mutate(depth_diff = depth_dem_mean - depth_m) |>
  ggplot(aes(x = depth_diff)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "tomato") +
  facet_wrap(~ survey_abbrev, scales = "free_y", ncol = 2) +
  labs(x = "DEM mean − original depth_m (m)", y = "Count")

# What proportion of blocks match the original stratum when using the DEM depth mean
# Not sure if either of these two plots are at all helpful; the boxplot probably shows this the best.
survey_bounds_lu <- strata_ranges |>
  group_by(survey_abbrev) |>
  summarise(lo = min(strata_min), hi = max(strata_max), .groups = "drop")

unassigned <- depth_summaries_sf |>
  st_drop_geometry() |>
  select(survey_abbrev, block_id, strata_depth_label, depth_dem_mean) |>
  anti_join(dem_stratum_lu, by = c("survey_abbrev", "block_id")) |>
  left_join(survey_bounds_lu, by = "survey_abbrev") |>
  mutate(dem_stratum = case_when(
    is.na(depth_dem_mean) ~ "No DEM data",
    depth_dem_mean < lo   ~ "Outside lower",
    TRUE                  ~ "Outside upper"
  )) |>
  select(survey_abbrev, block_id, strata_depth_label, dem_stratum)

dem_strata_levels <- c("No DEM data", "Outside lower", levels(strata_ranges$dem_stratum), "Outside upper")

# Option 1
bind_rows(dem_stratum_lu, unassigned) |>
  mutate(dem_stratum = factor(dem_stratum, levels = dem_strata_levels)) |>
  count(survey_abbrev, strata_depth_label, dem_stratum) |>
  ggplot(aes(x = strata_depth_label, y = dem_stratum, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), colour = "white", size = 3) +
  scale_fill_viridis_c() +
  facet_wrap(~ survey_abbrev, scales = "free") +
  labs(x = "Original stratum", y = "DEM-implied stratum", fill = "n blocks") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Option 2
bind_rows(dem_stratum_lu, unassigned) |>
    mutate(dem_stratum = factor(dem_stratum, levels = dem_strata_levels)) |>
    count(survey_abbrev, strata_depth_label, dem_stratum) |>
    group_by(survey_abbrev, strata_depth_label) |>
    mutate(prop = n / sum(n)) |>
    ggplot(aes(x = strata_depth_label, y = prop, fill = dem_stratum)) +
    geom_col() +
    scale_fill_manual(
      values = c(
        "Outside lower" = "#d73027",
        "40 - 71"       = "#4575b4",
        "71 - 100"      = "#74add1",
        "20 - 71"       = "#4575b4",
        "71 - 151"      = "#74add1",
        "151 - 260"     = "#abd9e9",
        "Outside upper" = "#f46d43",
        "No DEM data"   = "grey70"
      )
    ) +
    facet_wrap(~ survey_abbrev, scales = "free_x", ncol = 2) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Original stratum", y = "Proportion of blocks", fill = "DEM-implied stratum")

# Option 3. Proportion of blocks where the DEM depth matches the assigned stratum
bind_rows(dem_stratum_lu, unassigned) |>
  mutate(in_stratum = dem_stratum == strata_depth_label) |>
  count(survey_abbrev, strata_depth_label, in_stratum) |>
  group_by(survey_abbrev, strata_depth_label) |>
  mutate(prop = n / sum(n)) |>
  filter(in_stratum) |>
  ggplot(aes(x = strata_depth_label, y = prop)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "tomato") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ survey_abbrev, scales = "free_x", ncol = 2) +
  labs(x = "Original stratum", y = "% blocks with DEM depth in designated stratum")

# ------
tmap_mode("view")
  tm_shape(depth_summaries_sf) +
    tm_polygons(
      fill = "strata_depth_label",
      fill_alpha = 0.1,
      popup.vars = c("block_id", "strata_depth_label", "depth_dem_mean", "depth_dem_sd", "n_pixels", "depth_dem_min", "depth_dem_max")
    )

ssid_labels <- c("HBLL OUT N" = 22, "HBLL OUT S" = 36, "HBLL INS N" = 39, "HBLL INS S" = 40)

depth_summaries_sf |>
  st_drop_geometry() |>
  mutate(dem_minus_og = depth_dem_mean - depth_m) |>
  ggplot() +
  aes(depth_m, depth_dem_mean, colour = factor(survey_series_id, levels = ssid_labels, labels = names(ssid_labels))) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_brewer(palette = "Dark2") +
  # scale_y_continuous(trans = ggsidekick::fourth_root_power_trans()) +
  labs(x = "Grid depth_m", y = "DEM mean depth", colour = "Survey") +
  facet_wrap(~ survey_series_id, ncol = 2)

depth_summaries_sf |>
  st_drop_geometry() |>
  mutate(
    in_stratum = case_when(
      is.na(min_depth_m) | is.na(max_depth_m) ~ NA,
      stringr::str_detect(depth_operator, ">= MIN_DEPTH_M and < MAX_DEPTH_M")  ~ depth_dem_mean >= min_depth_m & depth_dem_mean < max_depth_m,
      stringr::str_detect(depth_operator, ">= MIN_DEPTH_M and <= MAX_DEPTH_M") ~ depth_dem_mean >= min_depth_m & depth_dem_mean <= max_depth_m,
      stringr::str_detect(depth_operator, "> MIN_DEPTH_M and < MAX_DEPTH_M")   ~ depth_dem_mean > min_depth_m  & depth_dem_mean < max_depth_m,
      stringr::str_detect(depth_operator, "> MIN_DEPTH_M and <= MAX_DEPTH_M")  ~ depth_dem_mean > min_depth_m  & depth_dem_mean <= max_depth_m,
      stringr::str_detect(depth_operator, ">= MIN_DEPTH_M") ~ depth_dem_mean >= min_depth_m,
      stringr::str_detect(depth_operator, "> MIN_DEPTH_M")  ~ depth_dem_mean > min_depth_m,
      TRUE ~ NA
    )
  ) |>
  ggplot() +
  aes(depth_m, depth_dem_mean, colour = in_stratum) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
    na.value = "grey60",
    labels = c("TRUE" = "In stratum", "FALSE" = "Out of stratum", "NA" = "No bounds")
  ) +
  # scale_x_continuous(trans = ggsidekick::fourth_root_power_trans()) +
  # scale_y_continuous(trans = ggsidekick::fourth_root_power_trans()) +
  labs(x = "Grid depth_m", y = "DEM mean depth", colour = NULL) +
  facet_wrap(~ survey_series_id, ncol = 2)
