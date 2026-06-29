source(here::here("R", "00-setup.R"))

sample_files <- list.files(file.path(sample_dir, "yelloweye-rockfish"), full.names = TRUE)
ye_samps <- purrr::map_dfr(sample_files, readRDS)

# Double check that allocations look correct
ye_samps |> filter(replicate == 1) |>
  filter(sim_mpa_trend == 1.009) |>
  filter(year %in% 2025:2026) |>
  group_by(survey_abbrev, plan) |>
  summarise(n = n(), .groups = "drop")


ye_samps |> filter(replicate == 1) |>
  filter(sim_mpa_trend == 1.009) |>
  filter(year %in% 2025:2026) |>
  group_by(survey_abbrev, plan, restricted) |>
  summarise(n = n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = restricted, values_from = n, values_fill = 0) |>
  mutate(prop_inside = `1` / (`0` + `1`))

distinct(ye_samps, plan)

test <- ye_samps |>
  filter(replicate == 1) |>
  # filter(plan == "status quo") |>
  # visual check to make sure historical bootstrapping looks like allocation-based sampling
  filter(plan %in% c("historical survey-year bootstrap - no MPA every 2nd survey", "MPAs every 4 years")) |> # uses allocations - NOT bootstrapped
  filter(sim_mpa_trend == 1.009) |>
  filter(year %in% 2025:2030) |>
  XY_to_sf() |>
  rotate_sf()
ggplot() +
  geom_sf(data = pacea::bc_coast |> rotate_sf(), fill = "grey94", colour = "grey90") +
  geom_sf(data = display_mpa |> rotate_sf(), fill = "#0072B2", colour = NA, alpha = 0.3) +
  geom_sf(data = test, aes(shape = factor(restricted))) +
  scale_shape_manual(name = "Restricted", values = c(`0` = 21, `1` = 19)) +
  facet_grid(plan ~ year) +
  theme(axis.text = element_blank()) +
  gfplot::coord_sf_auto(test)
