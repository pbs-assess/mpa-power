library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

source(here::here("R", "00-overlay-functions.R"))
source(here::here("R", "00-utils.R"))

overlay_dir <- here::here("data-generated", "overlays")
dir.create(overlay_dir, recursive = TRUE, showWarnings = FALSE)

overlay_fig_dir <- here::here("figures", "spatial-overlays")
dir.create(overlay_fig_dir, recursive = TRUE, showWarnings = FALSE)


# Load spatial data for plotting and overlay analysis
simple_coast <- pacea::bc_coast |>
  st_transform(crs = 3005) |>
  st_simplify(dTolerance = 100)

ne_coast <- rnaturalearth::ne_states(country = c("canada", "united states of america"), returnclass = "sf") |>
  st_transform(crs = 3005) |>
  filter(name %in% c("British Columbia", "Washington", "Alaska"))

simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-analytical-mpa.rds")) |>
  mutate(site = gsub("_", " ", common_site_name_site_profile), map_id = map_label)

display_mpa <- simple_mpa |> st_simplify(dTolerance = 100)

subregion_uid_lu <- simple_mpa |>
  st_drop_geometry() |>
  distinct(uid, subregion)

ecp_clean <- readRDS(here::here("data-generated", "ecp-clean.rds"))

n_mpa_sites <- length(unique(simple_mpa$site))
n_mpa_zones <- length(unique(simple_mpa$uid))

# Multi-species dive data
# ------------------------------------------------------------------------------
msd0 <- readRDS(here::here("data-generated", "msd-catch.rds")) |>
  mutate(species_science_name = gsub("nereocystis luetkeanus", "nereocystis luetkeana", species_science_name))
msd <- distinct(msd0, trip_id, transect_site, .keep_all = TRUE)

msd_lines <- mapply(
  create_linestring,
  msd$lon_start, msd$lat_start,
  msd$lon_end, msd$lat_end,
  SIMPLIFY = FALSE
)

msd_sf <- st_sfc(msd_lines, crs = 4326) %>%
  st_sf(msd, geometry = .)

msd_mpa_sf <- msd_sf |>
  st_transform(crs = st_crs(simple_mpa)) |>
  st_join(simple_mpa, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
msd_mpa_df <- st_drop_geometry(msd_mpa_sf)

plot_overlays(msd_mpa_df, msd_mpa_sf |> st_centroid(), title = "MSD dive locations",
              geom_size = 0.2,
              legend_title = "Number of MSD dives",
              breaks = c(-0.1, 0, 5, 20, 50, Inf),
              labels = c("0", "1-5", "6-20", "21-50", "50+"))
ggsave(file.path(overlay_fig_dir, "heatmap-msd-mpa-sf.png"), width = 4.5, height = 5.3)

msd_site_zone_counts <- msd_mpa_df |>
  filter(in_mpa == 1) |>
  reframe(n_sites = n_distinct(site), n_zones = n_distinct(uid)) |>
  mutate(n_sites_text = paste0(n_sites, "/", n_mpa_sites),
         n_zones_text = paste0(n_zones, "/", n_mpa_zones))

msd_mpas_ecps <- msd0 |>
  left_join(msd_mpa_df |> select(in_mpa, uid, site, transect_site, year)) |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  select(in_mpa, site, uid, year, transect_site, species_code, species, species_science_name, is_ecp)

msd_zones_with_ecps <- msd_mpas_ecps |>
  filter(is_ecp, in_mpa == 1) |>
  # distinct(site, uid, is_ecp) |>
  distinct(site, uid, .keep_all = TRUE) |>
  summarise(n_zones_with_ecps = n()) |>
  mutate(n_zones_ecp_text = paste0(n_zones_with_ecps, "/", n_mpa_zones))

msd_sites_with_ecps <- msd_mpas_ecps |>
  filter(is_ecp, in_mpa == 1) |>
  distinct(site, is_ecp) |>
  summarise(n_sites = n()) |>
  mutate(n_sites_ecp_text = paste0(n_sites, "/", n_mpa_sites))

msd_longest_timeseries_site <- msd_mpa_df |>
  filter(in_mpa == 1) |>
  group_by(site) |>
  summarise(longest_timeseries = max(year) - min(year),
            n_years = length(unique(year)),
            .groups = "drop") |>
  arrange(desc(longest_timeseries)) |>
  slice(1) |>
  mutate(text = paste0(n_years, " years (", site, ")"))

msd_summary_table <- msd_mpa_df |>
  mutate(survey_name = "Multispecies Benthic Invertebrate Dive Survey") |>
  reframe(
    `Survey name` = "Multispecies Benthic Invertebrate Dive Survey",
    `Years evaluated here` = paste(min(year), max(year), sep = "-"),
    `Total transects` = dplyr::n_distinct(transect_site),
    `years_not_in_mpa` = paste(setdiff(min(year):max(year), unique(year[in_mpa == 1])), collapse = ", "),
    `Years with surveys inside MPA boundaries` =
      paste0(
        n_distinct(year[in_mpa == 1]), " of ", max(year) - min(year) + 1,
        " (no surveys inside MPAs in ", `years_not_in_mpa`, ")"
      ),
    `Transects inside MPA boundaries` = dplyr::n_distinct(transect_site[in_mpa == 1]),
    `Transects outside MPA boundaries` = dplyr::n_distinct(transect_site[in_mpa == 0]),
    `% of effort inside MPA boundaries` = round(100 * `Transects inside MPA boundaries` / `Total transects`, 1),
    n_sites_sampled = msd_site_zone_counts$n_sites,
    `Number of MPAN sites sampled` = msd_site_zone_counts$n_sites_text,
    # `Number of MPAN zones sampled` = msd_site_zone_counts$n_zones_text,
    `Number of MPAN sites with E-CP sightings` = msd_sites_with_ecps$n_sites_ecp_text,
    # `Number of MPAN zones with E-CP sightings` = msd_zones_with_ecps$n_zones_ecp_text,
    `Longest timeseries at any single site` = msd_longest_timeseries_site$text
  ) |>
  mutate(across(everything(), as.character)) |>
  tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value")
saveRDS(msd_summary_table, file.path(overlay_dir, "msd-summary-table.rds"))


# MSD site summary table --------------------
# Site ID
# Common Site Name
# Total transects
# First survey
# Last survey
# No. years
# Unique species
# Mean species richness
msd_species_list <- #msd_raw |>
  msd0 |>
  left_join(msd_mpa_df |> select(in_mpa, uid, site, transect_site, year), by = c("transect_site", "year")) |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  filter(in_mpa == 1) |>  # Only keep observations inside MPAs
  group_by(site) |>
  summarise(
    `Unique species` = n_distinct(species),
    unique_species = paste(unique(species), collapse = ", "),
    .groups = "drop"
  )

# msd_raw |>
msd0 |>
  # select(algae_code:species) |>
  filter(stringr::str_detect(species_science_name, "zostera")) |>
  distinct() |>
  glimpse()

ecp_clean |>
  filter(stringr::str_detect(scientific_name, "zostera")) |>
  distinct() |>
  glimpse()

msd_site_richness <- readRDS(here::here("data-generated", "msd-transect-richness.rds")) |>
  left_join(msd_mpa_df |> select(in_mpa, uid, site, transect_site, year)) |>
  filter(in_mpa == 1) |>  # Only transects inside MPAs
  # group_by(site, uid) |>
  group_by(site) |>
  summarise(`Mean species richness` = mean(richness), .groups = "drop")

msd_site_summary_table <- msd_mpa_df |>
  filter(in_mpa == 1) |>
  # group_by(site, uid) |>
  group_by(site) |>
  summarise(
    # `Site ID` = unique(uid),
    `Common Site Name` = unique(site),
    `Total transects` = length(unique(transect_site)),
    `First survey` = min(year),
    `Last survey` = max(year),
    `No. years` = n_distinct(year),
    `Years in MPA` = paste(unique(year[in_mpa == 1]), collapse = ", "),
    .groups = "drop"
  ) |>
  left_join(msd_site_richness, by = c("site")) |>
  left_join(msd_species_list, by = c("site")) |>
  arrange(desc(`Total transects`))

saveRDS(msd_site_summary_table, file.path(overlay_dir, "msd-site-summary-table.rds"))

# MSD site summary table --------------------
msd_species_summary_table <-
  msd0 |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  distinct(year, transect_site, species_science_name, species_common_name, species_desc, is_ecp, species_group, catch_count,) |>
  arrange(is_ecp) |>
  left_join(msd_mpa_df |> select(in_mpa, uid, site, transect_site, year), by = c("transect_site", "year")) |>
  # mutate(in_mpa = ifelse(is.na(uid), 0, 1)) |>  # Only keep observations inside MPAs
  mutate(Species = stringr::str_to_sentence(species_science_name)) |>
  filter(!is.na(is_ecp) | stringr::str_detect(species_science_name, "pterygophora")) |>  # Only keep E-CP species
  group_by(Species, species_group, species_desc) |>
  summarise(
    `Sightings inside MPAs` = sum(catch_count[in_mpa == 1]),
    `Sightings outside MPAs` = sum(catch_count[in_mpa == 0]),
    `Proportion inside MPAs` = round(`Sightings inside MPAs` / (`Sightings inside MPAs` + `Sightings outside MPAs`), 2),
    .groups = "drop"
  ) |>
  select(Species, `Common name` = species_desc, `Sightings inside MPAs`, `Sightings outside MPAs`, `Proportion inside MPAs`, species_group) |>
  arrange(species_group)
saveRDS(msd_species_summary_table, file.path(overlay_dir, "msd-species-summary-table.rds"))

# Site-level E-CP sightings (for narrative: "most productive sites")
msd_site_ecp_sightings <-
  msd0 |>
  left_join(msd_mpa_df |> select(in_mpa, uid, site, transect_site, year), by = c("transect_site", "year")) |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  filter(is_ecp == TRUE, in_mpa == 1) |>
  group_by(site) |>
  summarise(total_ecp_sightings = sum(catch_count), .groups = "drop") |>
  arrange(desc(total_ecp_sightings))
saveRDS(msd_site_ecp_sightings, file.path(overlay_dir, "msd-site-ecp-sightings.rds"))
