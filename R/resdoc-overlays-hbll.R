library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

source(here::here("R", "00-overlay-functions.R"))
source(here::here("R", "00-utils.R"))
source(here::here("R", "00-setup.R"))

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

ecp_clean0 <- readRDS(here::here("data-generated", "ecp-clean.rds"))

fish_groups <- c(
    "pleuronectidae(righteye flounders)",
    "osmeridae(smelts)",
    "clupeidae(herrings)",
    "ammodytidae(sand lances)",
    "hexagrammidae(greenlings)",
    "anoplopomatidae(sablefishes)",
    "anarhichadidae(wolffishes)",
    "myctophidae(lanternfishes)",
    "bathylagidae(deepsea smelts)",
    "oncorhynchus(pacific salmon and native trout)",
    "salmoninae(pacific salmon and native trout)",
    "scombridae(mackerels and tunas)",
    "molidae(molas)",
    "sebastes",                      # rockfish genus
    "sebastolobus(thornyheads)",
    "gadidae",                       # cods
    "merlucciidae(merluccid hakes)",
    "acipenseridae",                 # sturgeons
    "embiotocidae(surfperches)",
    "hexanchidae(cow sharks)",
    "somniosidae",                   # sleeper sharks
    "squalidae(dogfish sharks)",
    "cetorhinidae(basking shark)",
    "carcharhinidae(requiem sharks)",
    "lamnidae(mackerel sharks)",
    "rajidae(skates)",
    "bathyraja"                      # skate genus
  )

ecp_clean <- ecp_clean0 |> mutate(is_fish = parent_taxonomic_unit %in% fish_groups)

n_fish_ecps <- sum(ecp_clean$is_fish)
n_non_fish_ecps <- sum(!ecp_clean$is_fish)

n_mpa_sites <- length(unique(simple_mpa$site))
n_mpa_zones <- length(unique(simple_mpa$uid))

# HBLL data
# ------------------------------------------------------------------------------
# All groundfish data
sp <- sp_to_hyphens("yelloweye rockfish")
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

sp_dat <- sp_dat0 |>
  filter(survey_series_id.x %in% c(hbll_ssids, syn_ssids)) |>
  mutate(fo_id = case_when(
    survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S") ~ "StARGF_01",
    survey_abbrev %in% c("HBLL INS N") ~ "StARGF_03",
    survey_abbrev %in% c("SYN HS") ~ "StARGF_02",
    survey_abbrev %in% c("SYN WCVI") ~ "StARGF_04",
    survey_abbrev %in% c("SYN QCS") ~ "StARGF_07",
    survey_abbrev %in% c("SYN WCHG") ~ "StARGF_08"
  )) |>
  select(fo_id, ssid = survey_series_id.x, fe = fishing_event_id, trip_id,
    survey_abbrev, year, species = species_common_name, doorspread_m,
    latitude_start = latitude, longitude_start = longitude,
    latitude_end, longitude_end) |>
  filter(survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S", "HBLL INS N"))

gf_line_geoms <- mapply(
  create_linestring,
  sp_dat$longitude_start, sp_dat$latitude_start,
  sp_dat$longitude_end, sp_dat$latitude_end,
  SIMPLIFY = FALSE
)

gf_line_sf <- st_sfc(gf_line_geoms, crs = 4326) %>%
  st_sf(sp_dat, geometry = .)

# HBLL transects
# -----------------------------------------------------------------------------
hbll_transects <- gf_line_sf |>
  filter(ssid %in% hbll_ssids) |>
  filter(survey_abbrev != "HBLL INS S") |>
  st_transform(crs = st_crs(simple_mpa))

if (!file.exists(file.path(overlay_dir, "hbll-mpa-sf.rds"))) {
hbll_mpa_sf <- st_join(hbll_transects, simple_mpa, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
saveRDS(hbll_mpa_sf, file.path(overlay_dir, "hbll-mpa-sf.rds"))
} else {
  hbll_mpa_sf <- readRDS(file.path(overlay_dir, "hbll-mpa-sf.rds"))
}

hbll_mpa_df <- st_drop_geometry(hbll_mpa_sf) |>
  # Sort before distinct to ensure deterministic row selection across systems
  # Need to filter out duplicate fishing event ids that I think arise because of
  # transect overlap with multiple zones are harder to deal with - makes sense
  # that we are looking at site level for the CSAS
  arrange(ssid, year, fe) |>
  distinct(ssid, fe, year, .keep_all = TRUE)

ye_data |> summarise(n_distinct(fishing_event_id))

hbll_mpa_df |> summarise(n_distinct(fe))

# Take peak at duplicate fishing event IDs
# test <- hbll_mpa_df |>
#   group_by(fe) |>
#   filter(n()>1)
# test |> glimpse()

# HBLL ECP data
# -----------------------------------------------------------------------------
gf_ecp_data0 <- readRDS(file.path("data-generated", "hbll-ecp-species.rds"))

gf_ecp_data <- gf_ecp_data0 |>
  filter(survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")) |>
  left_join(ecp_clean |> select(species_code, parent_taxonomic_unit, is_ecp, is_fish))

gf_ecp_site_summary <- gf_ecp_data |>
  distinct(survey_group, species_common_name, species_science_name, in_mpa, site) |>
  group_by(survey_group, site) |>
  summarise(n_species = n(), .groups = "drop")
gf_ecp_site_summary #|>
# saveRDS(file.path("data-generated", "ecp-site-summary.rds"))

gf_ecp_summary <- gf_ecp_data |>
  distinct(survey_group, species_common_name, species_science_name, in_mpa, is_fish, is_ecp) |>
  group_by(survey_group, in_mpa) |>
  summarise(n_species = n(), n_fish_ecps = sum(is_fish), n_non_fish_ecps = sum(!is_fish))
gf_ecp_summary

gf_ecp_data |>
  distinct(survey_group, fishing_event_id, species_common_name, species_science_name, in_mpa, .keep_all = TRUE) |>
  group_by(survey_group, species_common_name, species_science_name, in_mpa) |>
  summarise(encounters = n(), .groups = "drop")

# Complete to one row per (event, species) with zero catch where species was not observed
hbll_event_rows <- gf_ecp_data |>
  distinct(fishing_event_id, survey_series_id, survey_abbrev, year, month, day,
    latitude, longitude, latitude_end, longitude_end, survey_group, in_mpa, uid, site)
hbll_species_lu <- gf_ecp_data |> distinct(species_common_name, species_science_name)

gf_ecp_data_complete <- tidyr::expand_grid(hbll_event_rows, hbll_species_lu) |>
  left_join(
    gf_ecp_data |> select(fishing_event_id, species_common_name, in_mpa, uid, site, catch_count, catch_weight, present),
    by = c("fishing_event_id", "species_common_name", "in_mpa", "uid", "site")
  ) |>
  mutate(
    catch_count = tidyr::replace_na(catch_count, 0),
    catch_weight = tidyr::replace_na(catch_weight, 0),
    present = tidyr::replace_na(present, 0)
  )
saveRDS(gf_ecp_data_complete, file.path(overlay_dir, "gf-ecp-data-complete.rds"))

hbll_spp_encounter_rate <- gf_ecp_data_complete |>
  group_by(survey_group, species_common_name, species_science_name) |>
  summarise(encounters = sum(present), .groups = "drop",
            n_sets = n(),
            pos_sets = encounters / n_sets) |>
  arrange(desc(pos_sets))
saveRDS(hbll_spp_encounter_rate, file.path(overlay_dir, "hbll-spp-encounter-rate.rds"))

# Plot overlays -----------------------------------------------------------------
plot_overlays(hbll_mpa_df, hbll_mpa_sf, title = "HBLL set locations",
              geom_size = 0.5,
              legend_title = "Number of HBLL survey sets",
              breaks = c(-0.1, 0, 10, 50, 100, Inf),
              labels = c("0", "1-10", "11-50", "51-100", "100+"))
ggsave(file.path(overlay_fig_dir, "heatmap-hbll-mpa-sf.png"), width = 4.5, height = 5.3)

# Survey summary tables --------------------------------------------------------

hbll_site_zone_counts <- hbll_mpa_df |>
  filter(in_mpa == 1) |>
  reframe(n_sites = n_distinct(site), n_zones = n_distinct(uid)) |>
  mutate(n_sites_text = paste0(n_sites, "/", n_mpa_sites),
         n_zones_text = paste0(n_zones, "/", n_mpa_zones))

hbll_zones_with_ecps <- gf_ecp_data_complete |>
  filter(in_mpa == 1) |>
  distinct(site, uid) |>
  summarise(n_zones = n()) |>
  mutate(n_zones_text = paste0(n_zones, "/", n_mpa_zones))

hbll_sites_with_ecps <- gf_ecp_data_complete |>
  filter(in_mpa == 1) |>
  distinct(site) |>
  summarise(n_sites = n()) |>
  mutate(n_sites_text = paste0(n_sites, "/", n_mpa_sites))

hbll_longest_timeseries <- hbll_mpa_df |>
  filter(in_mpa == 1) |>
  # group_by(site, uid) |>
  group_by(site) |>
  summarise(longest_timeseries = max(year) - min(year),
            n_years = length(unique(year[in_mpa == 1])),
            .groups = "drop") |>
  arrange(desc(longest_timeseries)) |>
  mutate(text = paste0(n_years, " years (", site, ")"))

hbll_summary_table <- hbll_mpa_df |>
  mutate(survey_name = "Hard Bottom Longline Survey") |>
  reframe(
    `Survey name` = "Hard Bottom Longline Survey",
    `Years evaluated here` = paste(min(year), max(year), sep = "-"),
    `Survey schedule` = "Biennial (North and South regions in alternating years)",
    `Total sets` = paste0(n_distinct(fe), " fishing events"),
    n_sets = n_distinct(fe),
    `years_not_in_mpa` = paste(setdiff(min(year):max(year), unique(year[in_mpa == 1])), collapse = ", "),
    `Years with surveys inside MPA boundaries` =
      paste0(
        n_distinct(year[in_mpa == 1]), " of ", max(year) - min(year) + 1,
        " (no surveys inside MPAs in ", `years_not_in_mpa`, ")"
      ),
    `Sets inside MPA boundaries` = sum(in_mpa),
    `Sets outside MPA boundaries` = sum(in_mpa == 0),
    `% of effort inside MPA boundaries` = round(100 * `Sets inside MPA boundaries` / n_sets, 1),
    n_sites_sampled = hbll_sites_with_ecps$n_sites,
    `Number of MPAN sites sampled` = hbll_site_zone_counts$n_sites_text,
    # `Number of MPAN zones sampled` = hbll_site_zone_counts$n_zones_text,
    `Number of MPAN sites with E-CP sightings` = hbll_sites_with_ecps$n_sites_text,
    # `Number of MPAN zones with E-CP sightings` = hbll_zones_with_ecps$n_zones_text,
    `Longest timeseries at any single site` = hbll_longest_timeseries[1, ]$text
  ) |>
  mutate(across(everything(), as.character)) |>
  tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value")
saveRDS(hbll_summary_table, file.path(overlay_dir, "hbll-summary-table.rds"))


# HBLL site summary table --------------------
# Site ID
# Common Site Name
# Total transects
# First survey
# Last survey
# No. years
# Unique species
# Mean species richness

hbll_site_ecp_richness <- gf_ecp_data_complete |>
  filter(in_mpa == 1) |>
  filter(present == 1) |>
  group_by(site) |>
  summarise(`E-CP species count` = n_distinct(species_science_name), .groups = "drop")


hbll_zone_ecp_richness <- gf_ecp_data_complete |>
  filter(in_mpa == 1) |>
  filter(present == 1) |>
  group_by(site, uid) |>
  summarise(`E-CP species count` = n_distinct(species_science_name), .groups = "drop")


hbll_site_summary_table <-
  hbll_mpa_df |>
  filter(in_mpa == 1) |>
  # group_by(site, uid) |>
  group_by(site) |>
  summarise(
    # `Site ID` = unique(uid),
    `Common Site Name` = unique(site),
    `Total transects` = length(unique(fe)),
    `First survey` = min(year),
    `Last survey` = max(year),
    `No. years` = n_distinct(year),
    `Years in MPA` = paste(min(year), max(year), sep = "-"),
    .groups = "drop"
  ) |>
  left_join(hbll_site_ecp_richness, by = "site") |>
  # left_join(hbll_zone_ecp_richness, by = c("site", "uid")) |>
  arrange(desc(`Total transects`))

saveRDS(hbll_site_summary_table, file.path(overlay_dir, "hbll-site-summary-table.rds"))
