# Appendix figures for the survey spatial overlay analysis
# 1) Large table summarising number of transects in each MPA site (zone) by subregion
# 2) Map of survey transects by survey type and subregion/MPA site (zone) (make heat maps)

source(here::here("R", "00-setup.R"))
source(here::here("R", "00-utils.R"))

library(dplyr)
library(purrr)
library(tidyr)
library(sf)
library(officer)
library(flextable)
library(patchwork)
library(ggrepel)

survey_group_lu <- tribble(
  ~survey_series_id, ~survey_group,
  1, "Synoptic",
  3, "Synoptic",
  4, "Synoptic",
  16, "Synoptic",
  22, "HBLL",
  35, "Sablefish",
  36, "HBLL",
  39, "HBLL",
  40, "HBLL"
)

create_linestring <- function(lon_start, lat_start, lon_end, lat_end) {
  st_linestring(matrix(c(lon_start, lat_start, lon_end, lat_end),
                        ncol = 2, byrow = TRUE))
}


create_site_label_ranges <- function(data, site_col = site,
                                    label_col = map_label) {

# Format a vector of numbers as ranges
format_label_ranges <- function(labels) {
  # Sort labels
  labels_sorted <- sort(labels)

  # Find breaks (where diff > 1)
  breaks <- c(0, which(diff(labels_sorted) > 1), length(labels_sorted))

  # Create range strings
  ranges <- map_chr(seq_len(length(breaks) - 1), function(i) {
    start_idx <- breaks[i] + 1
    end_idx <- breaks[i + 1]

    start_val <- labels_sorted[start_idx]
    end_val <- labels_sorted[end_idx]

    if (start_val == end_val) {
      as.character(start_val)
    } else {
      paste0(start_val, "-", end_val)
    }
  })

  # Join with comma-space
  paste(ranges, collapse = ",")
}

  data |>
    arrange(subregion, {{ site_col }}, {{ label_col }}) |>
    group_by({{ site_col }}, subregion) |>
    summarise(
      zone_label = format_label_ranges({{ label_col }}),
      .groups = "drop"
    )
}


# Prepare directories
# -------------------
table_dir <- here::here("tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

overlay_dir <- here::here("data-generated", "spatial", "overlays")
dir.create(overlay_dir, recursive = TRUE, showWarnings = FALSE)

overlay_fig_dir <- here::here("figures", "spatial-overlays")
dir.create(overlay_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Load spatial data for plotting and overlay analysis
simple_coast <- pacea::bc_coast |>
  st_transform(crs = 3005) |>
  st_simplify(dTolerance = 100)

simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-analytical-mpa.rds")) |>
  mutate(site = gsub("_", " ", common_site_name_site_profile), map_id = map_label)

mpa_sf <- simple_mpa |> select(uid:area_km2, shape_area:map_id)

subregion_uid_lu <- mpa_sf |>
  st_drop_geometry() |>
  distinct(uid, subregion)

zone_labels <- mpa_sf |>
  st_drop_geometry() |>
  create_site_label_ranges(
    site_col = site,
    label_col = map_label
  )

all_subregion_sites <- mpa_sf |>
  st_drop_geometry() |>
  distinct(site, subregion) |>
  left_join(zone_labels, by = c("site", "subregion"))

display_mpa <- mpa_sf |> st_simplify(dTolerance = 100) |>
  left_join(zone_labels, by = c("site", "subregion"))
mpa_count <- unique(mpa_sf$uid) |> length()

subregions <- readRDS(here::here("data-generated", "spatial", "subregion-masks.rds"))
display_subregions <- subregions |> st_simplify(dTolerance = 100)

sp <- sp_to_hyphens("yelloweye rockfish")
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

# -----------------------------------------------------------------------------
# Groundfish surveys
# -----------------------------------------------------------------------------
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
    latitude_end, longitude_end)

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
  st_transform(crs = st_crs(mpa_sf))

if (!file.exists(file.path(overlay_dir, "hbll-mpa-sf.rds"))) {
hbll_mpa_sf <- st_join(hbll_transects, mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
saveRDS(hbll_mpa_sf, file.path(overlay_dir, "hbll-mpa-sf.rds"))
} else {
  hbll_mpa_sf <- readRDS(file.path(overlay_dir, "hbll-mpa-sf.rds"))
}

hbll_mpa_df <- st_drop_geometry(hbll_mpa_sf)

# Synoptic bottom trawl polygons (swath transects)
# -----------------------------------------------------------------------------
if (!file.exists(file.path(overlay_dir, "syn-mpa-sf.rds"))) {
  syn_transects <- gf_line_sf |> filter(ssid %in% syn_ssids)

  mean_doorspread <- mean(syn_transects$doorspread_m, na.rm = TRUE)
  doorspreads <- syn_transects |>
    mutate(doorspread_m = ifelse(is.na(doorspread_m), mean_doorspread, doorspread_m)) |>
    pull(doorspread_m)

  syn_polygons <- syn_transects |>
    st_transform(crs = 32609) |>
    st_buffer(dist = doorspreads) # Buffer to each doorspread (use mean if missing)

  # Create synoptic MPA overlay
  syn_mpa_sf <- syn_polygons |>
    st_transform(crs = st_crs(mpa_sf)) |>
    st_join(mpa_sf, join = st_intersects) |>
    mutate(in_mpa = ifelse(is.na(uid), 0, 1))

  saveRDS(syn_mpa_sf, file.path(overlay_dir, "syn-mpa-sf.rds"))
} else {
  syn_mpa_sf <- readRDS(file.path(overlay_dir, "syn-mpa-sf.rds"))
}
syn_mpa_df <- st_drop_geometry(syn_mpa_sf)

# How much fishing pressure in how many MPAs? - e.g., low medium high;
# Doesn't want to mix trawl and longline fishing pressure
# How many MPAs
# Tomorrow we will talk about the different zones situation to prep the simulation scenarios


# Multi-species dive data
# ------------------------------------------------------------------------------
msd0 <- readRDS(here::here("data-generated", "msd-catch.rds"))
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
  st_transform(crs = st_crs(mpa_sf)) |>
  st_join(mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
msd_mpa_df <- st_drop_geometry(msd_mpa_sf)


# hbll_overlay_sf <- left_join(display_mpa, hbll_overlay_df, by = "uid") |> mutate(survey = "HBLL")
# syn_overlay_sf <- left_join(display_mpa, syn_overlay_df, by = "uid") |> mutate(survey = "SYN")
# msd_overlay_sf <- left_join(display_mpa, msd_overlay_df, by = "uid") |> mutate(survey = "MSD")

# hbll_overlay_plot <-
# ggplot() +
#   geom_sf(data = simple_coast |> rotate_a(), fill = "grey95") +
#   geom_sf(data = hbll_overlay_sf |> rotate_a(), aes(fill = n_transects)) +
#   # geom_sf(data = hbll_mpa_sf |> rotate_a() |> st_centroid(),
#   #   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   scale_fill_viridis_c(na.value = "grey60", name = "Number of survey sets") +
#   theme(legend.position = "top") +
#   guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 10,
#     display = "rectangles")) +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0) +
#   ggtitle("HBLL Survey Sets")

# syn_overlay_plot <-
# ggplot() +
#   geom_sf(data = simple_coast |> rotate_a(), fill = "grey95") +
#   geom_sf(data = syn_overlay_sf |> rotate_a(), aes(fill = n_transects)) +
#   # geom_sf(data = syn_mpa_sf |> rotate_a() |> st_centroid(),
#   #   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   scale_fill_viridis_c(na.value = "grey60", name = "Number of survey sets") +
#   # theme(legend.position = "inside",
#   #       legend.position.inside = c(0.1, 0.1)) +
#   theme(legend.position = "top",
#         axis.text.y = element_blank()) +
#   guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 10,
#     display = "rectangles")) +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a()) +
#   ggtitle("Synoptic Survey Sets")

# msd_overlay_plot <-
# ggplot() +
#   geom_sf(data = simple_coast |> rotate_a(), fill = "grey95") +
#   geom_sf(data = msd_overlay_sf |> rotate_a(), aes(fill = n_transects)) +
#   # geom_sf(data = msd_mpa_sf |> rotate_a() |> st_centroid(),
#   #   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   scale_fill_viridis_c(na.value = "grey60", name = "Number of transects") +
#   theme(legend.position = "top",
#         axis.text.y = element_blank()) +
#   guides(fill = guide_coloursteps(title.position = "top", title.hjust = 0.5, barwidth = 10,
#     display = "rectangles")) +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0) +
#   ggtitle("Multispecies Dive Survey Sets")

# hbll_overlay_plot + syn_overlay_plot + msd_overlay_plot

# (hbll_overlay_plot +
#   geom_sf(data = hbll_mpa_sf |> rotate_a() |> st_centroid(),
#     shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0)
# ) +
# (syn_overlay_plot +
#   geom_sf(data = syn_mpa_sf |> rotate_a() |> st_centroid(),
#    shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0)
# ) +
# (msd_overlay_plot +
#   geom_sf(data = msd_mpa_sf |> rotate_a() |> st_centroid(),
#   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0)
# )

# IPHC spatial overlay ---------------------------------------------------------
# iphc <- gfdata::load_iphc_dat(species = "yelloweye rockfish")
# iphc_sf <- st_as_sf(iphc, coords = c("longitude", "latitude"), crs = 4326)
# iphc_mpa_sf <- iphc_sf |>
#   st_transform(crs = st_crs(mpa_sf)) |>
#   st_join(mpa_sf, join = st_intersects) |>
#   mutate(in_mpa = ifelse(is.na(uid), 0, 1))
# iphc_mpa_df <- st_drop_geometry(iphc_mpa_sf)

# # eDNA data --------------------------------------------------------------------
# # Direct from Emily
# edna0 <- readxl::read_excel(here::here("data-raw", "overlay", "SEAC_master_eDNA_spatial_file_updated_Dec_2024.xlsx")) |>
#   janitor::clean_names()

# edna_sf <- st_as_sf(edna0, coords = c("lon_dd", "lat_dd"), crs = 4326) |>
#   mutate(survey = "eDNA")
# edna_mpa_sf <- edna_sf |>
#   st_transform(crs = st_crs(mpa_sf)) |>
#   st_join(mpa_sf, join = st_intersects) |>
#   mutate(in_mpa = ifelse(is.na(uid), 0, 1))
# edna_mpa_df <- st_drop_geometry(edna_mpa_sf)

# ROV data ---------------------------------------------------------------------
# Emily requested from Rob Skelly from ROV database
rov_sf <- st_read(here::here("data-raw", "overlay", "rov_dives_0", "c4e42f2a4df449ec9b70897ee43e5415.shp")) |>
  janitor::clean_names()

rov_mpa_sf <- rov_sf |>
  st_transform(crs = st_crs(mpa_sf)) |>
  st_join(mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1)) |>
  mutate(year = lubridate::year(lubridate::ymd_hms(start_time)))

rov_mpa_df <- st_drop_geometry(rov_mpa_sf)

# Get ECP data -----------------------------------------------------------------

# Only works on DFO network
if (!file.exists(file.path("data-raw", "species-table.rds"))) {
  species_table <- run_sql(
    db = "GFBioSQL",
    query = "
    SELECT *
    FROM SPECIES
    "
  )
  saveRDS(species_table, file = "data-raw/species-table.rds")
} else {
  species_table <- readRDS(file = "data-raw/species-table.rds")
}

if (!file.exists(file.path("data-generated", "ecp-clean.rds"))) {
  ecps <- readxl::read_excel("data-raw/gale CPs.xlsx") |>
    janitor::clean_names() |>
    mutate(common_name = tolower(common_name),
          scientific_name = tolower(scientific_name)) |>
    # Taxonomic corrections - align with species table from GFBio
    mutate(scientific_name = gsub("nereocystis leutkeana", "nereocystis luetkeanus", scientific_name),
          scientific_name = gsub("macrocystis sp.", "macrocystis integrifolia", scientific_name),
          scientific_name = gsub("raja binoculata", "beringraja binoculata", scientific_name),
          scientific_name = gsub("cymatogaster aggregate", "cymatogaster aggregata", scientific_name),
          scientific_name = gsub("ammodytes hexapterus", "ammodytes personatus", scientific_name),
          scientific_name = gsub("salvelinus malma lordi", "salvelinus malma", scientific_name),
          scientific_name = gsub("theragra chalcogramma", "gadus chalcogrammus", scientific_name),
    )

  # Note that pandalus borealis is now split into pandalus eous and pandalus jordani
  # Groundfish surveys do not identify tunicates to species level (i.e., no data on invasive tunicates)
  species_to_add <- tribble(
    ~common_name,             ~scientific_name,
    "pink shrimp (spiny)",    "pandalus eous",
    "pink shrimp (smooth)",   "pandalus jordani",
    "giant red sea cucumber", "apostichopus californicus"
  )

  surf_eel_grass <- tribble(
  ~common_name, ~scientific_name, ~species_science_name, ~species_code, ~species_common_name, ~parent_taxonomic_unit, ~taxonomic_rank, ~itis_tsn,
 "surfgrass", "phyllospadix sp.", "phyllospadix", "PH", "surfgrass", "zosteraceae", "genus", 39070,
  "eelgrass", "zostera sp.", "zostera", "ZO", "eel grass", "zosteraceae", "genus", 39073
  )

  ecps <- bind_rows(ecps, species_to_add)

  # For each survey, I need to extract all of the catch data for each species in
  # the ecps list

  # Step 1 - get the species table from GFBio so that I can make sure that the
  # species names match the CP list

  # species_table <- run_sql(
  #   database = "GFBioSQL",
  #   query = " SELECT * FROM SPECIES"
  # )
  # saveRDS(species_table, file = "data-raw/species-table.rds")

  species_table0 <- readRDS(file = "data-raw/species-table.rds") |>
    janitor::clean_names() |>
    as_tibble() |>
    mutate(
      species_science_name = tolower(species_science_name),
      species_common_name = tolower(species_common_name),
      parent_taxonomic_unit = tolower(parent_taxonomic_unit)
    )
  species_table <- species_table0 |>
    select(species_code, species_common_name, species_science_name, parent_taxonomic_unit, taxonomic_rank, itis_tsn)

  ecp_clean <- left_join(ecps, species_table, by = c("scientific_name" = "species_science_name"),
  relationship = "many-to-many") |> # many to many because of life stage of chinook, and resident/offshore/etc, orcas
    distinct(scientific_name, .keep_all = TRUE) |>
    filter(scientific_name != "phyllospadix sp.") |> # cleaned this up
    bind_rows(surf_eel_grass) |>
    mutate(is_ecp = TRUE)

  # Double check that these species make sense to not be in GFBio
  anti_join(ecp_clean, species_table, by = c("scientific_name" = "species_science_name")) |>
    select(common_name, scientific_name)

  saveRDS(ecp_clean, file = "data-generated/ecp-clean.rds")
} else {
  ecp_clean <- readRDS(file.path("data-generated", "ecp-clean.rds"))
}

# Species not queried from ECP list
# # A tibble: 7 × 2
#   common_name                  scientific_name
#   <chr>                        <chr>
# 1 spiny/northern pink shrimp   pandalus borealis --> name change
# 2 neocalanus copepods          neocalanus sp.
# 3 other crustacean zooplankton other crustacean zooplankton
# 4 littorina snail              littorina sp.
# 5 non-crustacean zooplankton   non-crustacean zooplankton
# 6 phytoplankton                phytoplankton
# 7 surfgrass                    phyllospadix sp.

# Get GFBio data for ECPs-------------------------------------------------------
if (!file.exists(file.path("data-generated", "hbll-ecp-species.rds"))) {
  if (!file.exists(file.path("data-raw", "mpan-csas-gf-ecp-data.rds"))) {
  ssid_lu0 <- get_ssids()
  ssid_lu <- ssid_lu0 |> janitor::clean_names() |>
    select(-survey_abbrev)
  survey_sets <- get_all_survey_sets(species = ecps_clean$species_common_name,
                      ssid =c(1, 3, 4, 6, 16, 22, 35, 36, 39, 40, 41, 42, 43, 46))

  gf_ecp_data0 <- left_join(survey_sets, ssid_lu, by = "survey_series_id") |>
    saveRDS(gf_ecp_data0, file.path("data-raw", "mpan-csas-gf-ecp-data.rds"))
  } else {
    gf_ecp_data0 <- readRDS(file.path("data-raw", "mpan-csas-gf-ecp-data.rds"))
  }

  ssid_lu <- gf_ecp_data0 |>
    distinct(survey_series_id, survey_abbrev)

  in_mpa_lu <- bind_rows(syn_mpa_df, hbll_mpa_df) |>
    distinct(survey_series_id = ssid, survey_abbrev, fishing_event_id = fe, in_mpa, uid, site)

  gf_ecp_data <- gf_ecp_data0 |>
    select(survey_series_id, species_common_name, species_science_name,
      fishing_event_id,
      survey_abbrev, year, month, day, latitude, longitude, latitude_end, longitude_end,
      catch_count, catch_weight) |>
    left_join(survey_group_lu, by = "survey_series_id") |>
    mutate(present = if_else(catch_count > 0 | catch_weight > 0, 1, 0)) |>
    filter(survey_abbrev %in% c(#"SYN HS", "SYN QCS", "SYN WCHG",
    "HBLL OUT N", "HBLL OUT S", "HBLL INS N")) |>
    tidyr::drop_na(species_common_name) |>
    filter(present == 1) |>
    left_join(in_mpa_lu, by = c("survey_series_id", "survey_abbrev", "fishing_event_id"), relationship = "many-to-many") |>
    drop_na(in_mpa)
    saveRDS(gf_ecp_data, file.path("data-generated", "hbll-ecp-species.rds"))
  } else {
    gf_ecp_data <- readRDS(file.path("data-generated", "hbll-ecp-species.rds"))
}

gf_ecp_summary <- gf_ecp_data |>
  distinct(survey_group, species_common_name, species_science_name, in_mpa) |>
  group_by(survey_group, in_mpa) |>
  summarise(n_species = n())
gf_ecp_summary

gf_ecp_site_summary <- gf_ecp_data |>
  distinct(survey_group, species_common_name, species_science_name, in_mpa, site) |>
  group_by(survey_group, site) |>
  summarise(n_species = n(), .groups = "drop")
gf_ecp_site_summary #|>
# saveRDS(file.path("data-generated", "ecp-site-summary.rds"))

gf_ecp_data |>
  distinct(survey_group, species_common_name, species_science_name, in_mpa) |>
  group_by(survey_group, in_mpa) |>
  summarise(n_species = n())

gf_ecp_data |>
  distinct(survey_group, fishing_event_id, species_common_name, species_science_name, in_mpa, .keep_all = TRUE) |>
  group_by(survey_group, species_common_name, species_science_name, in_mpa) |>
  summarise(encounters = n(), .groups = "drop")

# Complete to one row per (event, species) with zero catch where species was not observed
hbll_event_rows <- gf_ecp_data |>
  distinct(fishing_event_id, survey_series_id, survey_abbrev, year, month, day,
    latitude, longitude, latitude_end, longitude_end, survey_group, in_mpa, uid, site)
hbll_species_lu <- gf_ecp_data |> distinct(species_common_name, species_science_name)

gf_ecp_data_complete <- expand_grid(hbll_event_rows, hbll_species_lu) |>
  left_join(
    gf_ecp_data |> select(fishing_event_id, species_common_name, in_mpa, uid, site, catch_count, catch_weight, present),
    by = c("fishing_event_id", "species_common_name", "in_mpa", "uid", "site")
  ) |>
  mutate(
    catch_weight = replace_na(catch_weight, 0),
    present = replace_na(present, 0)
  )

hbll_spp_encounter_rate <- gf_ecp_data_complete |>
  group_by(survey_group, species_common_name, species_science_name) |>
  summarise(encounters = sum(present), .groups = "drop",
            n_sets = n(),
            pos_sets = encounters / n_sets) |>
  arrange(desc(pos_sets))
saveRDS(hbll_spp_encounter_rate, file.path("data-generated", "hbll-spp-encounter-rate.rds"))



# -----------------------------------------------------------------------------

# Heatmap overlays -------------------------------------------------------------
library(rnaturalearth)

ne_coast <- ne_states(country = c("canada", "united states of america"), returnclass = "sf") |>
  # filter(region == "West") |>
  st_transform(crs = 3005)


plot_overlays <- function(survey_mpa_df, survey_mpa_sf, title, legend_title,
                          geom_size = waiver(),
                          return_data = FALSE, breaks = NULL, labels = NULL) {

  subregions <- unique(survey_mpa_df$subregion)
  factor_subregions <- factor(subregions, levels = c("HG", "NC", "CC", "NVI"))
  num_subregions <- length(sort(factor_subregions))
  message(paste0("Number of subregions: ", num_subregions, " (", paste(sort(factor_subregions), collapse = ", "), ")"))

  mpa_summary_sf <- survey_mpa_df |>
    group_by(uid) |>
    summarise(n = n(), .groups = "drop") |>
    left_join(x = display_mpa, y = _, by = "uid")

  if (return_data) {
    return(mpa_summary_sf)
  }

  north_arrow <- function() {
    ggspatial::annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(0.4, "cm"),
    width = unit(0.4, "cm"),
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.2, "cm"),
    style = ggspatial::north_arrow_orienteering(fill = c("black", "black"),
    text_size = 5)
  )}

  map_theme <- function() {
    theme(panel.background = element_rect(fill = ocean_colour),
        legend.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        legend.margin = margin(t = 0, r = 0, b = 0.1, l = 0),
        legend.box.spacing = unit(0.1, "cm"),
        legend.position = "top",
        legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.spacing.x = unit(0.15, "cm"),
        legend.text.align = 0.5,
        legend.text = element_text(size = 7))
  }

  # Bin data if breaks provided
  use_bins <- !is.null(breaks)

  if (use_bins) {
    mpa_summary_sf <- mpa_summary_sf |>
      mutate(n = ifelse(is.na(n), 0, n)) |>
      mutate(n_binned = cut(n, breaks = breaks, labels = labels,
                            include.lowest = TRUE, right = TRUE))
  }

  #zero_colour <- "white"
  zero_colour <- "#feedde"
  ocean_colour <- "#b8d4ed"
  land_colour <- "grey70"

    # Grey for "0", viridis for sample counts
    colours <- if (labels[1] == "0") {
      c(zero_colour, viridis::viridis(length(labels) - 1))
    } else {
      colours <- c(zero_colour, "#2E7D32", "#FDD835", "#FB8C00", "#D32F2F")
    }

    # colours <- c(zero_colour, RColorBrewer::brewer.pal(length(labels) - 1, "YlOrRd"))
    # colours <- c(zero_colour, "#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3")
    colours <- c(zero_colour, "#fdbe85", "#fd8d3c", "#d94701", "#a63603")

  tag_theme <- theme(
    plot.tag = element_text(size = 9),
    plot.tag.position = c(0.02, 0.95)
  )

  p1 <- ggplot(data = survey_mpa_sf |> rotate_sf()) +
      geom_sf(data = ne_coast |> rotate_sf(), fill = land_colour, linewidth = 0.08) +
      geom_sf(data = mpa_summary_sf |> rotate_sf(), fill = "white",
              colour = "grey30", linewidth = 0.1) +
      geom_sf(size = geom_size) +
      north_arrow() +
      guides(fill = "none") +
      map_theme() +
      gfplot::coord_sf_auto(mpa_sf |> rotate_sf(), buffer = 0) +
      ggtitle(title) +
      theme(plot.title = element_text(vjust = -7, size = 9)) +
      tag_theme

  p2 <- ggplot() +
    # geom_sf(data = ne_coast |> rotate_sf(), fill = "grey90") +
    geom_sf(data = ne_coast |> rotate_sf(), fill = land_colour, linewidth = 0.08) +
    geom_sf(data = mpa_summary_sf |> rotate_sf(),
            aes(fill = if (use_bins) n_binned else n),
            colour = "grey30", linewidth = 0.1) +
    north_arrow() +
    map_theme() +
    theme(axis.text.y = element_blank()) +
    gfplot::coord_sf_auto(mpa_sf |> rotate_sf(), buffer = 0) +
    guides(fill = guide_legend(
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE                             # Keep keys in order
    )) +
    tag_theme

  if (use_bins) {
    p2 <- p2 + scale_fill_manual(values = colours, na.value = "grey60", name = legend_title)
  } else {
    p2 <- p2 + scale_fill_viridis_c(na.value = "grey60", name = legend_title)
  }
  (p1 + p2) +
    plot_annotation(tag_levels = "a", tag_prefix = "", tag_suffix = ")")
}

# plot_overlays(syn_mpa_df, return_data = TRUE) |>
#   pull(n) |> quantile(probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)#range(na.rm = TRUE)
# plot_overlays(hbll_mpa_df, return_data = TRUE) |>
#   pull(n) |> quantile(probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)#range(na.rm = TRUE)
# plot_overlays(iphc_mpa_df, return_data = TRUE) |>
#   pull(n) |> quantile(probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)#range(na.rm = TRUE)
# plot_overlays(msd_mpa_df, return_data = TRUE) |>
#   pull(n) |> quantile(probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)#range(na.rm = TRUE)
# plot_overlays(rov_mpa_df, return_data = TRUE) |>
#   pull(n) |> quantile(probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)#range(na.rm = TRUE)

plot_overlays(hbll_mpa_df, hbll_mpa_sf, title = "HBLL set locations",
              geom_size = 0.5,
              legend_title = "Number of HBLL survey sets",
              breaks = c(-0.1, 0, 10, 50, 100, Inf),
              labels = c("0", "1-10", "11-50", "51-100", "100+"))
ggsave(file.path(overlay_fig_dir, "heatmap-hbll-mpa-sf.png"), width = 4.5, height = 5.3)
plot_overlays(syn_mpa_df, syn_mpa_sf, title = "SYN set locations",
              geom_size = 0.2,
              legend_title = "Number of SYN survey sets",
              breaks = c(-0.1, 0, 20, 50, 100, Inf),
              labels = c("0", "1-20", "21-50", "51-100", "100+"))
ggsave(file.path(overlay_fig_dir, "heatmap-syn-mpa-sf.png"), width = 4.5, height = 5.3)
plot_overlays(iphc_mpa_df, iphc_mpa_sf, title = "IPHC set locations",
              geom_size = 0.2,
              legend_title = "Number of IPHC survey sets",
              breaks = c(-0.1, 0, 10, 50, 100, Inf),
              labels = c("0", "1-10", "11-50", "51-100", "100+"))
ggsave(file.path(overlay_fig_dir, "heatmap-iphc-mpa-sf.png"), width = 4.5, height = 5.3)
plot_overlays(msd_mpa_df, msd_mpa_sf |> st_centroid(), title = "MSD dive locations",
              geom_size = 0.2,
              legend_title = "Number of MSD dives",
              breaks = c(-0.1, 0, 5, 20, 50, Inf),
              labels = c("0", "1-5", "6-20", "21-50", "50+"))
ggsave(file.path(overlay_fig_dir, "heatmap-msd-mpa-sf.png"), width = 4.5, height = 5.3)
# plot_overlays(edna_mpa_df, "Number of eDNA samples") # No NVI subregion in eDNA data
plot_overlays(rov_mpa_df, rov_mpa_sf, title = "ROV dive locations",
              geom_size = 0.8,
              legend_title = "Number of ROV dives",
              breaks = c(-0.1, 0, 10, 30, 100, Inf),
              labels = c("0", "1-10", "11-30", "31-100", "100+"))
ggsave(file.path(overlay_fig_dir, "heatmap-rov-mpa-sf.png"), width = 4.5, height = 5.3)

################################################################################
# Survey summary tables --------------------------------------------------------

# HBLL summary tables
hbll_sites_with_ecps <- gf_ecp_data_complete |>
  filter(in_mpa == 1) |>
  distinct(site) |>
  summarise(n_sites = n())

hbll_longest_timeseries <- hbll_mpa_df |>
  filter(in_mpa == 1) |>
  group_by(site, uid) |>
  summarise(longest_timeseries = max(year) - min(year),
            n_years = length(unique(year[in_mpa == 1])),
            .groups = "drop") |>
  arrange(desc(longest_timeseries)) |>
  mutate(text = paste0(n_years, " years (", site, ")"))

hbll_summary_table <- hbll_mpa_df |>
  mutate(survey_name = "Hard Bottom Longline Survey") |>
  reframe(
    `Survey name` = "Hard Bottom Longline Survey",
    `Years conducted` = paste(min(year), max(year), sep = "-"),
    `Survey schedule` = "Biennial (North and South regions in alternating years)",
    `Total transects` = paste0(n(), " fishing events"),
    `years_not_in_mpa` = paste(setdiff(min(year):max(year), unique(year[in_mpa == 1])), collapse = ", "),
    `Years with surveys inside MPA boundaries` =
      paste0(
        n_distinct(year[in_mpa == 1]), " of ", max(year) - min(year) + 1,
        " (no surveys inside MPAs in ", `years_not_in_mpa`, ")"
      ),
    `Transects inside MPA boundaries` = sum(in_mpa),
    `Transects outside MPA boundaries` = sum(in_mpa == 0),
    `% of effort inside MPA boundaries` = round(100 * `Transects inside MPA boundaries` / `Total transects`, 1),
    `Number of MPAN sites sampled` = n_distinct(uid, na.rm = TRUE),
    `Number of MPAs with E-CP sightings` = hbll_sites_with_ecps$n_sites,
    `Longest timeseries at any single site` = hbll_longest_timeseries[1, ]$text
  ) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value")
saveRDS(hbll_summary_table, file.path("data-generated", "hbll-summary-table.rds"))


# HBLL site summary table --------------------
# Site ID
# Common Site Name
# Total transects
# First survey
# Last survey
# No. years
# Unique species
# Mean species richness

hbll_ecp_richness <- gf_ecp_data_complete |>
  filter(in_mpa == 1) |>
  filter(present == 1) |>
  group_by(site, uid) |>
  summarise(`E-CP species count` = n_distinct(species_science_name), .groups = "drop")


hbll_site_summary_table <-
  hbll_mpa_df |>
  filter(in_mpa == 1) |>
  group_by(site, uid) |>
  summarise(
    `Site ID` = unique(uid),
    `Common Site Name` = unique(site),
    `Total transects` = length(unique(fe)),
    `First survey` = min(year),
    `Last survey` = max(year),
    `No. years` = n_distinct(year),
    `Years in MPA` = paste(min(year), max(year), sep = "-"),
    .groups = "drop"
  ) |>
  left_join(hbll_ecp_richness, by = c("site", "uid")) |>
  arrange(desc(`Total transects`))

saveRDS(hbll_site_summary_table, file.path("data-generated", "hbll-site-summary-table.rds"))




# [TO BE COMPLETED] Insert species sightings table once data are extracted from the
# spatial analysis — equivalent to Table 5.4.3 in the dive survey case study.
# Recommended columns: Latin name | Common name | MPAN status (E-CP / other) |
# Count inside MPAs | Count outside MPAs | Proportion inside MPAs |
# Proportion with biological samples. Note: figure described as
# 'total counts of each species sighted within NSB MPA sites' is already in the
# document and should be retained here.
#
# For comparability per unit effort, catch is scaled by effort (number of sets)
# in each zone (CPUE = catch per set). Encounter rates use zone-specific
# denominators (sets inside vs outside MPAs).

hbll_ecp_encounter_cpue <-
  gf_ecp_data_complete |>
  group_by(species_science_name, species_common_name) |>
  mutate(species_common_name = gsub("north pacific", "pacific", species_common_name)) |>
  summarise(
    n_sets_inside = n_distinct(.data$fishing_event_id[.data$in_mpa == 1]),
    n_sets_outside = n_distinct(.data$fishing_event_id[.data$in_mpa == 0]),
    encounters_inside = sum(present[in_mpa == 1]),
    encounters_outside = sum(present[in_mpa == 0]),
    encounter_rate_inside = encounters_inside / n_sets_inside,
    encounter_rate_outside = encounters_outside / n_sets_outside,
    total_catch_inside = sum(catch_count[in_mpa == 1], na.rm = TRUE),
    total_catch_outside = sum(catch_count[in_mpa == 0], na.rm = TRUE),
    cpue_inside = total_catch_inside / n_sets_inside,
    cpue_outside = total_catch_outside / n_sets_outside,
    .groups = "drop"
  ) |>
  arrange(desc(encounter_rate_inside))

saveRDS(hbll_ecp_encounter_cpue, file.path("data-generated", "hbll-ecp-encounter-cpue.rds"))
# Long data for tigures: encounter rate and CPUE by zone
hbll_ecp_encounter_long <-
  hbll_ecp_encounter_cpue |>
  select(species_common_name, encounter_rate_inside, encounter_rate_outside) |>
  tidyr::pivot_longer(
    cols = c(encounter_rate_inside, encounter_rate_outside),
    names_to = "zone", values_to = "val"
  ) |>
  mutate(zone = if_else(zone == "encounter_rate_inside", "Inside MPAs", "Outside MPAs"))
# saveRDS(hbll_ecp_encounter_long, file.path("data-generated", "hbll-ecp-encounter-long.rds"))

hbll_ecp_cpue_long <-
  hbll_ecp_encounter_cpue |>
  select(species_common_name, cpue_inside, cpue_outside) |>
  tidyr::pivot_longer(
    cols = c(cpue_inside, cpue_outside),
    names_to = "zone", values_to = "val"
  ) |>
  mutate(zone = if_else(zone == "cpue_inside", "Inside MPAs", "Outside MPAs"))

# Bar chart of CPUE inside vs outside MPAs ------------------------------------
hbll_ecp_cpue_bar <-
  hbll_ecp_cpue_long |>
  mutate(
    species_common_name = stringr::str_to_title(species_common_name),
    species_common_name = gsub(" North Pacific", "", species_common_name)
  ) |>
  group_by(species_common_name) |>
  mutate(max_cpue = max(val, na.rm = TRUE)) |>
  ungroup() |>
  mutate(
    species_common_name = forcats::fct_reorder(species_common_name, max_cpue)
  ) |>
  ggplot(aes(x = val, y = species_common_name, fill = zone)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(
    values = c(
      "Inside MPAs" = "#0072B2",  # blue, colourblind friendly
      "Outside MPAs" = "#D55E00"  # reddish orange, colourblind friendly
    )
  ) +
  labs(
    x = "CPUE (catch per set)",
    y = "",
    fill = "",
    title = "HBLL E-CP species: CPUE inside vs outside MPAs"
  ) +
  gfplot::theme_pbs() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(hjust = 1.03)
  )
hbll_ecp_cpue_bar


# Tigure (table + figure): tile heatmap with values, matching make_mssm_tigure style.
# df must have species_common_name, zone (x-axis), val (fill and text).
# show_species_labels: if FALSE, omit y-axis (species) labels for use as right panel.
make_hbll_ecp_tigure <- function(df, fill_limits = c(0, NA), padding = 0.5, digits = 2L,
                                show_species_labels = TRUE) {
  df$species_common_name <- stringr::str_to_title(df$species_common_name)
  df$species_common_name <- gsub(" North Pacific", "", df$species_common_name)
  sp <- df |>
    filter(zone == "Inside MPAs") |>
    arrange(val) |>
    pull(species_common_name) |>
    unique()
  if (length(sp) == 0) sp <- unique(df$species_common_name)
  df$species_common_name <- factor(df$species_common_name, levels = sp)
  df$zone <- factor(df$zone, levels = c("Inside MPAs", "Outside MPAs"))
  df$txt <- round(df$val, digits)
  g <- df |>
    ggplot(aes(x = zone, y = species_common_name)) +
    geom_tile(aes(fill = val), colour = "white") +
    geom_text(aes(label = txt), size = ggplot2::rel(3), hjust = 0.5, vjust = 0.5) +
    scale_fill_viridis_c(
      limits = fill_limits, begin = 0.15, end = 1, alpha = 0.9, option = "D", direction = 1
    ) +
    xlab("") +
    ylab("") +
    coord_cartesian(
      expand = FALSE,
      xlim = range(as.numeric(df$zone)) + c(-padding, padding),
      ylim = range(as.numeric(df$species_common_name)) + c(-padding - 0.5, padding + 0.5),
      clip = "off"
    ) +
    gfplot::theme_pbs() +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"))) +
    theme(
      plot.background = element_rect(fill = NA),
      plot.margin = margin(t = 0, r = -2, b = 0, l = -2),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(colour = "grey10"),
      axis.text.y = if (show_species_labels) element_text(hjust = 1.03) else element_blank(),
      legend.margin = margin(t = 0, r = 0.1, b = -15, l = 0),
      legend.position = "top",
      legend.box = "horizontal"
    ) +
    scale_x_discrete(position = "top")
  g
}

hbll_ecp_encounter_tigure <- make_hbll_ecp_tigure(
  hbll_ecp_encounter_long, fill_limits = c(0, 1), digits = 2L, show_species_labels = TRUE
) + labs(fill = "Encounter rate")

hbll_ecp_cpue_tigure <- make_hbll_ecp_tigure(
  hbll_ecp_cpue_long, fill_limits = c(0, NA), digits = 2L, show_species_labels = FALSE
) + labs(fill = "CPUE (catch per set)")

hbll_ecp_tigure_combined <-
  hbll_ecp_encounter_tigure + hbll_ecp_cpue_tigure +
  patchwork::plot_annotation(
    title = "HBLL E-CP species: inside vs outside MPAs (per unit effort)",
    theme = theme(plot.title = element_text(size = 11))
  )
hbll_ecp_tigure_combined


# MSD summary table ------------------------------------------------------------
# msd_raw <- readRDS(here::here("data-generated", "multi-species-data-no-abalone-cleaned.rds"))
msd_catch <- readRDS(here::here("data-generated", "msd-catch.rds"))
msd_mpas_ecps <- #msd_raw |>
  msd_catch |>
  left_join(msd_mpa_df |> select(uid, site, transect_site, year)) |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  select(site, uid, year, transect_site, species_code, species, species_science_name, is_ecp)

msd_sites_with_ecps <- msd_mpas_ecps |>
  filter(is_ecp) |>
  distinct(site, uid, is_ecp) |>
  drop_na(uid) |>
  summarise(n_sites = n())

msd_longest_timeseries <- msd_mpa_df |>
  filter(in_mpa == 1) |>
  group_by(site, uid) |>
  summarise(longest_timeseries = max(year) - min(year),
            n_years = length(unique(year[in_mpa == 1])),
            .groups = "drop") |>
  arrange(desc(longest_timeseries)) |>
  slice(1) |>
  mutate(text = paste0(n_years, " years (", site, ")"))

msd_summary_table <-
msd_mpa_df |>
  mutate(survey_name = "Multispecies Benthic Invertebrate Dive Survey") |>
  reframe(
    `Survey name` = "Multispecies Benthic Invertebrate Dive Survey",
    `Years conducted` = paste(min(year), max(year), sep = "-"),
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
    `Number of MPAN sites sampled` = n_distinct(uid, na.rm = TRUE),
    `Number of MPAs with E-CP sightings` = msd_sites_with_ecps$n_sites,
    `Longest timeseries at any single site` = msd_longest_timeseries$text
  ) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value")
saveRDS(msd_summary_table, file.path("data-generated", "msd-summary-table.rds"))


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
  msd_catch |>
  left_join(msd_mpa_df |> select(uid, site, transect_site, year), by = c("transect_site", "year")) |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  filter(!is.na(uid)) |>  # Only keep observations inside MPAs
  group_by(site, uid) |>
  summarise(
    `Unique species` = n_distinct(species),
    unique_species = paste(unique(species), collapse = ", "),
    .groups = "drop"
  )

# msd_raw |>
msd_catch |>
  # select(algae_code:species) |>
  filter(stringr::str_detect(species_science_name, "zostera")) |>
  distinct() |>
  glimpse()

ecp_clean |>
  filter(stringr::str_detect(scientific_name, "zostera")) |>
  distinct() |>
  glimpse()

msd_richness <- readRDS(here::here("data-generated", "msd-transect-richness.rds")) |>
  left_join(msd_mpa_df |> select(uid, site, transect_site, year)) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1)) |>
  filter(in_mpa == 1) |>  # Only transects inside MPAs
  group_by(site, uid) |>
  summarise(`Mean species richness` = mean(richness), .groups = "drop")

msd_site_summary_table <-
msd_mpa_df |>
  filter(in_mpa == 1) |>
  group_by(site, uid) |>
  summarise(
    `Site ID` = unique(uid),
    `Common Site Name` = unique(site),
    `Total transects` = length(unique(transect_site)),
    `First survey` = min(year),
    `Last survey` = max(year),
    `No. years` = n_distinct(year),
    `Years in MPA` = paste(unique(year[in_mpa == 1]), collapse = ", "),
    .groups = "drop"
  ) |>
  left_join(msd_richness, by = c("site", "uid")) |>
  left_join(msd_species_list, by = c("site", "uid")) |>
  arrange(desc(`Total transects`))

saveRDS(msd_site_summary_table, file.path("data-generated", "msd-site-summary-table.rds"))

# MSD site summary table --------------------
msd_species_summary_table <-
  msd_catch |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  distinct(year, transect_site, species_science_name, species_common_name, species_desc, is_ecp, species_group, catch_count) |>
  arrange(is_ecp) |>
  left_join(msd_mpa_df |> select(uid, site, transect_site, year), by = c("transect_site", "year")) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1)) |>  # Only keep observations inside MPAs
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
saveRDS(msd_species_summary_table, file.path("data-generated", "msd-species-summary-table.rds"))

# Site-level E-CP sightings (for narrative: "most productive sites")
msd_site_ecp_sightings <-
  msd_catch |>
  left_join(msd_mpa_df |> select(uid, site, transect_site, year), by = c("transect_site", "year")) |>
  left_join(ecp_clean |> select(species_code, is_ecp)) |>
  filter(is_ecp == TRUE) |>
  group_by(site) |>
  summarise(total_ecp_sightings = sum(catch_count), .groups = "drop") |>
  arrange(desc(total_ecp_sightings))
saveRDS(msd_site_ecp_sightings, file.path("data-generated", "msd-site-ecp-sightings.rds"))

# ------------------------------------------------------------------------------


################################################################################
################################################################################




metric_labels <- c(
  "monitoring_program" = "Monitoring program name",
  "sampling_frequency" = "Sampling frequency",
  "year_range" = "Year range",
  "n_years" = "Number of years",
  "n_years_in_mpa" = "Number of years in MPA",
  "n_transects" = "Number of transects",
  "n_transects_in_mpa" = "Number of transects in MPA",
  "n_transects_outside_mpa" = "Number of transects outside MPA",
  "prop_in_mpa" = "Proportion of transects in MPA (%)",
  "n_mpas" = "Number of MPAs"
)

hbll_flextable <- hbll_mpa_df |>
  mutate(
    program = case_when(
      survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S") ~ "Hard Bottom Longline Outside",
      survey_abbrev %in% c("HBLL INS N", "HBLL INS S") ~ "Hard Bottom Longline Inside",
      TRUE ~ survey_abbrev
    ),
    region = case_when(
      survey_abbrev %in% c("HBLL OUT N", "HBLL INS N") ~ "North",
      survey_abbrev %in% c("HBLL OUT S", "HBLL INS S") ~ "South",
      TRUE ~ survey_abbrev
    ),
    program_region = paste(program, region, sep = " ")
  ) |>
  select(-fo_id) |>
  group_by(program, region, program_region) |>
  summarise(
    program_region = unique(program_region),
    sampling_frequency = "biennial",
    year_range = paste(min(year), max(year), sep = "-"),
    n_years = n_distinct(year),
    n_years_in_mpa = n_distinct(year[in_mpa == 1]),
    n_transects = n(),
    n_transects_in_mpa = sum(in_mpa),
    prop_in_mpa = round(100 * n_transects_in_mpa / n_transects, 1),
    n_mpas = n_distinct(uid, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(across(-c(program, region, program_region), as.character)) |>
  pivot_longer(cols = -program_region, names_to = "metric", values_to = "value") |>
  pivot_wider(
    names_from = program_region,
    values_from = value,
    names_vary = "slowest"
  ) |>
  select(metric, `Hard Bottom Longline Outside North`, `Hard Bottom Longline Outside South`, `Hard Bottom Longline Inside North`) |>
  mutate(metric = metric_labels[metric]) |>
  drop_na(metric) |>
  flextable() |>
  add_header_row(
    values = c("Monitoring program name", "Hard Bottom Longline Outside", "Hard Bottom Longline Inside"),
    colwidths = c(1, 2, 1),
    top = TRUE
  ) |>
  set_header_labels(
    metric = "Survey region",
    `Hard Bottom Longline Outside North` = "North",
    `Hard Bottom Longline Outside South` = "South",
    `Hard Bottom Longline Inside North` = "North"
  ) |>
  align(align = "right", part = "all") |>
  autofit()

hbll_flextable
saveRDS(hbll_flextable, file.path(table_dir, "hbll-summary-flextable.rds"))

# HBLL aggregate (all surveys combined) ----------------------------------------
hbll_aggregate_flextable <- hbll_mpa_df |>
  select(-fo_id) |>
  summarise(
    monitoring_program = "Hard Bottom Longline",
    sampling_frequency = "annual",
    year_range = paste(min(year), max(year), sep = "-"),
    n_years = n_distinct(year),
    n_years_in_mpa = n_distinct(year[in_mpa == 1]),
    n_transects = n(),
    n_transects_in_mpa = sum(in_mpa),
    prop_in_mpa = round(100 * n_transects_in_mpa / n_transects, 1),
    n_mpas = n_distinct(uid, na.rm = TRUE)
  ) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(everything(), names_to = "metric", values_to = "value") |>
  mutate(metric = metric_labels[metric]) |>
  drop_na(metric) |>
  flextable(col_keys = c("metric", "value")) |>
  set_header_labels(metric = "Metric", value = "Hard Bottom Longline (all surveys)") |>
  align(align = "right", part = "all") |>
  autofit()

saveRDS(hbll_aggregate_flextable, file.path(table_dir, "hbll-summary-aggregate-flextable.rds"))

hbll_mpa_df |>
  filter(in_mpa == 1) |>
  distinct(site) |>
  pull(site) |>
  length()

hbll_mpa_df |>
  filter(in_mpa == 1) |>
  group_by(site) |>
  summarise(n_years = n_distinct(year)) |>
  arrange(desc(n_years))

hbll_mpa_df |>
  filter(in_mpa == 1) |>
  group_by(site, uid) |>
  summarise(
    n_transects = length(unique(fe)),
    first_year = min(year),
    last_year = max(year),
    n_years = n_distinct(year)
  ) |>
  arrange(desc(n_transects)) |>
  left_join(gf_ecp_site_summary |> filter(survey_group == "HBLL")) |>
  left_join(subregion_uid_lu, by = "uid") |>
  select(-survey_group) |>
flextable() |>
  set_header_labels(
    uid = "Site UID",site = "Site Name",
    subregion = "Subregion",
    n_transects = "Total fishing events",
    first_year = "First year",
    last_year = "Last year",
    n_years = "Number of years within-MPA events",
    n_species = "E-CP species count") |>
  colformat_int(j = c("first_year", "last_year", "n_transects", "n_years", "n_species"), big.mark = "") |>
  align(align = "right", part = "all") |>
  autofit() |>
saveRDS(file.path(table_dir, "hbll-site-summary.rds"))




# syn summary table ------------------------------------------------------------
# syn_flextable <- syn_mpa_df |>
#   mutate(survey_name = case_when(
#     survey_abbrev == "SYN HS" ~ "Hecate Strait",
#     survey_abbrev == "SYN QCS" ~ "Queen Charlotte Sound",
#     survey_abbrev == "SYN WCHG" ~ "West Coast Haida Gwaii",
#     TRUE ~ survey_abbrev
#   )) |>
#   select(-fo_id) |>
#   group_by(survey_name) |>
#   summarise(
#     sampling_frequency = "biennial",
#     year_range = paste(min(year), max(year), sep = "-"),
#     n_years = n_distinct(year),
#     n_years_in_mpa = n_distinct(year[in_mpa == 1]),
#     n_transects = n(),
#     n_transects_in_mpa = sum(in_mpa),
#     prop_in_mpa = round(100 * n_transects_in_mpa / n_transects, 1),
#     n_mpas = n_distinct(uid, na.rm = TRUE)
#   ) |>
#   mutate(across(everything(), as.character)) |>
#   pivot_longer(cols = -c(survey_name), names_to = "metric", values_to = "value") |>
#   pivot_wider(
#     names_from = survey_name,
#     values_from = value,
#     names_vary = "slowest"
#   ) |>
#   select(metric, `West Coast Haida Gwaii`, `Hecate Strait`, `Queen Charlotte Sound`) |>
#   mutate(metric = metric_labels[metric]) |>
#   flextable() |>
#   add_header_row(
#     values = c("Monitoring program name", "Synoptic Bottom Trawl Survey"),
#     colwidths = c(1, 3),
#     top = TRUE
#   ) |>
#   set_header_labels(metric = "Survey region") |>
#   align(i = 1, align = "center", part = "header") |>  # Center top header row
#   align(i = 2, align = "right", part = "header") |>  # Right align column names
#   align(j = 2, align = "center", part = "header") |>  # Right align column names
#   align(j = 2:4, align = "right", part = "body") |>    # Right align data columns
#   autofit()
# syn_flextable
# saveRDS(syn_flextable, file.path(table_dir, "syn-summary-flextable.rds"))

# # SYN aggregate (all surveys combined) -----------------------------------------
# syn_aggregate_flextable <- syn_mpa_df |>
#   select(-fo_id) |>
#   summarise(
#     monitoring_program = "Synoptic Bottom Trawl Survey",
#     sampling_frequency = "biennial",
#     year_range = paste(min(year), max(year), sep = "-"),
#     n_years = n_distinct(year),
#     n_years_in_mpa = n_distinct(year[in_mpa == 1]),
#     n_transects = n(),
#     n_transects_in_mpa = sum(in_mpa),
#     prop_in_mpa = round(100 * n_transects_in_mpa / n_transects, 1),
#     n_mpas = n_distinct(uid, na.rm = TRUE)
#   ) |>
#   mutate(across(everything(), as.character)) |>
#   pivot_longer(everything(), names_to = "metric", values_to = "value") |>
#   mutate(metric = metric_labels[metric]) |>
#   drop_na(metric) |>
#   flextable(col_keys = c("metric", "value")) |>
#   set_header_labels(metric = "Metric", value = "Synoptic Bottom Trawl Survey (all regions)") |>
#   align(align = "right", part = "all") |>
#   autofit()

# saveRDS(syn_aggregate_flextable, file.path(table_dir, "syn-summary-aggregate-flextable.rds"))

# syn_mpa_df |>
#   filter(in_mpa == 1) |>
#   group_by(site, uid) |>
#   summarise(
#     n_transects = length(unique(fe)),
#     first_year = min(year),
#     last_year = max(year),
#     n_years = n_distinct(year)
#   ) |>
#   arrange(desc(n_transects)) |>
#   left_join(gf_ecp_site_summary |> filter(survey_group == "Synoptic")) |>
#   left_join(subregion_uid_lu, by = "uid") |>
#   select(-survey_group) |>
# flextable() |>
#   set_header_labels(
#     uid = "Site UID",site = "Site Name",
#     subregion = "Subregion",
#     n_transects = "Total fishing events",
#     first_year = "First year",
#     last_year = "Last year",
#     n_years = "Number of years within-MPA events",
#     n_species = "E-CP species count") |>
#   colformat_int(j = c("first_year", "last_year", "n_transects", "n_years", "n_species"), big.mark = "") |>
#   align(align = "right", part = "all") |>
#   autofit() |>
# saveRDS(file.path(table_dir, "syn-site-summary.rds"))


# IPHC summary table ------------------------------------------------------------
# iphc_flextable <- iphc_mpa_df |>
#   mutate(survey_name = "IPHC") |>
#   reframe(
#     survey_name = "IPHC",
#     sampling_frequency = "annual",
#     year_range = paste(min(year), max(year), sep = "-"),
#     n_years = n_distinct(year),
#     n_years_in_mpa = n_distinct(year[in_mpa == 1]),
#     n_transects = n(),
#     n_transects_in_mpa = sum(in_mpa),
#     prop_in_mpa = round(100 * n_transects_in_mpa / n_transects, 1),
#     n_mpas = n_distinct(uid, na.rm = TRUE)
#   ) |>
#   mutate(across(everything(), as.character)) |>
#   pivot_longer(cols = -c(survey_name), names_to = "metric", values_to = "value") |>
#   pivot_wider(
#     names_from = survey_name,
#     values_from = value,
#     names_vary = "slowest"
#   ) |>
#   select(metric, IPHC) |>
#   mutate(metric = metric_labels[metric]) |>
#   rename(`International Pacific Halibut Commission Setline Survey` = IPHC) |>
#   flextable() |>
#   set_header_labels(metric = "Monitoring program name") |>
#   align(align = "right", part = "all") |>
#   autofit()

# saveRDS(iphc_flextable, file.path(table_dir, "iphc-summary-flextable.rds"))




# ROV summary table ------------------------------------------------------------
rov_flextable <- rov_mpa_df |>
  mutate(survey_name = "ROV") |>
  reframe(
    survey_name = "ROV",
    sampling_frequency = "annual",
    year_range = paste(min(year), max(year), sep = "-"),
    n_years = n_distinct(year),
    n_years_in_mpa = n_distinct(year[in_mpa == 1]),
    n_transects = n(),
    n_transects_in_mpa = sum(in_mpa),
    prop_in_mpa = round(100 * n_transects_in_mpa / n_transects, 1),
    n_mpas = n_distinct(uid, na.rm = TRUE),
  ) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(cols = -c(survey_name), names_to = "metric", values_to = "value") |>
  pivot_wider(
    names_from = survey_name,
    values_from = value,
    names_vary = "slowest"
  ) |>
  select(metric, ROV) |>
  mutate(metric = metric_labels[metric]) |>
  flextable() |>
  set_header_labels(metric = "") |>
  align(align = "right", part = "all") |>
  autofit()

saveRDS(rov_flextable, file.path(table_dir, "rov-summary-flextable.rds"))
