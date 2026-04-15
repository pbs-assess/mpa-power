library(dplyr)

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

  ssid_lu <- gf_ecp_data0 |>
    distinct(survey_series_id, survey_abbrev)

  # in_mpa_lu <- bind_rows(syn_mpa_df, hbll_mpa_df) |>
  hbll_mpa_df <- readRDS(file.path("data-generated", "overlays", "hbll-mpa-sf.rds")) |> sf::st_drop_geometry()
  in_mpa_lu <- hbll_mpa_df |>
    distinct(survey_series_id = ssid, survey_abbrev, fishing_event_id = fe, in_mpa, uid, site)

  gf_ecp_data <- gf_ecp_data0 |>
    select(survey_series_id, species_common_name, species_science_name, species_code,
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
    tidyr::drop_na(in_mpa)
    saveRDS(gf_ecp_data, file.path("data-generated", "hbll-ecp-species.rds"))
  } else {
    gf_ecp_data <- readRDS(file.path("data-generated", "hbll-ecp-species.rds"))
}