library(dplyr)

# See https://gitlab.com/dfo-msea/query-dfo-databases/-/blob/master/DataSources/sql/get-multispecies-records.sql?ref_type=heads

#sf_server <- "WDC-SQL2016-P\\SIOSP01"
#db <- "SFBioSQL"

db_connection <- function(server, database, driver = "SQL Server", ...) {
  DBI::dbConnect(
    odbc::odbc(),
    Driver   = driver,
    Server   = server,
    Database = database,
    ...
  )
}


run_sql <- function(database, query, driver = "SQL Server", ...) {
  # Map server keys to connection info
  key <- tolower(database)
  db_list <- list(
    sfbio = list(server = "WDC-SQL2016-P\\SIOSP01", database = "SFBioSQL"),
    gfbio = list(server = "DFBCV9TWVASP003\\SQL2016STD", database = "GFBioSQL")
  )

  if (!key %in% names(db_list)) stop("Invalid database. Use 'sfbio' or 'gfbio'.")

  sql <- paste(query, collapse = "\n")

  con <- db_connection(db_list[[key]]$server, db_list[[key]]$database, driver = driver, ...)
  on.exit(suppressWarnings(DBI::dbDisconnect(con)), add = TRUE)

  tryCatch(
    DBI::dbGetQuery(con, sql),
    error = function(e) {
      stop(sprintf(
        "Query failed on %s/%s: %s\nSQL (first 200 chars): %s",
        db_list[[key]]$server, db_list[[key]]$database,
        conditionMessage(e), substr(sql, 1, 200)
      ), call. = FALSE)
    }
  )
}

d <- run_sql(database = "sfbio",
             query =
            "WITH hdr AS (
               SELECT *
                 FROM vw_ms1_headers
             ),
             -- Header × Algae observations, mapped to SPECIES_CODE via the algae lookup view
             hdr_algae AS (
               SELECT
               h.TRIP_ID,
               h.TransectSite,
               a.ALGAE_CODE,
               alg_lu.SPECIES_CODE,
               'algae' AS source
               FROM hdr AS h
               LEFT JOIN vw_ms3_algae AS a
               ON a.TRIP_ID = h.TRIP_ID
               AND a.TransectSite = h.TransectSite
               LEFT JOIN vw_luALGAE_SPECIES_ms AS alg_lu
               ON alg_lu.ALGAE_CODE = a.ALGAE_CODE
             ),
             -- Header × SF observations (already have SPECIES_CODE)
             hdr_sf AS (
               SELECT
               h.TRIP_ID,
               h.TransectSite,
               CAST(NULL AS VARCHAR(50)) AS ALGAE_CODE,
               sf.SPECIES_CODE,
               'sf' AS source
               FROM hdr AS h
               LEFT JOIN vw_ms3_sf AS sf
               ON sf.TRIP_ID = h.TRIP_ID
               AND sf.TransectSite = h.TransectSite
             ),
             -- Stack both sources
             hdr_obs AS (
               SELECT TRIP_ID, TransectSite, ALGAE_CODE, SPECIES_CODE, source FROM hdr_algae
               UNION ALL
               SELECT TRIP_ID, TransectSite, ALGAE_CODE, SPECIES_CODE, source FROM hdr_sf
             )
             SELECT
             h.*,
             o.source,
             o.ALGAE_CODE,
             o.SPECIES_CODE as obs_SPECIES_CODE,
             sp.*   -- ← all columns from SPECIES
             FROM hdr AS h
             LEFT JOIN hdr_obs AS o
             ON o.TRIP_ID = h.TRIP_ID
             AND o.TransectSite = h.TransectSite
             LEFT JOIN SPECIES AS sp
             ON sp.SPECIES_CODE = o.SPECIES_CODE
             -- Optional exclusions that KEEP headers with no observations:
               -- WHERE (o.SPECIES_CODE IS NULL OR o.SPECIES_CODE <> '11L')
             ORDER BY
             h.TRIP_ID, h.TransectSite, o.source, o.SPECIES_CODE;")

d |>
  janitor::clean_names() |>
  as_tibble() |>
saveRDS(d, file = here::here("data-raw", "multi-species-data-no-abalone.rds"))
# readRDS(d, file = here::here("data-raw", "multi-species-data-no-abalone.rds"))

d0 <- readRDS(here::here("data-raw", "multi-species-data-no-abalone.rds"))
d <- d0 |>
  mutate(across(c(species_common_name, species_science_name, species_desc,
    taxonomic_rank, parent_taxonomic_unit), tolower)
  ) |>
  tidyr::drop_na(species_code) |>
  mutate(species_group = ifelse(is.na(algae_code), "invert", "algae"),
         species = ifelse(is.na(species_common_name), species_desc, species_common_name))
saveRDS(d, file = here::here("data-generated", "multi-species-data-no-abalone-cleaned.rds"))

msd_header <- d |>
  distinct(
    trip_id, survey, year, month, day, transect_site,
    lat_start, lon_start, lat_end, lon_end)

# Algae seem to be presence absence only
msd_algae <- d |>
  filter(source == "algae") |>
  group_by(trip_id, year, month, day, transect_site,
    species, species_code, species_science_name, species_desc, itis_tsn
  ) |>
  summarise(catch_count = n(), .groups = "drop") |>
  right_join(msd_header)

# d_with_cukes <- run_sql(database = "sfbio",
#              query = "select * from vw_ms2b_summarycounts_noab") |>
#   as_tibble()
# saveRDS(d_with_cukes, "multi-species-data-no-abalone-with-cukes-wide.rds")
d_with_cukes <- readRDS(here::here("data-raw", "multi-species-data-no-abalone-with-cukes-wide.rds")) |>
  as_tibble() |>
  tidyr::pivot_longer(cols = -c(TRIP_ID, TransectSite, Quadrat), names_to = "species_code", values_to = "count") |>
  janitor::clean_names()

wide_species_lu <- readRDS(here::here("data-raw", "species-table.rds")) |>
  as_tibble() |>
  janitor::clean_names() |>
  # select(species_code, species_science_name, species_common_name) |>
  mutate(across(c(species_science_name, species_common_name), tolower))

msd_inverts <- d_with_cukes |>
  group_by(trip_id, transect_site, species_code) |>
  summarise(catch_count = sum(count), .groups = "drop") |>
  left_join(wide_species_lu) |>
  right_join(msd_header)

# aggregated by transect
msd_catch <- bind_rows(
  msd_algae |> mutate(species_group = "algae"),
  msd_inverts |> mutate(species_group = "invert")
) |>
  mutate(across(c(species_desc, parent_taxonomic_unit), tolower)) |>
  filter(!is.na(species_code))


# test <- msd_catch |> filter(is.na(species_code)) |> pull(transect_site)
# test2 <- filter(d, transect_site %in% test) |> filter(!is.na(species_code)) |>
#   distinct(transect_site) |>
#   pull(transect_site)
# setdiff(test, test2)

msd_catch |>
  saveRDS(here::here("data-generated", "msd-catch.rds"))

# # Algae seem to be presence absence only
# msd_counts <- d |>
#   group_by(trip_id, year, month, day, transect_site,
#     species, species_code, species_science_name, species_desc, itis_tsn
#   ) |>
#   summarise(catch_count = n(), .groups = "drop")



msd_transect_richness <- msd_catch |>
  group_by(trip_id, year, month, day, transect_site) |>
  summarise(richness = n_distinct(species_code),
   u_species = paste(unique(species), collapse = ", "), .groups = "drop")
saveRDS(msd_transect_richness, here::here("data-generated", "msd-transect-richness.rds"))
# msd_catch <- left_join(msd_counts, msd_header) |>
#   left_join(msd_transect_richness)

# msd_catch |>
#   saveRDS(here::here("data-generated", "msd-catch.rds"))


# msd <- readRDS(here::here("data-generated", "msd-catch.rds"))
# # A few urchins aren't identified to species
# filter(msd, species_science_name == "strongylocentrotidae")




# msd_catch |>
#   filter(species_code %in% unique(test$species_code)) |>
#   select(trip_id, year, transect_site, species, species_code, catch_count) |>
#   left_join(test) |>
#   mutate(diff = count_by_transect - catch_count) |> View()


# filter(d, transect_site == 6125) |>
# filter(source != "algae") |>
# select(transect_site, quadrat, species_code, species, comments_inverts) |>
# View()
