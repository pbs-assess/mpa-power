# =============================================================================
# Fit conditioning models
# =============================================================================

# Load setup and functions
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

library(tidyr)
library(patchwork)
library(digest)

# Setup directories
fit_dir <- here::here("data-generated", "fits")
fig_dir <- here::here("draft-figures", "diagnostics")
sim_cache <- here::here("data-generated", "sim-cache")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sim_cache, recursive = TRUE, showWarnings = FALSE)

# Zero counts
get_zero_counts <- function(fit, historical, seed = NULL, zero_text = FALSE) {
  df <- plot_d(fit, historical, type = "ts", seed = seed, check_sanity = TRUE) |>
    group_by(d) |>
    summarise(species = unique(species_common_name),
      survey_abbrev = unique(survey_abbrev),
      mean_catch_count = mean(catch_count),
      zero_count = sum(catch_count == 0),
      sanity = unique(sanity)) |>
    select(species, survey_abbrev, d, mean_catch_count, zero_count, sanity)

  if (zero_text) {
    df |>
      select(d, mean_catch_count, zero_count) |>
      pivot_wider(names_from = d,
                  values_from = c(mean_catch_count, zero_count),
                  names_sep = "_") |>
      mutate(
        mean_catch_text = paste0("Mean catch ", round(mean_catch_count_historical, 1), " / ", round(mean_catch_count_simulated, 1)),
        zero_text = paste0("Zeroes: ", zero_count_historical, " / ", zero_count_simulated)
      ) |>
      select(mean_catch_text, zero_text)
  } else {
    return(df)
  }
}

clean_family_name <- function(fit) {
    if (!is.null(fit$family$delta)) {
      # Delta model - sdmTMB style: "Delta family1(link1), family2(link2)"
      fam1 <- fit$family$family[[1]]
      fam2 <- fit$family$family[[2]]
      link1 <- fit$family$link[[1]]
      link2 <- fit$family$link[[2]]
      paste0("Delta ", fam1, "(", link1, "), ", fam2, "(", link2, ")")
    } else {
      # Regular model - sdmTMB style: "family(link)"
      fam <- fit$family$family
      link <- fit$family$link
      paste0(fam, "(", link, ")")
    }
  }

mk_title <- function(fit, check_sanity = TRUE) {
  survey_name <- unique(fit$data$survey_abbrev)
  family_name <- clean_family_name(fit)
  formula_str <- deparse(fit$formula[[1]])
  spatiotemporal <- fit$spatiotemporal
  species <- unique(fit$data$species_common_name)
  sanity_str <- ifelse(check_sanity, summarise_sanity(fit), "")
  aic <- round(AIC(fit))

  paste0(survey_name, " - ", species, " - ", sanity_str,
        "\n-   f: ", formula_str,
        "\n-   st: ", spatiotemporal,
        "\n-   ", family_name,
        "\n-   aic: ", aic)
}

make_timeseries <- function(fit, historical, nsim = 1, seed = NULL,
  append_years = TRUE, model = NA, cache_dir = NULL) {

  # Simple file-based caching
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

    # Create cache key from key parameters
    cache_key <- digest::digest(list(
      model_hash = digest::digest(list(fit$model$par, AIC(fit))),
      nsim = nsim,
      seed = seed,
      model = model
    ))

    cache_file <- file.path(cache_dir, paste0("timeseries_", cache_key, ".rds"))

    if (file.exists(cache_file)) {
      message("Loading cached timeseries data")
      return(readRDS(cache_file))
    }
  }

  .survey_abbrev <- unique(fit$data$survey_abbrev)
  sim_mat <- simulate(fit, nsim = nsim, type = "mle-eb", seed = seed, model = model)

  sim_df <- as.data.frame(sim_mat) |>
    mutate(row_id = row_number()) |>
    pivot_longer(cols = -row_id, names_to = "replicate", values_to = "catch_count") |>
    mutate(replicate = as.numeric(gsub("V", "", replicate)),
           d = "simulated",
           survey_abbrev = .survey_abbrev) |>
    arrange(replicate)

  # Process historical data for this survey
  hist_data <- historical |>
    filter(survey_abbrev == .survey_abbrev) |>
    mutate(d = "historical",
           row_id = row_number(),
           replicate = 0)

  sim_data <- left_join(sim_df, hist_data |> select(-catch_count, -replicate, -d))

  if ("Beta" %in% family(fit)$family) {
    sim_data <- sim_data |>
      mutate(catch_count = catch_count * hook_count)
  }

  if (append_years) {
    sim_data <- sim_data |> mutate(year = year - min(year) + 2 + max(year))
  }

  # Combine
  result <- bind_rows(sim_data, hist_data)

  # Save to cache if caching enabled
  if (!is.null(cache_dir)) {
    saveRDS(result, cache_file)
    message("Saved timeseries data to cache")
  }

  result
}


# Create dynamic title from model components
plot_d <- function(fit, historical, nsim = 1, type = "sim",
                   append_years = TRUE, seed = NULL, check_sanity = TRUE, cache_dir = NULL) {
  type <- match.arg(type, c("ts", "sim", "cumul-rank", "dharma"))

  comb <- make_timeseries(fit, historical = historical, nsim = nsim, append_years = append_years, seed = seed, cache_dir = cache_dir) |>
    mutate(aic = round(AIC(fit)))
  if (check_sanity) comb$sanity <- summarise_sanity(fit)

  if (type == "ts") {
    return(comb)
  }

  zero_text <- get_zero_counts(fit, historical, seed = seed, zero_text = TRUE)

  title <- mk_title(fit, check_sanity)

  if (type == "sim") {
    p <- comb |>
      ggplot() +
      geom_jitter(aes(x = year, y = catch_count, color = d), shape = 21) +
      geom_text(data = zero_text, aes(x = Inf, y = Inf, label = mean_catch_text),
        hjust = 1.1, vjust = -2) +
      geom_text(data = zero_text, aes(x = Inf, y = Inf, label = zero_text),
        hjust = 1.1, vjust = -4) +
      coord_cartesian(clip = "off") +
      scale_colour_manual(values = c("steelblue", "orange")) +
      labs(title = title)
    return(p)
  }

  if (type == "cumul-rank") {
    plot_data <- comb |>
      group_by(d, replicate) |>
      arrange(catch_count) |>
      mutate(cumul_catch = cumsum(catch_count),
             rank_obs = rank(catch_count))
    p <- plot_data |>
      ggplot(aes(x = rank_obs, y = cumul_catch)) +
      # geom_line(aes(color = d, group = interaction(d, replicate))) +
      geom_line(data = filter(plot_data, d == "historical"),
                colour = "steelblue", linewidth = 1.2) +
      geom_line(data = filter(plot_data, d == "simulated"),
                colour = "orange", alpha = 0.5) +
      geom_text(data = zero_text, aes(x = Inf, y = Inf, label = mean_catch_text),
        hjust = 1.1, vjust = -2) +
      geom_text(data = zero_text, aes(x = Inf, y = Inf, label = zero_text),
        hjust = 1.1, vjust = -4) +
      coord_cartesian(clip = "off") +
      scale_colour_manual(values = c("steelblue", "orange")) +
      # geom_line(data = filter(plot_data, d == "simulated"),
      #           colour = "black", alpha = 0.3) +
      # geom_line(data = filter(plot_data, d == "historical"),
      #           colour = "red", linewidth = 1.2) +
      labs(x = "Rank (ordered observations)", y = "Cumulative catch", title = title)
    return(p)
  }

  if (type == "dharma") {
    s_nb2 <- simulate(fit, nsim = 100, type = "mle-mvn")
    r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
    plot(r_nb2)
    rect(par("usr")[1], par("usr")[4] + 0.1 * diff(par("usr")[3:4]),  # bottom edge higher
        par("usr")[2], par("usr")[4] + 0.18 * diff(par("usr")[3:4]),  # top edge higher
        col = "white", border = NA, xpd = TRUE)
    title(main = title, line = 2, cex.main = 1.2, font = 2)
  }
}

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
# Grids
hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") #|>
  # depth not currently used in sim
  # mutate(depth_mean = mean(log(depth_m), na.rm = TRUE),
  #        depth_sd = sd(log(depth_m), na.rm = TRUE),
  #        depth_scaled = (log(depth_m) - depth_mean[1]) / depth_sd[1])

hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S")

hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
comm_ll_activity_status <- readRDS(here::here("data-generated", "spatial", "comm-ll-draft-activity-status.rds"))
mpa_shape_simplified <- comm_ll_activity_status |> st_simplify(dTolerance = 100)

# HBLL species data
sp <- sp_to_hyphens("yelloweye rockfish")
sp <- sp_to_hyphens("north pacific spiny dogfish")
sp <- sp_to_hyphens("lingcod")
sp <- sp_to_hyphens("quillback rockfish")

devtools::load_all("~/R_DFO/sdmTMB")

source(here::here("R", "caching-functions.R"))

sp_list <- c("yelloweye rockfish", "north pacific spiny dogfish", "lingcod", "quillback rockfish")
for (sp in sp_to_hyphens(sp_list)) {
message(paste0("Simulation diagnostics for ", sp))
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
  prep_hbll_data(bait_counts = bait_counts) |>
  mutate(
    obs_id = factor(row_number()),
    catch_prop = catch_count / hook_count,
    log_depth = log(depth_m)
  )

historical <- sp_dat |>
  mutate(x = X * 1000, y = Y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 3156) |>
  st_join(comm_ll_activity_status |> st_transform(crs = 3156), join = st_within) |>
  mutate(activity_status_label = if_else(is.na(activity_status_label), "outside", activity_status_label)) |>
  mutate(restricted = ifelse(activity_status_label == "outside", 0, 1)) |>
  st_join(hbll_grid_poly |> select(block_id, grouping_code) |> st_transform(crs = 3156), join = st_within) |>
  st_drop_geometry() |>
  select(ssid, survey_abbrev, year, species_common_name, fishing_event_id, latitude, longitude, X, Y,
    block_id, fe_grouping_code = grouping_code.x, grouping_code = grouping_code.y,
    depth_m, offset, hook_count,
    catch_count, restricted) |>
  mutate(after = 0) |>
  left_join(hbll_allocations,
    by = c("survey_abbrev", "grouping_code", "ssid" = "survey_series_id"))

# sp_dat |> pull(depth_m) |> summary()

# combined <- st_intersection(
#     st_as_sfc(st_bbox(hbll_grid_poly |> st_transform(crs = st_crs(mpa_shape_simplified)))),
#     st_as_sfc(st_bbox(mpa_shape_simplified))
#   ) |>
#   st_as_sf()

# plot_limits_combined <- get_plot_limits(combined, buffer = 1000)


# Set up fit directory for caching
fit_dir <- here::here("data-generated", "fits")
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)

# Prepare data and meshes
d_IN <- sp_dat |> filter(survey_abbrev == "HBLL INS N")
d_IN$weights <- d_IN$hook_count / mean(d_IN$hook_count)
mesh_IN <- local(make_mesh(d_IN, xy_cols = c("X", "Y"), cutoff = 10))
d_OS <- sp_dat |> filter(survey_abbrev == "HBLL OUT S")
d_OS$weights <- d_OS$hook_count / mean(d_OS$hook_count)
mesh_OS <- local(make_mesh(d_OS, xy_cols = c("X", "Y"), cutoff = 10))
d_ON <- sp_dat |> filter(survey_abbrev == "HBLL OUT N")
d_ON$weights <- d_ON$hook_count / mean(d_ON$hook_count)
mesh_ON <- local(make_mesh(d_ON, xy_cols = c("X", "Y"), cutoff = 10))

check_cache <- TRUE
silent <- TRUE
link <- "cloglog"
# link <- "logit"

# Beta binomial ----------------------------------------------------------------
fit_ON_bb <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-beta-binomial-off-", link),
  fit_dir = fit_dir,
  data = d_ON,
  formula = catch_prop ~ 0 + fyear,
  mesh = mesh_ON,
  family = betabinomial(link = "cloglog"),
  spatial = "off",
  spatiotemporal = "off",
  anisotropy = FALSE,
  weights = d_ON$hook_count,
  check_cache = check_cache,
  silent = silent
)

fit_OS_bb <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-beta-binomial-off-", link),
  fit_dir = fit_dir,
  data = d_OS,
  formula = catch_prop ~ 0 + fyear,
  mesh = mesh_OS,
  family = betabinomial(link = "cloglog"),
  spatial = "off",
  spatiotemporal = "off",
  anisotropy = FALSE,
  weights = d_OS$hook_count,
  check_cache = check_cache,
  silent = silent
)

fit_IN_bb <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-beta-binomial-off-", link),
  fit_dir = fit_dir,
  data = d_IN,
  formula = catch_prop ~ 0 + fyear,
  mesh = mesh_IN,
  family = betabinomial(link = "cloglog"),
  spatial = "off",
  spatiotemporal = "off",
  anisotropy = FALSE,
  weights = d_IN$hook_count,
  check_cache = check_cache,
  silent = silent
)

fit_ON_bb_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-beta-binomial-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_ON_bb,
  formula = catch_prop ~ 0 + fyear,
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  silent = silent,
  check_cache = check_cache
)

fit_OS_bb_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-beta-binomial-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_OS_bb,
  formula = catch_prop ~ 0 + fyear,
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  silent = silent,
  check_cache = check_cache
)

fit_IN_bb_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-beta-binomial-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_IN_bb,
  formula = catch_prop ~ 0 + fyear,
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  silent = silent,
  check_cache = check_cache
)

fit_ON_bb_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-beta-binomial-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_ON_bb,
  formula = catch_prop ~ 0 + fyear,
  spatial = "on",
  spatiotemporal = "ar1",
  time = "year",
  anisotropy = FALSE,
  silent = silent,
  check_cache = check_cache
)

fit_OS_bb_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-beta-binomial-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_OS_bb,
  formula = catch_prop ~ 0 + fyear,
  spatial = "on",
  spatiotemporal = "ar1",
  time = "year",
  anisotropy = FALSE,
  silent = FALSE,
  check_cache = check_cache
)

fit_IN_bb_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-beta-binomial-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_IN_bb,
  formula = catch_prop ~ 0 + fyear,
  spatial = "on",
  spatiotemporal = "ar1",
  time = "year",
  anisotropy = FALSE,
  silent = FALSE,
  check_cache = check_cache
)


# Binomial ------------------------------------------------------------------
fit_ON_bin <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-bin-st-off-", link),
  fit_dir = fit_dir,
  data = d_ON,
  formula = catch_prop ~ 0 + fyear + (1 | obs_id),
  mesh = mesh_ON,
  family = binomial(link = "cloglog"),
  spatial = "on",
  anisotropy = FALSE,
  spatiotemporal = "off",
  weights = d_ON$hook_count,
  check_cache = check_cache,
  silent = silent
)

fit_OS_bin <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-bin-st-off-", link),
  fit_dir = fit_dir,
  data = d_OS,
  formula = catch_prop ~ 0 + fyear + (1 | obs_id),
  mesh = mesh_OS,
  family = binomial(link = "cloglog"),
  spatial = "on",
  anisotropy = FALSE,
  spatiotemporal = "off",
  weights = d_OS$hook_count,
  check_cache = check_cache,
  silent = silent
)

fit_IN_bin <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-bin-st-off-", link),
  fit_dir = fit_dir,
  data = d_IN,
  formula = catch_prop ~ 0 + fyear + (1 | obs_id),
  mesh = mesh_IN,
  family = binomial(link = "cloglog"),
  spatial = "on",
  anisotropy = FALSE,
  spatiotemporal = "off",
  weights = d_IN$hook_count,
  check_cache = check_cache,
  silent = silent
)

# IID -------------
fit_ON_bin_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-binomial-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_ON_bin,
  formula = catch_prop ~ 0 + fyear,
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  share = FALSE,
  check_cache = check_cache,
  silent = silent
)

fit_OS_bin_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-binomial-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_OS_bin,
  formula = catch_prop ~ 0 + fyear,
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  share = FALSE,
  check_cache = check_cache,
  silent = silent
)

fit_IN_bin_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-binomial-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_IN_bin,
  formula = catch_prop ~ 0 + fyear,
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  share = FALSE,
  check_cache = check_cache,
  silent = silent
)

# AR1 -------------
fit_ON_bin_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-binomial-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_ON_bin,  # Reference to original model,
  formula = catch_prop ~ 0 + fyear,
  spatiotemporal = "ar1",
  anisotropy = FALSE,
  time = "year",
  extra_time = sdmTMB:::find_missing_time(fit_ON_bin$data$year),
  share = FALSE,
  check_cache = check_cache,
  silent = silent
)

fit_OS_bin_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-binomial-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_OS_bin,  # Reference to original model,
  formula = catch_prop ~ 0 + fyear,
  spatiotemporal = "ar1",
  anisotropy = FALSE,
  time = "year",
  extra_time = sdmTMB:::find_missing_time(fit_OS_bin$data$year),
  share = FALSE,
  check_cache = check_cache,
  silent = silent
)

fit_IN_bin_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-binomial-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_IN_bin,  # Reference to original model
  formula = catch_prop ~ 0 + fyear,
  spatiotemporal = "ar1",
  anisotropy = FALSE,
  time = "year",
  extra_time = sdmTMB:::find_missing_time(fit_IN_bin$data$year),
  share = FALSE,
  check_cache = check_cache,
  silent = silent
)

# Delta beta ------------------------------------------------------------------
fit_ON_db <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-delta-beta-off-", link),
  fit_dir = fit_dir,
  formula = catch_prop ~ 0 + fyear,
  data = d_ON,
  mesh = mesh_ON,
  family = delta_beta(link1 = "cloglog", link2 = "cloglog"),
  spatial = "on",
  spatiotemporal = "off",
  weights = d_ON$weights
)

fit_OS_db <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-delta-beta-off-", link),
  fit_dir = fit_dir,
  formula = catch_prop ~ 0 + fyear,
  data = d_OS,
  mesh = mesh_OS,
  family = delta_beta(link1 = "cloglog", link2 = "cloglog"),
  spatial = "on",
  spatiotemporal = "off",
  time = "year",
  weights = d_OS$weights
)

fit_IN_db <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-delta-beta-off-", link),
  fit_dir = fit_dir,
  formula = catch_prop ~ 0 + fyear,
  data = d_IN,
  mesh = mesh_IN,
  family = delta_beta(link1 = "cloglog", link2 = "cloglog"),
  spatial = "on",
  spatiotemporal = "off",
  weights = d_IN$weights
)

# IID -------------
fit_ON_db_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-delta-beta-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_ON_db,
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  share = TRUE,
  check_cache = check_cache
)

fit_OS_db_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-delta-beta-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_OS_db,
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  share = TRUE,
  check_cache = check_cache
)

fit_IN_db_iid <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-delta-beta-iid-", link),
  fit_dir = fit_dir,
  update_from = fit_IN_db,
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  share = TRUE,
  check_cache = check_cache
)

# AR1 -------------
fit_IN_db_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-INS-N-delta-beta-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_IN_db,
  spatiotemporal = "ar1",
  time = "year",
  extra_time = sdmTMB:::find_missing_time(fit_IN_db$data$year),
  anisotropy = FALSE,
  share = TRUE,
  check_cache = check_cache
)

fit_ON_db_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-N-delta-beta-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_ON_db,
  spatiotemporal = "ar1",
  time = "year",
  extra_time = sdmTMB:::find_missing_time(fit_ON_db$data$year),
  anisotropy = FALSE,
  share = TRUE,
  check_cache = check_cache
)

fit_OS_db_ar1 <- fit_cached_sdmTMB(
  model_tag = paste0(sp, "-HBLL-OUT-S-delta-beta-ar1-", link),
  fit_dir = fit_dir,
  update_from = fit_OS_db,
  spatiotemporal = "ar1",
  time = "year",
  extra_time = sdmTMB:::find_missing_time(fit_OS_db$data$year),
  anisotropy = FALSE,
  share = TRUE,
  check_cache = check_cache
)

# check_sanity_fits(list(fit_IN_bin_iid, fit_ON_bin_iid, fit_OS_bin_iid))
# check_sanity_fits(list(fit_IN_bin_ar1, fit_ON_bin_ar1, fit_OS_bin_ar1))

# AIC(fit_IN_bin, fit_IN_bin_iid, fit_IN_bin_ar1)
# AIC(fit_ON_bin, fit_ON_bin_iid, fit_ON_bin_ar1)
# AIC(fit_OS_bin, fit_OS_bin_iid, fit_OS_bin_ar1)

seed <- 42

# bind_rows(get_zero_counts(fit_IN_bin, historical, seed = seed),
#   get_zero_counts(fit_ON_bin, historical, seed = seed),
#   get_zero_counts(fit_OS_bin, historical, seed)) |>
#   mutate(model = "1 | obs_id")

# bind_rows(get_zero_counts(fit_IN_bin_iid, historical, seed = seed),
#   get_zero_counts(fit_ON_bin_iid, historical, seed = seed),
#   get_zero_counts(fit_OS_bin_iid, historical, seed = seed)) |>
#   mutate(model = "iid")

# Create plotting function for 3x3 grid
create_model_plots <- function(models, historical, plot_type = "sim", filename_suffix,
  nsim = NULL, save_plot = TRUE, width = 15.3, height = 13, cache_dir = NULL) {

  # Build arguments for plot_d
  plot_args <- list(historical = historical, type = plot_type)
  if (!is.null(cache_dir)) {
    plot_args$cache_dir <- cache_dir
  }
  if (!is.null(nsim)) {
    plot_args$nsim <- nsim
  }

  plots <- purrr::map(models, ~do.call(plot_d, c(list(.x), plot_args)))

  # Arrange in 3x3 grid
  combined_plot <- wrap_plots(plots, ncol = 3, guides = "collect")

  # Save plot
  if (save_plot) {
  ggsave(here::here(fig_dir, paste0(sp, "-", filename_suffix, "-", plot_type, ".pdf")),
    plot = combined_plot, width = width, height = height)
  }
  return(combined_plot)
}

# Define model sets
bb_models <- list(fit_IN_bb, fit_ON_bb, fit_OS_bb,
                  fit_IN_bb_iid, fit_ON_bb_iid, fit_OS_bb_iid,
                  fit_IN_bb_ar1, fit_ON_bb_ar1, fit_OS_bb_ar1)

bin_models <- list(fit_IN_bin, fit_ON_bin, fit_OS_bin,
                   fit_IN_bin_iid, fit_ON_bin_iid, fit_OS_bin_iid,
                   fit_IN_bin_ar1, fit_ON_bin_ar1, fit_OS_bin_ar1)

db_models <- list(fit_IN_db, fit_ON_db, fit_OS_db,
                  fit_IN_db_iid, fit_ON_db_iid, fit_OS_db_iid,
                  fit_IN_db_ar1, fit_ON_db_ar1, fit_OS_db_ar1)

# AIC summary:
# model_names <- c("IN", "ON", "OS", "IN_iid", "ON_iid", "OS_iid", "IN_ar1", "ON_ar1", "OS_ar1")
# # Combine into a data frame
# aic_summary <- data.frame(
#   model = model_names,
#   beta_binomial = purrr::map_dbl(bb_models, AIC),
#   binomial = purrr::map_dbl(bin_models, AIC),
#   delta_beta = purrr::map_dbl(db_models, AIC)
# )
# aic_summary

# Time series simulation plots
create_model_plots(bb_models, historical, plot_type = "sim", filename_suffix = "beta-binomial", cache_dir = sim_cache)
create_model_plots(bin_models, historical, plot_type = "sim", filename_suffix = "binomial", cache_dir = sim_cache)
create_model_plots(db_models, historical, plot_type = "sim", filename_suffix = "delta-beta", cache_dir = sim_cache)

create_model_plots(
  list(fit_IN_bin_ar1, fit_ON_bin_ar1, fit_OS_bin_ar1,
    fit_IN_bb_ar1, fit_ON_bb_ar1, fit_OS_bb_ar1,
    fit_IN_db_ar1, fit_ON_db_ar1, fit_OS_db_ar1),
  historical, nsim = 1, plot_type = "sim",
  filename_suffix = "ts-family-comparison", cache_dir = sim_cache)

create_model_plots(
  list(fit_IN_bin_ar1, fit_ON_bin_ar1, fit_OS_bin_ar1,
    fit_IN_bb_ar1, fit_ON_bb_ar1, fit_OS_bb_ar1,
    fit_IN_db_ar1, fit_ON_db_ar1, fit_OS_db_ar1),
  historical, nsim = 10, plot_type = "cumul-rank",
  filename_suffix = "cr-family-comparison", cache_dir = sim_cache)


create_model_plots(bb_models, historical, nsim = 10, plot_type = "cumul-rank",
  filename_suffix = "beta-binomial", cache_dir = sim_cache)
create_model_plots(bin_models, historical, nsim = 10, plot_type = "cumul-rank",
  filename_suffix = "binomial", cache_dir = sim_cache)
create_model_plots(db_models, historical, nsim = 10, plot_type = "cumul-rank",
  filename_suffix = "delta-beta", cache_dir = sim_cache)

}

# DHARMa plots don't work with current functions and patchwork.
meep()# DHARMa plots - 2x3 layout

par(mfrow = c(2, 3))
plot_d(fit_IN_bin, historical, type = "dharma")
plot_d(fit_ON_bin, historical, type = "dharma")
plot_d(fit_OS_bin, historical, type = "dharma")
plot_d(fit_IN_bin_iid, historical, type = "dharma")
plot_d(fit_ON_bin_iid, historical, type = "dharma")
plot_d(fit_OS_bin_iid, historical, type = "dharma")
plot_d(fit_IN_bin_ar1, historical, type = "dharma")
plot_d(fit_ON_bin_ar1, historical, type = "dharma")
plot_d(fit_OS_bin_ar1, historical, type = "dharma")

# TODO: clean up this code (mostly scratch now) and maybe add the nb2 and delta_nb2 models
# to the diagnostics supplement above for completeness.

# Fit conditioning models
# ------------------------------------------------------------
fit_all_surveys <- function(
  .dat = sp_dat,
  .formula = catch_count ~ 0 + fyear,
  .family = nbinom2(link = "log"),
  .species = sp,
  .tag = "nb2",
  .check_cache = TRUE,
  .mesh_cutoff = 10,
  .fit_dir = fit_dir,
  .spatiotemporal = list("OS" = "iid", "ON" = "iid", "IN" = "off"),
  .anisotropy = list("OS" = FALSE, "ON" = FALSE, "IN" = FALSE),
  .offset = NULL,
  .weights = NULL,
  .use_extra_time = FALSE,
  .surveys = c("OS", "ON", "IN") # not sure if this feature is working
) {
  stopifnot(length(.spatiotemporal) == length(.surveys))

  if ("OS" %in% .surveys) {
    fit_OS <- fit_hbll(dat = .dat,
      survey_type = "HBLL OUT S",
      formula = .formula,
      family = .family,
      species = .species,
      spatiotemporal = .spatiotemporal[["OS"]],
      use_extra_time = .use_extra_time,
      time = "year",
      fit_dir = .fit_dir,
      tag = .tag,
      check_cache = .check_cache,
      mesh_cutoff = .mesh_cutoff,
      offset = .offset,
      weights = .weights,
      anisotropy = .anisotropy[["OS"]]
    )
  }

  if ("ON" %in% .surveys) {
  fit_ON <- fit_hbll(dat = .dat,
    survey_type = "HBLL OUT N",
    formula = .formula,
    family = .family,
    species = .species,
    spatiotemporal = .spatiotemporal[["ON"]],
    use_extra_time = .use_extra_time,
    time = "year",
    fit_dir = .fit_dir,
    tag = .tag,
    check_cache = .check_cache,
    mesh_cutoff = .mesh_cutoff,
    offset = .offset,
    weights = .weights,
    anisotropy = .anisotropy[["ON"]]
    )
  }

  if ("IN" %in% .surveys) {
    fit_IN <- fit_hbll(dat = .dat, # didn't converge with spatiotemporal = "iid"
    survey_type = "HBLL INS N",
    formula = .formula,
    family = .family,
    species = .species,
    spatiotemporal = .spatiotemporal[["IN"]],
    use_extra_time = .use_extra_time,
    time = "year",
    fit_dir = .fit_dir,
    tag = .tag,
    check_cache = .check_cache,
    mesh_cutoff = .mesh_cutoff,
    offset = .offset,
    weights = .weights,
    anisotropy = .anisotropy[["IN"]]
    )
  }

  list(fit_OS = fit_OS, fit_ON = fit_ON, fit_IN = fit_IN)
}

check_sanity_fits <- function(fits) {
  lapply(fits, function(fit) {
    s <- sanity(fit)
    unlist(s)
  })
}
stop()
# Binomial ------------------------------------------------------------------

fits_binomial <- fit_all_surveys(
  .formula = catch_prop ~ 0 + fyear + (1 | obs_id),
  .family = binomial(link = "logit"),
  .spatiotemporal = list("OS" = "off", "ON" = "off", "IN" = "off"),
  .anisotropy = list("OS" = FALSE, "ON" = FALSE, "IN" = FALSE),
  .tag = "binomial",
  .check_cache = TRUE,
  .offset = NULL,
  .weights = "hook_count"
)
check_sanity_fits(fits_binomial)
fit_OS <- fits_binomial$fit_OS; fit_ON <- fits_binomial$fit_ON; fit_IN <- fits_binomial$fit_IN


# NB2 ------------------------------------------------------------------------
fits_nb2 <- fit_all_surveys(
  .formula = catch_count ~ 0 + fyear,
  .family = nbinom2(link = "log"),
  .spatiotemporal = list("OS" = "iid", "ON" = "iid", "IN" = "off"),
  .tag = "nb2",
  .check_cache = TRUE
)
fit_OS <- fits_nb2$fit_OS; fit_ON <- fits_nb2$fit_ON; fit_IN <- fits_nb2$fit_IN
check_sanity_fits(fits_nb2)
stop()

# Depth
fits_depth_nb2 <- fit_all_surveys(
  .formula = catch_count ~ 0 + fyear + poly(log_depth, 2),
  .family = nbinom2(link = "log"),
  .tag = "depth-nb2",
  .check_cache = TRUE
)
check_sanity_fits(fits_depth_nb2)
fit_OS <- fits_depth_nb2$fit_OS; fit_ON <- fits_depth_nb2$fit_ON; fit_IN <- fits_depth_nb2$fit_IN

# Delta NB2 ------------------------------------------------------------------
fits_delta_nb2 <- fit_all_surveys(
  .formula = catch_count ~ 0 + fyear,
  .family = delta_truncated_nbinom2(link1 = "logit", link2 = "log"),
  .tag = "delta-nb2",
  .check_cache = TRUE,
  .spatiotemporal = list("IN" = "off", "ON" = "iid", "OS" = "off")
)
check_sanity_fits(fits_delta_nb2)
fit_OS <- fits_delta_nb2$fit_OS; fit_ON <- fits_delta_nb2$fit_ON; fit_IN <- fits_delta_nb2$fit_IN

# Depth
fits_depth_delta_nb2 <- fit_all_surveys(
  .formula = catch_count ~ 0 + fyear + poly(log_depth, 2),
  .family = delta_truncated_nbinom2(link1 = "logit", link2 = "log"),
  .tag = "depth-delta-nb2",
  .check_cache = TRUE,
  .spatiotemporal = list("IN" = "off", "ON" = "iid", "OS" = "off")
)
check_sanity_fits(fits_depth_delta_nb2)
# fit_OS <- fits_depth_delta_nb2$fit_OS; fit_ON <- fits_depth_delta_nb2$fit_ON; fit_IN <- fits_depth_delta_nb2$fit_IN

meep()

stop()
# TODO: add sanity checks
# TODO: evaluate and compare conditioning models: see - https://github.com/mis-assess/shrimp_surveydesign_csas/blob/794abdf0d4657dff5ed3316fe876b58afab0dd83/Reproducible_Examples/coastwide-density.R#L157

OS <- list(
  "binomial" = fits_binomial$fit_OS,
  "nb2" = fits_nb2$fit_OS,
  "depth-nb2" = fits_depth_nb2$fit_OS,
  "delta-nb2" = fits_delta_nb2$fit_OS,
  "depth-delta-nb2" = fits_depth_delta_nb2$fit_OS
)
ON <- list(
  "binomial" = fits_binomial$fit_ON,
  "nb2" = fits_nb2$fit_ON,
  "depth-nb2" = fits_depth_nb2$fit_ON,
  "delta-nb2" = fits_delta_nb2$fit_ON,
  "depth-delta-nb2" = fits_depth_delta_nb2$fit_ON
)
IN <- list(
  "binomial" = fits_binomial$fit_IN,
  "nb2" = fits_nb2$fit_IN,
  "depth-nb2" = fits_depth_nb2$fit_IN,
  "delta-nb2" = fits_delta_nb2$fit_IN,
  "depth-delta-nb2" = fits_depth_delta_nb2$fit_IN
)

# Create AIC comparison tables
data.frame(
  model = c("binomial", "nb2", "depth-nb2", "delta-nb2", "depth-delta-nb2"),
  AIC = sapply(OS, AIC)
)
data.frame(
  model = c("nb2", "depth-nb2", "delta-nb2", "depth-delta-nb2"),
  AIC = sapply(OS, AIC)
)
data.frame(
  model = c("binomial", "nb2", "depth-nb2", "delta-nb2", "depth-delta-nb2"),
  AIC = sapply(ON, AIC)
)
data.frame(
  model = c("nb2", "depth-nb2", "delta-nb2", "depth-delta-nb2"),
  AIC = sapply(ON, AIC)
)
data.frame(
  model = c("nb2", "depth-nb2", "delta-nb2", "depth-delta-nb2"),
  AIC = sapply(IN, AIC)
)


get_preds <- function(fits) {
  purrr::imap_dfr(fits, function(fit, model_name) {
  fam <- family(fit)$family
  message("Fitted family for ", model_name, ": ", paste(fam, collapse = ", "))
  pred <- predict(fit, type = "link", se_fit = FALSE, re_form = NULL)
  if (length(fam) == 1) {
    pred$is_delta <- FALSE
    pred$catch_est <- exp(pred$est)
  } else {
    pred$is_delta <- TRUE
    pred$encounter_est <- plogis(pred$est1)
    pred$positive_est <- exp(pred$est2)
    pred$catch_est <- pred$encounter_est * pred$positive_est * exp(pred$offset)
  }
  pred |>
    mutate(model = model_name) |>
    select(model, ssid, survey_abbrev, fishing_event_id,
      year, X, Y, depth_m,
      catch_count, hook_count, offset,
      is_delta,
      catch_est,
      any_of(c("est", "est_non_rf", "est_rf", "omega_s", "epsilon_st",
      "est1", "est2", "omega_s1", "omega_s2"))
    )
  })
}
test_OS <- get_preds(OS) |> left_join(historical)
test_ON <- get_preds(ON) |> left_join(historical)
test_IN <- get_preds(IN) |> left_join(historical)

library(ggsidekick)
test_ON |>
  filter(model == "nb2") |>
  ggplot() +
  geom_point(aes(x = catch_est, y = catch_count + 1, colour = X), size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_colour_viridis_c() +
  facet_wrap(~ model, scales = "free_y") +
  # facet_wrap(~ pfma, scales = "free_y") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  ggtitle(unique(test$survey_abbrev))

test |> glimpse()
test |> filter(is_delta == TRUE) |> glimpse()

# Compare conditioning model predictions vs observations
get_preds <- function(
  surveys_to_check = list(
    list(fit = fit_IN, data = sp_dat |> filter(stringr::str_detect(survey_abbrev, "INS N")), name = "INS N"),
    list(fit = fit_ON, data = sp_dat |> filter(stringr::str_detect(survey_abbrev, "OUT N")), name = "OUT N"),
    list(fit = fit_OS, data = sp_dat |> filter(stringr::str_detect(survey_abbrev, "OUT S")), name = "OUT S")
  )
) {

  map_dfr(surveys_to_check, function(survey) {
    pred <- predict(survey$fit, model = NA,)
    data.frame(
      survey = survey$name,
      observed = survey$data$catch_count,
      fitted_mu = exp(pred$est),  # Convert from link to response scale
      fitted_eta = pred$est,
      residual = survey$data$catch_count - pred$observed#exp(pred$est)
    )
  })
}
p <- get_preds()
fit_OS_nb2 <- readRDS(here::here(fit_dir, "yelloweye-rockfish-HBLL-OUT-S.rds"))
fit_OS_nb2_mix <- readRDS(here::here(fit_dir, "yelloweye-rockfish-HBLL-OUT-S-nb2mixed.rds"))
p <- get_preds(
  surveys_to_check = list(
    list(fit = fit_OS_nb2, data = sp_dat |> filter(survey_abbrev == "HBLL OUT S"), name = "OUT S nb2"),
    list(fit = fit_OS_nb2_mix, data = sp_dat |> filter(survey_abbrev == "HBLL OUT S"), name = "OUT S nb2mixed")
  )
)

d <- fit_OS$data

p <- predict(fit_OS, nsim = 100)
p_exp <- exp(p)
test <- fit_OS$data |>
  select(survey_abbrev, year, fishing_event_id, catch_count, hook_count)
test$mean_mu <- rowMeans(p_exp)
test$sd_mu <- apply(p_exp, 1, sd)

test |>
  ggplot() +
  aes(x = mean_mu, y = catch_count) +
  geom_point() +
  facet_wrap(~survey_abbrev, scales = "free")

test <- p |>
  select(survey_abbrev, year, fishing_event_id, X, Y, catch_count, hook_count,
         est1:omega_s2) |>
  mutate(est_mu = plogis(est1) * exp(est2),
         residual = catch_count - est_mu,
         present = ifelse(est_mu > 0, 1, 0))
test |>
  summarise(
    mean_observed = mean(catch_count),
    mean_fitted = mean(est_mu),
    ratio = mean_fitted / mean_observed,
    var_observed = var(catch_count),
    var_fitted = var(est_mu),
    zero_count_obs = sum(catch_count == 0),
    zero_count_fitted = sum(est_mu == 0)
  )

# Look at mean levels and overdispersion by survey
p |>
  group_by(survey) |>
  summarise(
    mean_observed = mean(observed),
    mean_fitted = mean(fitted_mu),
    ratio = mean_fitted / mean_observed,
    var_observed = var(observed),
    var_fitted = var(fitted_mu),
    zero_count_obs = sum(observed == 0),
    zero_count_fitted = sum(fitted_mu < 0.5)  # Approximation
  ) |> as.data.frame()

# Number of zeroes observed over time doesn't really seem to change?
historical |>
  group_by(survey_abbrev, year) |>
  summarise(zero_prop = mean(catch_count == 0)) |>
  ggplot(aes(year, zero_prop, color = survey_abbrev)) + geom_line()

sp_dat |>
  mutate(depth_bin = cut(depth_m, breaks = 10)) |>
  group_by(depth_bin) |>
  summarise(zero_prop = mean(catch_count == 0),
            mean_catch = mean(catch_count))

sp_dat |>
  ggplot() +
  aes(longitude, latitude, color = depth_m, shape = presence) +
  geom_point() +
  scale_shape_manual(values = c(21, 19)) +
  viridis::scale_colour_viridis() +
  facet_wrap(~survey_abbrev, scales = "free") +
  labs(title = "Spatial distribution of zeros")

sp_dat |>
  mutate(depth_bin = cut(depth_m, breaks = 10)) |>
  group_by(depth_bin, survey_abbrev) |>
  summarise(prop_zero = sum(catch_count == 0) / n()) |>
  ggplot() +
  aes(depth_bin, prop_zero) +
  geom_point() +
  facet_wrap(~ survey_abbrev, scales = "free")

sp_dat |>
  ggplot() +
  aes(depth_m, fill = presence) +
  geom_histogram() +
  facet_wrap(~ survey_abbrev, scales = "free")

sp_dat |>
  ggplot() +
  aes(x = depth_m, y = as.numeric(catch_count == 0)) +
  geom_smooth() +
  facet_wrap(~survey_abbrev, scales = "free") +
  labs(title = "Zero catch counts by depth")


# More likely to get zero catch counts when the hook count is at the
# lower end of the range
# HBLL INS seems to have more consistent hook counts over time so this doesn't
# seem to be as much of a problem
historical |>
  mutate(hook_bins = cut_number(hook_count, n = 9)) |>
  group_by(survey_abbrev, hook_bins) |>
  summarise(zero_prop = mean(catch_count == 0), n = n()) |>
  mutate(bin_id = as.numeric(hook_bins)) |>
  as.data.frame() |>
  ggplot() +
  aes(x = bin_id, y = zero_prop, color = survey_abbrev) +
  geom_point() +
  stat_smooth(method = "lm", aes(group = survey_abbrev), se = FALSE)

historical |>
  ggplot(aes(x = hook_count)) +
  geom_histogram(binwidth = 10) +  # Added binwidth to avoid the warning
  geom_vline(
    data = . %>% group_by(survey_abbrev) %>% summarise(mean_hooks = mean(hook_count)),
    aes(xintercept = mean_hooks),
    color = "red"
  ) +
  facet_wrap(~ survey_abbrev, scales = "free")

library(DHARMa)
# fit <- fit_OS
model_name <- "HBLL-OUT-N-nb2"
model_name <- "HBLL-OUT-S-nb2"
model_name <- "HBLL-INS-N-nb2"
model_name <- "HBLL-OUT-N-binomial"
model_name <- "HBLL-OUT-S-binomial"
model_name <- "HBLL-INS-N-binomial"
# model_name <- "HBLL-OUT-S-delta-nb2"
# model_name <- "HBLL-OUT-N-delta-nb2"
# model_name <- "HBLL-OUT-S-depth-tweedie"
# model_name <- "HBLL-OUT-S-poisson"
fit <- readRDS(here::here(fit_dir, paste0("yelloweye-rockfish-", model_name, ".rds")))
s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn")
r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
# dev.set(2)
plot(r_nb2, title = model_name)

# Function to create dharma residuals
create_dharma_residuals <- function(model_name, fit_dir, return_sims = FALSE) {
  fit <- readRDS(here::here(fit_dir, paste0("yelloweye-rockfish-", model_name, ".rds")))
  s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn", mle_mvn_samples = "multiple")
  if (return_sims) {
    return(s_nb2)
  } else {
    dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
  }
}

# Create residuals for all models
model_names <- c("HBLL-OUT-N-binomial", "HBLL-OUT-S-binomial", "HBLL-INS-N-binomial")
model_names <- c("HBLL-OUT-N-nb2", "HBLL-OUT-S-nb2", "HBLL-INS-N-nb2")
model_names <- c("HBLL-OUT-N-depth-nb2", "HBLL-OUT-S-depth-nb2", "HBLL-INS-N-depth-nb2")
sim_list <- map(model_names, ~create_dharma_residuals(.x, fit_dir, return_sims = TRUE)) %>%
  purrr::set_names(model_names)
residuals_list <- map(model_names, ~create_dharma_residuals(.x, fit_dir)) %>%
  purrr::set_names(model_names)

# Residuals:
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

for(i in seq_along(model_names)) {
  res_obj <- residuals_list[[i]]

  # Plot 1: QQ plot
  plotQQunif(res_obj, main = paste("QQ Plot -", model_names[i]))

  # Plot 2: Residuals vs predicted
  plotResiduals(res_obj) #main = paste("Residuals -", model_names[i]))
}
par(mfrow = c(1, 1))

# Look at cumulative catch by rank order using simulate.sdmTMB():
survey_names <- c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")  # e.g., c("HBLL INS N", "HBLL OUT N", "HBLL OUT S")
all_surveys_comb <- map_dfr(1:length(sim_list), function(i) {
  .sa <- survey_names[i]
  sim_matrix <- sim_list[[i]]

  # Process simulation data
  sim_data <- as.data.frame(sim_matrix) |>
    mutate(row_id = row_number()) |>
    pivot_longer(cols = -row_id, names_to = "replicate", values_to = "observed") |>
    mutate(replicate = as.numeric(gsub("V", "", replicate))) |>
    mutate(d = "simulated", survey_abbrev = .sa) |>
    arrange(replicate)

  # Process historical data for this survey
  hist_data <- historical |>
    filter(survey_abbrev == .sa) |>
    mutate(d = "historical",
           row_id = row_number(),
           replicate = 0)

  sim_data <- left_join(sim_data, hist_data |> select(-observed, -replicate, -d))

  # Combine
  bind_rows(sim_data, hist_data)
})
glimpse(all_surveys_comb)
# p_depth_nb2 <-
# p_nb2 <-
  ggplot(data = all_surveys_comb |> filter(replicate %in% c(1:5))) +
  facet_wrap(survey_abbrev ~ replicate, scales = "free", ncol = 5) +
  scale_colour_manual(values = c("orange")) +
  geom_jitter(mapping = aes(x = year, y = observed, color = d), shape = 21) +
  geom_jitter(data = all_surveys_comb |> filter(replicate %in% c(0)) |> select(-replicate),
    mapping = aes(x = year, y = observed), colour = "steelblue", shape = 21, alpha = 0.5) +
  # ggtitle("Depth-NB2")
  ggtitle("NB2")

# Pre-calculate the cumulative catch for all data
plot_data <- all_surveys_comb |>
  group_by(d, survey_abbrev, replicate) |>
  arrange(observed) |>
  mutate(cumul_catch = cumsum(observed),
          rank_obs = rank(observed))

plot_data |>
    ggplot(aes(x = rank_obs, y = cumul_catch)) +
    geom_line(data = filter(plot_data, d == "simulated"),
              aes(group = interaction(d, replicate)),
              colour = "black", alpha = 0.3) +
    geom_line(data = filter(plot_data, d == "historical"),
              aes(group = interaction(d, replicate)),
              colour = "red", linewidth = 1.2) +
    facet_wrap(~ survey_abbrev, scales = "free") +
    labs(x = "Rank (ordered observations)", y = "Cumulative catch", title = "simulate.sdmTMB")

mod1 <- simulate(
  fit_OS,
  nsim = 100L,
  seed = sample.int(1e+06, 1L),
  type = "mle-mvn", #c("mle-eb", "mle-mvn"),
  model = 1, #c(NA, 1, 2),
  newdata = NULL,
  re_form = NULL,
  mle_mvn_samples = "multiple", #c("single", "multiple"),
  mcmc_samples = NULL,
  return_tmb_report = FALSE,
  observation_error = TRUE,
  size = NULL,
  silent = FALSE
)
test <-mod1 |>
  as.data.frame() |>
    mutate(row_id = row_number()) |>
    pivot_longer(cols = -row_id, names_to = "replicate", values_to = "observed") |>
    mutate(replicate = as.numeric(gsub("V", "", replicate))) |>
    mutate(d = "simulated", survey_abbrev = "HBLL OUT S") |>
    arrange(replicate) |>
    left_join(fit_OS$data |> mutate(row_id = row_number())) |>
    bind_rows(fit_OS$data |> mutate(d = "historical", row_id = row_number(), replicate = 0))

test |>
  filter(replicate %in% 0:9) |>
  ggplot() +
  aes(x = year, y = observed, color = d) +
  geom_point() +
  facet_wrap(~ replicate, scales = "free")


testZeroInflation(r_nb2)

fit <- readRDS(here::here(fit_dir, "yelloweye-rockfish-HBLL-OUT-S-depth.rds"))
s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn")
r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
dev.set(3)
plot(r_nb2, title = "depth")




fit <- readRDS(here::here(fit_dir, "yelloweye-rockfish-HBLL-OUT-S-nb2mixed.rds"))
s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn")
r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
dev.set(3)
plot(r_nb2, title = "depth")

# DHARMa::testResiduals(r_nb2)
# sp_r <- DHARMa::recalculateResiduals(s_nb2, group = fit$data$fyear)
# DHARMa::testSpatialAutocorrelation(r_nb2,
#   x = fit$data$X,
#   y = fit$data$Y)
# DHARMa::testZeroInflation(r_nb2)
