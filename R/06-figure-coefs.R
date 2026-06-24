source(here::here("R", "00-setup.R"))

fit_dir <- here::here("data-generated", "fits")
target_terms <- c("restricted", "phi", "sigma_E", "sigma_O", "range")

ar1_params <- readRDS(ar1_parameters_file)

split_species_label <- function(x) {
  dplyr::if_else(
    x == "Pacific Spiny Dogfish",
    "Pacific Spiny\nDogfish",
    stringr::str_replace(x, "^(\\S+)\\s+", "\\1\n")
  )
}

log_phi_only <- function(x, term) {
  phi_rows <- term == "phi" & !is.na(x)
  x[phi_rows] <- log(x[phi_rows])
  x
}

parse_fit_filename <- function(path) {
  file_name <- basename(path)
  fit_parts <- stringr::str_match(
    file_name,
    "^(.*?)-(HBLL-(?:OUT-N|OUT-S|INS-N))-.*\\.rds$"
  )

  if (is.na(fit_parts[1, 1])) {
    stop("Could not parse species/survey from filename: ", file_name)
  }

  tibble(
    fit_path = path,
    fit_file = file_name,
    species = sp_from_hyphens(fit_parts[1, 2]),
    survey = gsub("-", " ", fit_parts[1, 3], fixed = TRUE)
  )
}

extract_fit_terms <- function(path, terms = target_terms) {
  fit <- readRDS(path)
  fit_meta <- parse_fit_filename(path)

  fit_terms <- bind_rows(
    sdmTMB::tidy(fit, conf.int = TRUE),
    sdmTMB::tidy(fit, effects = "ran_pars", conf.int = TRUE)
  ) |>
    filter(term %in% terms) |>
    select(term, estimate, std.error, conf.low, conf.high)

  tibble(term = terms) |>
    left_join(fit_terms, by = "term") |>
    mutate(
      component = if_else(term == "restricted", "fixed", "ran_pars"),
      .before = estimate
    ) |>
    (\(x) bind_cols(fit_meta, x))()
}

fit_coef_df <- list.files(
  fit_dir,
  pattern = "\\.rds$",
  full.names = TRUE
) |>
  sort() |>
  purrr::map_dfr(extract_fit_terms)

ar1_coef_df <- ar1_params |>
  transmute(
    species,
    survey = survey_abbrev,
    rho_V,
    sigma_V
  ) |>
  tidyr::pivot_longer(
    cols = c(rho_V, sigma_V),
    names_to = "term",
    values_to = "estimate"
  ) |>
  mutate(
    fit_path = NA_character_,
    fit_file = NA_character_,
    component = "ar1_params",
    std.error = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    .before = estimate
  ) |>
  select(fit_path, fit_file, species, survey, term, component, estimate, std.error, conf.low, conf.high)

fit_coef_df <- bind_rows(fit_coef_df, ar1_coef_df)

term_lu <- tibble(
  term = c("restricted", "phi", "sigma_E", "sigma_O", "range", "sigma_V", "rho_V"),
  term_label = c(
    "Inside MPAN",
    "log(Precision)",
    "Spatiotemporal\nSD",
    "Spatial SD",
    "Range",
    "Interannual SD",
    "Interannual\ncorrelation"
  )
)

fit_coef_summary <- fit_coef_df |>
  filter(species != "pacific halibut") |>
  left_join(term_lu, by = "term") |>
  mutate(
    survey = factor(
      survey,
      levels = c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")
    ),
    estimate = log_phi_only(estimate, term),
    conf.low = log_phi_only(conf.low, term),
    conf.high = log_phi_only(conf.high, term),
    term_label = factor(term_label, levels = term_lu$term_label)
  ) |>
  select(species, survey, term, term_label, estimate, std.error, conf.low, conf.high) |>
  mutate(survey = forcats::fct_rev(survey))

fit_coef_summary |>
  mutate(
    species = stringr::str_to_title(species),
    species = dplyr::recode(species, "North Pacific Spiny Dogfish" = "Pacific Spiny Dogfish"),
    species = split_species_label(species)
  ) |>
  ggplot(aes(survey, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline(
    data = tibble(
      term_label = c("Inside MPAN", "Interannual\ncorrelation"),
      yintercept = 0
    ),
    aes(yintercept = yintercept),
    linetype = "dashed", colour = "grey50",
    inherit.aes = FALSE
  ) +
  geom_linerange(
    data = \(x) dplyr::filter(x, !is.na(estimate) & !is.na(conf.low) & !is.na(conf.high)),
    linewidth = 0.4,
    colour = "grey30"
  ) +
  geom_point(
    data = \(x) dplyr::filter(x, !is.na(estimate) & !is.na(conf.low) & !is.na(conf.high)),
    pch = 21,
    colour = "grey30"
  ) +
  geom_point(
    data = \(x) dplyr::filter(x, !is.na(estimate) & (is.na(conf.low) | is.na(conf.high))),
    pch = 21,
    colour = "grey30"
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = scales::label_number(trim = TRUE)
  ) +
  coord_flip() +
  facet_grid(species~term_label, scales = "free_x") +
  xlab("Survey") + ylab("Estimate")
ggsave("figures/conditioning-coefs.png", width = 8, height = 5.3)
