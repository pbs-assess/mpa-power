library(dplyr)
library(sdmTMB)
library(ggplot2)

source(here::here("R", "00-utils.R"))

SURVEY_ABBREV <- "HBLL OUT N"

for (SPECIES in c("lingcod", "yelloweye rockfish", "quillback rockfish")) {

  if (SPECIES == "lingcod") fit <- readRDS(file.path(fit_dir, "lingcod-HBLL-OUT-N-betabinomial-restricted-depth-on-iid-13b5638708a67147.rds"))
  if (SPECIES == "yelloweye rockfish") fit <- readRDS(file.path(fit_dir, "yelloweye-rockfish-HBLL-OUT-N-betabinomial-restricted-depth-on-iid-5d7c2f4e2b7a5a0c.rds"))
  if (SPECIES == "quillback rockfish") fit <- readRDS(file.path(fit_dir, "quillback-rockfish-HBLL-OUT-N-betabinomial-restricted-depth-on-iid-4ee1e6841eb55b7f.rds"))
  print(fit) # don't crash!

  one_sample_posterior <- function(object, use_names = TRUE) {
    tmp <- object$tmb_obj$env$MC(n = 1L, keep = TRUE, antithetic = FALSE)
    re_samp <- as.vector(attr(tmp, "samples"))
    lp <- object$tmb_obj$env$last.par.best
    p <- numeric(length(lp))
    fe <- object$tmb_obj$env$lfixed()
    re <- object$tmb_obj$env$lrandom()
    p[re] <- re_samp
    p[fe] <- lp[fe]
    if (use_names) names(p) <- names(lp)
    p
  }

  set.seed(1)
  osp <- one_sample_posterior(fit)
  omega_draw <- osp[grepl("omega_s", names(osp))]

  b <- tidy(fit, "ran_pars")
  bf <- tidy(fit)

  intercept <- bf$estimate[bf$term == "fyear2023"]
  restricted <- bf$estimate[bf$term == "restricted"]
  sigma_E <- b$estimate[b$term == "sigma_E"]
  phi <- b$estimate[b$term == "phi"]
  range_val <- b$estimate[b$term == "range"]
  b_log_depth <- bf$estimate[bf$term == "log_depth"]
  b_log_depth2 <- bf$estimate[bf$term == "I(log_depth^2)"]

  dat <- fit$data
  lastdat <- filter(dat, year == 2023)
  mesh <- make_mesh(lastdat, c("X", "Y"), mesh = fit$spde$mesh)

  out_list <- purrr::map(1:5, \(i) {
    sim_dat <- sdmTMB::simulate_new(
      formula = ~ 1 + restricted + log_depth + I(log_depth^2),
      data = lastdat,
      mesh = mesh,
      family = fit$family,
      time = "year",
      sigma_E = sigma_E,
      phi = phi,
      range = range_val,
      fixed_re = list(omega_s = matrix(omega_draw), epsilon_st = NULL, zeta_s = NULL),
      B = c(intercept, restricted, b_log_depth, b_log_depth2),
      weights = lastdat$hook_count,
      seed = i * 2
    )
  })

  out <- bind_rows(out_list, .id = "iteration")
  out$species <- SPECIES
  out$survey_abbrev <- SURVEY_ABBREV

  out2 <- bind_rows(out,
    mutate(lastdat, iteration = "0") |> select(observed = catch_count, X, Y, iteration)
  )

  n_sims <- length(unique(out$iteration))
  sim_labels <- paste0("Simulation\nreplicate ", seq_len(n_sims))
  iteration_levels <- c("0", as.character(seq_len(n_sims)))
  iteration_labels <- c("Observed", sim_labels)
  out2$iteration <- factor(out2$iteration, levels = iteration_levels, labels = iteration_labels)

  out_sf <- XY_to_sf(out2)

  ggplot() +
    geom_sf(data = pacea::bc_coast, fill = "grey90", colour = "grey90") +
    geom_sf(data = out_sf,aes(size = observed, colour = observed), alpha = 0.8) +
    gfplot::coord_sf_auto(out_sf) +
    scale_colour_viridis_c(name = "Count") +
    scale_size_continuous(name = "Count") +
    gfplot::theme_pbs() +
    facet_wrap(~iteration) +
    theme(axis.title = element_blank()) +
    ggtitle(stringr::str_to_title(unique(out$species))) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3))
  ggsave(paste0("figures/turing-map-", unique(out$species), "-", gsub(" ", "-", unique(out$survey_abbrev)), ".png"), width = 10, height = 8)

  rm(fit)
  gc()
}
