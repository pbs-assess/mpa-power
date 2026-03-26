library(dplyr)
library(sdmTMB)
library(ggplot2)

source(here::here("R", "00-utils.R"))

fit <- readRDS("data-generated/fits/north-pacific-spiny-dogfish-HBLL-OUT-N-betabinomial-on-iid-e40c7b759e26ff69.rds")
fit <- readRDS("data-generated/fits/lingcod-HBLL-OUT-N-betabinomial-on-iid-2a49c4ed06e10dc5.rds")
fit

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

dat <- fit$data
lastdat <- filter(dat, year == 2023)
mesh <- make_mesh(lastdat, c("X", "Y"), mesh = fit$spde$mesh)

out_list <- purrr::map(1:5, \(i) {
  sim_dat <- sdmTMB::simulate_new(
    formula = ~ 1 + restricted,
    data = lastdat,
    mesh = mesh,
    family = fit$family,
    time = "year",
    sigma_E = sigma_E,
    phi = phi,
    range = range_val,
    fixed_re = list(omega_s = matrix(omega_draw), epsilon_st = NULL, zeta_s = NULL),
    B = c(intercept, restricted),
    weights = lastdat$hook_count,
    seed = i * 2
  )
})

out <- bind_rows(out_list, .id = "iteration")

out2 <- bind_rows(out,
  mutate(lastdat, iteration = "0") |> select(observed = catch_count, X, Y, iteration)
  )

out_sf <- XY_to_sf(out2)

ggplot() +
  geom_sf(data = pacea::bc_coast, fill = "grey90", colour = "grey90") +
  geom_sf(data = out_sf,aes(size = observed, colour = observed), alpha = 0.8) +
  gfplot::coord_sf_auto(out_sf) +
  scale_colour_viridis_c() +
  gfplot::theme_pbs() +
  facet_wrap(~iteration) +
  theme(axis.title = element_blank())

