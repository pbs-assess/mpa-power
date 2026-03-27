library(dplyr)
library(sdmTMB)
library(ggplot2)

source(here::here("R", "00-utils.R"))

fig_dir <- here::here("figures")

fit <- readRDS("data-generated/fits/north-pacific-spiny-dogfish-HBLL-OUT-N-betabinomial-on-iid-e40c7b759e26ff69.rds")
fit <- readRDS("data-generated/fits/lingcod-HBLL-OUT-N-betabinomial-on-iid-2a49c4ed06e10dc5.rds")
fit <- readRDS("data-generated/fits/yelloweye-rockfish-HBLL-OUT-N-betabinomial-on-iid-144f4b8c390be8df.rds")
fit <- readRDS("data-generated/fits/quillback-rockfish-HBLL-OUT-N-betabinomial-on-iid-211d46c156192c75.rds")
print(fit)

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

dat <- fit$data |>
  mutate(last_sampled_year = max(year), dataset = "Historical") |>
  select(species = species_common_name, survey_abbrev, X, Y, restricted, year,
         catch_count, hook_count, last_sampled_year, dataset)

future_data <- purrr::map_dfr(seq(2, 25, by = 2), function(new_year) {
  sampled_year <- sample(unique(dat$year), size = 1)

  dat |>
    filter(year == sampled_year) |>
    mutate(year = last_sampled_year + new_year)
})

glimpse(future_data)

mesh <- make_mesh(future_data, c("X", "Y"), mesh = fit$spde$mesh)

out <- purrr::map_dfr(1:6, \(i) {
  sim_dat <- sdmTMB::simulate_new(
    formula = ~ 1 + restricted,
    data = future_data,
    mesh = mesh,
    family = fit$family,
    time = "year",
    sigma_E = sigma_E,
    phi = phi,
    range = range_val,
    fixed_re = list(omega_s = matrix(omega_draw), epsilon_st = NULL, zeta_s = NULL),
    B = c(intercept, restricted),
    weights = future_data$hook_count,
    seed = i * 2
  )
  bind_rows(dat, sim_dat |> mutate(catch_count = observed, dataset = "Simulated")) |>
    mutate(iteration = i)
})

out

ggplot(out) +
  geom_point(aes(x = year, y = catch_count, colour = dataset), shape = 21) +
  gfplot::theme_pbs() +
  scale_colour_manual(values = c("Historical" = "dodgerblue", "Simulated" = "orange")) +
  facet_wrap(~iteration, ncol = 3, labeller = labeller(iteration = ~paste("Rep", .))) +
  theme(axis.title = element_blank()) +
  labs(x = "Year", y = "Catch count", colour = "") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9)) +
  ggtitle(paste0(stringr::str_to_title(unique(dat$species)), " | ", unique(dat$survey_abbrev)))

ggsave(file.path(fig_dir, paste0("sim-posterior-", unique(dat$species), "-", unique(dat$survey_abbrev), ".png")),
  width = 7, height = 5.2)
