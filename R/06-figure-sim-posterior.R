library(dplyr)
library(sdmTMB)
library(ggplot2)

source(here::here("R", "00-setup.R"))
source(here::here("R", "00-utils.R"))

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

# fit_file <- list.files(fit_dir, pattern = "^yelloweye-rockfish.*OUT-N.*rds$", full.names = TRUE)
# fit_file <- list.files(fit_dir, pattern = "^lingcod.*OUT-N.*rds$", full.names = TRUE)
# fit_file <- list.files(fit_dir, pattern = "^quillback-rockfish.*OUT-N.*rds$", full.names = TRUE)
# fit_file <- list.files(fit_dir, pattern = "^north-pacific-spiny-dogfish.*OUT-N.*rds$", full.names = TRUE)
fit_file <- list.files(fit_dir, pattern = "^north-pacific-spiny-dogfish.*OUT-S.*rds$", full.names = TRUE)
# Note: dogfish seems to have less variation than expected in the simulated outputs
# I can't remember what it looked like without depth.
fit <- readRDS(fit_file)
fit # I think it was segfaulting if I didn't print the model before working with it

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

dat <- fit$data |>
  mutate(last_sampled_year = max(year), dataset = "Historical") |>
  select(species = species_common_name, survey_abbrev, X, Y, restricted, log_depth,
         year, catch_count, hook_count, last_sampled_year, dataset)

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
    formula = ~ 1 + restricted + log_depth + I(log_depth^2),
    data = future_data,
    mesh = mesh,
    family = fit$family,
    time = "year",
    sigma_E = sigma_E,
    phi = phi,
    range = range_val,
    fixed_re = list(omega_s = matrix(omega_draw), epsilon_st = NULL, zeta_s = NULL),
    B = c(intercept, restricted, b_log_depth, b_log_depth2),
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
