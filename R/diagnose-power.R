power_df <- readRDS("data-generated/power-results-df.rds")
theta <- readRDS("data-generated/fit_characteristics.rds")

unique(power_df$eval_year)
unique(power_df$sim_mpa_trend)

x <- filter(power_df, eval_year == 2034, sampling_plan == "status quo", sim_mpa_trend == 1.009)

x <- left_join(x, theta)

library(ggplot2)

ggplot(x, aes(encountered_rate_avg, power_signed, colour = species)) + geom_point()
ggplot(x, aes(encountered_count_per_year, power_signed, colour = species)) + geom_point()
ggplot(x, aes(sigma_O, power_signed, colour = species)) + geom_point()
ggplot(x, aes(sigma_E, power_signed, colour = species)) + geom_point()
ggplot(x, aes(convergence_rate, power_signed, colour = species)) + geom_point()

ggplot(power_df, aes(convergence_rate, species)) + geom_point(position = position_jitter(height = 0.1))
