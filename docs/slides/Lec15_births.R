library(tidyverse)

births = read_csv("data/births_usa_1969.csv") |>
  mutate(
    has_leap_year = sum(month == 2 & day == 29),
    .by = year
  ) |> 
  filter(!has_leap_year) |>
  summarize(
    births = mean(births),
    .by = c(day_of_year)
  ) |>
  mutate(
    scaled_log_births = c(scale(log(births)))
  )

births = read_csv("data/births_usa_1969.csv") |>
  mutate(
    date = paste(year,month,day,sep="-") |> ymd(),
    scaled_log_births = c(scale(log(births)))
  )

library(tsibble)
library(feasts)

tsbirths = births |>
  as_tsibble(index=date)

feasts::gg_tsdisplay(tsbirths, births, plot_type = "season")


model = rstan::stan_model(file = "Lec15_births_gp.stan")

fit = rstan::sampling(
  model,
  data = list(
    N = nrow(births),
    y = births$scaled_log_births,
    x = births$id / max(births$id),
    M_f1 = 20,
    M_f2 = 20,
    period = 365.25 / max(births$id)
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4
)

pars = fit |>
  tidybayes::gather_draws(`sigma_.*`,`lengthscale_.*`,`sigma`, regex = TRUE)

ggplot(pars, aes(x=.iteration, y=.value, color=as_factor(.chain))) +
  geom_line() +
  facet_wrap(~.variable, scale="free_y") +
  guides(color="none")

ggplot(pars, aes(x=.value, y=as_factor(.variable))) +
  tidybayes::stat_halfeye(.width = c(.90, .5)) +
  labs(y="Parameter", x="value")

ggplot(pars, aes(x=.value, fill=.variable)) +
  geom_density() +
  facet_wrap(~.variable, scale="free") +
  guides(fill="none")


f = tidybayes::spread_draws(fit, f1[id], f2[id]) |>
  ungroup() |>
  full_join(births, by="id")

f |>
  summarize(
    post_mean = mean(f1), .by = c(date, .chain)
  ) |>
  ggplot() +
    geom_line(aes(x=date, y=post_mean, color=as_factor(.chain)), linewidth=1.5) +
    geom_point(data=births, aes(y=scaled_log_births, x=date), color="black", alpha=0.1)

f |>
  summarize(
    post_mean = mean(f2), 
    avg_births = mean(scaled_log_births),
    .by = c(day_of_year, .chain)
  ) |>
  ggplot() +
    geom_line(aes(x=day_of_year, y=post_mean, color=as_factor(.chain)), linewidth=1.5) +
    geom_point(aes(y=avg_births, x=day_of_year), color="black", alpha=0.1)


f |>
  summarize(
    post_mean = mean(f1+f2), .by = c(date, .chain)
  ) |>
  ggplot() +
  geom_line(aes(x=date, y=post_mean, color=as_factor(.chain)), linewidth=1.5) +
  geom_point(data=births, aes(y=scaled_log_births, x=date), color="black", alpha=0.1)
