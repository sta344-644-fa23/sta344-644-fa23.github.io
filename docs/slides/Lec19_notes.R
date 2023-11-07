library(sf)
library(tidyverse)

nc = st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE) %>% 
  select(-(AREA:CNTY_ID), -(FIPS:CRESS_ID)) %>%
  st_transform(st_crs("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

A = st_touches(nc, sparse=FALSE)
listW = spdep::mat2listw(A)

sar_model = "
data {
  int<lower=0> N;
  vector[N] y;
  matrix[N,N] A;
}
transformed data {
  vector[N] nb = A * rep_vector(1, N);
  matrix[N,N] D = diag_matrix(nb);
  matrix[N,N] D_inv = diag_matrix(1.0/nb);
  matrix[N,N] I = diag_matrix(rep_vector(1, N));
}
parameters {
  vector[N] w_s;
  real beta;
  real<lower=0> sigma2;
  real<lower=0> sigma2_w;
  real<lower=0,upper=1> phi;
}
model {
  matrix[N,N] C = I - phi * D_inv * A;
  matrix[N,N] Sigma_inv = C' * D * C / sigma2;  
  
  w_s ~ multi_normal_prec(rep_vector(0,N), Sigma_inv);

  beta ~ normal(0,10);
  sigma2 ~ cauchy(0,5);
  sigma2_w ~ cauchy(0,5);

  y ~ normal(beta + w_s, sigma2_w);
}
"

write_file(sar_model, "sar_model.stan")

car_stan = cmdstanr::cmdstan_model("sar_model.stan")

data = list(
  N = nrow(nc),
  y = 1000 * nc$SID74 / nc$BIR74,
  A = A
)
car_fit = car_stan$sample(
  data = data, chains=4, parallel_chains=4, seed = 12345, 
  iter_sampling = 2000, iter_warmup = 1000
)

car_post = tidybayes::gather_draws(
  car_fit, beta, sigma2, sigma2_w, phi
)
  
car_post |>
  ggplot(aes(x=.iteration, y=.value, color=as.factor(.chain))) +
    geom_line(alpha=0.5) +
    facet_wrap(~.variable, scale="free_y")

car_post |>
  ggplot(aes(x=.value, fill=as.factor(.chain))) +
  geom_density(alpha=0.25) +
  facet_wrap(~.variable, scale="free")

car_post_pred = tidybayes::spread_draws(
  car_fit, beta, w_s[i]
) |>
  mutate(y_hat = beta + w_s) |>
  group_by(.chain, i) |>
  summarize(
    y_hat_mean = mean(y_hat),
    y_hat_q025 = quantile(y_hat, 0.025),
    y_hat_q975 = quantile(y_hat, 0.975),
    .groups = "drop"
  )

nc_post = left_join(
  nc |> mutate(i = seq_len(n())),
  car_post_pred,
  by = "i"
) |>
  filter(.chain == 1) |>
  mutate(
    resid = 1000*SID74/BIR74 - y_hat_mean
  )

(
  ggplot(nc_post, aes(fill=y_hat_mean)) +
    geom_sf()
) / (
  ggplot(nc_post, aes(fill=resid)) +
    geom_sf()
)

spdep::moran.test(nc_post$resid, listW)


ggplot(nc_post, aes(x=1000*SID74/BIR74, y=y_hat_mean)) +
  geom_abline(intercept=0, slope=1, color="grey") +
  geom_point() +
  geom_errorbar(aes(ymin=y_hat_q025, ymax=y_hat_q975), alpha=0.25) +
  coord_fixed()
