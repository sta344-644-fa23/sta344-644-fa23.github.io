library(tidyverse)
library(patchwork)

ggplot2::theme_set(theme_minimal())






model = "model{
  y ~ dmnorm(rep(0,N), inverse(Sigma))

  for (i in 1:(length(y)-1)) {
    for (j in (i+1):length(y)) {
      k1[i,j] <- sigma2[1] * exp(- pow(l[1] * d[i,j],2))
      k2[i,j] <- sigma2[2] * exp(- pow(l[2] * d[i,j],2) - 2 * pow(l[3] * sin(pi*d[i,j] / per), 2))
      
      Sigma[i,j] <- k1[i,j] + k2[i,j]
      Sigma[j,i] <- Sigma[i,j]
    }
  }

  for (i in 1:length(y)) {
    Sigma[i,i] <- sigma2[1] + sigma2[2] + sigma2[3]
  }  

  for(i in 1:3){
    #sigma2[i] ~ dnorm(0, 1) T(0,)
    #l[i] ~ dnorm(0, 1) T(0,)
    sigma2[i] ~ dt(0, 2.5, 1) T(0,)
    l[i] ~ dt(0, 2.5, 1) T(0,)
  }
}"

dir.create("Lec15_results/", showWarnings = FALSE, recursive = TRUE)

flexiblas::flexiblas_load_backend("OPENBLAS-THREADS") |>
  flexiblas::flexiblas_switch()
#flexiblas::flexiblas_set_num_threads(24)

if (!file.exists("Lec15_results/init_model.rds")) {
  m = rjags::jags.model(
    textConnection(model), 
    data = list(
      y = births$scaled_log_births,
      d = dist(births$day_of_year / max(births$day_of_year)) %>% as.matrix(),
      per = 7 / max(births$day_of_year),
      pi = pi,
      N = nrow(births)
    ),
    n.adapt=5000,
    n.chains = 1
  )
  saveRDS(m, file="Lec15_results/init_model.rds")
} else {
  m = readRDS("Lec15_results/init_model.rds")
}

if (!file.exists("Lec15_results/fit_model.rds")) {
  gp_coda = rjags::coda.samples(
    m, variable.names=c("sigma2", "l"),
    n.iter=5000,
    thin=5
  )
  saveRDS(gp_coda, file="Lec15_results/fit_model.rds")
} else {
  gp_coda = readRDS("Lec15_results/fit_model.rds")
}


r = tidybayes::gather_draws(gp_coda, sigma2[i], l[i]) |>
  ungroup() |>
  mutate(
    var = glue::glue("{.variable}[{i}]")
  )

ggplot(r, aes(x=.iteration, y=.value, color=as_factor(.chain))) +
  geom_line() +
  facet_wrap(~var, scales = "free_y") +
  guides(fill = "none")

ggplot(r, aes(x=.value, fill=as_factor(.chain))) +
  geom_density() +
  facet_wrap(~var, scales = "free") +
  guides(fill = "none")

cov_C1 = function(d, sigma2, l, per) {
 sigma2[1] * exp(- (l[1] * d)^2)
}

cov_C2 = function(d, sigma2, l, per) {
  sigma2[2] * exp(- (l[2] * d)^2 - 2 * (l[3] * sin(pi*d / per))^2)
}

cov_full = function(d, sigma2, l, per) {
  cov_C1(d, sigma2, l, per) +
  cov_C2(d, sigma2, l, per) +
  ifelse(abs(d)<1e-6, sigma2[3] + 1e-6, 0)
}



y = births$scaled_log_births
x = births$day_of_year / max(births$day_of_year)
x_pred = runif(1000)

post = tidybayes::gather_draws(gp_coda, sigma2[i], l[i]) |>
  ungroup() |>
  summarize(post_mean = mean(.value), .by=c(i,.variable)) |>
  pivot_wider(id_cols = i, names_from = .variable, values_from = post_mean)

full = dukestm::cond_predict(
  y=y, x=x, x_pred=x_pred, cov_full, sigma2 = post$sigma2, l=post$l, per=7 / max(births$day_of_year)
)

comp_C1 = dukestm::cond_predict(
  y=y, x=x, x_pred=x_pred, sigma2 = post$sigma2, l=post$l, per=7 / max(births$day_of_year),
  cov_f_o  = cov_full,
  cov_f_p  = cov_C1,
  cov_f_po = cov_C1
)

comp_C2 = dukestm::cond_predict(
  y=y, x=x, x_pred=x_pred, sigma2 = post$sigma2, l=post$l, per=7 / max(births$day_of_year),
  cov_f_o  = cov_full,
  cov_f_p  = cov_C2,
  cov_f_po = cov_C2
)

pred = tibble(
  x_pred,
  day_of_year = x_pred * max(births$day_of_year),
  y_hat = apply(full,1,mean),
  Cov1 = apply(comp_C1, 1, mean),
  Cov2 = apply(comp_C2, 1, mean)
)


( 
  ggplot(pred, aes(x=day_of_year)) +
    geom_line(aes(y=y_hat), color='blue', linewidth=1.1) +
    geom_line(data=births, aes(y=scaled_log_births), alpha=0.3)
) / (
  ( 
    ggplot(pred, aes(x=day_of_year)) +
      geom_line(aes(y=Cov1), color='red') +
      geom_line(data=births, aes(y=scaled_log_births), alpha=0.3) 
  ) + (
    ggplot(pred, aes(x=day_of_year)) +
      geom_line(aes(y=Cov2), color='green') #+
      #geom_line(data=births, aes(y=scaled_log_births), alpha=0.3) 
  )
)








