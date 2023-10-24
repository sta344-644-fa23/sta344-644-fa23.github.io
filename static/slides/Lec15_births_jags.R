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
    sigma2[i] ~ dt(0, 2.5, 1) T(0,)
  }
  for(i in 1:3)
    l[i] ~ dt(0, 2.5, 1) T(0,)
  }

  alpha ~ dt(0, 2.5, 1) T(0,)
}"

dir.create("results/", showWarnings = FALSE, recursive = TRUE)

if (!file.exists("results/init_model.rds")) {
  m = rjags::jags.model(
    textConnection(model), 
    data = list(
      y = births$scaled_log_births,
      d = dist(births$day_of_year / max(births$day_of_year)) %>% as.matrix(),
      per = 7 / max(births$day_of_year),
      pi = pi,
      N = nrow(births)
    ),
    n.adapt=1000
  )
  saveRDS(m, file="results/init_model.rds")
} else {
  m = readRDS("results/init_model.rds")
}