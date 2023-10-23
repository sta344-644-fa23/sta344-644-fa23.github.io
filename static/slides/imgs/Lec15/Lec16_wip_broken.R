# Mauna Loa Example

## Atmospheric CO$_2$

library(tidyverse)

loa = data_frame(
  y = co2 |> dukestm::strip_attrs(),
  x = time(co2) |> dukestm::strip_attrs()
)

noaa = readr::read_csv("data/noaa_mauna_loa.csv") |>
  transmute(x = year+(month-1)/12, y=average) |>
  filter(x > max(loa$x))

rbind(
  loa |> mutate(Source="Scripps (co2 in R)"),
  noaa |> mutate(Source="NOAA")
) |>
  ggplot(aes(x=x,y=y, color=Source)) + 
  geom_line() 


ml_model = "model{
  y ~ dmnorm(mu, inverse(Sigma))

  for (i in 1:(length(y)-1)) {
    for (j in (i+1):length(y)) {
      k1[i,j] <- pow(theta[1],2) * exp(-0.5 * pow( d[i,j] / theta[2], 2))
      k2[i,j] <- pow(theta[3],2) * exp(-0.5 * pow(d[i,j] / theta[4],2) - 2 * pow(sin(pi*d[i,j] / per) / theta[5], 2))
      k3[i,j] <- pow(theta[6],2) * pow(1 + 0.5 * pow(d[i,j] / theta[7], 2)/theta[8], -theta[8])
      k4[i,j] <- pow(theta[9],2) * exp(-0.5 * pow(d[i,j] / theta[10],2))
      
      Sigma[i,j] <- k1[i,j] + k2[i,j] + k3[i,j] + k4[i,j]
      Sigma[j,i] <- Sigma[i,j]
    }
  }

  for (i in 1:length(y)) {
    Sigma[i,i] <- pow(theta[1],2) + pow(theta[3],2) + pow(theta[6],2) + pow(theta[9],2) + pow(theta[11],2)
  }  

  for(i in 1:11){
    theta[i] ~ dt(0, 2.5, 1) T(0,)
  }
}"

dir.create("results", showWarnings = FALSE)

if (!file.exists("results/init_model.rds")) {
  m = rjags::jags.model(
    textConnection(ml_model), 
    data = list(
      mu = mean(loa$y) * rep(1, nrow(loa)),
      y = loa$y,
      #x = loa$x,
      d = dist(loa$x) |> as.matrix(),
      per = 1,
      pi = pi
    ),
    inits = list(
      theta = c(66,67,2.4,90,1.3,0.66,1.2,0.78,0.18,1.6,0.19)
    ),
    n.adapt=10000
  )
  saveRDS(m, file="results/init_model.rds")
} else {
  m = readRDS("results/init_model.rds")
}

if (!file.exists("results/fit_model.Rdata")) {
  gp_coda = coda.samples(
    m, variable.names=c("sigma2", "l", "alpha"),
    n.iter=10000,
    thin=10
  )
  save(gp_coda, file="results/fit_model.Rdata")
} else {
  load("results/fit_model.Rdata")
}


sigma2 = get_coda_parameter(gp_coda, "sigma")
l      = get_coda_parameter(gp_coda, "^l")
alpha  = get_coda_parameter(gp_coda, "alpha")



## Diagnostics

r = tidybayes::gather_samples(gp_coda, sigma2[i], l[i], alpha) |>
  ungroup() |>
  mutate(term = fixterm(term, i))

r |>  
  mutate(term = forcats::as_factor(term)) |>
  group_by(term) |>
  slice(seq(1,n(),length.out = 300)) |>
  ggplot(aes(x=.iteration, y=estimate)) +
  geom_line() +
  facet_wrap(~term, scale="free_y", ncol=5)



## Fit Components

```{r}
#| echo: false
if (!file.exists("results/ml_gp_fit_comp.Rdata"))
{
  x_pred = loa$x
  y_fit_comp = lapply(1:4, function(i) matrix(NA, nrow=n_post_samp, ncol=length(x_pred)))
  
  mu = rep(mean(y), length(x))
  mu_pred = rep(0, length(x_pred))
  
  dist_o = rdist(x)
  dist_p = rdist(x_pred)
  dist_op = rdist(x, x_pred)
  
  p = progress_estimated(n_post_samp)
  for(i in 1:nrow(sigma2))
  {
    cov_o  = calc_cov(dist_o, i)$S
    inv = solve(cov_o)
    
    for(j in 1:4)
    {
      cov_name = paste0("S",j)
      cov_p  = calc_cov(dist_p, i)[[cov_name]]
      cov_op = calc_cov(dist_op, i)[[cov_name]]
      
      cond_cov = cov_p - t(cov_op) %*% inv %*% cov_op
      diag(cond_cov) = diag(cond_cov) + 1e-4
      cond_mu  = mu_pred + t(cov_op) %*% inv %*% (y - mu)
      
      y_fit_comp[[j]][i,] = cond_mu + t(chol(cond_cov)) %*% matrix(rnorm(length(x_pred)), ncol=1)
    }
    p$tick()$print()
  }
  
  save(y_fit_comp, file="results/ml_gp_fit_comp.Rdata")
} else {
  load(file="results/ml_gp_fit_comp.Rdata")
}
```

```{r}
#| echo: false
lapply(
  1:4, 
  function(i) 
    post_summary(y_fit_comp[[i]]) |> 
    mutate(cov=paste0("Sigma_",i), x=loa[["x"]])
) |> 
  bind_rows() |>
  ggplot(aes(x=x, y=post_mean, col=cov, fill=cov)) +
  geom_ribbon(aes(ymin=post_lower, ymax=post_upper, color=NULL), alpha=0.2) +
  geom_line() +
  facet_wrap(~cov, scale="free_y")
```




## Forecasting

```{r}
#| echo: false
y = loa$y
x = loa$x

n_post_samp = 1000

calc_cov = function(d,i)
{
  S1 = sigma2[i,1] * exp(- (l[i,1] * d)^2)
  S2 = sigma2[i,2] * exp(- (l[i,2] * d)^2 - 2 * (l[i,3] * sin(pi*d))^2)
  S3 = sigma2[i,3] * (1+(l[i,4] * d)^2/alpha[i,1])^-alpha[i,1]
  S4 = sigma2[i,4] * exp(- (l[i,5] * d)^2)
  S5 = ifelse(abs(d)<1e-4, sigma2[i,5] + 1e-6, 0)
  S = S1 + S2 + S3 + S4 + S5
  
  list(S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S=S)
}
```

```{r}
#| echo: false
if (!file.exists("results/ml_gp_fore.Rdata"))
{
  x_pred = noaa$x
  y_fore = matrix(NA, nrow=n_post_samp, ncol=length(x_pred))
  colnames(y_fore) = paste0("Y_pred[", 1:length(x_pred), "]")
  
  mu = rep(mean(y), length(x))
  mu_pred = rep(mean(y), length(x_pred))
  
  dist_o = fields::rdist(x)
  dist_p = fields::rdist(x_pred)
  dist_op = fields::rdist(x, x_pred)
  
  for(i in 1:nrow(sigma2))
  {
    cov_o  = calc_cov(dist_o, i)$S
    cov_p  = calc_cov(dist_p, i)$S
    cov_op = calc_cov(dist_op, i)$S
    
    inv = solve(cov_o, cov_op)
    cond_cov = cov_p - t(cov_op) %*% inv
    diag(cond_cov) = diag(cond_cov) + 1e-5
    cond_mu  = mu_pred + t(inv) %*% (y - mu)
    
    y_fore[i,] = cond_mu + t(chol(cond_cov)) %*% matrix(rnorm(length(x_pred)), ncol=1)
  }
  
  save(y_fore, file="results/ml_gp_fore.Rdata")
} else {
  load(file="results/ml_gp_fore.Rdata")
}
```

```{r}
#| echo: false
summ_80 = post_summary(y_fore, ci_width = 0.8) |> mutate(x = noaa[["x"]])
summ_95 = post_summary(y_fore, ci_width = 0.95) |> mutate(x = noaa[["x"]])

ggplot(summ_95, aes(x=x, y=post_mean)) +
  geom_ribbon(data=summ_95, aes(ymin=post_lower, ymax=post_upper), fill="#D3D2D6") +
  geom_ribbon(data=summ_80, aes(ymin=post_lower, ymax=post_upper), fill="#A3A5C3") +
  geom_line(color="#1007FC") +
  geom_line(data=loa, aes(x=x,y=y), alpha=0.5) +
  geom_line(data=noaa, aes(x=x,y=y), alpha=0.5, col='red')
```

## Forecasting (zoom)

```{r}
#| echo: false
ggplot(summ_95, aes(x=x, y=post_mean)) +
  geom_ribbon(data=summ_95, aes(ymin=post_lower, ymax=post_upper), fill="#D3D2D6") +
  geom_ribbon(data=summ_80, aes(ymin=post_lower, ymax=post_upper), fill="#A3A5C3") +
  geom_line(color="#1007FC", size=0.8) +
  geom_line(data=loa, aes(x=x,y=y), alpha=0.5, size=0.8) +
  geom_line(data=noaa, aes(x=x,y=y), alpha=0.5, col='red', size=0.8) +
  xlim(1998, 2017) + ylim(360,410)
```

## Forecasting ARIMA (auto)


```{r}
#| echo: false
m5 = forecast::auto.arima(co2)
m5_fore = forecast::forecast(m5,h = nrow(noaa))

autoplot(m5_fore) +
  geom_line(data=noaa, aes(x=x,y=y), alpha=0.5, col='red', size=0.8)
```



## Forecasting RMSE

```{r}
#| echo: false
noaa_pred_arima = data_frame(
  y_hat = m5_fore$mean |> strip_attrs(),
  x = time(m5_fore$mean) |> strip_attrs()
) |> 
  mutate(y = noaa[["y"]])

noaa_pred_gp = summ_95 |>
  select(y_hat = post_mean, x) |>
  mutate(y = noaa[["y"]])

filters = c(2003, 2008, 2013, 2018)

noaa_rmse = function(year, d) 
{
  d |> 
    filter(x < year) |> 
    {(.$y - .$y_hat)^2} |> 
    mean() |> 
    sqrt() |> 
    round(3)
}

data_frame(
  "dates"        = c("Jan 1998 - Jan 2003", "Jan 1998 - Jan 2008", "Jan 1998 - Jan 2013", "Jan 1998 - Mar 2017"),
  "RMSE (arima)" = purrr::map_dbl(filters, noaa_rmse, d=noaa_pred_arima),
  "RMSE (gp)"    = purrr::map_dbl(filters, noaa_rmse, d=noaa_pred_gp)
) |>
  knitr::kable()
```




## Forecasting Components

```{r}
#| echo: false
if (!file.exists("results/ml_gp_fore_comp.Rdata"))
{
  x_pred = noaa$x
  y_comp = lapply(1:4, function(i) matrix(NA, nrow=n_post_samp, ncol=length(x_pred)))
  
  mu = rep(mean(y), length(x))
  mu_pred = rep(0, length(x_pred))
  
  dist_o = rdist(x)
  dist_p = rdist(x_pred)
  dist_op = rdist(x, x_pred)
  
  p = progress_estimated(n_post_samp)
  for(i in 1:nrow(sigma2))
  {
    cov_o  = calc_cov(dist_o, i)$S
    inv = solve(cov_o)
    
    for(j in 1:4)
    {
      cov_name = paste0("S",j)
      cov_p  = calc_cov(dist_p, i)[[cov_name]]
      cov_op = calc_cov(dist_op, i)[[cov_name]]
      
      cond_cov = cov_p - t(cov_op) %*% inv %*% cov_op
      diag(cond_cov) = diag(cond_cov) + 1e-4
      cond_mu  = mu_pred + t(cov_op) %*% inv %*% (y - mu)
      
      y_comp[[j]][i,] = cond_mu + t(chol(cond_cov)) %*% matrix(rnorm(length(x_pred)), ncol=1)
    }
    p$tick()$print()
  }
  
  save(y_comp, file="results/ml_gp_fore_comp.Rdata")
} else {
  load(file="results/ml_gp_fore_comp.Rdata")
}
```

```{r}
#| echo: false
lapply(
  1:4, 
  function(i) 
    post_summary(y_comp[[i]]) |> 
    mutate(cov=paste0("Sigma_",i), x=noaa[["x"]])
) |> 
  bind_rows() |>
  ggplot(aes(x=x, y=post_mean, col=cov, fill=cov)) +
  geom_ribbon(aes(ymin=post_lower, ymax=post_upper, color=NULL), alpha=0.2) +
  geom_line() +
  facet_wrap(~cov, scale="free_y")
```

