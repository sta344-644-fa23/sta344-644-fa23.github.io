

# Back to our example 

## RMSE

```{r}
rmse = df_pred %>%
  group_by(.iteration) %>%
  summarize(rmse = sqrt( sum( (y - y_pred)^2 ) / n())) %>%
  pull(rmse)

length(rmse)

head(rmse)

mean(rmse)

modelr::rmse(l, data = d)
```

## RMSE ($\mu$)

```{r}
rmse = df_pred %>%
  group_by(.iteration) %>%
  summarize(rmse = sqrt( sum( (y - mu)^2 ) / n())) %>%
  pull(rmse)

length(rmse)

head(rmse)

mean(rmse)

modelr::rmse(l, data = d)
```

## CRPS

```{r}
crps = df_pred %>% 
  group_by(i) %>% 
  summarise(crps = calc_crps(y_pred, y)) %>%
  pull(crps)

length(crps)

head(crps)

mean(crps)
```


## Empirical Coverage



```{r}
df_cover = df_pred %>% 
  group_by(x,y) %>% 
  tidybayes::mean_hdi(y_pred, .prob = c(0.5,0.9,0.95)) 

df_cover %>%
  mutate(contains = y >= .lower & y <= .upper) %>%
  group_by(prob=.width) %>%
  summarize(emp_cov = sum(contains)/n())
```

## Posterior predictive distribution ($y_{pred}$)

```{r echo=TRUE, fig.height=4}
df_pred %>% ungroup() %>%
  ggplot(aes(x=x)) + 
  tidybayes::stat_lineribbon(aes(y=y_pred), alpha=0.5) +
  geom_point(data=d, aes(y=y))
```


# Compared to what?

## Polynomial fit

```{r}
ggplot(d, aes(x=x,y=y)) + 
  geom_line() + 
  geom_smooth(
    method='lm', color="blue", se = FALSE, 
    formula = y~poly(x,5,simple=TRUE)
  )
```

## Model

```{r}
l_p = lm(y~poly(x,5,simple=TRUE), data=d)
summary(l_p)
```

## JAGS Model

```{r}
poly_model = 
  "model{
  # Likelihood
  for(i in 1:length(y)){
    y[i] ~ dnorm(mu[i], tau)
    y_pred[i] ~ dnorm(mu[i], tau)
    mu[i] = beta[1]        + beta[2]*x[i]   + beta[3]*x[i]^2 +
            beta[4]*x[i]^3 + beta[5]*x[i]^4 + beta[6]*x[i]^5
  }

  # Prior for beta
  for(j in 1:6){
    beta[j] ~ dnorm(0,1/1000)
  }

  # Prior for sigma / tau2
  tau ~ dgamma(1, 1)
  sigma2 = 1/tau
}"
```

```{r echo=FALSE, message=FALSE}
n_burn = 1000; n_iter = 5000

m = rjags::jags.model(
  textConnection(poly_model), data=d, 
  quiet=TRUE, n.chains = 1
) 
update(m, n.iter=n_burn, progress.bar="none")

df_poly = rjags::coda.samples(
  m, variable.names=c("beta","sigma2","mu","y_pred","y","x"), 
  n.iter=n_iter, progress.bar="none"
) %>%
  tidybayes::spread_draws(y_pred[i], y[i], x[i], mu[i]) %>%
  mutate(resid = y - mu)
```

## Posterior Predictive Distribution

```{r}
df_poly %>% ungroup() %>%
  ggplot(aes(x=x)) + 
  tidybayes::stat_lineribbon(aes(y=y_pred), alpha=0.5) +
  geom_point(data=d, aes(y=y))
```

## Comparing Results

```{r echo=FALSE}
rmse = rbind(
  mutate(df_pred, model="y~x"),
  mutate(df_poly, model="y~poly(x,5)")
) %>% 
  group_by(model,.iteration) %>%
  transmute(
    y_pred = (y-y_pred)^2  %>% mean() %>% sqrt(),
    mu     = (y-mu    )^2  %>% mean() %>% sqrt()
  ) %>% 
  tidyr::gather(parameter, rmse, -.iteration, -model) %>%
  group_by(parameter, model) %>%
  summarize(rmse = mean(rmse))
```

```{r echo=FALSE}
crps = rbind(
  mutate(df_pred, model="y~x"),
  mutate(df_poly, model="y~poly(x,5)")
) %>% 
  group_by(model,i) %>% 
  summarise(crps = calc_crps(y_pred, y)) %>%
  select(model,i,crps) %>%
  group_by(model) %>%
  summarize(crps = mean(crps)) %>%
  mutate(parameter = "y_pred")
```


```{r echo=FALSE}
empc = rbind(
  mutate(df_pred, model="y~x"),
  mutate(df_poly, model="y~poly(x,5)")
) %>% 
  group_by(x,y,model) %>% 
  tidybayes::mean_hdi(y_pred, .width = c(0.9)) %>%
  mutate(contains = y >= .lower & y <= .upper) %>%
  group_by(model) %>%
  summarize("emp_cov (90%)" = sum(contains)/n()) %>%
  mutate(parameter = "y_pred")
```

```{r echo=FALSE}
full_join(rmse, crps) %>% 
  full_join(empc) %>%
  tidyr::gather(metric, value, -parameter, -model) %>%
  tidyr::spread(model, value) %>%
  select(2,1,4,3) %>% 
  mutate(metric = factor(metric, levels = c("rmse", "crps", "emp_cov (90%)"))) %>%
  arrange(metric) %>%
  knitr::kable(digits = 3)
```