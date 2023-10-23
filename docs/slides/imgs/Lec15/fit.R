library(tidyverse)
library(sta344)

ggplot2::theme_set(ggplot2::theme_bw())

get_coda_parameter = function(coda, pattern) {
  w = coda[[1]] %>% colnames() %>% stringr::str_detect(pattern)
  coda[[1]][,w,drop=FALSE]
}

post_summary = function(m, ci_width=0.95) {
  d = data_frame(
    post_mean  = apply(m, 2, mean),
    post_med   = apply(m, 2, median),
    post_lower = apply(m, 2, quantile, probs=(1-ci_width)/2),
    post_upper = apply(m, 2, quantile, probs=1 - (1-ci_width)/2)
  )
  
  if (!is.null(colnames(m)))
    d = d %>% mutate(param = colnames(m)) %>% select(param,post_mean:post_upper)
  
  d
}


loa = data_frame(
  y = co2 %>% strip_attrs(),
  x = time(co2) %>% strip_attrs()
)

noaa = readr::read_csv("noaa_mauna_loa.csv") %>%
  transmute(x = year+(month-1)/12, y=average) %>%
  filter(x > max(loa$x))

load("fit_model.Rdata")
load("ml_gp_fit_comp.Rdata")

sigma2 = get_coda_parameter(gp_coda, "sigma")
l      = get_coda_parameter(gp_coda, "^l")
alpha  = get_coda_parameter(gp_coda, "alpha")

y = loa$y
x = loa$x

n_post_samp = 1000

calc_cov = function(d,i) {
  S1 = sigma2[i,1] * exp(- (l[i,1] * d)^2)
  S2 = sigma2[i,2] * exp(- (l[i,2] * d)^2 - 2 * (l[i,3] * sin(pi*d))^2)
  S3 = sigma2[i,3] * (1+(l[i,4] * d)^2/alpha[i,1])^-alpha[i,1]
  S4 = sigma2[i,4] * exp(- (l[i,5] * d)^2)
  S5 = ifelse(abs(d)<1e-4, sigma2[i,5] + 1e-6, 0)
  S = S1 + S2 + S3 + S4 + S5
  
  list(S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S=S)
}

x_pred = c(loa$x |> jitter(), noaa$x)
y_fore = matrix(NA, nrow=n_post_samp, ncol=length(x_pred))
colnames(y_fore) = paste0("Y_pred[", 1:length(x_pred), "]")

mu = rep(mean(y), length(x))
mu_pred = rep(mean(y), length(x_pred))

dist_o = fields::rdist(x)
dist_p = fields::rdist(x_pred)
dist_op = fields::rdist(x, x_pred)

for(i in 1:nrow(sigma2)) {
  cov_o  = calc_cov(dist_o, i)$S
  cov_p  = calc_cov(dist_p, i)$S
  cov_op = calc_cov(dist_op, i)$S
  
  inv = solve(cov_o, cov_op)
  cond_cov = cov_p - t(cov_op) %*% inv
  diag(cond_cov) = diag(cond_cov) + 1e-5
  cond_mu  = mu_pred + t(inv) %*% (y - mu)
  
  y_fore[i,] = cond_mu + t(chol(cond_cov)) %*% matrix(rnorm(length(x_pred)), ncol=1)
}

summ_80 = post_summary(y_fore, ci_width = 0.8) %>% mutate(x = x_pred)
summ_95 = post_summary(y_fore, ci_width = 0.95) %>% mutate(x = x_pred)

summ_80_fit = summ_80 |> filter(x < 1998)
summ_95_fit = summ_95 |> filter(x < 1998)


g1 = ggplot(summ_95, aes(x=x, y=post_mean)) +
  geom_ribbon(data=summ_95, aes(ymin=post_lower, ymax=post_upper), fill="#D3D2D6") +
  geom_ribbon(data=summ_80, aes(ymin=post_lower, ymax=post_upper), fill="#A3A5C3") +
  geom_line(color="#1007FC") +
  geom_line(data=loa, aes(x=x,y=y), alpha=0.5) +
  geom_line(data=noaa, aes(x=x,y=y), alpha=0.5, col='red')

ggsave("fit_forecast.png", g1, width = 9, height = 6)

g2 = ggplot(summ_95_fit, aes(x=x, y=post_mean)) +
  geom_ribbon(data=summ_95_fit, aes(ymin=post_lower, ymax=post_upper), fill="#D3D2D6") +
  geom_ribbon(data=summ_80_fit, aes(ymin=post_lower, ymax=post_upper), fill="#A3A5C3") +
  geom_line(color="#1007FC") +
  geom_line(data=loa, aes(x=x,y=y), alpha=0.5)

ggsave("fit.png", g2, width = 9, height = 6)
