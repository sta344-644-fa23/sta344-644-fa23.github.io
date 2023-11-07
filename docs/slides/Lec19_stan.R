# https://github.com/mbjoseph/CARstan
#
# https://discourse.mc-stan.org/t/correcting-spatial-autocorrelation-using-sar-different-results-with-brms-vs-spatialreg-packages/23269/3

## ?brms:::posterior_predict_gaussian_lagsar

niter <- 1E4   # definitely overkill, but good for comparison
nchains <- 4
  
car_m = rstan::stan_model(
  model_code = "
    data {
      int<lower=1> n;
      int<lower=1> p;
      matrix[n, p] X;
      vector[n] y;
      matrix<lower=0, upper=1>[n, n] W;
      matrix[n, n] D;
    }
    parameters {
      vector[p] beta;
      real<lower=0> tau;
      real<lower=0, upper=1> alpha;
    }
    transformed parameters {
      vector[n] y_cond = X * beta + (D - alpha * W) * (y - X * beta);
    }
    model {
      y ~ multi_normal_prec(X * beta, tau * (D - alpha * W));
      beta ~ normal(0, 1);
      tau ~ gamma(2, 2);
    }
  "
)

X = model.matrix(~1, data = nc)

d = list(
  n = nrow(X),         # number of observations
  p = ncol(X),         # number of coefficients
  X = X,               # design matrix
  y = 1000*nc$SID74/nc$BIR74,               # observed number of cases
  W = A*1,
  D = diag(rowSums(A))
)   

s = rstan::sampling(car_m, data = d, cores = 4)

gather_draws(s, yp[i]) |>
  left_join(
    nc |>
      as_tibble() |>
      transmute(
        rate = 1000*SID74/BIR74,
        i = row_number()
      )
  ) |>
  filter(.iteration %in% 1:2) |>
  ggplot(aes(x=rate, y=.value, color=.iteration)) +
    geom_point() +
    facet_wrap(~.iteration) +
    geom_abline(slope=1, intercept=0, color="black", linetype="dashed")




# Hacky test

rlang::env_unlock(env = asNamespace('brms'))
rlang::env_binding_unlock(env = asNamespace('brms'))
patch = function (i, prep, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(prep, diag(nobs) - ac$lagsar[s] * ac$Msar)
    mu <- as.numeric( mu[s, ] )
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    brms:::rmulti_normal(1, mu = mu, Sigma = Sigma)
  }
  mu <- brms:::get_dpar(prep, "mu")
  sigma <- brms:::get_dpar(prep, "sigma")
  brms:::rblapply(seq_len(prep$ndraws), .predict)
}
assign('posterior_predict_gaussian_lagsar', patch, envir = asNamespace('brms'))
rlang::env_binding_lock(env = asNamespace('brms'))
rlang::env_lock(asNamespace('brms'))

predict(b_sar)

