library(tidyverse)

d = readRDS("data/avg_temp_df.rds") %>%
  slice(1:(52*3)) %>%
  mutate(week = as.numeric(date - date[1])/7)

print.spLM = function(x, ...) {
  print("spLM model object")
  invisible(x)
}

#fit_spbayes = function(formula, data, coords, )

chains = 4
verbose = FALSE
n_batch = 200
batch_len = 100
n_samples = n_batch * batch_len
burnin = n_samples*0.5 + 1
thin = 10
cov_model = "gaussian"

coords = cbind(d$week, 0)
n = nrow(coords)

get_mcmc = function(
    m,
    burnin_frac = 0.5,
    thin = 1, verbose = FALSE
) {
  n_samples = nrow(m$model$p.theta.samples)
  burnin = floor(n_samples*burnin_frac + 1)
  
  rec = spBayes::spRecover(
    m$model, 
    start = burnin, thin = thin, 
    get.w = FALSE, 
    verbose = verbose
  )
  
  posterior::bind_draws(
    posterior::as_draws_matrix(m$model$p.theta.samples) %>%
      posterior::subset_draws(iteration = burnin:n_samples) %>%
      posterior::thin_draws(thin = thin),
    posterior::as_draws_matrix(rec$p.beta.recover.samples),
    along = "variable"
  )
}


fit_model = function(
    formula,
    data,
    coords,
    starting,
    priors,
    tuning = list(
      "phi"=1, "sigma.sq"=1, "tau.sq"=1
    ),
    n_batch = 200,
    batch_len = 100,
    cov_model = "gaussian",
    accept_rate = 0.43,
    verbose = FALSE
) {
  args = as.list(environment())
  
  list(
    model = spBayes::spLM(
      formula, 
      data = data, 
      coords = coords,
      starting = starting,
      priors = priors,
      tuning = tuning,
      amcmc = list(
        "n.batch"=n_batch, "batch.length"=batch_len, "accept.rate"=accept_rate
      ),
      cov.model = cov_model,
      verbose = verbose
    ),
    args = args
  )
}



gplm = function(
  formula,
  data,
  coords,
  chains = 4,
  ...,
  thin = 1
) {
  l = lapply(
    seq_len(chains),
    function(i) {
      message("Fitting chain ", i)
      fit_model(formula=formula, data=data, coords=coords, ...)
    }
  )
  
  mcmc = lapply(
    l, get_mcmc, thin=thin
  ) %>%
    purrr::reduce(posterior::bind_draws, along="chain")
  
  structure(
    list(
      models = l,
      mcmc = mcmc
    ),
    class = "gplm_fit"
  )
}

z = gplm(
  avg_temp~1, data = d, coords = cbind(d$week, 0),
  starting=list(
    "phi"=sqrt(3)/4, "sigma.sq"=1, "tau.sq"=1
  ),
  tuning=list(
    "phi"=1, "sigma.sq"=1, "tau.sq"=1
  ),
  priors=list(
    "phi.unif"=c(sqrt(3)/52, sqrt(3)/1),
    "sigma.sq.ig"=c(2, 1),
    "tau.sq.ig"=c(2, 1)
  ),
  thin=10
)

plot.gplm_fit = function(x, combo = c("dens", "trace"), ...) {
  bayesplot::mcmc_combo(x$mcmc, combo = combo, ...)  
}

summary.gplm_fit = function(object, ...) {
  summary(object$mcmc)
}

predict.gplm_fit = function(
  object,
  newdata,
  coords,
  burnin_frac = 0.5,
  thin = 1,
  verbose = FALSE
) {
  lapply(
    object$models,
    function(m) {
      n_samples = nrow(m$model$p.theta.samples)
      burnin = floor(n_samples*burnin_frac + 1)
      n_y = nrow(newdata)
      
      covars = model.matrix(
        m$args$formula %>% update(NULL ~ .),
        newdata
      )
      
      p = spBayes::spPredict(
        m$model, 
        pred.coords = coords, 
        pred.covars = covars,
        start = burnin, thin = thin,
        verbose = FALSE
      )
      post = p$p.y.predictive.samples
      rownames(post) = stringr::str_glue("y[{seq_len(n_y)}]")
      
      t(post) %>%
        posterior::as_draws()
    }
  ) %>%
    purrr::reduce(posterior::bind_draws, along="chain")
}


newdata = tibble (
  week = seq(0,3.5*52)
) 

coords = cbind(
  newdata$week %>% jitter(),
  0 
)

pred = predict(z, newdata=newdata, coords = coords, thin=10)


saveRDS(z, "Lec15_temp_m.rds")
saveRDS(pred, "Lec15_temp_pred.rds")