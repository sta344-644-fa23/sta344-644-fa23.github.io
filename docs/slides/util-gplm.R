plot.gplm_fit = function(x, combo = c("dens", "trace"), ...) {
  bayesplot::mcmc_combo(x$mcmc, combo = combo, ...)
}

print.gplm_fit = function(x, ...) {
  nchains = posterior::nchains(m$mcmc)
  ndraws = posterior::ndraws(m$mcmc)
  nvar = posterior::nvariables(m$mcmc)
  
  cat( stringr::str_glue(
    "# A gplm model (spBayes spLM) with {nchains} chains, {nvar} variables, and {ndraws} iterations."
  ), "\n" )
  summary(x)
}

summary.gplm_fit = function(object, ...) {
  print(summary(object$mcmc))
}

print.spLM = function(x, ...) {
  print("spLM model object")
  invisible(x)
}
