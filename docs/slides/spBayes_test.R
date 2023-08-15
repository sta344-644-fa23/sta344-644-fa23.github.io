d = readRDS("data/frn_example.rds")


#fit_spbayes = function(formula, data, coords, )

n_batch = 200
batch_len = 200
n_samples = n_batch * batch_len
burnin = n_samples*0.5 + 1

coords = cbind(d$day, 0)
n = nrow(coords)

m = spBayes::spLM(
  pm25~day+I(day^2), data=d, coords=coords,
  starting=list(
    "beta"=rep(0,3), "phi"=sqrt(3)/7, "sigma.sq"=1, "tau.sq"=1
  ),
  tuning=list(
    "phi"=1, "sigma.sq"=0.5, "tau.sq"=0.5
  ),
  priors=list(
    "phi.unif"=c(sqrt(3)/100, sqrt(3)/1),
    "sigma.sq.ig"=c(2, 6),
    "tau.sq.ig"=c(2, 12)
  ),
  amcmc=list("n.batch"=n_batch, "batch.length"=batch_len, "accept.rate"=0.43),
  cov.model="gaussian", verbose=TRUE, n.report=10
)
  
plot(m$p.theta.samples)

#m_rec = spBayes::spRecover(m, start = burnin, thin = 50, get.w = FALSE)

#plot(m_rec$p.beta.recover.samples)


df_p = tibble (
  day = seq(0,365) + rnorm(366,mean=0, sd=0.01)
) %>%
  model.matrix(~day+I(day^2), data=.)
p = spBayes::spPredict(
  m, pred.coords = cbind(df_p$day, 0), 
  pred.covars = df_p, 
  start = burnin, thin=20
)


plot(df_p$day, rowMeans(p$p.y.predictive.samples), type='l', ylim=c(0,15), col='blue')
lines(d$day, d$pm25, col='red')
lines(df_p$day, apply(p$p.y.predictive.samples,1, quantile, probs=0.025))
lines(df_p$day, apply(p$p.y.predictive.samples,1, quantile, probs=0.975))
 



