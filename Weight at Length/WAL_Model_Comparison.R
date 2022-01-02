library(loo)
library(rstan)

# load models
load("StanFitPredlog_2018_03_31.RData")
load("StanFitLatRElog_short_2018_03_31.RData")
load("StanFitCARlog_short_2018_04_02.RData")

loglike1 <- extract_log_lik(StanFitPred)
loglike2 <- extract_log_lik(StanFitLatRElog)
loglike3 <- extract_log_lik(StanFitLatCAR)

loo1 <- loo(loglike1)
loo2 <- loo(loglike2)
loo3 <- loo(loglike3)

compare(loo1, loo2)
compare(loo2, loo3)

waic1 <- waic(loglike1)
waic2 <- waic(loglike2)
waic3 <- waic(loglike3)

