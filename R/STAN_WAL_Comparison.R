library(loo)
library(rstan)

# load models
load("WAL_Stan_models_2018_03_31.RData")

library(rstan)
library(loo)

log_lik1 <- extract_log_lik(mod_list[[1]])
log_lik2 <- extract_log_lik(mod_list[[2]])
log_lik3 <- extract_log_lik(mod_list[[3]])


loo1 <- loo(log_lik1)
loo2 <- loo(log_lik2)
loo3 <- loo(log_lik3)


diff1 <- compare(loo1, loo2)
diff2 <- compare(loo2, loo3)
diff3 <- compare(loo1, loo3)

diff1
diff2
diff3

waic(log_lik1)
waic(log_lik2)
waic(log_lik3)
