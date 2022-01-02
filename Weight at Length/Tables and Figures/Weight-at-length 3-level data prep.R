###################################################################################################
############################# Von Bertalanffy Sheepshead Plots ################ ###################
###################################################################################################

# Load the data
setwd("~/Dropbox/2017_Sheepshead_Growth_Variation/Sheepshead R Work/Weight-at-length")

# Load previous parameter estimates
prev_coefs <- read.csv("WAL_Prev_Param_Estimates.csv")

load("C:/Users/Punt Lab/Dropbox/2017_Sheepshead_Growth_Variation/Sheepshead R Work/Weight-at-length/JAGS_weight_at_length_lognormal_3_level_no_dic_2017_09_30.RData")
source("Weight_at_length_data_prep.R")
setwd(paste(getwd(),"/Figures", sep = ""))


# Get the JAGs parameters
library(runjags)
library(coda)
summ = summary(runJagsOut)
# 
# 
a.est <- summ[substr(rownames(summ), 1, 5) == "a.est",c("Mean", "Median", "Lower95", "Upper95", "SSeff","psrf")]
b.est <- summ[substr(rownames(summ), 1, 5) == "b.est",c("Mean", "Median", "Lower95", "Upper95", "SSeff","psrf")]
# 
mu.a.est <- summ[substr(rownames(summ), 1, 9) == "mu.a.est[",c("Mean", "Median", "Lower95", "Upper95", "SSeff","psrf")]
mu.b.est <- summ[substr(rownames(summ), 1, 9) == "mu.b.est[",c("Mean", "Median", "Lower95", "Upper95", "SSeff","psrf")]
# 
grand.a.est <- summ[substr(rownames(summ), 1, 11) == "grand.a.est",c("Mean", "Median", "Lower95", "Upper95", "SSeff","psrf")]
grand.b.est <- summ[substr(rownames(summ), 1, 11) == "grand.b.est",c("Mean", "Median", "Lower95", "Upper95","SSeff", "psrf")]

state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS"), state_no = c(1:length(unique(dat$state))))

# Print results
results <- signif(rbind(grand.a.est, grand.b.est, mu.a.est, mu.b.est, data.frame(a.est)[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16),], data.frame(b.est)[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16),] ),3)
results$param <- rownames(results)
results$sex <- c(NA,NA,rep(c("f","m"),18))
results$region <- c(rep(NA,6), as.character(rep(state.id$state, each = 2, time = 2)))
results <- results[,c("param","sex","region","Mean", "Median", "Lower95", "Upper95", "SSeff","psrf")]
write.csv(results, file = "Weight-at-length parameters.csv")


# Create the MCMC for each parameter
mcmc.jags <- as.mcmc(runJagsOut)
mcmc.jags <- data.frame((mcmc.jags ))

#mcmc.jags <- data.frame(rbind(mcmc.jags[[1]], mcmc.jags[[2]], mcmc.jags[[3]] ))
a.mcmc <- mcmc.jags[,substr(colnames(mcmc.jags), 1, 5) == "a.est"]
b.mcmc <- mcmc.jags[,substr(colnames(mcmc.jags), 1, 5) == "b.est"]

mu.a.mcmc <- mcmc.jags[,substr(colnames(mcmc.jags), 1, 8) == "mu.a.est"]
mu.b.mcmc <- mcmc.jags[,substr(colnames(mcmc.jags), 1, 8) == "mu.b.est"]

grand.a.mcmc <- mcmc.jags[,substr(colnames(mcmc.jags), 1, 11) == "grand.a.est"]
grand.b.mcmc <- mcmc.jags[,substr(colnames(mcmc.jags), 1, 11) == "grand.b.est"]

grand.a.raw.mcmc <- mcmc.jags[,substr(colnames(mcmc.jags), 1, 15) == "grand.a.est.raw"]
# Create table of results



