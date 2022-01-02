#############################################################################
############################# Von Bertalanffy Sheepshead  ###################
#############################################################################

# clear the workspace
rm(list = ls(all = T))

# SET THE WORKING DIRECTORY TO THE LOCATION OF THIS FILE
# Session > Set Working Directory > To Source File Location

# load packages
library(rstan)
library(coda)
library(tidyr)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

##### RE AD AND PREPARE DATA #####
source("VGBM Data Prep/VGBM_Data_Prep.R")

# Make neighbor matrix
Neighbor_Mat <- matrix(0, nrow = 11, ncol = 11)
for(i in 1:ncol(Neighbor_Mat)){
  Neighbor_Mat[i,i+1] <- 1
  Neighbor_Mat[i+1,i] <- 1
}

# put into a list to pass to BUGS
Predictors <- matrix(c(rep(1, length(dat$sex)), dat$sex, (dat$Lat - min(dat$Lat))), nrow = nrow(dat))

dataList = list(y = dat$FL_mm, region = dat$state_no, age = dat$fractional_age, n = nrow(dat), Jtrunc = 10, J = 11, selectivity = (dat$sex.selectivity), W = Neighbor_Mat, D = diag(rowSums(Neighbor_Mat)), Pred = Predictors, n_pred = ncol(Predictors), W_n = sum(Neighbor_Mat)/2, state = factor(dat$state, levels = c("TX", "LA", "MS", "AL", "FL_gulf", "FL_atlantic", "SC", "NC", "VA_ocean", "VA_bay")))

# With longitude
Predictors <- matrix(c(rep(1, length(dat$sex)), dat$sex, (dat$Lat - min(dat$Lat)), (dat$Lon - min(dat$Lon))), nrow = nrow(dat))
dataListLong <- dataList
dataListLong$Pred <- Predictors
dataListLong$n_pred <- ncol(Predictors) 

##### MCMC DIMENSIONS #####
ni = 5000
burn = 5000
thin = 5
nChains = 1

##### RUN THE BASE MODEL IN STAN WITHOUT TN #####
StanFitBase <- stan('STAN Files/vbgf_base.stan', data = dataList, iter = 3000, chains = nChains, verbose = FALSE, warmup = 1000)

# Check to see if there are structure in the residuals
mean_params <- summary(StanFitBase)$summary[,1]
dataList$predicted <- exp(mean_params[1]) * (1 - exp(-exp(mean_params[2]) * (dataList$age - mean_params[3])))
plot(x = dataList$age, y = dataList$y, ylim = c(0, 600), col = "grey", pch = 16)
points(x = dataList$age, y = dataList$predicted)
dataList$residual <- dataList$predicted - dataList$y
boxplot(dataList$residual ~ dataList$Pred[,2])
boxplot(dataList$residual ~ dataList$Pred[,3])

plot(x = dataList$predicted, y = dataList$residual)


##### RUN THE BASE MODEL IN STAN WITH PREDICTORS #####
StanFitPred <- stan('STAN Files/vbgf_base_w_predictors.stan', data = dataList, iter = 3000, chains = nChains, verbose = FALSE, warmup = 1000)

# Check to see if there are spatial structure in the residuals
mean_params <- summary(StanFitPred)$summary[,1]
dataList$predicted <- exp(dataList$Pred %*% mean_params[1:3]) * (1 - exp(-exp(dataList$Pred %*% mean_params[4:6]) * (dataList$age - (dataList$Pred %*% mean_params[7:9]))))
plot(x = dataList$age, y = dataList$y, ylim = c(0, 600), col = "grey", pch = 16)
points(x = dataList$age, y = dataList$predicted)
dataList$residual <- dataList$predicted - dataList$y
boxplot(dataList$residual ~ dataList$state)
boxplot(dataList$residual ~ dat$Lon, xlab = "Lon")

# Judging from the residuals it seems that there is still some variation, no explained by the model

##### RUN THE BASE MODEL IN STAN WITH LAT AND LONG #####
StanFitLatLon <- stan('STAN Files/vbgf_base_w_predictors.stan', data = dataListLong, iter = 3000, chains = 1, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14))

# Check to see if there are spatial structure in the residuals
mean_params <- summary(StanFitLatLon)$summary[,1]
dataListLong$predicted <- exp(dataListLong$Pred %*% mean_params[1:4]) * (1 - exp(-exp(dataListLong$Pred %*% mean_params[5:8]) * (dataListLong$age - (dataListLong$Pred %*% mean_params[9:12]))))
plot(x = dataListLong$age, y = dataListLong$y, col = "grey", pch = 16, ylim = c(0, 600))
points(x = dataListLong$age, y = dataListLong$predicted)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)


##### RUN THE MODEL IN STAN WITH LAT AND LONG AND RANDOM EFFECTS #####
StanFitLatRE <- stan('STAN Files/vbgf_non_centered_v2.stan', data = dataList, iter = 3000, chains = 1, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14))

StanFitLatLonRE <- stan('STAN Files/vbgf_non_centered_v2.stan', data = dataListLong, iter = 3000, chains = 1, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14))

mean_params <- summary(StanFitLatLon)$summary[,1]
dataList$predicted <- exp(dataListLong$Pred %*% mean_params[1:4]) * (1 - exp(-exp(dataListLong$Pred %*% mean_params[5:8]) * (dataListLong$age - (dataListLong$Pred %*% mean_params[9:12]))))
plot(x = dataListLong$age, y = dataListLong$y, col = "grey", pch = 16, ylim = c(0, 600))
points(x = dataListLong$age, y = dataListLong$predicted)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)

##### RUN THE MODEL IN STAN WITH LAT AND LONG AND SPARSE CAR RANDOM EFFECTS #####
StanFitLatCAR <- stan('STAN Files/CAR Models/vbgf_car_sparse_v2.stan', data = dataList, iter = 3000, chains = 1, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14))

StanFitLatLonCAR <- stan('STAN Files/CAR Models/vbgf_car_sparse_v2.stan', data = dataListLong, iter = 3000, chains = 1, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14))

##### RUN THE MODEL IN STAN WITH LAT AND LONG AND NON-CENTERED CAR RANDOM EFFECTS #####
StanFitLatCARnc <- stan('STAN Files/CAR Models/vbgf_car_non_centered_v3.stan', data = dataList, iter = 3000, chains = 1, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14))

StanFitLatLonCARnc <- stan('STAN Files/CAR Models/vbgf_car_non_centered_v3.stan', data = dataListLong, iter = 3000, chains = 1, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14))


##### SAVE MODEL LIST #####
mod_list <- list(StanFitBase, StanFitPred, StanFitLatLon, StanFitLatRE, StanFitLatLonRE, StanFitLatCAR, StanFitLatLonCAR, StanFitLatCARnc, StanFitLatLonCARnc)


file_name <- paste("vbgf_model_progression",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(mod_list, file = file_name)


