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

dataList = list(log_y = log(dat$FL_mm), region = dat$state_no, age = dat$fractional_age, n = nrow(dat), Jtrunc = 10, J = 11, selectivity = (dat$sex.selectivity), W = Neighbor_Mat, D = diag(rowSums(Neighbor_Mat)), Pred = Predictors, n_pred = ncol(Predictors), W_n = sum(Neighbor_Mat)/2, state = factor(dat$state, levels = c("TX", "LA", "MS", "AL", "FL_gulf", "FL_atlantic", "SC", "NC", "VA_ocean", "VA_bay")))

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
#StanFitBase <- stan('STAN Files/vbgf_base.stan', data = dataList, iter = 3000, chains = nChains, verbose = FALSE, warmup = 1000)

##### RUN THE BASE MODEL IN STAN WITH LAT AND LONG #####
StanFitLat <- stan('STAN Files/vbgf_base_w_predictors.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))
loo_1 <- loo(extract_log_lik(StanFitLat))

##### RUN THE MODEL IN STAN WITH LAT AND LONG AND RANDOM EFFECTS #####
StanFitLatRE <- stan('STAN Files/vbgf_non_centered_v2.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))

##### RUN THE MODEL IN STAN WITH LAT AND LONG AND SPARSE CAR RANDOM EFFECTS #####
StanFitLatCAR <- stan('STAN Files/CAR Models/vbgf_car_non_centered_v3.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))

mod_list <- list(StanFitLat, StanFitLatRE, StanFitLatCAR)


file_name <- paste("3_vbgf_model_progression_",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(mod_list, file = file_name)
  

