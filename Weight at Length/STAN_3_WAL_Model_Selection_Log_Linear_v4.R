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
source("Weight_at_length_data_prep_v2_clean_v2.R")

# Make neighbor matrix
Neighbor_Mat <- matrix(0, nrow = 9, ncol = 9)
for(i in 1:ncol(Neighbor_Mat)){
  Neighbor_Mat[i,i+1] <- 1
  Neighbor_Mat[i+1,i] <- 1
}

# put into a list to pass to BUGS
  Predictors <- matrix(c(rep(1, length(dat$sex)), dat$sex, (dat$Lat - min(dat$Lat))), nrow = nrow(dat))

dataList = list(log_weight = log(dat$wet_weight_grams), weight = (dat$wet_weight_grams), region = dat$state_no, log_length = log(dat$FL_mm), length = (dat$FL_mm), n = nrow(dat), Jtrunc = 8, J = 9, W = Neighbor_Mat, D = diag(rowSums(Neighbor_Mat)), Pred = Predictors, n_pred = ncol(Predictors), W_n = sum(Neighbor_Mat)/2, state = factor(dat$state, levels = c("TX", "LA", "MS", "AL", "FL_gulf", "FL_atlantic", "SC", "NC", "VA_ocean", "VA_bay")))


##### MCMC DIMENSIONS #####
ni = 5000
burn = 5000
thin = 5
nChains = 1   


##### RUN THE BASE MODEL IN STAN WITH PREDICTORS #####
StanFitPred <- stan('STAN Files/wal_base_w_predictors_log_linear.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))
file_name <- paste("StanFitPredlog_",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(StanFitPred, file = file_name)


##### RUN THE MODEL IN STAN WITH LAT AND LONG AND RANDOM EFFECTS #####
dataListRELog = dataList
dataListRELog$J = 8
dataListRELog$region = dat$state_no_check

StanFitLatRElog <- stan('STAN Files/wal_non_centered_log_linear.stan', data = dataListRELog, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))
file_name <- paste("StanFitLatRElog_short_",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(StanFitLatRElog, file = file_name)

##### RUN THE MODEL IN STAN WITH LAT AND LONG AND SPARSE CAR RANDOM EFFECTS #####
StanFitLatCAR <- stan('STAN Files/CAR Models/wal_car_non_centered_log_linear_v2.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))
file_name <- paste("StanFitCARlog_short_",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(StanFitLatCAR, file = file_name)


 ##### SAVE MODEL LIST #####
mod_list <- list(StanFitPred, StanFitLatRElog, StanFitLatCAR)


file_name <- paste("wal_3_model_progression_log_linear",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(mod_list, file = file_name)


