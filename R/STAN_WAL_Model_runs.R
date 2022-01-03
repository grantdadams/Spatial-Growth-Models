#############################################################################
############################# Weight-at-length Sheepshead  ##################
#############################################################################

# clear the workspace
rm(list = ls(all = T))

# load packages
library(rstan)
library(coda)
library(tidyr)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

##### READ AND PREPARE DATA #####
source("Data/Data prep/WAL_Data_prep.R")

# Make neighbor matrix
neighbor_mat <- matrix(0, nrow = 9, ncol = 9)
for(i in 1:ncol(neighbor_mat)){
  neighbor_mat[i,i+1] <- 1
  neighbor_mat[i+1,i] <- 1
}

# put into a list to pass to Stand
design_mat <- matrix(c(rep(1, length(dat$sex)), dat$sex, (dat$Lat - min(dat$Lat))), nrow = nrow(dat))

dataList = list(
  n_i = nrow(dat), # Sample size
  n_r = 9, # Number of regions
  n_pred = ncol(design_mat), # Number of predictors
  log_length_i = log(dat$FL_mm), # Vector of log fork length of fish i
  log_weight_i = log(dat$wet_weight_grams), # Vector of age of fish i
  design_mat = design_mat, # n_i * n_pred design matrix of linear predictors
  r_i = dat$state_no, # Integer vector of region of fish i
  W = neighbor_mat, # Weights matrix where W_ij = 1 when regions i and j are neighbors
  D = diag(rowSums(neighbor_mat)) # Matrix where diagonal is number of neighbors of region i
)

##### MCMC DIMENSIONS #####
ni = 5000
burn = 1000
nChains = 2   


##### RUN THE BASE MODEL IN STAN WITH PREDICTORS #####
StanFitPred <- stan('Stan models/WAL_Model_1.stan', data = dataList, iter = ni, chains = nChains, cores = 2, verbose = TRUE, warmup = burn, control = list(max_treedepth = 14, adapt_delta = 0.9))

##### RUN THE MODEL IN STAN WITH LAT AND RANDOM EFFECTS #####
StanFitLatRE <- stan('Stan models/WAL_Model_2.stan', data = dataList, iter = ni, chains = nChains, cores = 2, verbose = TRUE, warmup = burn, control = list(max_treedepth = 14, adapt_delta = 0.9))

##### RUN THE MODEL IN STAN WITH LAT AND SPARSE CAR RANDOM EFFECTS #####
StanFitLatCAR <- stan('Stan models/WAL_Model_3.stan', data = dataList, iter = ni, chains = nChains, cores = 2, verbose = TRUE, warmup = burn, control = list(max_treedepth = 14, adapt_delta = 0.9))


 ##### SAVE MODEL LIST #####
mod_list <- list(StanFitPred, StanFitLatRE, StanFitLatCAR)
file_name <- paste("Stan models/Models/WAL_Stan_models_",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(mod_list, file = file_name)


