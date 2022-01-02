#############################################################################
############################# Von Bertalanffy Sheepshead  ###################
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

##### RE AD AND PREPARE DATA #####
source("Data/Data prep/VBGM_Data_prep.R")

# Make neighbor matrix
neighbor_mat <- matrix(0, nrow = 11, ncol = 11)
for(i in 1:ncol(neighbor_mat)){
  neighbor_mat[i,i+1] <- 1
  neighbor_mat[i+1,i] <- 1
} # Will get an error for +1 but that is OK

# put into a list to pass to Stan
design_mat <- matrix(c(rep(1, length(dat$sex)), dat$sex, (dat$Lat - min(dat$Lat))), nrow = nrow(dat))

dataList = list(
  n_i = nrow(dat), # Sample size
  J = 11, # Number of regions
  n_pred = ncol(design_mat), # Number of predictors
  log_length_i = log(dat$FL_mm), # Vector of log fork length of fish i
  age_i = dat$fractional_age, # Vector of age of fish i
  design_mat = design_mat, # n_i * n_pred design matrix of linear predictors
  r_i = dat$state_no, # Integer vector of region of fish i
  W = neighbor_mat, # Weights matrix where W_ij = 1 when regions i and j are neighbors
  D = diag(rowSums(neighbor_mat)) # Matrix where diagonal is number of neighbors of region i
)


##### MCMC DIMENSIONS #####
ni = 5000
burn = 5000
thin = 5
nChains = 1

##### RUN THE BASE MODEL IN STAN WITH LAT #####
StanFitLat <- stan('Stan models/VBGF_Model_1.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))
loo_1 <- loo(extract_log_lik(StanFitLat))

##### RUN THE MODEL IN STAN WITH LAT AND RANDOM EFFECTS #####
StanFitLatRE <- stan('Stan models/VBGF_Model_2.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))

##### RUN THE MODEL IN STAN WITH LAT AND SPARSE CAR RANDOM EFFECTS #####
StanFitLatCAR <- stan('Stan models/VBGF_Model_3.stan', data = dataList, iter = 5000, chains = 2, cores = 2, verbose = FALSE, warmup = 1000, control = list(max_treedepth = 14, adapt_delta = 0.9))

##### SAVE #####
mod_list <- list(StanFitLat, StanFitLatRE, StanFitLatCAR)
file_name <- paste("Stan models/Models/VBGF_Stan_models",gsub("-","_",as.Date(trunc(Sys.time(),"day"))),".RData", sep = "")
save(mod_list, file = file_name)


