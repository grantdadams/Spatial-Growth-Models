##### READ AND PREPARE DATA #####
##### RE AD AND PREPARE DATA #####
source("Weight_at_length_data_prep_v2_clean.R")

# Load data
load("StanFitLatRE_2018_02_19.RData")

# Get parameters
model <- StanFitLatRElog
library(rstan)
######################################################
# Get Linf, K, and t0
draws <- as.data.frame(model)
parameters_to_save <- c(paste("Betas_log_a[",1:3,"]", sep = ""), paste("Betas_log_b[",1:3,"]", sep = ""), "sigma_a", "sigma_b", "sigma", paste("log_a_re[",c(1:4,6:9), "]", sep = ""), paste("log_b_re[",c(1:4,6:9), "]", sep = ""))
draws <- draws[,which(colnames(draws) %in% parameters_to_save)]

# Get average lat
# put into a list to pass to BUGS
Predictors <- matrix(c(rep(1, length(dat$sex)), dat$sex, (dat$Lat - min(dat$Lat))), nrow = nrow(dat))

dataList = list(weight = dat$wet_weight_grams, region = dat$state_no, length = dat$FL_mm, n = nrow(dat), Jtrunc = 8, J = 9, Pred = Predictors, n_pred = ncol(Predictors), state = factor(dat$state, levels = c("TX", "LA", "MS", "AL", "FL_gulf", "FL_atlantic", "SC", "NC", "VA_ocean", "VA_bay")))

test <- c()

# Assess model fit using sum of squares
for(i in 1:nrow(draws)){
  # Get parameters
  # Get parameters
  #log_a_re <- as.numeric(draws[i,paste("log_a_re[",dataList$region,"]", sep = "")]) # Subset the mcmc chains
  betas_a_re <- as.numeric(draws[i,paste("Betas_log_a[",1:3,"]", sep = "")]) # Subset the mcmc chains
  a.sub <- exp(as.matrix(dataList$Pred) %*% as.vector(betas_a_re)) # + log_a_re)
  
  # log_b_re <- as.numeric(draws[i,paste("log_b_re[",dataList$region,"]", sep = "")]) # Subset the mcmc chains
  betas_b_re <- as.numeric(draws[i,paste("Betas_log_b[",1:3,"]", sep = "")]) # Subset the mcmc chains
  b.sub <- exp(as.matrix(dataList$Pred) %*% as.vector(betas_b_re)) # + log_b_re)
  
  # Calculate predicted and residuals
  pred.weight <- a.sub * dataList$length ^ b.sub 
  resid <- (dataList$weight) - pred.weight
  sum_sq <- resid^2
  
  #### PREDICTED BIT ####
  # generate new data and comput fit
  sigma <- as.numeric(draws[i,"sigma"])
  weight.new <- rnorm(n = length(dataList$weight), mean = (pred.weight), sd = sigma)
  resid.new <- (weight.new) - pred.weight
  sum_sq.new <- resid.new^2
  
  #### CALCULATE SUMMARY STATS ####
  fit <- sum(sum_sq)
  fit.new <- sum(sum_sq.new)
  
  test[i] <- ifelse(fit.new - fit > 0, 1, 0)
  print(i)
}

mean(test)

