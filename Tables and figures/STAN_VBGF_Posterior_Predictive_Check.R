##### READ AND PREPARE DATA #####
source("VGBM Data Prep/VGBM_Data_Prep.R")
prev_coefs <- read.csv("VBGF_Coefs_Previous_Studies.csv")

# Load model
load("3_vbgf_model_progression_2018_04_01.RData")

# Get parameters
model <- mod_list[[2]]

######################################################
# Get Linf, K, and t0
draws <- as.data.frame(model)
# Get average lat
# put into a list to pass to BUGS
Predictors <- matrix(c(rep(1, length(dat$sex)), dat$sex, (dat$Lat - min(dat$Lat))), nrow = nrow(dat))

dataList = list(log_y = log(dat$FL_mm), region = dat$state_no, age = dat$fractional_age, n = nrow(dat), Jtrunc = 10, J = 11, selectivity = (dat$sex.selectivity), Pred = Predictors, n_pred = ncol(Predictors), state = factor(dat$state, levels = c("TX", "LA", "MS", "AL", "FL_gulf", "FL_atlantic", "SC", "NC", "VA_ocean", "VA_bay")))

test <- c()

# Assess model fit using sum of squares
for(i in 1:nrow(draws)){
  # Get parameters
  log_linf_re <- as.numeric(draws[i,paste("log_linf_re[",dataList$region,"]", sep = "")]) # Subset the mcmc chains
  betas_linf_re <- as.numeric(draws[i,paste("Betas_log_linf[",1:3,"]", sep = "")]) # Subset the mcmc chains
  linf.sub <- exp(as.matrix(cbind(dataList$Pred[,1], dataList$Pred[,2], dataList$Pred[,3])) %*% as.vector(betas_linf_re) + log_linf_re)
  
  log_k_re <- as.numeric(draws[i,paste("log_k_re[",dataList$region,"]", sep = "")]) # Subset the mcmc chains
  betas_k_re <- as.numeric(draws[i,paste("Betas_log_k[",1:3,"]", sep = "")]) # Subset the mcmc chains
  k.sub <- exp(as.matrix(cbind(dataList$Pred[,1], dataList$Pred[  ,2], dataList$Pred[,3])) %*% as.vector(betas_k_re) + log_k_re)
  
  log_t0_re <- as.numeric(draws[i,paste("t0_re[",dataList$region,"]", sep = "")]) # Subset the mcmc chains
  betas_t0_re <- as.numeric(draws[i,paste("Betas_t0[",1:3,"]", sep = "")]) # Subset the mcmc chains
  t0.sub <- as.matrix(cbind(dataList$Pred[,1], dataList$Pred[,2], dataList$Pred[,3])) %*% as.vector(betas_t0_re) + log_t0_re  
  
  # Calculate predicted and residuals
  pred.length <- log(linf.sub * (1 - exp (-k.sub * (dataList$age - t0.sub))))
  resid <- (dataList$log_y) - pred.length
  sum_sq <- resid^2
  
  #### PREDICTED BIT ####
  # generate new data and comput fit
  sigma <- as.numeric(draws[i,"sigma"])
  y.new <- rnorm(n = length(dataList$age), mean = (pred.length), sd = sigma)
  resid.new <- (y.new) - pred.length
  sum_sq.new <- resid.new^2
  
  #### CALCULATE SUMMARY STATS ####
  fit <- sum(sum_sq)
  fit.new <- sum(sum_sq.new)
  
  test[i] <- ifelse(fit.new - fit > 0, 1, 0)
  print(i)
}

mean(test)


for(i in unique(dat$state_no)){
  subs <- which(dat$state_no == i)
  plot(y = exp(dataList$log_y[subs]), x = (dataList$age[subs]) )
  points(y = exp(y.new[subs]), x = (dataList$age[subs]), col = 2 )
  legend("topright", legend = i)
}
