 # Compare models using looAIC
library(loo)

# Load models
load("vbgf_model_progression2018_02_03.RData")


# Model 2
mean_params <- summary(mod_list[[2]])[[1]][,1]
dataList$predicted <- exp(dataList$Pred %*% mean_params[1:3]) * (1 - exp(-exp(dataList$Pred %*% mean_params[4:6]) * (dataList$age - (dataList$Pred %*% mean_params[7:9]))))
dataList$post_predicted <- rlnorm( length(dataList$predicted) , log(dataList$predicted), mean_params["sigma"])
plot(x = dataList$age, y = dataList$y, col = "grey", pch = 16, ylim = c(0, 650))
points(x = dataList$age, y = dataList$post_predicted)
points(x = dataList$age, y = dataList$predicted, col = 2, pch = 16)
dataList$residual <- dataList$predicted - dataList$y
boxplot(dataList$residual ~ dataList$state)

# Model 3
mean_params <- summary(mod_list[[3]])[[1]][,1]
dataListLong$predicted <- exp(dataListLong$Pred %*% mean_params[1:4]) * (1 - exp(-exp(dataListLong$Pred %*% mean_params[5:8]) * (dataListLong$age - (dataListLong$Pred %*% mean_params[9:12]))))
dataListLong$post_predicted <- rlnorm( length(dataListLong$predicted) , log(dataListLong$predicted), mean_params["sigma"])
plot(x = dataListLong$age, y = dataListLong$y, col = "grey", pch = 16, ylim = c(0, 650))
points(x = dataListLong$age, y = dataListLong$post_predicted)
points(x = dataListLong$age, y = dataListLong$predicted, col = 2, pch = 16)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)


# Model 4
mean_params <- summary(mod_list[[4]])[[1]][,1]
dataList$predicted <- exp(dataList$Pred %*% mean_params[1:3] + mean_params[paste("log_linf_re[",dataList$region,"]", sep = "")]) * (1 - exp(-exp(dataList$Pred %*% mean_params[4:6] + mean_params[paste("log_k_re[",dataList$region,"]", sep = "")]) * (dataList$age - (dataList$Pred %*% mean_params[7:9] + mean_params[paste("t0_re[",dataList$region,"]", sep = "")]))))
dataList$post_predicted <- rlnorm( length(dataList$predicted) , log(dataList$predicted), mean_params["sigma"])
plot(x = dataList$age, y = dataList$y, col = "grey", pch = 16, ylim = c(0, 650))
points(x = dataList$age, y = dataList$post_predicted)
points(x = dataList$age, y = dataList$predicted, col = 2, pch = 16)
dataList$residual <- dataList$predicted - dataList$y
boxplot(dataList$residual ~ dataList$state)

# Model 5
mean_params <- summary(mod_list[[5]])[[1]][,1]
dataListLong$predicted <- exp(dataListLong$Pred %*% mean_params[1:4] + mean_params[paste("log_linf_re[",dataListLong$region,"]", sep = "")]) * (1 - exp(-exp(dataListLong$Pred %*% mean_params[5:8] + mean_params[paste("log_k_re[",dataListLong$region,"]", sep = "")]) * (dataListLong$age - (dataListLong$Pred %*% mean_params[9:12] + mean_params[paste("t0_re[",dataListLong$region,"]", sep = "")]))))
dataListLong$post_predicted <- rlnorm( length(dataListLong$predicted) , log(dataListLong$predicted), mean_params["sigma"])
plot(x = dataListLong$age, y = dataListLong$y, col = "grey", pch = 16, ylim = c(0, 650))
points(x = dataListLong$age, y = dataListLong$post_predicted)
points(x = dataListLong$age, y = dataListLong$predicted, col = 2, pch = 16)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)

# Model 6
mean_params <- summary(mod_list[[6]])[[1]][,1]
dataList$predicted <- exp(dataList$Pred %*% mean_params[paste("Betas_log_linf[",1:3,"]", sep = "")] + mean_params[paste("log_linf_re[",dataList$region,"]", sep = "")]) * (1 - exp(-exp(dataList$Pred %*% mean_params[paste("Betas_log_k[",1:3,"]", sep = "")]  + mean_params[paste("log_k_re[",dataList$region,"]", sep = "")]) * (dataList$age - (dataList$Pred %*% mean_params[paste("Betas_t0[",1:3,"]", sep = "")] + mean_params[paste("log_t0_re[",dataList$region,"]", sep = "")]))))
dataList$post_predicted <- rlnorm( length(dataList$predicted) , log(dataList$predicted), mean_params["sigma"])
plot(x = dataList$age, y = dataList$y, col = "grey", pch = 16, ylim = c(0, 650))
points(x = dataList$age, y = dataList$post_predicted)
points(x = dataList$age, y = dataList$predicted, col = 2, pch = 16)
dataList$residual <- dataList$predicted - dataList$y
boxplot(dataList$residual ~ dataList$state)

# Model 7
mean_params <- summary(mod_list[[7]])[[1]][,1]
dataListLong$predicted <- exp(dataListLong$Pred %*% mean_params[paste("Betas_log_linf[",1:4,"]", sep = "")] + mean_params[paste("log_linf_re[",dataListLong$region,"]", sep = "")]) * (1 - exp(-exp(dataListLong$Pred %*% mean_params[paste("Betas_log_k[",1:4,"]", sep = "")]  + mean_params[paste("log_k_re[",dataListLong$region,"]", sep = "")]) * (dataListLong$age - (dataListLong$Pred %*% mean_params[paste("Betas_t0[",1:4,"]", sep = "")] + mean_params[paste("log_t0_re[",dataListLong$region,"]", sep = "")]))))
dataListLong$post_predicted <- rlnorm( length(dataListLong$predicted) , log(dataListLong$predicted), mean_params["sigma"])
plot(x = dataListLong$age, y = dataListLong$y, col = "grey", pch = 16, ylim = c(0, 700))
points(x = dataListLong$age, y = dataListLong$post_predicted)
points(x = dataListLong$age, y = dataListLong$predicted, col = 2, pch = 16)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)

# Model 8
mean_params <- summary(mod_list[[8]])[[1]][,1]
dataList$predicted <- exp(dataList$Pred %*% mean_params[paste("Betas_log_linf[",1:3,"]", sep = "")] + mean_params[paste("log_linf_re[",dataList$region,"]", sep = "")]) * (1 - exp(-exp(dataList$Pred %*% mean_params[paste("Betas_log_k[",1:3,"]", sep = "")]  + mean_params[paste("log_k_re[",dataList$region,"]", sep = "")]) * (dataList$age - (dataList$Pred %*% mean_params[paste("Betas_t0[",1:3,"]", sep = "")] + mean_params[paste("t0_re[",dataList$region,"]", sep = "")]))))
dataList$post_predicted <- rlnorm( length(dataList$predicted) , log(dataList$predicted), mean_params["sigma"])
plot(x = dataList$age, y = dataList$y, col = "grey", pch = 16, ylim = c(0, 650))
points(x = dataList$age, y = dataList$post_predicted)
points(x = dataList$age, y = dataList$predicted, col = 2, pch = 16)
dataList$residual <- dataList$predicted - dataList$y
boxplot(dataList$residual ~ dataList$state)

# Model 9
mean_params <- summary(mod_list[[9]])[[1]][,1]
dataListLong$predicted <- exp(dataListLong$Pred %*% mean_params[paste("Betas_log_linf[",1:4,"]", sep = "")] + mean_params[paste("log_linf_re[",dataListLong$region,"]", sep = "")]) * (1 - exp(-exp(dataListLong$Pred %*% mean_params[paste("Betas_log_k[",1:4,"]", sep = "")]  + mean_params[paste("log_k_re[",dataListLong$region,"]", sep = "")]) * (dataListLong$age - (dataListLong$Pred %*% mean_params[paste("Betas_t0[",1:4,"]", sep = "")] + mean_params[paste("t0_re[",dataListLong$region,"]", sep = "")]))))
dataListLong$post_predicted <- rlnorm( length(dataListLong$predicted) , log(dataListLong$predicted), mean_params["sigma"])
plot(x = dataListLong$age, y = dataListLong$y, col = "grey", pch = 16, ylim = c(0, 700))
points(x = dataListLong$age, y = dataListLong$post_predicted)
points(x = dataListLong$age, y = dataListLong$predicted, col = 2, pch = 16)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)
