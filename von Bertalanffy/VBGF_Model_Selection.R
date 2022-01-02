 # Compare models using looAIC
library(loo)

# Load models
load("vbgf_model_progression2018_02_03.RData")

# Model 4
mean_params <- summary(mod_list[[4]])[[1]][,1]
dataList$predicted <- exp(dataList$Pred %*% mean_params[1:3] + mean_params[paste("log_linf_re[",dataList$region,"]", sep = "")]) * (1 - exp(-exp(dataList$Pred %*% mean_params[4:6] + mean_params[paste("log_k_re[",dataList$region,"]", sep = "")]) * (dataList$age - (dataList$Pred %*% mean_params[7:9] + mean_params[paste("t0_re[",dataList$region,"]", sep = "")]))))
plot(x = dataList$age, y = dataList$y, col = "grey", pch = 16, ylim = c(0, 600))
points(x = dataList$age, y = dataList$predicted)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)

# Model 5
mean_params <- summary(mod_list[[4]])[[1]][,1]
dataList$predicted <- exp(dataList$Pred %*% mean_params[1:3] + mean_params[paste("log_linf_re[",dataList$region,"]", sep = "")]) * (1 - exp(-exp(dataList$Pred %*% mean_params[4:6] + mean_params[paste("log_k_re[",dataList$region,"]", sep = "")]) * (dataList$age - (dataList$Pred %*% mean_params[7:9] + mean_params[paste("t0_re[",dataList$region,"]", sep = "")]))))
plot(x = dataList$age, y = dataList$y, col = "grey", pch = 16, ylim = c(0, 600))
points(x = dataList$age, y = dataList$predicted)
dataListLong$residual <- dataListLong$predicted - dataListLong$y
boxplot(dataListLong$residual ~ dataListLong$state)
