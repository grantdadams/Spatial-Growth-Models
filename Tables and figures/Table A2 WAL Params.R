library(rstan)

##### RE AD AND PREPARE DATA #####
source("Data/Data prep/WAL_Data_prep.R")
################################################
# reorder state variables
state_no_df <- subset(dat, !duplicated(state_no_check))
state_no_df <- state_no_df[,c("state_no", "state_no_check", "state")]
state_no_df <- state_no_df[order(state_no_df$state_no),]
state_no_df$state_no_2 <- 1:8

# Load data
load("Stan models/Models/WAL_Stan_models_2018_02_19.RData")

# Get parameters
model <- StanFitLatRElog
mod_summ <- data.frame(summary(model)[[1]])
mod_summ$parameter <- rownames(mod_summ)

# Subset what we want
mod_summ <- mod_summ[,c("parameter","mean","X50.", "X2.5.", "X97.5.", "n_eff" , "Rhat")]
parameters_to_save <- c(paste("Betas_log_a[",1:3,"]", sep = ""), paste("Betas_log_b[",1:3,"]", sep = ""), "sigma_a", "sigma_b", "sigma", paste("log_a_re[",c(state_no_df$state_no_check), "]", sep = ""), paste("log_b_re[",c(state_no_df$state_no_check), "]", sep = ""))
mod_summ_sub <- mod_summ[which(mod_summ$parameter %in% parameters_to_save),]
mod_summ_sub$sex <- "NA"
mod_summ_sub$state_no_check <- "NA"
mod_summ_sub <- mod_summ_sub[match(parameters_to_save, mod_summ_sub$parameter),]
mod_summ_sub$order <- 1:nrow(mod_summ_sub)

order_ind <- nrow(mod_summ_sub) + 1


######################################################
# Get Linf, K, and t0
draws <- as.data.frame(model)
# Get average lat
lat_lon_df <- subset(dat, !duplicated(Lat))
lat_lon_df <- lat_lon_df[, c("state_no_check", "Lat")]
lat_lon_df$Lat <- (lat_lon_df$Lat - min(lat_lon_df$Lat))




######################################################
# a
a_names <- paste("log_a_re[", 1:11,"]", sep = "")
a_rows <- which(mod_summ$parameter %in% a_names)
a_df <- mod_summ[a_rows,]
m_a_df <- a_df
f_a_df <- a_df
m_a_df$sex <- "m"
f_a_df$sex <- "f"
m_a_df$order <- "NA"
f_a_df$order <- "NA"

for(i in c(state_no_df$state_no_check)){
  log_a_re <- draws[,paste("log_a_re[",i,"]", sep = "")]
  betas_a_re <- draws[,paste("Betas_log_a[",1:3,"]", sep = "")]
  
  F_a <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, 1, lat_lon_df[which(lat_lon_df$state_no_check == i),2])) + log_a_re)
  M_a <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, 2, lat_lon_df[which(lat_lon_df$state_no_check == i),2])) + log_a_re)
  
  # Get mean
  f_a_df[i,"mean"] <- mean(F_a)
  m_a_df[i,"mean"] <- mean(M_a)
  
  # Get mean
  f_a_df[i,3:5] <- as.numeric(quantile(F_a, probs = c(.5, 0.025, 0.975)))
  m_a_df[i,3:5] <- as.numeric(quantile(M_a, probs = c(.5, 0.025, 0.975)))
  
  f_a_df$state_no_check[i] <- 
  m_a_df$state_no_check[i] <- i
  
  f_a_df$order[i] <- order_ind
  order_ind <- order_ind + 1
  m_a_df$order[i] <- order_ind
  order_ind <- order_ind + 1
}

# # Remove greorgia
# f_a_df <- f_a_df[-c(5, 10, 11),]
# m_a_df <- m_a_df[-c(5, 10, 11),]

###############################################
# K
b_names <- paste("log_b_re[", 1:11,"]", sep = "")
b_rows <- which(mod_summ$parameter %in% b_names)
b_df <- mod_summ[b_rows,]
m_b_df <- b_df
f_b_df <- b_df
m_b_df$sex <- "m"
f_b_df$sex <- "f"
f_b_df$order <- "NA"
m_b_df$order <- "NA"

for(i in c(state_no_df$state_no_check)){
  log_b_re <- draws[,paste("log_b_re[",i,"]", sep = "")]
  betas_b_re <- draws[,paste("Betas_log_b[",1:3,"]", sep = "")]
  
  F_b <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, 1, lat_lon_df[which(lat_lon_df$state_no_check == i),2])) + log_b_re)
  M_b <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, 2, lat_lon_df[which(lat_lon_df$state_no_check == i),2])) + log_b_re)
  
  # Get mean
  f_b_df[i,"mean"] <- mean(F_b)
  m_b_df[i,"mean"] <- mean(M_b)
  
  # Get mean
  f_b_df[i,3:5] <- as.numeric(quantile(F_b, probs = c(.5, 0.025, 0.975)))
  m_b_df[i,3:5] <- as.numeric(quantile(M_b, probs = c(.5, 0.025, 0.975)))
  
  f_b_df$state_no_check[i] <- i
  m_b_df$state_no_check[i] <- i
  
  f_b_df$order[i] <- order_ind
  order_ind <- order_ind + 1
  m_b_df$order[i] <- order_ind
  order_ind <- order_ind + 1
}

# # Remove greorgia
# f_b_df <- f_b_df[-c(5, 10, 11),]
# m_b_df <- m_b_df[-c(5, 10, 11),]

##############################################################################
# Treatment for saving
# Combine
derived_params <- rbind(f_a_df, m_a_df, f_b_df, m_b_df)
mod_summ_sub <- rbind(mod_summ_sub, derived_params)

# Round values
#mod_summ_sub$Rhat <- round(mod_summ_sub$Rhat, 2)
#mod_summ_sub$n_eff <- round(mod_summ_sub$n_eff, 0)
#mod_summ_sub[,2:5] <- round(mod_summ_sub[,2:5], 3)

mod_summ_sub <- mod_summ_sub[,c(1,9,8,2,3,4,5,6,7,10)]

library(dplyr)
# mod_summ_sub <- mod_summ_sub %>% arrange(order) 
mod_summ_sub <- mod_summ_sub[order(as.numeric(mod_summ_sub$order)),]

write.csv(mod_summ_sub, file = "Tables and Figures/Table_A2_Parameter_Summ.csv")




