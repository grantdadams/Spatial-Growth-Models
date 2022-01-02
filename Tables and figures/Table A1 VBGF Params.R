library(rstan)

##### RE AD AND PREPARE DATA #####
source("VGBM Data Prep/VGBM_Data_Prep.R")

# Load data
load("3_vbgf_model_progression_2018_04_01.RData")

# Get parameters
model <- mod_list[[2]]
mod_summ <- data.frame(summary(model)[[1]])
mod_summ$parameter <- rownames(mod_summ)

# Subset what we want
mod_summ <- mod_summ[,c("parameter","mean","X50.", "X2.5.", "X97.5.", "n_eff" , "Rhat")]
parameters_to_save <- c(paste("log_linf_re[",c(1:4,6:11), "]", sep = ""), paste("log_k_re[",c(1:4,6:11), "]", sep = ""), paste("t0_re[",c(1:4,6:11), "]", sep = ""), paste("Betas_log_linf[",1:3,"]", sep = ""), paste("Betas_log_k[",1:3,"]", sep = ""), paste("Betas_t0[",1:3,"]", sep = ""), "sigma_linf", "sigma_k", "sigma_t0", "sigma")
mod_summ_sub <- mod_summ[which(mod_summ$parameter %in% parameters_to_save),]
mod_summ_sub$sex <- "NA"
mod_summ_sub$state_no <- "NA"
mod_summ_sub <- mod_summ_sub[match(parameters_to_save, mod_summ_sub$parameter),]
mod_summ_sub$order <- 1:nrow(mod_summ_sub)

order_ind <- nrow(mod_summ_sub) + 1


######################################################
# Get Linf, K, and t0
draws <- as.data.frame(model)
# Get average lat
lat_lon_df <- subset(dat, !duplicated(Lat))
lat_lon_df <- lat_lon_df[, c("state_no", "Lat")]
lat_lon_df$Lat <- (lat_lon_df$Lat - min(lat_lon_df$Lat))

######################################################
# Linf
linf_names <- paste("log_linf_re[", 1:11,"]", sep = "")
linf_rows <- which(mod_summ$parameter %in% linf_names)
Linf_df <- mod_summ[linf_rows,]
m_linf_df <- Linf_df
f_linf_df <- Linf_df
m_linf_df$sex <- "m"
f_linf_df$sex <- "f"
m_linf_df$order <- "NA"
f_linf_df$order <- "NA"

for(i in c(1:4,6:11)){
  log_linf_re <- draws[,paste("log_linf_re[",i,"]", sep = "")]
  betas_linf_re <- draws[,paste("Betas_log_linf[",1:3,"]", sep = "")]
  
  F_Linf <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, 0, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_linf_re)
  M_Linf <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, 1, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_linf_re)
  
  # Get mean
  f_linf_df[i,"mean"] <- mean(F_Linf)
  m_linf_df[i,"mean"] <- mean(M_Linf)
  
  # Get mean
  f_linf_df[i,3:5] <- as.numeric(quantile(F_Linf, probs = c(.5, 0.025, 0.975)))
  m_linf_df[i,3:5] <- as.numeric(quantile(M_Linf, probs = c(.5, 0.025, 0.975)))
  
  f_linf_df$state_no[i] <- i
  m_linf_df$state_no[i] <- i
  
  f_linf_df$order[i] <- order_ind
  order_ind <- order_ind + 1
  m_linf_df$order[i] <- order_ind
  order_ind <- order_ind + 1
}

# Remove greorgia
f_linf_df <- f_linf_df[-5,]
m_linf_df <- m_linf_df[-5,]

###############################################
# K
k_names <- paste("log_k_re[", 1:11,"]", sep = "")
k_rows <- which(mod_summ$parameter %in% k_names)
k_df <- mod_summ[k_rows,]
m_k_df <- k_df
f_k_df <- k_df
m_k_df$sex <- "m"
f_k_df$sex <- "f"
f_k_df$order <- "NA"
m_k_df$order <- "NA"

for(i in c(1:4,6:11)){
  log_k_re <- draws[,paste("log_k_re[",i,"]", sep = "")]
  betas_k_re <- draws[,paste("Betas_log_k[",1:3,"]", sep = "")]
  
  F_k <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, 0, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_k_re)
  M_k <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, 1, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_k_re)
  
  # Get mean
  f_k_df[i,"mean"] <- mean(F_k)
  m_k_df[i,"mean"] <- mean(M_k)
  
  # Get mean
  f_k_df[i,3:5] <- as.numeric(quantile(F_k, probs = c(.5, 0.025, 0.975)))
  m_k_df[i,3:5] <- as.numeric(quantile(M_k, probs = c(.5, 0.025, 0.975)))
  
  f_k_df$state_no[i] <- i
  m_k_df$state_no[i] <- i
  
  f_k_df$order[i] <- order_ind
  order_ind <- order_ind + 1
  m_k_df$order[i] <- order_ind
  order_ind <- order_ind + 1
}

# Remove greorgia
f_k_df <- f_k_df[-5,]
m_k_df <- m_k_df[-5,]

###############################################################################
# t0
t0_names <- paste("t0_re[", 1:11,"]", sep = "")
t0_rows <- which(mod_summ$parameter %in% t0_names)
t0_df <- mod_summ[t0_rows,]
m_t0_df <- t0_df
f_t0_df <- t0_df
m_t0_df$sex <- "m"
f_t0_df$sex <- "f"
f_t0_df$order <- "NA"
m_t0_df$order <- "NA"

for(i in c(1:4,6:11)){
  log_t0_re <- draws[,paste("t0_re[",i,"]", sep = "")]
  betas_t0_re <- draws[,paste("Betas_t0[",1:3,"]", sep = "")]
  
  F_t0 <- (as.matrix(betas_t0_re) %*% as.vector(c(1, 0, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_t0_re)
  M_t0 <- (as.matrix(betas_t0_re) %*% as.vector(c(1, 1, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_t0_re)
  
  # Get mean
  f_t0_df[i,"mean"] <- mean(F_t0)
  m_t0_df[i,"mean"] <- mean(M_t0)
  
  # Get mean
  f_t0_df[i,3:5] <- as.numeric(quantile(F_t0, probs = c(.5, 0.025, 0.975)))
  m_t0_df[i,3:5] <- as.numeric(quantile(M_t0, probs = c(.5, 0.025, 0.975)))
  
  f_t0_df$state_no[i] <- i
  m_t0_df$state_no[i] <- i
  
  f_t0_df$order[i] <- order_ind
  order_ind <- order_ind + 1
  m_t0_df$order[i] <- order_ind
  order_ind <- order_ind + 1
}

# Remove greorgia
f_t0_df <- f_t0_df[-5,]
m_t0_df <- m_t0_df[-5,]

###############################################################################
# W
k_names <- paste("log_k_re[", 1:11,"]", sep = "")
linf_names <- paste("log_linf_re[", 1:11,"]", sep = "")

w_df <- data.frame(parameter = paste("w[",1:11,"]", sep = ""), mean = rep(NA, 11), x50. = rep(NA,11), x2.5. = rep(NA, 11), x97.5. = rep(NA, 11), n_eff = rep(NA, 11), Rhat = rep(NA, 11))
m_w_df <- w_df
f_w_df <- w_df
m_w_df$sex <- "m"
f_w_df$sex <- "f"
m_w_df$order <- "NA"
f_w_df$order <- "NA"


for(i in c(1:4,6:11)){
  # Get linf
  log_linf_re <- draws[,paste("log_linf_re[",i,"]", sep = "")]
  betas_linf_re <- draws[,paste("Betas_log_linf[",1:3,"]", sep = "")]
  
  F_Linf <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, 0, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_linf_re)
  M_Linf <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, 1, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_linf_re)
  
  # Get k
  log_k_re <- draws[,paste("log_k_re[",i,"]", sep = "")]
  betas_k_re <- draws[,paste("Betas_log_k[",1:3,"]", sep = "")]
  
  F_k <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, 0, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_k_re)
  M_k <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, 1, lat_lon_df[which(lat_lon_df$state_no == i),2])) + log_k_re)
  
  # Calculate w
  F_w <- F_Linf * F_k
  M_w <- M_Linf * M_k
  
  # Get mean
  f_w_df[i,"mean"] <- mean(F_w)
  m_w_df[i,"mean"] <- mean(M_w)
  
  # Get mean
  f_w_df[i,3:5] <- as.numeric(quantile(F_w, probs = c(.5, 0.025, 0.975)))
  m_w_df[i,3:5] <- as.numeric(quantile(M_w, probs = c(.5, 0.025, 0.975)))
  
  f_w_df$state_no[i] <- i
  m_w_df$state_no[i] <- i
  
  f_w_df$order[i] <- order_ind
  order_ind <- order_ind + 1
  m_w_df$order[i] <- order_ind
  order_ind <- order_ind + 1
}

# Remove greorgia
f_w_df <- f_w_df[-5,]
m_w_df <- m_w_df[-5,]

colnames(f_w_df) <- colnames(m_linf_df)
colnames(m_w_df) <- colnames(m_linf_df)

# Combine
derived_params <- rbind(f_linf_df, m_linf_df, f_k_df, m_k_df, f_t0_df, m_t0_df, f_w_df, m_w_df)
mod_summ_sub <- rbind(mod_summ_sub, derived_params)

##############################################################################
# Treatment for saving

# Round values
mod_summ_sub$n_eff <- round(mod_summ_sub$n_eff, 0)
mod_summ_sub[,2:5] <- round(mod_summ_sub[,2:5], 3)

mod_summ_sub <- mod_summ_sub[,c(1,9,8,2,3,4,5,6,7,10)]

library(dplyr)
# mod_summ_sub <- mod_summ_sub %>% arrange(order) 
mod_summ_sub <- mod_summ_sub[order(as.numeric(mod_summ_sub$order)),]

write.csv(mod_summ_sub, file = "Tables and Figures/Table_A1_Parameter_Summ.csv")




