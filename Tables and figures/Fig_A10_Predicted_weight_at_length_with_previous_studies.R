# clear the workspace
rm(list = ls(all = T))

library(rstan)

##### RE AD AND PREPARE DATA #####
source("Data/Data prep/WAL_Data_prep.R")

# Load data
load("Stan models/Models/WAL_Stan_models_2018_02_19.RData")

# Get parameters
model <- mod_list[[2]]
mod_summ <- data.frame(summary(model)[[1]])
mod_summ$parameter <- rownames(mod_summ)

######################################################
# reorder state variables
state_no_df <- subset(dat, !duplicated(state_no_check))
state_no_df <- state_no_df[,c("state_no", "state_no_check", "state.x")]
state_no_df <- state_no_df[order(state_no_df$state_no),]
state_no_df$state_no_2 <- 1:8

# Get a, K, and t0
draws <- as.data.frame(model)
# Get average lat
lat_lon_df <- subset(dat, !duplicated(Lat))
lat_lon_df <- lat_lon_df[, c("state_no_check", "Lat")]
lat_lon_df$Lat <- (lat_lon_df$Lat - min(lat_lon_df$Lat))

state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS"), state_no_check = state_no_df$state_no_check, state_no = state_no_df$state_no)

parameters_to_save <- c(paste("B_log_a[",1:3,"]", sep = ""), paste("B_log_b[",1:3,"]", sep = ""), "sigma_a", "sigma_b", "sigma", paste("log_a_re[",c(1:8), "]", sep = ""), paste("log_b_re[",c(1:8), "]", sep = ""))
draws <- draws[,which(colnames(draws) %in% parameters_to_save)]
# This code creates a to illustrate the 95% and 80% credible intervals of parameter estimates from a hierarchical weight-at-length power equation following Midway et al. 2015
library(matrixStats)


################
# Midway plot ##
################
cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")

tiff(file="Tables and Figures/Fig_A10_Predicted_weight_at_length_w_previous_estimates.tiff" , height= 77.27, pointsize=12,  width=170 , res=300  , units = "mm", family = "serif")
par(mar=c(3.7 , 0.5 , .25 , .125) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:4), 1, 4, byrow = F ))


# Plot length-at-age
length.all <- c(200,300,400,500)

# Get vectors of all VBGF Parameters
library(reshape)
a.all <- c()
b.all <- c()
t0.all <- c()

for(j in 1:8){
  for ( i in 1:2){
    state_ind <- c(state_no_df$state_no_check)[j] # Indexing variable
    
    log_a_re <- draws[,paste("log_a_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_a_re <- draws[,paste("B_log_a[",1:3,"]", sep = "")] # Subset the mcmc chains
    a.sub <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, i - 1, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_a_re)
    a.all <- c(a.all, a.sub)
    
    log_b_re <- draws[,paste("log_b_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_b_re <- draws[,paste("B_log_b[",1:3,"]", sep = "")] # Subset the mcmc chains
    b.sub <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, i - 1, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_b_re)
    b.all <- c(b.all, b.sub)
  }
}


for(m in length.all){
  
  # Calculate max and min length
  pred.weight.all <- a.all * m ^ b.all # Predict weight-at-length
  quantiles <- quantile(pred.weight.all,c(0.005, 0.995), na.rm = T)/1000
  xlim.low <- min(quantiles)
  xlim.high <- max(quantiles)
  
  plot(NA,NA, main = NA, xlab =paste("Weight (kg) at ", m," mm FL", sep = ""), yaxt = "n", ylab= NA, xlim = c(xlim.low, xlim.high), ylim = c(-1.5,6.5), cex.lab = 1.25)
  sex <- c("Female","Male")
  for ( j in 1:length(unique(dat$state_no))){
    for ( i in 1:2){
      state_ind <- c(state_no_df$state_no_check)[j] # Indexing variable
      
      # Get parameters
      log_a_re <- draws[,paste("log_a_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
      betas_a_re <- draws[,paste("B_log_a[",1:3,"]", sep = "")] # Subset the mcmc chains
      a.sub <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, i - 1, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_a_re)
      
      log_b_re <- draws[,paste("log_b_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
      betas_b_re <- draws[,paste("B_log_b[",1:3,"]", sep = "")] # Subset the mcmc chains
      b.sub <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, i - 1, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_b_re)
      
      pred.weight <- a.sub * m ^ b.sub 
      quantiles_025_975 <- quantile(pred.weight,c(0.025, 0.975), na.rm = T)/1000
      quantiles_10_90 <- quantile(pred.weight,c(0.1, 0.9), na.rm = T)/1000
      
      
      spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
      
      # Plot the credible intervals
      segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
      segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
      
      # Plot the median estimate
      med <- median(pred.weight)/1000
      segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
      
      
      # Label lines
      if(i==1){
        text(mean(pred.weight)/1000, spacing.ind + 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
      }
      
      # Plot the other previous studies params params
      if(length(which(prev_coefs$state_no_check == j)) > 0){
        prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no == j),]
        if(length(which(prev_coefs_sub$sex == i)) > 0){
          prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$sex == i),]
          for(k in 1:nrow(prev_coefs_sub_sex)){
            param.est <-0.0000547 *  (m ^ 2.997)/1000
            points(param.est,spacing.ind, pch = 21, bg = cols.points[i])
          }
        }
        # Non sex
        if(length(which(is.na(prev_coefs_sub$sex)==T)) > 0 & i == 2){
          prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$sex)),]
          for(k in 1:nrow(prev_coefs_sub_sex)){
            param.est <-(prev_coefs_sub_sex$a[k] *  m ^ prev_coefs_sub_sex$b[k])/1000
            points(param.est,spacing.ind + 0.175, pch = 21, bg = "grey")
          }
        }
      }
      
    }
  }
  legend("topleft", c(paste("(", letters, ")", sep = ""))[which(length.all == m)], cex = 1, bty ="n", adj = 1)
  if(m == length.all[1]){
    legend(y = 6.3, x = 0.193, c(sex,"Combined","Previous"), pch = c(NA,NA,NA,21), lty = c(1,1,1,NA), lwd=c(3,3,3,NA), col = c(cols,"grey", 1), cex = 0.8, bty ="n")
  }
}

dev.off()
