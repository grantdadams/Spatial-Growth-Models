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

state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS"), state_no_check = state_no_df$state_no_check, state_no = state_no_df$state_no)

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

tiff(file="Tables and Figures/Fig_8_State_sex_specific_plots.tiff" , height= 103, pointsize=11,  width=85 , res=300  , units = "mm", family = "serif")
par(mar=c(0.25 , 0.25 , 0 , 0) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:15), 5, 3, byrow = T ), width = c(0.2,1,1), height = c(1,1,1,1,.4))
loop.ind <- c(NA,1,2,NA,3,4,NA,5,6,NA,7,8,NA,NA,NA) 
state.ind <- c(NA,1,2,NA,3,4,NA,5,6,NA,7,8,NA,NA,NA) 

# Put mcmc of parameters into a list to index later


for ( j in 1:15){
  if(is.na(loop.ind[j]) == T){plot.new()}
  else{
    for ( i in 1:2){
      # Subject dat
      dat.sub <- dat[which(dat$state_no_check == state.ind[j] & dat$sex == i - 1),]

      if (i == 1){
        plot(y = dat.sub$wet_weight_grams/1000 , x = dat.sub$FL_mm, col = cols.points[i], ylab = "Weight (kg)", xlab = "Fork length (mm)", cex = .5, cex.lab = 1.25, pch = 20, ylim = c(min(dat$wet_weight_grams/1000), max(dat$wet_weight_grams/1000)), xlim = c(min(dat$FL_mm), max(dat$FL_mm)), xaxt = "n" , yaxt = "n")
      } else {
        points(y = dat.sub$wet_weight_grams/1000 , x = dat.sub$FL_mm, col = cols.points[i], cex = .5, pch = 20)
      } 
      
      # Label plots
      if(i==1){
        legend(y = 7.2, x = 60, c(paste("(", letters, ") ", state.id[which(state.id$state_no_check == loop.ind[j]),1], sep = ""))[loop.ind[j]], cex = 1, bty ="n")
      }
    }
    
    
    # Plot lines
    for (i in 1:2){
      lengths <- seq(from = 0, to = 600, by = 1)
      
      ind <- loop.ind[j] + length(unique(dat$state_no_check)) * (i-1)# Indexing variable
      # Get parameters
      log_a_re <- draws[,paste("log_a_re[",state.ind[j],"]", sep = "")] # Subset the mcmc chains
      betas_a_re <- draws[,paste("B_log_a[",1:3,"]", sep = "")] # Subset the mcmc chains
      a.param.sub <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no_check == state.ind[j]),2])) + log_a_re)
      
      log_b_re <- draws[,paste("log_b_re[",state.ind[j],"]", sep = "")] # Subset the mcmc chains
      betas_b_re <- draws[,paste("B_log_b[",1:3,"]", sep = "")] # Subset the mcmc chains
      b.param.sub <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no_check == state.ind[j]),2])) + log_b_re)
      
      weights.mat <- matrix(NA, ncol = length(lengths), nrow = length(a.param.sub))
      
      for (m in 1: length(lengths)){
        weights.mat[,m] <- (a.param.sub *  (lengths[m] ^ b.param.sub))/1000
      }
      
      median.lengths <- apply(weights.mat, 2, median)
      lines(lengths, median.lengths, col = cols[i], lty = 1, lwd = 2)
    }
    if(i==1){
      
      ###################################### FLAG FOR LATER
      legend("topleft", paste(state.id[which(state.id$state_no_check == loop.ind[j]),1],sep=" ") ,  bty = "n")
    }
    
  }
  if(loop.ind[j] %in% 1){
    legend("bottomright", c("Female","Male"), lwd=3, col = cols, cex = 0.8, bty ="n")
  }
  
  if (loop.ind[j] %in% c(7,8)){axis(side = 1);mtext(side = 1, "Fork length (mm)", line = 2, cex = .75)}
  if (loop.ind[j] %in% c(1,3,5,7)){axis(side = 2)}
  if (loop.ind[j] %in% c(5)){axis(side = 2);mtext(side = 2, "Weight (kg)", line = 2, cex = 0.75, adj = 2.05)}
}


dev.off()
