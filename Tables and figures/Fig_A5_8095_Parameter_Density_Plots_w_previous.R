
# This code creates a to illustrate the 90% and 60% credible intervals of parameter estimates from a hierarchical von Bertalanffy growth function (VBGF) following Midway et al. 2015

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
lat_lon_df <- subset(dat, !duplicated(Lat))
lat_lon_df <- lat_lon_df[, c("state_no", "Lat")]
lat_lon_df$Lat <- (lat_lon_df$Lat - min(lat_lon_df$Lat))

state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS","LA","TX"), state_no = c(1:length(unique(dat$state))))


################
# Midway plot ##
################
cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")

tiff(file="Tables and Figures/Fig_A5_VBGF midway Plots 8095_w_OMEGA_with_previous.tiff" , height= 77.27, pointsize=12,  width=170 , res=300  , units = "mm", family = "serif")
par(mar=c(3.7 , 0.5 , .25 , .125) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:4), 1, 4, byrow = F ))


# Plot Loo density
plot(NA,NA, main = NA, xlab =expression(italic("L"[paste(infinity,"r,s",sep="")])~(mm)), yaxt = "n", ylab= NA, xlim = c(350, 600), ylim = c(-3.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")
for ( j in 1:length(unique(dat$state))){
  for ( i in 1:2){
    state_ind <- c(1:4,6:11)[j] # Indexing variable
    
    # Get parameters
    log_linf_re <- draws[,paste("log_linf_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_linf_re <- draws[,paste("B_log_linf[",1:3,"]", sep = "")] # Subset the mcmc chains
    
    dat.sub <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_linf_re)
    quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
    quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
    
    spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
    
    # Plot the credible intervals
    segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
    segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
    
    # Plot the median estimate
    med <- median(dat.sub)
    segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
    
    # Plot the other state params
    if(length(which(prev_coefs$state_no == j)) > 0){
      prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no == j),]
      if(length(which(prev_coefs_sub$Sex == i)) > 0){
        prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$Sex == i),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$Linf[m]
          points(med,spacing.ind, pch = 21, cex = 1.18, bg = cols.points[i])
        }
      }
      if(length(which(is.na(prev_coefs_sub$Sex)==T)) > 0 & i == 2){
        prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$Sex)),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$Linf[m]
          points(med,spacing.ind + 0.175, pch = 21, cex = 1.18, bg = "grey")
        }
      }
    }
    
    # Label lines
    if(i==1){
      text(mean(dat.sub), spacing.ind+ 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
    }
    
  }
}
legend(x = 345, y = 6.3, c(sex,"Combined","Previous"), pch = c(NA,NA,NA,21), lty = c(1,1,1,NA), lwd=c(3,3,3,NA), col = c(cols,"grey", 1), cex = 0.8, bty ="n")
legend("topleft", c("(a)"), cex = 1, bty ="n", adj = 1)

#plot K density
plot(NA,NA, main = NA, xlab =expression(italic("k"["s,r"])~(yr^-1)), yaxt = "n", ylab= NA, xlim = c(.13, 0.6), ylim = c(-3.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")
for ( j in 1:length(unique(dat$state))){
  for ( i in 1:2){
    state_ind <- c(1:4,6:11)[j] # Indexing variable
    
    # Get parameters
    log_linf_re <- draws[,paste("log_k_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_linf_re <- draws[,paste("B_log_k[",1:3,"]", sep = "")] # Subset the mcmc chains
    
    dat.sub <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_linf_re)
    quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
    quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
    
    spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
    
    # Plot the credible intervals
    segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
    segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
    
    # Plot the median estimate
    med <- median(dat.sub)
    segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
    
    # Plot the other state params
    if(length(which(prev_coefs$state_no == j)) > 0){
      prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no == j),]
      if(length(which(prev_coefs_sub$Sex == i)) > 0){
        prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$Sex == i),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$K[m]
          points(med,spacing.ind, pch = 21, cex = 1.18, bg = cols.points[i])
        }
      }
      if(length(which(is.na(prev_coefs_sub$Sex)==T)) > 0 & i == 2){
        prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$Sex)),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$K[m]
          points(med,spacing.ind + 0.175, pch = 21, cex = 1.18, bg = "grey")
        }
      }
    }
    
    # Label lines
    if(i==1){
      text(mean(dat.sub), spacing.ind+ 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
    }
    
  }
}
#legend("topright", sex, lwd=3, col = cols, cex = 0.75, bty ="n")
legend("topleft", c("(b)"), cex = 1, bty ="n", adj = 1)

#plot t0 density
plot(NA,NA, main = NA, xlab =expression(italic(t[paste(0,"s,r",sep="")])~(yr)), yaxt = "n", ylab= NA, xlim = c(-3,0), ylim = c(-3.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")
for ( j in 1:length(unique(dat$state))){
  for ( i in 1:2){
    state_ind <- c(1:4,6:11)[j] # Indexing variable
    
    # Get parameters
    log_linf_re <- draws[,paste("t0_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_linf_re <- draws[,paste("B_t0[",1:3,"]", sep = "")] # Subset the mcmc chains
    
    dat.sub <- (as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_linf_re)
    quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
    quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
    
    spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
    
    # Plot the credible intervals
    segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
    segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
    
    # Plot the median estimate
    med <- median(dat.sub)
    segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
    
    # Plot the other state params
    if(length(which(prev_coefs$state_no == j)) > 0){
      prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no == j),]
      if(length(which(prev_coefs_sub$Sex == i)) > 0){
        prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$Sex == i),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$t0[m]
          points(med,spacing.ind, pch = 21, cex = 1.18, bg = cols.points[i])
        }
      }
      if(length(which(is.na(prev_coefs_sub$Sex)==T)) > 0 & i == 2){
        prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$Sex)),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$t0[m]
          points(med,spacing.ind + 0.175, pch = 21, cex = 1.18, bg = "grey")
        }
      }
    }
    
    # Label lines
    if(i==1){
      text(mean(dat.sub), spacing.ind+ 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
    }
    
  }
}
#legend("topleft", sex, lwd=3, col = cols, cex = 0.75, bty ="n")
legend("topleft", c("(c)"), cex = 1, bty ="n", adj = 1)

#plot Gallucci OMEGA density
plot(NA,NA, main = NA, xlab =expression(omega["s,r"]~(mm~yr^-1)), yaxt = "n", ylab= NA, xlim = c(75,250), ylim = c(-3.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")
for ( j in 1:length(unique(dat$state))){
  for ( i in 1:2){
    
    state_ind <- c(1:4,6:11)[j] # Indexing variable
    
    # Get linf
    log_linf_re <- draws[,paste("log_linf_re[",state_ind,"]", sep = "")]
    betas_linf_re <- draws[,paste("B_log_linf[",1:3,"]", sep = "")]
    
    Linf <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_linf_re)

    # Get k
    log_k_re <- draws[,paste("log_k_re[",state_ind,"]", sep = "")]
    betas_k_re <- draws[,paste("B_log_k[",1:3,"]", sep = "")]
    
    k <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_k_re)
   
    
    # Calculate w
    dat.sub <- Linf * k
    
    quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
    quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
    
    spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
    
    # Plot the credible intervals
    segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
    segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
    
    # Plot the median estimate
    med <- median(dat.sub)
    segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
    
    # Plot the other state params
    if(length(which(prev_coefs$state_no == j)) > 0){
      prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no == j),]
      if(length(which(prev_coefs_sub$Sex == i)) > 0){
        prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$Sex == i),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$Linf[m] * prev_coefs_sub_sex$K[m]
          points(med,spacing.ind, pch = 21, cex = 1.18, bg = cols.points[i])
        }
      }
      if(length(which(is.na(prev_coefs_sub$Sex)==T)) > 0 & i == 2){
        prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$Sex)),]
        for(m in 1:nrow(prev_coefs_sub_sex)){
          med <- prev_coefs_sub_sex$Linf[m] * prev_coefs_sub_sex$K[m]
          points(med,spacing.ind + 0.175, pch = 21, cex = 1.18, bg = "grey")
        }
      }
    }
    
    # Label lines
    if(i==1){
      text(mean(dat.sub), spacing.ind+ 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
    }
    
  }
}
#legend("topright", sex, lwd=3, col = cols, cex = 0.75, bty ="n")
legend("topleft", c("(d)"), cex = 1, bty ="n", adj = 1)
dev.off()


