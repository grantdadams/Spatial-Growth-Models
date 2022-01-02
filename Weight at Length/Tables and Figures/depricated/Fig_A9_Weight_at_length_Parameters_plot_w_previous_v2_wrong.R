


library(rstan)

##### RE AD AND PREPARE DATA #####
source("Weight_at_length_data_prep.R")
prev_coefs <- read.csv("WAL_Prev_Param_Estimates.csv")

# Load data
load("StanFitLatRElog_short_2018_03_31.RData")

# Get parameters
model <- StanFitLatRElog
mod_summ <- data.frame(summary(model)[[1]])
mod_summ$parameter <- rownames(mod_summ)

######################################################
# Get a, K, and t0
draws <- as.data.frame(model)
# Get average lat
lat_lon_df <- subset(dat, !duplicated(Lat))
lat_lon_df <- lat_lon_df[, c("state_no_check", "Lat")]
lat_lon_df$Lat <- (lat_lon_df$Lat - min(lat_lon_df$Lat))

state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS"), state_no_check = c(1:length(unique(dat$state))), state_no = c(1:4,6:9))

prev_coefs <- merge(prev_coefs[,c(1,3,4,5,6,7)], state.id, by = "state_no")
prev_coefs$sex <- prev_coefs$sex - 1

parameters_to_save <- c(paste("Betas_log_a[",1:3,"]", sep = ""), paste("Betas_log_b[",1:3,"]", sep = ""), "sigma_a", "sigma_b", "sigma", paste("log_a_re[",c(1:8), "]", sep = ""), paste("log_b_re[",c(1:8), "]", sep = ""))
draws <- draws[,which(colnames(draws) %in% parameters_to_save)]


################
# Midway plot ##
################
cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")

tiff(file="Tables and Figures/Fig_A9_Weight_at_length_parameters_with_previous.tiff" , height= 77.27, pointsize=12,  width=170 , res=300  , units = "mm", family = "serif")
par(mar=c(3.7 , 0.5 , .25 , .125) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0.5, 0 , 0.5), cex = 0.75)
layout(matrix(c(1:4), 1, 4, byrow = F ))

# Put mcmc of parameters into a list to index later
param.names <- c(expression(italic("a"["r,s"])), expression(italic("b"["r,s"])))

plot(NA,NA, main = NA, xlab =param.names[1], yaxt = "n", ylab= NA, xlim = c(min(c(0.000015, prev_coefs$a)), max(c(0.000075, prev_coefs$a))), ylim = c(-1.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")
m = 1
for ( j in 1:length(unique(dat$state_no_check))){
  for ( i in 1:2){
    state_ind <- c(1:8)[j] # Indexing variable
    
    # Get parameters
    log_a_re <- draws[,paste("log_a_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_a_re <- draws[,paste("Betas_log_a[",1:3,"]", sep = "")] # Subset the mcmc chains
    
    dat.sub <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_a_re)
    
    
    
    quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
    quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
    
    spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
    
    # Plot the credible intervals
    segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
    segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
    
    # Plot the median estimate
    med <- median(dat.sub)
    segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
    
    
    # Label lines
    if(i==1){
      text(mean(dat.sub), spacing.ind + 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
    }
    
    # Plot the other previous studies params params
    if(length(which(prev_coefs$state_no_check == j)) > 0){
      prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no_check == j),]
      if(length(which(prev_coefs_sub$sex == i)) > 0){
        prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$sex == i),]
        for(k in 1:nrow(prev_coefs_sub_sex)){
          param.est <- ifelse(m == 1, prev_coefs_sub_sex$a[k], prev_coefs_sub_sex$b[k])
          points(param.est,spacing.ind, pch = 21, bg = cols.points[i])
        }
      }
      # Non sex
      if(length(which(is.na(prev_coefs_sub$sex)==T)) > 0 & i == 2){
        prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$sex)),]
        for(k in 1:nrow(prev_coefs_sub_sex)){
          param.est <- ifelse(m == 1, prev_coefs_sub_sex$a[k], prev_coefs_sub_sex$b[k])
          points(param.est,spacing.ind + 0.175, pch = 21, bg = "grey")
        }
      }
    }
    
  }
}
legend("topleft", c(sex), pch = c(NA,NA), lty = c(1,1), lwd=c(3,3), col = c(cols), cex = 0.8, bty ="n")




plot(NA,NA, main = NA, xlab =param.names[2], yaxt = "n", ylab= NA, xlim = c(min(c(2.85, prev_coefs$b)), max(c(3.1, prev_coefs$b))), ylim = c(-1.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")
m = 2
for ( j in 1:length(unique(dat$state_no_check))){
  for ( i in 1:2){
    state_ind <- c(1:8)[j] # Indexing variable
    
    # Get parameters
    log_b_re <- draws[,paste("log_b_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_b_re <- draws[,paste("Betas_log_b[",1:3,"]", sep = "")] # Subset the mcmc chains
    
    dat.sub <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_b_re)
    
    quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
    quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
    
    spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
    
    # Plot the credible intervals
    segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
    segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
    
    # Plot the median estimate
    med <- median(dat.sub)
    segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
    
    
    # Label lines
    if(i==1){
      text(mean(dat.sub), spacing.ind + 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
    }
    
    # Plot the other previous studies params params
    if(length(which(prev_coefs$state_no_check == j)) > 0){
      prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no_check == j),]
      if(length(which(prev_coefs_sub$sex == i)) > 0){
        prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$sex == i),]
        for(k in 1:nrow(prev_coefs_sub_sex)){
          param.est <- ifelse(m == 1, prev_coefs_sub_sex$a[k], prev_coefs_sub_sex$b[k])
          points(param.est,spacing.ind, pch = 21, bg = cols.points[i])
        }
      }
      # Non sex
      if(length(which(is.na(prev_coefs_sub$sex)==T)) > 0 & i == 2){
        prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$sex)),]
        for(k in 1:nrow(prev_coefs_sub_sex)){
          param.est <- ifelse(m == 1, prev_coefs_sub_sex$a[k], prev_coefs_sub_sex$b[k])
          points(param.est,spacing.ind + 0.175, pch = 21, bg = "grey")
        }
      }
    }
    
  }
}
dev.off()
