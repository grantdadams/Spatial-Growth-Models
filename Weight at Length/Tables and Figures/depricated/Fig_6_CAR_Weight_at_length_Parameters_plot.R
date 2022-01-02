library(rstan)

##### RE AD AND PREPARE DATA #####
source("Weight_at_length_data_prep.R")

# Load data
load("StanFitLatRE_2018_02_10.RData")

# Get parameters
model <- StanFitLatCAR
mod_summ <- data.frame(summary(model)[[1]])
mod_summ$parameter <- rownames(mod_summ)

######################################################
# Get a, K, and t0
draws <- as.data.frame(model)
# Get average lat
lat_lon_df <- subset(dat, !duplicated(Lat))
lat_lon_df <- lat_lon_df[, c("state_no", "Lat")]
lat_lon_df$Lat <- (lat_lon_df$Lat - min(lat_lon_df$Lat))

state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS"), state_no = c(1:length(unique(dat$state))))

parameters_to_save <- c(paste("Betas_log_a[",1:3,"]", sep = ""), paste("Betas_log_b[",1:3,"]", sep = ""), "sigma_a", "sigma_b", "sigma", paste("log_a_re[",c(1:4,6:9), "]", sep = ""), paste("log_b_re[",c(1:4,6:9), "]", sep = ""))
draws <- draws[,which(colnames(draws) %in% parameters_to_save)]


################
# Midway plot ##
################
cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")

tiff(file="Tables and Figures/Fig_6_CAR_Weight_at_length_parameters.tiff" , height= 77.27, pointsize=12,  width=170 , res=300  , units = "mm", family = "serif")
par(mar=c(3.7 , 0.5 , .25 , .125) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0.5, 0 , 0.5), cex = 0.75)
layout(matrix(c(1:4), 1, 4, byrow = F ))

# Put mcmc of parameters into a list to index later
param.names <- c(expression(italic("a"["r,s"])), expression(italic("b"["r,s"])))

plot(NA,NA, main = NA, xlab =param.names[1], yaxt = "n", ylab= NA, xlim = c(0.000015, 0.000115), ylim = c(-1.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")

for ( j in 1:length(unique(dat$state_no))){
  for ( i in 1:2){
    state_ind <- c(1:4,6:11)[j] # Indexing variable
    
    # Get parameters
    log_a_re <- draws[,paste("log_a_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_a_re <- draws[,paste("Betas_log_a[",1:3,"]", sep = "")] # Subset the mcmc chains
    
    dat.sub <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_a_re)
    
    
    
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
    
  }
}
legend("topleft", c(sex), pch = c(NA,NA), lty = c(1,1), lwd=c(3,3), col = c(cols), cex = 0.8, bty ="n")




plot(NA,NA, main = NA, xlab =param.names[2], yaxt = "n", ylab= NA, xlim = c(2.75, 3.1), ylim = c(-1.5,6.5), cex.lab = 1.25)
sex <- c("Female","Male")

for ( j in 1:length(unique(dat$state_no))){
  for ( i in 1:2){
    state_ind <- c(1:4,6:9)[j] # Indexing variable
    
    # Get parameters
    log_b_re <- draws[,paste("log_b_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_b_re <- draws[,paste("Betas_log_b[",1:3,"]", sep = "")] # Subset the mcmc chains
    
    dat.sub <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_b_re)
    
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
    
  }
}
dev.off()
