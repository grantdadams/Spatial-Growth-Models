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



tiff(file="Tables and Figures/Fig_A6_Predicted_Length_at_Age__with_previous_estimates.tiff" , height= 77.27, pointsize=12,  width=170 , res=300  , units = "mm", family = "serif")
par(mar=c(3.7 , 0.5 , .25 , .125) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:4), 1, 4, byrow = F ))


# Plot length-at-age
age.all <- c(1,5,10,15)

# Get vectors of all VBGF Parameters
linf.all <- c()
k.all <- c()
t0.all <- c()

for(j in 1:10){
  for ( i in 1:2){
    state_ind <- c(1:4,6:11)[j] # Indexing variable
    
    log_linf_re <- draws[,paste("log_linf_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_linf_re <- draws[,paste("Betas_log_linf[",1:3,"]", sep = "")] # Subset the mcmc chains
    linf.sub <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_linf_re)
    linf.all <- c(linf.all, linf.sub)
    
    log_k_re <- draws[,paste("log_k_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_k_re <- draws[,paste("Betas_log_k[",1:3,"]", sep = "")] # Subset the mcmc chains
    k.sub <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_k_re)
    k.all <- c(k.all, k.sub)
    
    log_t0_re <- draws[,paste("t0_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_t0_re <- draws[,paste("Betas_t0[",1:3,"]", sep = "")] # Subset the mcmc chains
    t0.sub <- (as.matrix(betas_t0_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_t0_re)
    t0.all <- c(t0.all, t0.sub)
  }
}


for(m in age.all){
  
  # Calculate max and min length
  pred.length.all <- linf.all * (1 - exp (-k.all * (m - t0.all) )) # Predict length-at-age m
  pred_prev <- prev_coefs$Linf * (1 - exp (- prev_coefs$K * (m - prev_coefs$t0)))
  quantiles <- quantile(pred.length.all,c(0.005, 0.995))
  xlim.low <- min(quantiles, pred_prev)
  xlim.high <- max(quantiles, pred_prev)
  
  plot(NA,NA, main = NA, xlab =paste("Length-at-age-", m," (mm)", sep = ""), yaxt = "n", ylab= NA, xlim = c(xlim.low, xlim.high), ylim = c(-3.5,6.5), cex.lab = 1.25)
  sex <- c("Female","Male")
  for ( j in 1:length(unique(dat$state))){
    for ( i in 1:2){
      
      state_ind <- c(1:4,6:11)[j] # Indexing variable
      
      # Get parameters
      log_linf_re <- draws[,paste("log_linf_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
      betas_linf_re <- draws[,paste("Betas_log_linf[",1:3,"]", sep = "")] # Subset the mcmc chains
      linf.sub <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_linf_re)
      
      log_k_re <- draws[,paste("log_k_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
      betas_k_re <- draws[,paste("Betas_log_k[",1:3,"]", sep = "")] # Subset the mcmc chains
      k.sub <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_k_re)
      
      log_t0_re <- draws[,paste("t0_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
      betas_t0_re <- draws[,paste("Betas_t0[",1:3,"]", sep = "")] # Subset the mcmc chains
      t0.sub <- (as.matrix(betas_t0_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_t0_re)
      
      pred.length <- linf.sub * (1 - exp (-k.sub * (m - t0.sub) ))
      quantiles_025_975 <- quantile(pred.length,c(0.025, 0.975))
      quantiles_10_90 <- quantile(pred.length,c(0.1, 0.9))
      
      spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
      
      # Plot the credible intervals
      segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[i])
      segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[i])
      
      # Plot the median estimate
      med <- median(pred.length)
      segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
      
      
      # Label lines
      if(i==1){
        text(mean(pred.length), spacing.ind + 0.3, paste(state.id[j,1],sep=" ") , cex = 0.75) 
      }
      
      # Plot the other state params
      if(length(which(prev_coefs$state_no == j)) > 0){
        prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no == j),]
        if(length(which(prev_coefs_sub$Sex == i)) > 0){
          prev_coefs_sub_sex <- prev_coefs_sub[which(prev_coefs_sub$Sex == i),]
          for(k in 1:nrow(prev_coefs_sub_sex)){
            pred.length <- prev_coefs_sub_sex$Linf[k] * (1 - exp (-prev_coefs_sub_sex$K[k] * (m - prev_coefs_sub_sex$t0[k]) ))
            points(pred.length,spacing.ind, pch = 21, cex = 1.18, bg = cols.points[i])
          }
        }
        if(length(which(is.na(prev_coefs_sub$Sex)==T)) > 0 & i == 2){
          prev_coefs_sub_sex <- prev_coefs_sub[which(is.na(prev_coefs_sub$Sex)),]
          for(k in 1:nrow(prev_coefs_sub_sex)){
            pred.length <- prev_coefs_sub_sex$Linf[k] * (1 - exp (-prev_coefs_sub_sex$K[k] * (m - prev_coefs_sub_sex$t0[k]) ))
            points(pred.length,spacing.ind + 0.175, pch = 21, cex = 1.18, bg = "grey")
          }
        }
      }
      
      
    }
  }
  legend("topleft", c("(a)", "(b)", "(c)", "(d)")[which(m == age.all)], cex = 1, bty ="n", adj = 1)
  
  if(m == 1){
    legend(y = 6.3, x = 140, c(sex,"Combined","Previous"), pch = c(NA,NA,NA,21), lty = c(1,1,1,NA), lwd=c(3,3,3,NA), col = c(cols,"grey", 1), cex = 0.8, bty ="n")
  }else{
    #legend("topleft", sex, lwd=3, col = cols, cex = 0.8, bty ="n")
  }
}

dev.off()
