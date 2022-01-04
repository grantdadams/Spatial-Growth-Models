##### READ AND PREPARE DATA #####
source("Data/Data prep/VBGM_Data_prep.R")
prev_coefs <- read.csv("Data/VBGF_Coefs_Previous_Studies.csv")

# Load model
load("VBGF_Stan_models_2018_04_01.RData")

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


library(matrixStats)


################
cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")

tiff(file="Tables and Figures/Fig_4_VBGM_state_sex_specific_plots.tiff" , height= 121.4, pointsize=11,  width=85 , res=300  , units = "mm", family = "serif")
par(mar=c(0.25 , 0.25 , 0 , 0) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:18), 6, 3, byrow = T ), width = c(0.2,1,1), height = c(1,1,1,1,1,.4))
loop.ind <- c(NA,1,2,NA,3,4,NA,5,6,NA,7,8,NA,9,10,NA,NA,NA) 
state.ind <- c(NA,1,2,NA,3,4,NA,6,7,NA,8,9,NA,10,11,NA,NA,NA)  

for ( j in 1:18){
  if(j %in% c(1,4,7,10,13,16,17,18)){plot.new()}
  else{
    for ( i in 1:2){    
      dat.sub <- dat[which(dat$state_no == state.ind[j] & dat$sex == i-1),]
      
      if (i == 1){
        plot(y = dat.sub$FL_mm , x = dat.sub$fractional_age, col = cols.points[i], ylab = "Fork length (mm)", xlab = "Age (yr)", cex = .5, cex.lab = 1.25, pch = 20, ylim = c(min(dat$FL_mm), max(dat$FL_mm)), xlim = c(min(dat$fractional_age), max(dat$fractional_age)), xaxt = "n" , yaxt = "n")
      } else {
        points(y = dat.sub$FL_mm , x = dat.sub$fractional_age, col = cols.points[i], cex = .5, pch = 20)
      } 
      
      # Label plots
      if(i==1){
        legend(x = 18.5, y = 210, c(paste("(", letters, ") ", state.id[which(state.id$state_no == loop.ind[j]),1], sep = ""))[loop.ind[j]], cex = 1, bty ="n")
        
      }
      
      if(i == 1 & state.ind[j] == 1)
        legend(x = 22, y = 350, c("Female","Male"), lwd=3, col = cols, cex = 0.8, bty ="n")
      
    }
    
    
    # Plot lines
    for ( i in 1:2){
      ages <- seq(from = -5, to = 50, by = .2)
      
      # Get parameters
      log_linf_re <- draws[,paste("log_linf_re[",state.ind[j],"]", sep = "")] # Subset the mcmc chains
      betas_linf_re <- draws[,paste("Betas_log_linf[",1:3,"]", sep = "")] # Subset the mcmc chains
      linf.sub <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state.ind[j]),2])) + log_linf_re)
      
      log_k_re <- draws[,paste("log_k_re[",state.ind[j],"]", sep = "")] # Subset the mcmc chains
      betas_k_re <- draws[,paste("Betas_log_k[",1:3,"]", sep = "")] # Subset the mcmc chains
      k.sub <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state.ind[j]),2])) + log_k_re)
      
      log_t0_re <- draws[,paste("t0_re[",state.ind[j],"]", sep = "")] # Subset the mcmc chains
      betas_t0_re <- draws[,paste("Betas_t0[",1:3,"]", sep = "")] # Subset the mcmc chains
      t0.sub <- (as.matrix(betas_t0_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state.ind[j]),2])) + log_t0_re)
      
      lengths.mat <- matrix(NA, ncol = length(ages), nrow = length(t0.sub))
      
      for (m in 1: length(ages)){
        lengths.mat[,m] <- linf.sub * (1 - exp(-k.sub * (ages[m] - t0.sub)))
      }
      
      median.lengths <- apply(lengths.mat, 2, median)
      lines(ages, median.lengths, col = cols[i], lty = 1, lwd = 2)
    }
    
    if (loop.ind[j] %in% c(9,10)){axis(side = 1)
      mtext(side = 1, "Age (yr)", line = 2,cex=0.75)}
    if (loop.ind[j] %in% c(1,3,5,7,9)){axis(side = 2)}
    
    if (loop.ind[j] %in% c(5)){
      mtext(side = 2, "Fork Length (mm)", line = 2, cex = 0.75)}
  }
}


dev.off()
