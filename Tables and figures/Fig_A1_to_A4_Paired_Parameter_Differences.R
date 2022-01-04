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

cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")


# Get vectors of all VBGF Parameters
linf.all <- list(matrix(NA, nrow = nrow(draws), ncol = 10), matrix(NA, nrow = nrow(draws), ncol = 10))
k.all <- list(matrix(NA, nrow = nrow(draws), ncol = 10), matrix(NA, nrow = nrow(draws), ncol = 10))
t0.all <- list(matrix(NA, nrow = nrow(draws), ncol = 10), matrix(NA, nrow = nrow(draws), ncol = 10))
w0.all <- list(matrix(NA, nrow = nrow(draws), ncol = 10), matrix(NA, nrow = nrow(draws), ncol = 10))

for(j in 1:10){
  for ( i in 1:2){
    state_ind <- c(1:4,6:11)[j] # Indexing variable
    
    log_linf_re <- draws[,paste("log_linf_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_linf_re <- draws[,paste("B_log_linf[",1:3,"]", sep = "")] # Subset the mcmc chains
    linf.sub <- exp(as.matrix(betas_linf_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_linf_re)
    linf.all[[i]][,j] <- linf.sub
    
    log_k_re <- draws[,paste("log_k_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_k_re <- draws[,paste("B_log_k[",1:3,"]", sep = "")] # Subset the mcmc chains
    k.sub <- exp(as.matrix(betas_k_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_k_re)
    k.all[[i]][,j] <- k.sub
    
    log_t0_re <- draws[,paste("t0_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_t0_re <- draws[,paste("B_t0[",1:3,"]", sep = "")] # Subset the mcmc chains
    t0.sub <- (as.matrix(betas_t0_re) %*% as.vector(c(1, i-1, lat_lon_df[which(lat_lon_df$state_no == state_ind),2])) + log_t0_re)
    t0.all[[i]][,j] <- t0.sub
    
    w0.all[[i]][,j] <- k.sub * linf.sub
  }
}

params_list <- list(linf.all, k.all, t0.all, w0.all)
sex <- c("Female","Male")

########################################################
figure_names <- c("Tables and Figures/Fig_A1_Paired_Linf_Parameter_Differences.tiff", 
                  "Tables and Figures/Fig_A2_Paired_K_Parameter_Differences.tiff", 
                  "Tables and Figures/Fig_A3_Paired_t0_Parameter_Differences.tiff", 
                  "Tables and Figures/Fig_A4_Paired_w0_Parameter_Differences.tiff")

parameter_names <- list(
  expression(paste("Paired differences in ",italic("L"[paste(infinity,"r,s",sep="")])~(mm))),
  expression(paste("Paired differences in ",italic("k"["s,r"])~(yr^-1))),
  expression(paste("Paired differences in ",italic(t[paste(0,"s,r",sep="")])~(yr))),
  expression(paste("Paired differences in ",omega["s,r"]~(mm~yr^-1))))

# Set limits and parameter
xlim <- list(c(-70, 150),
             c(-.38, .1),
             c(-1.8, 1.8),
             c(-158, 71))

xaxis_tick <- list(seq(-50, 150, by = 50), 
                   seq(-.35, 0.1, by = .1),
                   seq(-1, 1, by = 1),
                   seq(-150, 50, by = 50))

xaxis_lab <- list(c(NA, 0, NA, 100, NA),
                  c(NA, -0.2, NA, 0, NA),
                  c(NA,0,1),
                  c(NA, -100, NA, 0, NA))


#### Plot the thing ####
for(k in 1:4){
  
  tiff(file = figure_names[k] , height= 140, pointsize=12,  width=170 , res=300  , units = "mm", family = "serif")
  par(mar=c(0 , 0 , 0 , 0), tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
  layout(matrix(c(1:121), 11, 11, byrow = F ), widths = c(0.1, rep(1, 10)))
  
  # Add left border
  for(j in 1:11){
    plot.new()
  }
  
  # Data frame of unique values
  combo_df <- expand.grid(c(1:10), c(1:10))
  combo_df <- combo_df[,c(2,1)]
  
  for( j in 1:nrow(combo_df)){
    
    
    # Dont plot if same parameter
    if(combo_df[j,1] == combo_df[j,2]) {plot.new()
      # Add text
      text(x = 0.5, y = .5, labels = as.character(state.id[combo_df[j,1], 1]), adj = c(0.5, 0.5), srt = 0, cex = 1.25)  }
    
    # Plot the paired differences
    else{
      
      if(combo_df[j,1] > combo_df[j,2]) {
        plot.new()
      } else{
        
        # Calculate paired differences
        diff_fem <- params_list[[k]][[1]][,combo_df[j,1]] -  params_list[[k]][[1]][,combo_df[j,2]]
        diff_male <- params_list[[k]][[2]][,combo_df[j,1]] -  params_list[[k]][[2]][,combo_df[j,2]]
        
        # Put into a list
        diff_list <- list(diff_fem, diff_male)
        
        plot(NA,NA, main = NA, xlab =NA, yaxt = "n",xaxt = "n", ylab= NA, xlim = c(xlim[[k]][1], xlim[[k]][2]), ylim = c(0,1), cex.lab = 1.25)
        
        for ( i in 1:2){
          
          quantiles_025_975 <- quantile(diff_list[[i]],c(0.025, 0.975))
          quantiles_10_90 <- quantile(diff_list[[i]],c(0.1, 0.9))
          
          # Color background if 95% CI cross 0
          if(prod(sign(sign(quantiles_025_975))) < 0){
            if(i == 2){
              rect(par("usr")[1], par("usr")[3], par("usr")[2], 0.5, col = "#E7E7E7", border = NA)
            }
            if(i == 1){
              rect(par("usr")[1], 0.5, par("usr")[2], par("usr")[4], col = "#E7E7E7", border = NA)
            }
          }
          
          # Color background if 0% CI cross 0
          if(prod(sign(sign(quantiles_10_90))) < 0){
            if(i == 2){
              rect(par("usr")[1], par("usr")[3], par("usr")[2], 0.5, col = "#AAAAAA", border = NA)
            }
            if(i == 1){
              rect(par("usr")[1], 0.5, par("usr")[2], par("usr")[4], col = "#AAAAAA", border = NA)
            }
          }
          
          spacing.ind <- 1 - (0.25 + (i-1) * 0.5) # Spacing index for plots
          
          # Plot the credible intervals
          segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 6 , col = cols.points[i], lend = 1)
          segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 9 , col = cols[i], lend = 1)
          
          # Plot the median estimate
          med <- median(diff_list[[i]])
          segments(med,spacing.ind-.07,med,spacing.ind+.07, lwd = 3 , col = 1, lend = 2)
          
          # Add zero line
          abline(v = 0, col = 1, lty = 2)
          box()
        }
        
        # Add x-axis to the bottom
        if(j %in% seq(from = 10, to = 90, by = 10)){
          #axis(side = 1, cex.axis = 1.2)
          axis(side = 1, at = xaxis_tick[[k]] , labels = xaxis_lab[[k]], cex.axis = 1.2)
        }
        
        # Add x-axis label
        if(j == 50){
          mtext(side = 1, parameter_names[[k]], cex = 1, line = 2.5)
        }
        
        # Add space at the bottom
        if(j %in% seq(from = 10, to = 100, by = 10)){
          plot.new()
        }
      }
    }
  }
  
  
  dev.off()
}
