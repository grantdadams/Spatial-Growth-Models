library(rstan)

##### RE AD AND PREPARE DATA #####
source("Data/Data prep/WAL_Data_prep.R")
prev_coefs <- read.csv("Data/WAL_Coefs_Previous_Studies.csv")

# Load data
load("Stan models/Models/WAL_Stan_models_2018_02_19.RData")

# Get parameters
model <- mod_list[[2]]
mod_summ <- data.frame(summary(model)[[1]])
mod_summ$parameter <- rownames(mod_summ)

######################################################
# reorder state variables
state_no_df <- subset(dat, !duplicated(state_no_check))
state_no_df <- state_no_df[,c("state_no", "state_no_check", "state")]
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

cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")


# Get vectors of all VBGF Parameters
a.all <- list(matrix(NA, nrow = nrow(draws), ncol = 8), matrix(NA, nrow = nrow(draws), ncol = 8))
b.all <- list(matrix(NA, nrow = nrow(draws), ncol = 8), matrix(NA, nrow = nrow(draws), ncol = 8))


for(j in 1:8){
  for ( i in 1:2){
    state_ind <- c(state_no_df$state_no_check)[j] # Indexing variable
    
    log_a_re <- draws[,paste("log_a_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_a_re <- draws[,paste("B_log_a[",1:3,"]", sep = "")] # Subset the mcmc chains
    a.sub <- exp(as.matrix(betas_a_re) %*% as.vector(c(1, i, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_a_re)
    a.all[[i]][,j] <- a.sub
    
    log_b_re <- draws[,paste("log_b_re[",state_ind,"]", sep = "")] # Subset the mcmc chains
    betas_b_re <- draws[,paste("B_log_b[",1:3,"]", sep = "")] # Subset the mcmc chains
    b.sub <- exp(as.matrix(betas_b_re) %*% as.vector(c(1, i, lat_lon_df[which(lat_lon_df$state_no_check == state_ind),2])) + log_b_re)
    b.all[[i]][,j] <- b.sub
  }
}

params_list <- list(a.all, b.all)
sex <- c("Female","Male")

########################################################
figure_names <- c("Tables and Figures/Fig_A7_Paired_a_Parameter_Differences.tiff", 
                  "Tables and Figures/Fig_A8_Paired_b_Parameter_Differences.tiff")

parameter_names <- list(
  expression(paste("Paired differences in ",italic("a"["s,r"]))),
  expression(paste("Paired differences in ",italic("b"["s,r"]))))

# Set limits and parameter
xlim <- list(c(-.000024, .000053),
             c(-0.19, 0.145))

xaxis_tick <- list(seq(-.00002, .00004, by = 0.00002), 
                   seq(-0.1, .2, by = .1))

xaxis_lab <- list(c(NA,0, NA, 4e-5),
                  c(NA, 0, 0.1, NA))


#### Plot the thing ####
for(k in 1:2){
  
  tiff(file = figure_names[k] , height= 140, pointsize=12,  width=170 , res=300  , units = "mm", family = "serif")
  par(mar=c(0 , 0 , 0 , 0), tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
  layout(matrix(c(1:81), 9, 9, byrow = F ), widths = c(0.1, rep(1, 8)))
  
  # Add left border
  for(j in 1:9){
    plot.new()
  }
  
  # Data frame of unique values
  combo_df <- expand.grid(c(1:8), c(1:8))
  combo_df <- combo_df[c(2,1)]
  
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
            if(i == 1){
              rect(par("usr")[1], par("usr")[3], par("usr")[2], 0.5, col = "#E1E1E1", border = NA)
            }
            if(i ==2){
              rect(par("usr")[1], 0.5, par("usr")[2], par("usr")[4], col = "#E1E1E1", border = NA)
            }
          }
          
          # Color background if 0% CI cross 0
          if(prod(sign(sign(quantiles_10_90))) < 0){
            if(i == 1){
              rect(par("usr")[1], par("usr")[3], par("usr")[2], 0.5, col = "#AAAAAA", border = NA)
            }
            if(i ==2){
              rect(par("usr")[1], 0.5, par("usr")[2], par("usr")[4], col = "#AAAAAA", border = NA)
            }
          }
          
          spacing.ind <- 0.25 + (i-1) * 0.5 # Spacing index for plots
          
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
        if(j %in% seq(from = 8, to = 64, by = 8)){
          #axis(side = 1, cex.axis = 1.2)
          axis(side = 1, at = xaxis_tick[[k]] , labels = xaxis_lab[[k]], cex.axis = 1.2)
        }
        
        # Add x-axis label
        if(j == 32){
          mtext(side = 1, parameter_names[[k]], cex = 1, line = 2.5)
        }
        
        # Add space at the bottom
        if(j %in% seq(from = 8, to = 64, by = 8)){
          plot.new()
        }
      }
    }
  }
  
  
  dev.off()
}
