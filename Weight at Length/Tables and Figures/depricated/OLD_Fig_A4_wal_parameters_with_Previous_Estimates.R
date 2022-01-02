setwd("C:/Users/w966213/Dropbox/Sheepshead/Sheepshead R Work/Weight-at-length/Figures")
source("Weight-at-length 3-level data prep.R")
setwd("C:/Users/w966213/Dropbox/Sheepshead/Sheepshead R Work/Weight-at-length/Figures")



# This code creates a to illustrate the 95% and 80% credible intervals of parameter estimates from a hierarchical weight-at-length power equation following Midway et al. 2015


################
# Midway plot ##
################
cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")

tiff(file="Weight_at_length_parameters_with_previous.tiff" , height= 5, pointsize=18,  width=11 , res=300  , units = "in", family = "serif")
par(mar=c(3.7 , 0.5 , .25 , .125) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0.5, 0 , 0.5), cex = 0.75)
layout(matrix(c(1,2,3,4), 1, 4, byrow = F ))

# Put mcmc of parameters into a list to index later
param.list <- list(a.mcmc, b.mcmc)
param.names <- c(expression(italic("a"["s,r"])), expression(italic("b"["s,r"])))

for(m in 1:length(param.list)){
  
  # Get xlim from quantiles of mcmc
  xlim.low <- ifelse(m==1, -0.00001, quantile(as.matrix(param.list[[m]]), c(0.001,.9975))[1])
  xlim.high <- quantile(as.matrix(param.list[[m]]), c(0.001,.9975))[2]
  
  plot(NA,NA, main = NA, xlab =param.names[m], yaxt = "n", ylab= NA, xlim = c(xlim.low, xlim.high), ylim = c(-3.5,6.5), cex.lab = 1.25)
  sex <- c("Female","Male")
  
  for ( j in 1:length(unique(dat$state_no))){
    for ( i in 1:2){
      ind <- j + length(unique(dat$state_no)) * (i-1)# Indexing variable
      
      dat.sub <- param.list[[m]][,ind] # Subset the JAGS mcmc chains
      quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
      quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
      
      spacing.ind <- 7-(j + (.35*(i-1))) # Spacing index for plots
      
      # Plot the credible intervals
      segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 4 , col = cols.points[i])
      segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 6 , col = cols[i])
      
      # Plot the median estimate
      med <- median(dat.sub)
      segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
      
      # Plot the other previous studies params params
      if(length(which(prev_coefs$state_no == j)) > 0){
        prev_coefs_sub <- prev_coefs[which(prev_coefs$state_no == j),]
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
      
      
      # Label lines
      if(i==1){
        text(mean(dat.sub), spacing.ind+.25, paste(state.id[which(state.id$state_no == j),1],sep=" ") , cex = 0.75) 
      }
      
    }
  }
  if(m == 1){
  legend("bottomright", c(sex,"Combined","Previous"), pch = c(NA,NA,NA,21), lty = c(1,1,1,NA), lwd=c(3,3,3,NA), col = c(cols,"grey", 1), cex = 0.75, bty ="n")
  }
}
dev.off()
