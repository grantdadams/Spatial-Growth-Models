library(rstan)

##### RE AD AND PREPARE DATA #####
source("Weight_at_length_data_prep.R")


tiff(file="von Bertalanffy/Figures/Histogram_of_weights2.tiff" , height= 10, pointsize=18,  width=7 , res=300  , units = "in", family = "serif")
par(mar=c(0.5 , 1.75 , .25 , .5) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:18), 6, 3, byrow = T ))
loop.ind <- c(NA,1,2,NA,3,4,NA,5,6,NA,7,8,NA,NA,NA) 
state.ind <- c(NA,1,2,NA,3,4,NA,6,7,NA,8,9,NA,NA,NA)  

# Put mcmc of parameters into a list to index later


for ( j in 1:15){
  if(is.na(loop.ind[j]) == T){plot.new()}

    dat.sub.females <- dat[which(dat$state_no == state.ind[j] & dat$sex == 1),]
    dat.sub.males <- dat[which(dat$state_no == state.ind[j] & dat$sex == 2),]
    
    female.hist <- hist(dat.sub.females$wet_weight_grams, breaks = seq(from = 0, to = 7000, by = 500), plot = F)
    male.hist <- hist(dat.sub.males$wet_weight_grams, breaks = seq(from = 0, to = 7000, by = 500), plot = F)
    
    ylim.hist <- c(0, max(female.hist$counts, male.hist$counts))
    
    
    # plot them
    plot(NA,NA, main=NA, xlim = c(0, 7000),las = 1, xlab = NA, ylim = ylim.hist, xaxt = "n")
    segments(female.hist$mids - (500/4), 0, female.hist$mids - (500/4), female.hist$counts, lwd = 2, col = cols[1])
    segments(male.hist$mids + (500/4), 0, male.hist$mids + (500/4), male.hist$counts, lwd = 2, col = cols[2])
    
    
    
    if(loop.ind[j] %in% c(1)){
      legend("topleft", c(paste(state.id$state[which(state.id$state_no == loop.ind[j])], sep = ""),"Females", "Males"), lwd = c(NA,2,2), bty = "n", col = c(NA,cols), cex = .75)
    }
    if(loop.ind[j] %in% c(2:10)){
      legend("topright", paste(state.id$state[which(state.id$state_no == loop.ind[j])], sep = ""), bty = "n", cex = .75)
    }
    
  }
  
  
  
  
  if (loop.ind[j] %in% c(7,8)){
    axis(side = 1)
    mtext(side = 1, "Weight (g)", line = 2,cex=0.75)}
  if (loop.ind[j] %in% c(1,3,5,7)){
    mtext(side = 2, "Count", line = 2.25, cex = 0.75)
  }
}

dev.off()

