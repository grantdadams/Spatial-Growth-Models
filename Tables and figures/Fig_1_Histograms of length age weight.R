


source("Data/Data prep/VBGM_Data_prep.R")
state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS","LA","TX"), state_no = c(1:length(unique(dat$state))))
library(FSA)

cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")

tiff(file="Fig1_Histogram_of_Lengths_BW.tiff" , height= 6, pointsize=14,  width= 14, res=300  , units = "in", family = "serif")
par(mar=c(0.5 , 1.75 , .25 , .5) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75, fg = "white", bg = NA, col.axis = "white", col.lab ="white")
layout(matrix(c(1:18), 3, 6, byrow = T ))

loop.ind <- c(NA,1,2,3,4,5,NA,6,7,8,9,10,NA,NA,NA)    
for ( j in 1:18){
  if(is.na(loop.ind[j]) == T){plot.new()}
  else{
    
    dat.sub.females <- dat[which(dat$state_no == loop.ind[j] & dat$sex == 1),]
    dat.sub.males <- dat[which(dat$state_no == loop.ind[j] & dat$sex == 2),]
    
    female.hist <- hist(dat.sub.females$FL_mm, breaks = seq(from = 0, to = 700, by = 25), plot = F)
    male.hist <- hist(dat.sub.males$FL_mm, breaks = seq(from = 0, to = 700, by = 25), plot = F)
    
    ylim.hist <- c(0, max(female.hist$counts, male.hist$counts))
    
    
    # plot them
    plot(NA,NA, main=NA, xlim = c(0, 700),las = 1, xlab = NA, ylim = ylim.hist, xaxt = "n")
    segments(female.hist$mids - (25/4), 0, female.hist$mids - (25/4), female.hist$counts, lwd = 2, col = cols[1])
    segments(female.hist$mids + (25/4), 0, female.hist$mids + (25/4), male.hist$counts, lwd = 2, col = cols[2])

    # Axis
    if(loop.ind[j] %in% c(6:10)){
      mtext(side = 1, line = 2, text = "Fork Length (mm)")
    }
    
    
    if(loop.ind[j] %in% c(1)){
      legend("topleft", c(paste(state.id$state[which(state.id$state_no == loop.ind[j])], sep = ""),"Females", "Males"), lwd = c(NA,2,2), bty = "n", col = c(NA,cols), cex = .75)
    }
    if(loop.ind[j] %in% c(2:10)){
      legend("topright", paste(state.id$state[which(state.id$state_no == loop.ind[j])], sep = ""), bty = "n", cex = 1)
    }
    
    
  }
  
  
  
  
  if (loop.ind[j] %in% c(6:10)){
    axis(side = 1)
    #mtext(side = 1, "Fork Length (mm)", line = 2,cex=0.75)
    }
  if (loop.ind[j] %in% c(1,6)){
    mtext(side = 2, "Count", line = 2.25, cex = 0.75)
    }
}

dev.off()




tiff(file="Figures/Histogram_of_ages.tiff" , height= 10, pointsize=18,  width=7 , res=300  , units = "in", family = "serif")
par(mar=c(0.5 , 1.85 , .25 , .5) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:18), 6, 3, byrow = T ))

loop.ind <- c(NA,1,2,NA,3,4,NA,5,6,NA,7,8,NA,9,10,NA,NA,NA)    
for ( j in 1:18){
  if(j %in% c(1,4,7,10,13,16,17,18)){plot.new()}
  else{
    
    dat.sub.females <- dat[which(dat$state_no == loop.ind[j] & dat$sex == 1),]
    dat.sub.males <- dat[which(dat$state_no == loop.ind[j] & dat$sex == 2),]
    
    female.hist <- hist(dat.sub.females$fractional_age, breaks = seq(from = 0, to = 42, by = 2), plot = F)
    male.hist <- hist(dat.sub.males$fractional_age, breaks = seq(from = 0, to = 42, by = 2), plot = F)
    
    ylim.hist <- c(0, max(female.hist$counts, male.hist$counts))
    
    
    # plot them
    plot(NA,NA, main=NA, xlim = c(0, 40),las = 1, xlab = NA, ylim = ylim.hist, xaxt = "n")
    segments(female.hist$mids - (2/4), 0, female.hist$mids - (2/4), female.hist$counts, lwd = 2, col = cols[1])
    segments(female.hist$mids + (2/4), 0, female.hist$mids + (2/4), male.hist$counts, lwd = 2, col = cols[2])
    

    if(loop.ind[j] %in% c(1)){
      legend("topright", c(paste(state.id$state[which(state.id$state_no == loop.ind[j])], sep = ""),"Females", "Males"), lwd = c(NA,2,2), bty = "n", col = c(NA,cols), cex = .75)
    }
    if(loop.ind[j] %in% c(2:10)){
      legend("topright", paste(state.id$state[which(state.id$state_no == loop.ind[j])], sep = ""), bty = "n", cex = .75)
    }

    
  }
  
  
  
  
  if (loop.ind[j] %in% c(9,10)){
    axis(side = 1)
    mtext(side = 1, "Fractional Age (yr)", line = 2,cex=0.75)}
  if (loop.ind[j] %in% c(1,3,5,7,9)){
    mtext(side = 2, "Count", line = 2.25, cex = 0.75)
  }
}

dev.off()


source("Data/Data prep/WAL_Data_prep.R")


tiff(file="von Bertalanffy/Figures/Histogram_of_weights2.tiff" , height= 10, pointsize=18,  width=7 , res=300  , units = "in", family = "serif")
par(mar=c(0.5 , 1.75 , .25 , .5) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)
layout(matrix(c(1:18), 6, 3, byrow = T ))

loop.ind <- c(NA,1,2,NA,3,4,NA,5,6,NA,7,8,NA,9,10,NA,NA,NA)    
for ( j in 1:18){
  if(j %in% c(1,4,7,10,13,14,15,16,17,18)){plot.new()}
  else{
    
    dat.sub.females <- dat[which(dat$state_no == loop.ind[j] & dat$sex == 1),]
    dat.sub.males <- dat[which(dat$state_no == loop.ind[j] & dat$sex == 2),]
    
    female.hist <- hist(dat.sub.females$wet_weight_grams, breaks = seq(from = 0, to = 7000, by = 500), plot = F)
    male.hist <- hist(dat.sub.males$wet_weight_grams, breaks = seq(from = 0, to = 7000, by = 500), plot = F)
    
    ylim.hist <- c(0, max(female.hist$counts, male.hist$counts))
    
    
    # plot them
    plot(NA,NA, main=NA, xlim = c(0, 7000),las = 1, xlab = NA, ylim = ylim.hist, xaxt = "n")
    segments(female.hist$mids - (500/4), 0, female.hist$mids - (500/4), female.hist$counts, lwd = 2, col = cols[1])
    segments(female.hist$mids + (500/4), 0, female.hist$mids + (500/4), male.hist$counts, lwd = 2, col = cols[2])
    
    
    
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

