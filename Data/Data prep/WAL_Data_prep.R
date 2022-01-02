
load("Data/Sheepshead_data_GOM.RData" )

# Female is 1, male is 2
dat$sex <- ifelse(dat$sex == "F", 0, ifelse(dat$sex == "M", 1, NA))
dat <- dat[which(!is.na(dat$sex)),]
dat = dat[which(!is.na(dat$wet_weight_grams)),]
dat = dat[which(dat$wet_weight_grams < 7500 & dat$wet_weight_grams > 0),]
dat = dat[which(!is.na(dat$state_no)),]
dat = dat[which(!is.na(dat$FL_mm)),]
dat = dat[which(!is.na(dat$fi_dependent)),]
dat$source <- ifelse(dat$fi_dependent == "FI", 0 , 1)


state.id <- data.frame(state = c("VA Ocean","VA Bay","NC","SC","FL Atlantic","FL Gulf","AL","MS","LA","TX"), state_no = c(1:4,6:11), state_no_check = c(1:10))

dat <- merge(dat, state.id, by = "state_no")


# Subset data outside prediction interval
mod_simp <- lm(log(dat$wet_weight_grams) ~ log(dat$FL_mm))
#mod_norm <- nls(LW.dat[,2] ~ a* LW.dat[,1] ^b, start = list(a = 0.00003, b = 3))

pred <- data.frame(predict(mod_simp, data.frame(log(dat$FL_mm)), interval = "prediction", level=.9999))

dat$wet_weight_grams <- ifelse(log(dat$wet_weight_grams) > pred$upr, NA, dat$wet_weight_grams)
dat$wet_weight_grams <- ifelse(log(dat$wet_weight_grams) < pred$lwr, NA, dat$wet_weight_grams)
dat = dat[which(!is.na(dat$wet_weight_grams)),]

plot(log(dat$FL_mm), log(dat$wet_weight_grams))


# Plot it
for ( j in 1:length(unique(dat$state_no))){
  dat.sub <- subset(dat, dat$state_no == j)
  if (j == 1){
    plot(y = dat.sub$wet_weight_grams , x = dat.sub$FL_mm, col = j, xlab = "Fork length (mm)", ylab = "Weight (g)", cex = .5, cex.lab = 1.25, pch = 20, ylim = c(min(dat$wet_weight_grams), max(dat$wet_weight_grams)), xlim = c(min(dat$FL_mm), max(dat$FL_mm)))
    
  }  else {
    if(j == 5){
      points(dat.sub$TL_mm , dat.sub$wet_weight_grams, col = j, cex = .5, pch = 20)
    } else{
      points(dat.sub$FL_mm , dat.sub$wet_weight_grams, col = j, cex = .5, pch = 20)
    }
  } 
}

legend("topleft", legend = c(unique(dat$state_no)), col = c(1:8), pch = 16, cex = 1.25, bty ="n", inset = -0.02)
legend("topright", c("Female","Male"), pch = 20, col = c(1:8), pt.cex = 2, cex = 1.25, bty ="n")


# Write table 4
filename <- "Tables and figures/Table_4_Weight_at_Length_data_summary.csv"
library(tidyr)
# Sample size
summ_stat <- data.frame(table(dat$state_no, dat$sex))
summ_stat <- tidyr::spread(summ_stat, Var2, Freq)
colnames(summ_stat) <- c("state_no", "n_f", "n_m")

# Add state names and reorder
state.id <- data.frame(state = c("VA_ocean","VA_bay","NC","SC", "GA","FL_atlantic","FL_gulf","AL","MS","LA","TX"), state_no = c(1:11))
summ_stat <- merge(summ_stat, state.id, by = "state_no" )
summ_stat <- summ_stat[order(summ_stat$state_no),]


