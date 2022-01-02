load("Data/Sheepshead_data_GOM.RData" )
head(dat); tail(dat)

dat$sex <- ifelse(dat$sex == "F", 0, ifelse(dat$sex == "M", 1, NA))
dat <- dat[which(!is.na(dat$sex)),]
dat = dat[which(!is.na(dat$fractional)),]
dat = dat[which(!is.na(dat$state_no)),]
dat = dat[which(!is.na(dat$FL_mm)),]
dat = dat[which(!is.na(dat$fi_dependent)),]
dat$source <- ifelse(dat$fi_dependent == "FI", 0 , 1)



#########################
# Min size in group - 1
#########################

min <- aggregate(dat[,c("state","FL_mm")], by = list(dat[,c("state","FL_mm")]$state), FUN = function(x) min(x, na.rm = T))
min <- min[,c(2:3)]
colnames(min)[2] <- "selectivity"
min$selectivity <- floor(min$selectivity - 1)


max <- aggregate(dat[,c("state","FL_mm")], by = list(dat$state), FUN = max, na.rm=TRUE)
max <- max[,c(2:3)]
colnames(max)[2] <- "max_selectivity"
max$max_selectivity <- ceiling(max$max_selectivity+1)

min <- merge(min,max,by=c("state"))
dat <- merge(dat, min, by = c("state"))


#########################
# Min size in group - 1 for sex
#########################

min <- aggregate(dat[,c("state","FL_mm","sex")], by = list(dat[,c("state","FL_mm","sex")]$state , dat[,c("state","FL_mm","sex")]$sex), FUN = function(x) min(x, na.rm = T))
min <- min[,c(3:5)]
colnames(min)[2] <- "sex.selectivity"
min$sex.selectivity <- floor(min$sex.selectivity - 1)


max <- aggregate(dat[,c("state","FL_mm","sex")], by = list(dat$state , dat$sex), FUN = max, na.rm=TRUE)
max <- max[,c(3:5)]
colnames(max)[2] <- "max_sex.selectivity"
max$max_sex.selectivity <- ceiling(max$max_sex.selectivity+1)

min <- merge(min,max,by=c("state","sex"))

dat <- merge(dat, min, by = c("state","sex"))


rm(min, max)


# Add spatial data
library(rgdal)
library(sp)
state.id <- data.frame(state = c("VA_ocean","VA_bay","NC","SC", "GA","FL_atlantic","FL_gulf","AL","MS","LA","TX"), state_no = c(1:11))
management_centroids <- rgdal::readOGR("Data/State boundaries", "coastal_management_centroids_w_ga")
management_centroids <- as.data.frame(management_centroids)[,c(16,15,7)]
management_centroids <- management_centroids[order(management_centroids$NAME),]
management_centroids$NAME <- sort(c("VA_ocean","VA_bay","NC","SC", "GA","FL_atlantic","FL_gulf","AL","MS","LA","TX"))
colnames(management_centroids) <- c("Lat","Lon","state")
management_centroids <- merge(management_centroids, state.id, by = "state")
management_centroids <- management_centroids[order(management_centroids$state_no),]
rm(state.id)

state_no_df <- subset(dat, !duplicated(state_no))
state_no_df <- state_no_df[,c("state_no", "Lat", "state")]
state_no_df <- state_no_df[order(state_no_df$state_no),]
write.csv(state_no_df, file = "Data/State_centroid_latitudes.csv")
