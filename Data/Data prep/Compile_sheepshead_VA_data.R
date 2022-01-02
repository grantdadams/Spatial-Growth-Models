library(readxl)
library(dplyr)



# Virginia
setwd("Data/VBGM Data Prep/CSV Data Files")
virginia.data <- c()
for (j in 1:(length(dir()))) {
  dat <- read.csv(dir()[j])
  colnames(dat) <- tolower(colnames(dat))
  dat <- dat[,which(colnames(dat)%in%c("year","month","day","weight","gram","fork","sl..mm.","total","sex","gear","type","area","geartype","otoage","state", "mi","prev.inc"))]
  virginia.data <- bind_rows(virginia.data, dat)
}

colnames(virginia.data)[which(colnames(virginia.data)=="mi")] <- "margin_code"
colnames(virginia.data)[which(colnames(virginia.data)=="sl..mm.")] <- "SL_mm"
colnames(virginia.data)[which(colnames(virginia.data)=="fork")] <- "FL_mm"
colnames(virginia.data)[which(colnames(virginia.data)=="total")] <- "TL_mm"
colnames(virginia.data)[which(colnames(virginia.data)=="gram")] <- "wet_weight_grams"
colnames(virginia.data)[which(colnames(virginia.data)=="otoage")] <- "no_rings"
virginia.data$state <- ifelse(virginia.data$area == "BAY", "VA_bay", "VA_ocean")
virginia.data$common_name <- "Sheepshead"
virginia.data <- virginia.data[which(!is.na(virginia.data$no_rings)),]
virginia.data$date <- as.Date(paste(virginia.data$year, virginia.data$month, virginia.data$day, sep="-"), "%Y-%m-%d")
virginia.data$sex <- as.character(virginia.data$sex)
virginia.data$fi_dependent <- ifelse(virginia.data$type == "FI", "FI", ifelse(virginia.data$type == "C", "COMMERCIAL", "RECREATIONAL"))

