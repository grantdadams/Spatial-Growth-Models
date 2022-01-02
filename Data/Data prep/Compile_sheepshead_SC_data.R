# Compile SC data
library(readxl)
library(dplyr)



# Virginia
sc.data <- read.csv("Data\\Raw\\South Carolina\\SCDNR_Sheepshead.csv")
sc.data <- sc.data[which(sc.data$AgeMethod=="otolith"),]
setwd("Data\\Raw\\South Carolina\\EXCEL_DATA\\")
sc.mar <- c()
for (j in 1:(length(dir()))) {
  dat <- as.data.frame(read_excel(dir()[j], sheet = 1))
  colnames(dat) <- gsub(" ","",colnames(dat))
  colnames(dat)[which(colnames(dat)=="Coll_Num")] <- "CollectionNumber"
  colnames(dat)[which(colnames(dat)=="Collection#")] <- "CollectionNumber"
  colnames(dat)[which(colnames(dat)=="Spec#")] <- "Spec_Num"
  colnames(dat)[which(colnames(dat)=="Specimen#")] <- "Spec_Num"
  colnames(dat)[which(colnames(dat)=="Spec_Num")] <- "Spec_Num"
  colnames(dat)[which(colnames(dat)=="#ofannuli")] <- "Annuli"
  dat$CollectionNumber <- ifelse(nchar(dat$CollectionNumber)==6, paste("20",dat$CollectionNumber, sep = ""),dat$CollectionNumber)
  dat <- dat[,which(colnames(dat)%in%c("CollectionNumber","Spec_Num","Annuli"))]
  sc.mar <- rbind(sc.mar, dat)
}
sc.mar[,1] <- as.character(sc.mar[,1])
sc.mar[,2] <- as.character(sc.mar[,2])
sc.data[,1] <- as.character(sc.data[,1])
sc.data[,2] <- as.character(sc.data[,2])

sc.data <- merge(sc.data, sc.mar, by = c("CollectionNumber","Spec_Num"))
sc.data$Annuli<- grepl("?","",sc.data$Annuli)
sc.data$Annuli<-as.numeric(sc.data$Annuli)

colnames(sc.data) <- tolower(colnames(sc.data))
colnames(sc.data)[which(colnames(sc.data)=="tl_mm")] <- "TL_mm"
colnames(sc.data)[which(colnames(sc.data)=="fishweight_g")] <- "wet_weight_grams"
colnames(sc.data)[which(colnames(sc.data)=="annuli")] <- "no_rings"
sc.data$state <- "NC"
sc.data$common_name <- "Sheepshead"
sc.data <- sc.data[which(!is.na(sc.data$no_rings)),]
sc.data$date <- as.Date(sc.data$datefull, "%d-%b-%y")
sc.data$sex <- ifelse(sc.data$sex == "FEMALE", "F" ,  ifelse(sc.data$sex == "MALE", "M", "U"))
sc.data$fi_dependent <- ifelse(sc.data$gear_description == 'GILL NET - MESH G - 3" STRETCH', "FI", ifelse(sc.data$type == "TRAMMEL NET - 200 YD X 8 FT - 14 & 2.5 IN STR MESH", "FI", "RECREATIONAL"))



