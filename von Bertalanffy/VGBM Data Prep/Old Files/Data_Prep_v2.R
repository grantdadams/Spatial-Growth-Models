#####################################################################################################
####################### Step 1: Data Prep and Consolidation for VBGM ################################
#####################################################################################################
# This file takes length-at-age data from state fisheries-independent sampling programs and the GulfFin fishery-dependent sampling data and compiles them together in one format.

# Following Vanderkooy 2009 the biological birthdata is assinged to April 1. 

# Set WD
setwd("C:/Users/w966213/Dropbox/Sheepshead/Sheepshead R Work/von Bertalanffy/VGBM Data Prep/CVS Data Files")

# Load Data
florida.data <- read.csv("Sheepshead_age_data_Florida.csv")
alabama.data <- read.csv("Sheepshead_age_data_Alabama.csv")
mississippi.data <- read.csv("Sheepshead_age_data_Mississippi.csv")
louisiana.data <- read.csv("Sheepshead_age_data_Louisiana.csv")
gulffin.data <- read.csv("Sheepshead_age_data_GulfFin.csv")
north.carolina.data <- read.csv("Sheepshead_age_data_North_Carolina.csv")


library(dplyr)
library(plyr)
library(readxl)


# Filter out species that are not Sheepshead
alabama.data <- alabama.data[which(alabama.data$spec=="shd"),]
louisiana.data <- louisiana.data[which(louisiana.data$Common_Name=="Sheepshead"),]
gulffin.data <- gulffin.data[which(gulffin.data$Common.Name=="SHEEPSHEAD"),]

###########################################
# Standardize column names for the datasets
###########################################
# Florida
colnames(florida.data) <- c("field_no" ,"ID","sex", "bay" , "SL_mm" , "FL_mm" , "TL_mm" , "wet_weight_grams" ,"no_rings" , "margin_code", "age_1jan", "date", "birthdate" , "inch_FL" ,"year","gear","program","coast","project","month","day")

florida.data$state <- ifelse(florida.data$coast == "Gulf", "FL_gulf", "FL_atlantic")
florida.data$common_name <- "Sheepshead"
florida.data <- florida.data[which(!is.na(florida.data$no_rings)),]
florida.data$date <- as.Date(florida.data$date, "%m/%d/%Y")
florida.data$birthdate <- as.Date(florida.data$birthdate, "%m/%d/%Y")
florida.data$sex <- as.character(florida.data$sex)
#florida.data$no_rings <- florida.data$no_rings -1
florida.data$fi_dependent <- "FI" # 1 indicates fishery dependent sampleing 0 indicates fishery-independent sampling
#florida.data$selectivity <- 0
florida.data <- florida.data[,which(colnames(florida.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


# Florida gulf-fin
florida.data.GSMFC <- gulffin.data[which(gulffin.data$State.Sampled=="FL"),]
colnames(florida.data.GSMFC) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code" , "Reader.Id.1" , "Access.Site.Id")
florida.data.GSMFC$wet_weight_grams <- florida.data.GSMFC$wet_wt_kg * 1000
florida.data.GSMFC <- florida.data.GSMFC[which(!is.na(florida.data.GSMFC$no_rings)),]
# Date function for format
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
florida.data.GSMFC$date <- as.Date(paste(florida.data.GSMFC$year, gsub(" ", "",capwords(tolower(florida.data.GSMFC$month))), florida.data.GSMFC$day, sep="-"), "%Y-%B-%d")
florida.data.GSMFC$sex <- as.character(florida.data.GSMFC$sex)
florida.data.GSMFC$sex <- ifelse(florida.data.GSMFC$sex == "FEMALE", "F" ,  ifelse(florida.data.GSMFC$sex == "MALE", "M", "U"))
florida.data.GSMFC$fi_dependent <- ifelse(florida.data.GSMFC$fishery == "COMMERCIAL","COMMERCIAL", "RECREATIONAL")
#florida.data.GSMFC$selectivity <- 304.8 # 12 inch size limit
florida.data.GSMFC$state <- "FL_gulf"
florida.data.GSMFC <- florida.data.GSMFC[,which(colnames(florida.data.GSMFC)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]



# Alabama
colnames(alabama.data) <- c("year","month","day", "field_no","common_name", "set","area","mesh","FL_mm","TL_mm", "wet_wt_kg", "sex", "gonad_weight_grams", "no_rings" ,"margin_code", "X")
alabama.data$state <- "AL"
alabama.data$common_name <- "Sheepshead"
alabama.data$wet_weight_grams <- as.numeric(alabama.data$wet_wt_kg) * 1000
alabama.data$no_rings <- as.character(alabama.data$no_rings)
alabama.data <- alabama.data[which(alabama.data$no_rings!="" & alabama.data$no_rings!="missing stones" ),]
alabama.data$no_rings <- as.numeric(alabama.data$no_rings)
alabama.data <- alabama.data[which(!is.na(alabama.data$no_rings)),]
alabama.data$date <- as.Date(paste(alabama.data$year, alabama.data$month, alabama.data$day, sep="-"), "%Y-%m-%d")
alabama.data$sex <- as.character(alabama.data$sex)
alabama.data$fi_dependent <- "FI" # 1 indicates fishery dependent sampleing 0 indicates fishery-independent sampling
#alabama.data$selectivity <- 0
alabama.data <- alabama.data[,which(colnames(alabama.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]

# Alabama GSMFC
alabama.data.GSMFC <- gulffin.data[which(gulffin.data$State.Sampled=="AL"),]
colnames(alabama.data.GSMFC) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code" , "Reader.Id.1" , "Access.Site.Id")
alabama.data.GSMFC$wet_weight_grams <- alabama.data.GSMFC$wet_wt_kg * 1000
alabama.data.GSMFC <- alabama.data.GSMFC[which(!is.na(alabama.data.GSMFC$no_rings)),]
alabama.data.GSMFC <- alabama.data.GSMFC[which(alabama.data.GSMFC$fishery!="FISHERY INDEPENDENT"),]
# Date function for format
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
alabama.data.GSMFC$date <- as.Date(paste(alabama.data.GSMFC$year, gsub(" ", "",capwords(tolower(alabama.data.GSMFC$month))), alabama.data.GSMFC$day, sep="-"), "%Y-%B-%d")
alabama.data.GSMFC$sex <- as.character(alabama.data.GSMFC$sex)
alabama.data.GSMFC$sex <- ifelse(alabama.data.GSMFC$sex == "FEMALE", "F" ,  ifelse(alabama.data.GSMFC$sex == "MALE", "M", "U"))
alabama.data.GSMFC$fi_dependent <- ifelse(alabama.data.GSMFC$fishery == "COMMERCIAL","COMMERCIAL", "RECREATIONAL")
alabama.data.GSMFC <- alabama.data.GSMFC[,which(colnames(alabama.data.GSMFC)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]

##############################
#alabama.data.GSMFC$selectivity <- 0
#alabama.data.GSMFC$selectivity[which(alabama.data.GSMFC$year > 2012)] <- 304.8

# Mississippi
colnames(mississippi.data) <- c("date","month" , "spec_no" , "TL_mm", "FL_mm","SL_mm","wet_weight_grams","sex","gross_maturity","biological_age","margin_code","gonad_weight_grams","gutted_weight_grams","GSI","histological_class","development_subclass")
mississippi.data$state <- "MS"
#mississippi.data$age_1jan <- ifelse(mississippi.data$month > 3 , mississippi.data$age, mississippi.data$age + 1)
mississippi.data$gear <- NA
mississippi.data$common_name <- "Sheepshead"
mississippi.data$no_rings <- mississippi.data$biological_age
mississippi.data$no_rings[which(mississippi.data$no_rings==-1)] = NA
mississippi.data$margin_code[which(mississippi.data$margin_code==-1)] = NA
mississippi.data <- mississippi.data[which(!is.na(mississippi.data$no_rings)),]

mississippi.data$date <- as.Date(mississippi.data$date, "%m/%d/%Y")
mississippi.data$sex <- ifelse(mississippi.data$sex == 1 , "M", "F")
mississippi.data$margin_code[which(mississippi.data$margin_code==3)] = 1# Convert margin codes
mississippi.data$margin_code[which(mississippi.data$margin_code==4)] = 2# Convert margin codes
mississippi.data$margin_code[which(mississippi.data$margin_code==5)] = 3# Convert margin codes
mississippi.data$margin_code[which(mississippi.data$margin_code==6)] = 4# Convert margin codes
mississippi.data$fi_dependent <- "FI"

#mississippi.data$selectivity <- 0
mississippi.data <- mississippi.data[,which(colnames(mississippi.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


### GSMFC MS
mississippi.data.GSMFC <- gulffin.data[which(gulffin.data$State.Sampled=="MS"),]
colnames(mississippi.data.GSMFC) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code" , "Reader.Id.1" , "Access.Site.Id")
mississippi.data.GSMFC$wet_weight_grams <- mississippi.data.GSMFC$wet_wt_kg * 1000
mississippi.data.GSMFC <- mississippi.data.GSMFC[which(!is.na(mississippi.data.GSMFC$no_rings)),]
# Date function for format
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
mississippi.data.GSMFC$date <- as.Date(paste(mississippi.data.GSMFC$year, gsub(" ", "",capwords(tolower(mississippi.data.GSMFC$month))), mississippi.data.GSMFC$day, sep="-"), "%Y-%B-%d")
mississippi.data.GSMFC$sex <- as.character(mississippi.data.GSMFC$sex)
mississippi.data.GSMFC$sex <- ifelse(mississippi.data.GSMFC$sex == "FEMALE", "F" ,  ifelse(mississippi.data.GSMFC$sex == "MALE", "M", "U"))
mississippi.data.GSMFC$fi_dependent <- ifelse(mississippi.data.GSMFC$fishery == "COMMERCIAL","COMMERCIAL", "RECREATIONAL")
#mississippi.data.GSMFC$selectivity <- 0
mississippi.data.GSMFC <- mississippi.data.GSMFC[,which(colnames(mississippi.data.GSMFC)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]

# Lousisiana
colnames(louisiana.data) <- c("month" , "day" , "year" , "Interview_id_Num" , "common_name" , "Age_Tag_1" , "FL_mm" , "TL_mm" , "wet_weight_grams" , "gutted_weight_grams" , "sex" , "Mode" , "Distance" , "Gear_Id" , "gear", "Otolith_Cond" , "Cut_Date" , "Section_Annuli_1" , "Edge_S1" , "Section_Annuli_2" , "Edge_S2" , "Fishery" , "Parish_Cnty_Name")
louisiana.data <- louisiana.data[which(louisiana.data$Section_Annuli_1==louisiana.data$Section_Annuli_2),]
louisiana.data$no_rings <- louisiana.data$Section_Annuli_1
louisiana.data$margin_code <- (louisiana.data$Edge_S1 + louisiana.data$Edge_S2)/2
louisiana.data <- louisiana.data[which(!is.na(louisiana.data$no_rings)),]
louisiana.data$state <- "LA"
louisiana.data$date <- as.Date(paste(louisiana.data$year, louisiana.data$month, louisiana.data$day, sep="-"), "%Y-%m-%d")
louisiana.data$sex <- ifelse(louisiana.data$sex == "3 : Female", "F" ,  ifelse(louisiana.data$sex == "2 : Male", "M", "U"))
louisiana.data$fi_dependent <- toupper(louisiana.data$Fishery)
#louisiana.data$selectivity <- 0
#louisiana.data$selectivity[which(louisiana.data$Fishery == "Commercial")] <- 254
louisiana.data <- louisiana.data[,which(colnames(louisiana.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


# Texas
texas.data <- gulffin.data[which(gulffin.data$State.Sampled=="TX"),]
colnames(texas.data) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code" , "Reader.Id.1" , "Access.Site.Id")
texas.data$wet_weight_grams <- texas.data$wet_wt_kg * 1000
texas.data <- texas.data[which(!is.na(texas.data$no_rings)),]
# Date function for format
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
texas.data$date <- as.Date(paste(texas.data$year, gsub(" ", "",capwords(tolower(texas.data$month))), texas.data$day, sep="-"), "%Y-%B-%d")
texas.data$sex <- as.character(texas.data$sex)
texas.data$sex <- ifelse(texas.data$sex == "FEMALE", "F" ,  ifelse(texas.data$sex == "MALE", "M", "U"))
texas.data$fi_dependent <- ifelse(texas.data$fishery == "COMMERCIAL","COMMERCIAL", "RECREATIONAL")
texas.data <- texas.data[,which(colnames(texas.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


# Virginia
setwd("C:/Users/w966213/Dropbox/Sheepshead/Data/Virginia/CSV_files")
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
virginia.data$sex <- ifelse(virginia.data$sex == 2, "F" ,  ifelse(virginia.data$sex == 1, "M", "U"))
virginia.data$fi_dependent <- ifelse(virginia.data$type == "FI", "FI", ifelse(virginia.data$type == "C", "COMMERCIAL", "RECREATIONAL"))
virginia.data <- virginia.data[,which(colnames(virginia.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


# SC
sc.data <- read.csv("C:/Users/w966213/Dropbox/Sheepshead/Data/South Carolina/SCDNR_Sheepshead.csv")
sc.data <- sc.data[which(sc.data$AgeMethod=="otolith"),]
sc.data$TL_mm <- as.numeric(sc.data$TL_mm)
setwd("C:/Users/w966213/Dropbox/Sheepshead/Data/South Carolina/EXCEL_DATA/")
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
sc.mar$no_rings <- as.numeric(sc.mar$Annuli)
sc.data[,1] <- as.character(sc.data[,1])
sc.data[,2] <- as.character(sc.data[,2])

sc.data <- merge(sc.data, sc.mar, by = c("CollectionNumber","Spec_Num"))

colnames(sc.data) <- tolower(colnames(sc.data))
colnames(sc.data)[which(colnames(sc.data)=="tl_mm")] <- "FL_mm"
colnames(sc.data)[which(colnames(sc.data)=="fishweight_g")] <- "wet_weight_grams"
sc.data$state <- "SC"
sc.data$common_name <- "Sheepshead"
sc.data <- sc.data[which(!is.na(sc.data$no_rings)),]
sc.data$date <- as.Date(sc.data$datefull, "%d-%b-%y")
sc.data$sex <- ifelse(sc.data$sex == 2, "F" ,  ifelse(sc.data$sex == 1, "M", "U"))
sc.data$gear_description <- gsub(" ","",sc.data$gear_description)
sc.data$gear_description <- gsub('"',"",sc.data$gear_description)
sc.data$fi_dependent <- ifelse(sc.data$gear_description %in% c("GILLNET-MESHG-3STRETCH", "TRAMMELNET-200YDX8FT-14&2.5INSTRMESH"), "FI", "RECREATIONAL")
sc.data <- sc.data[,which(colnames(sc.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]



## North Carolina

north.carolina.data$fi_dependent <- ifelse(north.carolina.data$FLEET == "independent ", "FI", ifelse(north.carolina.data$FLEET == "commercial", "COMMERCIAL", "RECREATIONAL"))
north.carolina.data$FL_mm <- north.carolina.data$LENGTH..FL..mm.
north.carolina.data$TL_mm <- north.carolina.data$LENGTH..TL..mm.
north.carolina.data$sex <- ifelse(north.carolina.data$SEX == "female", "F" ,  ifelse(north.carolina.data$SEX == "male", "M", "U"))
north.carolina.data$no_rings <- north.carolina.data$AGE
north.carolina.data$wet_weight_grams <- north.carolina.data$WEIGHT..kg.* 1000
north.carolina.data$date <- as.Date(paste(north.carolina.data$YEAR, north.carolina.data$MONTH, north.carolina.data$DAY, sep="-"), "%y-%m-%d")
north.carolina.data$month <- north.carolina.data$MONTH 
north.carolina.data$state <- "NC"
north.carolina.data <- north.carolina.data[,which(colnames(north.carolina.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month","biological_age_day"))]

setwd("C:/Users/w966213/Dropbox/Sheepshead/Sheepshead R Work/von Bertalanffy/VGBM Data Prep/CVS Data Files")




##########################
# Biological age asignment for the GULF
##########################
library(lubridate)

gulf.data <- rbind.fill(florida.data, florida.data.GSMFC,alabama.data,alabama.data.GSMFC,mississippi.data,mississippi.data.GSMFC, louisiana.data,texas.data)
gulf.data <- gulf.data[which(!is.na(gulf.data$margin_code)),]
gulf.data$no_rings[which(gulf.data$no_rings>138)] <- 14

# Bump the annuli. For specimens caught between March 1st and March 31st with a margin code of 1, a constant of 1 was subtracted from the annulus count. While for specimens caught between April 1st and April 31st with a margin code of 4, a constant of 1 was added to the annulus count. 
gulf.data$no_rings <- ifelse( gulf.data$month == 3 & gulf.data$margin_code == 1, gulf.data$no_rings - 1, gulf.data$no_rings  )
gulf.data$no_rings <- ifelse( gulf.data$month == 4 & gulf.data$margin_code >= 3.5, gulf.data$no_rings + 1, gulf.data$no_rings  )

gulf.data$biological_age_day <- gulf.data$no_rings + ((-121  + yday(gulf.data$date))/365)
  
  ifelse(gulf.data$date >= as.Date(paste(year(gulf.data$date), "04" , "01", sep = "-"), "%Y-%m-%d") , gulf.data$no_rings + ((-91 + yday(gulf.data$date))/365), (274 + (gulf.data$no_rings * 365  + yday(gulf.data$date)))/365  )

##########################
# Biological age asignment for the Atlantic
##########################
atlantic.data <- rbind.fill(virginia.data, sc.data, north.carolina.data)
atlantic.data$biological_age_day <-  atlantic.data$no_rings + ((-121  + yday(atlantic.data$date))/365)

#######################
# Compile all the data!
#######################
library(plyr)
gulf.data <- rbind.fill(gulf.data, atlantic.data)



gulf.data <- gulf.data[which(!is.na(gulf.data$no_rings)),]

########################
# SL, FL, TL Regressions
########################
fl.tl <- gulf.data[which(!is.na(gulf.data$FL_mm) & !is.na(gulf.data$TL_mm)),]
fl.tl.regression <- lm(FL_mm ~ TL_mm, data = fl.tl)
fl.tl.intercept <- fl.tl.regression$coefficients[1]
fl.tl.slope <- fl.tl.regression$coefficients[2]
gulf.data$FL_mm[which(is.na(gulf.data$FL_mm))] <- fl.tl.intercept + fl.tl.slope * gulf.data$TL_mm[which(is.na(gulf.data$FL_mm))]

fl.sl <- gulf.data[which(!is.na(gulf.data$FL_mm) & !is.na(gulf.data$SL_mm)),]
fl.sl.regression <- lm(FL_mm ~ SL_mm, data = fl.sl)
fl.sl.intercept <- fl.sl.regression$coefficients[1]
fl.sl.slope <- fl.sl.regression$coefficients[2]
gulf.data$FL_mm[which(is.na(gulf.data$FL_mm))] <- fl.sl.intercept + fl.sl.slope * gulf.data$SL_mm[which(is.na(gulf.data$FL_mm))]


gulf.data <- gulf.data[which(!is.na(gulf.data$FL_mm)),]


##########################
# Plots
##########################

# State Specific
for ( i in 1:length(unique(gulf.data$state))){
  state.vec <- c(unique(gulf.data$state))
  #gulf.data2 <- gulf.data[which(!is.na(gulf.data$biological_age_d)),]
  gulf.data.sub <- subset(gulf.data, gulf.data$state == state.vec[i])
  if ( i == 1){
    plot(gulf.data.sub$biological_age_day,gulf.data.sub$FL_mm, col = i, main = "Without change", xlab = "Age (yrs)", ylab = "Fork length (mm)")
  } else {
    points(gulf.data.sub$biological_age_day,gulf.data.sub$FL_mm , col = i)
  }
}
legend("bottomright", state.vec, pch = 16, col = c(1:length(state.vec)))



# M v F
plot(gulf.data$biological_age_day[which(gulf.data$sex=="F")], gulf.data$FL_mm[which(gulf.data$sex=="F")], col =1)
points(gulf.data$biological_age_day[which(gulf.data$sex=="M")], gulf.data$FL_mm[which(gulf.data$sex=="M")], col = 2)

# M v F and State Specific
for ( i in 1:length(unique(gulf.data$state))){
  state.vec <- c(unique(gulf.data$state))
  gulf.data2 <- gulf.data[which(!is.na(gulf.data$sex)),]
  gulf.data.sub.f <- subset(gulf.data, gulf.data$state == state.vec[i] & gulf.data$sex == "F")
  gulf.data.sub.m <- subset(gulf.data, gulf.data$state == state.vec[i] & gulf.data$sex == "M")
  if ( i == 1){
    plot(gulf.data.sub.f$biological_age_day,gulf.data.sub.f$FL_mm, col = i)
    points(gulf.data.sub.m$biological_age_day,gulf.data.sub.m$FL_mm, col = i , pch = 0)
  } else {
    points(gulf.data.sub.f$biological_age_day,gulf.data.sub.f$FL_mm , col = i)
    points(gulf.data.sub.m$biological_age_day,gulf.data.sub.m$FL_mm , col = i, pch = 0)
  }
}
legend("bottomright", c(state.vec,"Female","Male"), pch = c(rep(16,length(state.vec)),1,0), col = c(1:length(state.vec),1,1))



#########################
# Min size in group - 1
#########################

min <- aggregate(gulf.data[,c("state","FL_mm","fi_dependent")], by = list(gulf.data[,c("state","FL_mm","fi_dependent")]$state , gulf.data[,c("state","FL_mm","fi_dependent")]$fi_dependent), FUN = function(x) min(x))
min <- min[,c(3:5)]
colnames(min)[2] <- "selectivity"
min$selectivity <- floor(min$selectivity - 1)


max <- aggregate(gulf.data[,c("state","FL_mm","fi_dependent")], by = list(gulf.data$state , gulf.data$fi_dependent), FUN = max, na.rm=TRUE)
max <- max[,c(3:5)]
colnames(max)[2] <- "max_selectivity"
max$max_selectivity <- ceiling(max$max_selectivity+1)

min <- merge(min,max,by=c("state","fi_dependent"))

gulf.data <- merge(gulf.data, min, by = c("state","fi_dependent"))







#########################
# Save
#########################
dat <- gulf.data

state.id <- data.frame(state = c("VA_ocean","VA_bay","NC","SC","FL_atlantic","FL_gulf","AL","MS","LA","TX"), state_no = c(1:length(unique(dat$state))))
dat <- merge(dat, state.id , by = "state")

dat <- dat[which(!is.na(dat$biological_age_day)),]
dat <- dat[which(!is.na(dat$state_no)),]
dat <- dat[which(!is.na(dat$FL_mm)),]

setwd("C:/Users/w966213/Dropbox/Sheepshead/Sheepshead R Work/von Bertalanffy")
save(dat, file = "Sheepshead_data_GOM.RData" )

  
  