#####################################################################################################
####################### Step 1: Data Prep and Consolidation for VBGM ################################
#####################################################################################################
# This file takes length-at-age data from state fisheries-independent sampling programs and the GulfFin fishery-dependent sampling data and compiles them together in one format.

# Following Vanderkooy 2009 the biological birthdata is assinged to April 1. 

# Load Data
florida.data <- read.csv("Data/Data prep/CSV Data Files/Sheepshead_age_data_Florida.csv")
alabama.data <- read.csv("Data/Data prep/CSV Data Files/Sheepshead_age_data_Alabama.csv")
mississippi.data <- read.csv("Data/Data prep/CSV Data Files/Sheepshead_age_data_Mississippi.csv")
louisiana.data <- read.csv("Data/Data prep/CSV Data Files/Sheepshead_age_data_Louisiana.csv")
gulffin.data <- read.csv("Data/Data prep/CSV Data Files/Sheepshead_age_data_GulfFin.csv")
north.carolina.data <- read.csv("Data/Data prep/CSV Data Files/Sheepshead_age_data_North_Carolina.csv")


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
#florida.data <- florida.data[which(!is.na(florida.data$no_rings)),]
florida.data$date <- as.Date(florida.data$date, "%m/%d/%Y")
florida.data$birthdate <- as.Date(florida.data$birthdate, "%m/%d/%Y")
florida.data$sex <- as.character(florida.data$sex)
#florida.data$no_rings <- florida.data$no_rings -1
florida.data$fi_dependent <- "FI" # 1 indicates fishery dependent sampleing 0 indicates fishery-independent sampling
#florida.data$selectivity <- 0
florida.data <- florida.data[,which(colnames(florida.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


# Florida gulf-fin
florida.data.GSMFC <- gulffin.data[which(gulffin.data$State.Sampled=="FL"),]
colnames(florida.data.GSMFC) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code" , "Access.Site.Id")
florida.data.GSMFC$wet_weight_grams <- florida.data.GSMFC$wet_wt_kg * 1000
#florida.data.GSMFC <- florida.data.GSMFC[which(!is.na(florida.data.GSMFC$no_rings)),]
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
alabama.data$wet_weight_grams <- as.numeric(as.character((alabama.data$wet_wt_kg))) * 1000
alabama.data$no_rings <- as.character(alabama.data$no_rings)
alabama.data <- alabama.data[which(alabama.data$no_rings!="missing stones" ),]
alabama.data$no_rings[which(alabama.data$no_rings=="")] <- NA 
alabama.data$no_rings <- as.numeric(alabama.data$no_rings)
#alabama.data <- alabama.data[which(!is.na(alabama.data$no_rings)),]
alabama.data$date <- as.Date(paste(alabama.data$year, alabama.data$month, alabama.data$day, sep="-"), "%Y-%m-%d")
alabama.data$sex <- as.character(alabama.data$sex)
alabama.data$fi_dependent <- "FI" # 1 indicates fishery dependent sampleing 0 indicates fishery-independent sampling
#alabama.data$selectivity <- 0
alabama.data <- alabama.data[,which(colnames(alabama.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]

# Alabama GSMFC
alabama.data.GSMFC <- gulffin.data[which(gulffin.data$State.Sampled=="AL"),]
colnames(alabama.data.GSMFC) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code"  , "Access.Site.Id")
alabama.data.GSMFC$wet_weight_grams <- alabama.data.GSMFC$wet_wt_kg * 1000
#alabama.data.GSMFC <- alabama.data.GSMFC[which(!is.na(alabama.data.GSMFC$no_rings)),]
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
colnames(mississippi.data) <- c("date","month" , "spec_no" , "TL_mm", "FL_mm","SL_mm","wet_weight_grams","sex","gross_maturity","biological_age","margin_code","gonad_weight_grams","gutted_weight_grams","GSI","maturity", "source")
mississippi.data$state <- "MS"
#mississippi.data$age_1jan <- ifelse(mississippi.data$month > 3 , mississippi.data$age, mississippi.data$age + 1)
mississippi.data$gear <- NA
mississippi.data$common_name <- "Sheepshead"
mississippi.data$no_rings <- mississippi.data$biological_age
mississippi.data$no_rings[which(mississippi.data$no_rings==-1)] = NA
mississippi.data$margin_code[which(mississippi.data$margin_code==-1)] = NA
#mississippi.data <- mississippi.data[which(!is.na(mississippi.data$no_rings)),]

mississippi.data$date <- as.Date(mississippi.data$date, "%m/%d/%Y")
mississippi.data$sex <- ifelse(mississippi.data$sex == 1 , "M", "F")
mississippi.data$margin_code[which(mississippi.data$margin_code==3 & mississippi.data$source == "gcrl")] = 1# Convert margin codes
mississippi.data$margin_code[which(mississippi.data$margin_code==4 & mississippi.data$source == "gcrl")] = 2# Convert margin codes
mississippi.data$margin_code[which(mississippi.data$margin_code==5 & mississippi.data$source == "gcrl")] = 3# Convert margin codes
mississippi.data$margin_code[which(mississippi.data$margin_code==6 & mississippi.data$source == "gcrl")] = 4# Convert margin codes
mississippi.data$fi_dependent <- "FI"

#mississippi.data$selectivity <- 0
mississippi.data <- mississippi.data[,which(colnames(mississippi.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


### GSMFC MS
mississippi.data.GSMFC <- gulffin.data[which(gulffin.data$State.Sampled=="MS"),]
colnames(mississippi.data.GSMFC) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code"  , "Access.Site.Id")
mississippi.data.GSMFC$wet_weight_grams <- mississippi.data.GSMFC$wet_wt_kg * 1000
#mississippi.data.GSMFC <- mississippi.data.GSMFC[which(!is.na(mississippi.data.GSMFC$no_rings)),]
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
louisiana.data2 <- louisiana.data[which(louisiana.data$Section_Annuli_1!=louisiana.data$Section_Annuli_2),]
louisiana.data <- louisiana.data[which(louisiana.data$Section_Annuli_1==louisiana.data$Section_Annuli_2),]

par(mfrow = c(1,2))
hist(louisiana.data2$FL_mm, main = "Individuals Excluded", xlab = "Fork Length (mm)", xlim = c(200,700))
hist(louisiana.data$FL_mm, main = "Individuals Retained", xlab = "Fork Length (mm)", xlim = c(200,700))

louisiana.data$no_rings <- louisiana.data$Section_Annuli_1
louisiana.data$margin_code <- (louisiana.data$Edge_S1 + louisiana.data$Edge_S2)/2
#louisiana.data <- louisiana.data[which(!is.na(louisiana.data$no_rings)),]
louisiana.data$state <- "LA"
louisiana.data$date <- as.Date(paste(louisiana.data$year, louisiana.data$month, louisiana.data$day, sep="-"), "%Y-%m-%d")
louisiana.data$sex <- ifelse(louisiana.data$sex == "3 : Female", "F" ,  ifelse(louisiana.data$sex == "2 : Male", "M", "U"))
louisiana.data$fi_dependent <- toupper(louisiana.data$Fishery)
#louisiana.data$selectivity <- 0
#louisiana.data$selectivity[which(louisiana.data$Fishery == "Commercial")] <- 254
louisiana.data <- louisiana.data[,which(colnames(louisiana.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


# Texas
texas.data <- gulffin.data[which(gulffin.data$State.Sampled=="TX"),]
colnames(texas.data) <- c("year" , "Species.Itis" , "common_name" , "Sampler.Id" , "month" , "day" , "state" , "County.Sampled" , "State.Landed" , "County.Landed" , "Data.Source" , "fishery" , "Fishing.Mode.Code" , "Area.Name" , "Area.Code" , "Sub.Area.Name" , "Sub.Area.Code" , "Depth" , "Gear.Code" , "gear" , "Duplicate.Length" , "FL_mm" , "TL_mm" , "wet_wt_kg" , "sex" , "Age.Tag.Number" , "Lab.Catalog.Number" , "Age.Tissue.Type" , "Structure.Disp.Code" , "Type.Of.Processing" , "no_rings" , "margin_code" , "Access.Site.Id")
texas.data$wet_weight_grams <- texas.data$wet_wt_kg * 1000
#texas.data <- texas.data[which(!is.na(texas.data$no_rings)),]
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
virginia_dir <- "Data/Raw/Virginia/CSV_files"
virginia.data <- c()
for (j in 1:(length(dir(virginia_dir)))) {
  dat <- read.csv(paste0(virginia_dir,"/", dir(virginia_dir)[j]))
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
#virginia.data <- virginia.data[which(!is.na(virginia.data$no_rings)),]
virginia.data$date <- as.Date(paste(virginia.data$year, virginia.data$month, virginia.data$day, sep="-"), "%Y-%m-%d")
virginia.data$sex <- ifelse(virginia.data$sex == 2, "F" ,  ifelse(virginia.data$sex == 1, "M", "U"))
virginia.data$fi_dependent <- ifelse(virginia.data$type == "FI", "FI", ifelse(virginia.data$type == "C", "COMMERCIAL", "RECREATIONAL"))
virginia.data <- virginia.data[,which(colnames(virginia.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]


# SC
sc.data <- read.csv("Data/Raw/South Carolina/SCDNR_Sheepshead.csv")
sc.data <- sc.data[which(sc.data$AgeMethod=="otolith"),]
sc.data$TL_mm <- as.numeric(as.character(sc.data$TL_mm))
sc.data$FishWeight_g <- as.numeric(as.character(sc.data$FishWeight_g))
sc.data$no_rings <- as.numeric(as.character(sc.data$Age_Yrs))

colnames(sc.data) <- tolower(colnames(sc.data))
colnames(sc.data)[which(colnames(sc.data)=="tl_mm")] <- "TL_mm"
colnames(sc.data)[which(colnames(sc.data)=="fishweight_g")] <- "wet_weight_grams"
sc.data$state <- "SC"
sc.data$common_name <- "Sheepshead"
#sc.data <- sc.data[which(!is.na(sc.data$no_rings)),]
sc.data$date <- as.Date(sc.data$datefull, "%d-%b-%y")
sc.data$sex <- ifelse(sc.data$sex == 2, "F" ,  ifelse(sc.data$sex == 1, "M", "U"))
sc.data$gear_description <- gsub(" ","",sc.data$gear_description)
sc.data$gear_description <- gsub('"',"",sc.data$gear_description)
sc.data$fi_dependent <- ifelse(sc.data$gear_description %in% c("GILLNET-MESHG-3STRETCH", "TRAMMELNET-200YDX8FT-14&2.5INSTRMESH"), "FI", "RECREATIONAL")
sc.data <- sc.data[,which(colnames(sc.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month"))]



## North Carolina
north.carolina.data$fi_dependent <- ifelse(north.carolina.data$FLEET == "independent ", "FI", ifelse(north.carolina.data$FLEET == "commercial", "COMMERCIAL", "RECREATIONAL"))
north.carolina.data$FL_mm <- as.numeric(as.character(north.carolina.data$LENGTH..FL..mm.))
north.carolina.data$TL_mm <- as.numeric(as.character(north.carolina.data$LENGTH..TL..mm.))
north.carolina.data$sex <- ifelse(north.carolina.data$SEX == "female", "F" ,  ifelse(north.carolina.data$SEX == "male", "M", "U"))
north.carolina.data$no_rings <- north.carolina.data$AGE
north.carolina.data$wet_weight_grams <- north.carolina.data$WEIGHT..kg.* 1000
north.carolina.data$date <- as.Date(paste(north.carolina.data$YEAR, north.carolina.data$MONTH, north.carolina.data$DAY, sep="-"), "%y-%m-%d")
north.carolina.data$month <- north.carolina.data$MONTH 
north.carolina.data$state <- "NC"
north.carolina.data <- north.carolina.data[,which(colnames(north.carolina.data)%in%c("margin_code","no_rings","state","wet_weight_grams","TL_mm","FL_mm","sex","fi_dependent","SL_mm", "date","month","biological_age_day"))]




##########################
# Biological age asignment for the GULF
##########################
library(lubridate)

gulf.data <- rbind.fill(florida.data, florida.data.GSMFC,alabama.data,alabama.data.GSMFC,mississippi.data,mississippi.data.GSMFC, louisiana.data,texas.data)
#gulf.data <- gulf.data[which(!is.na(gulf.data$margin_code)),]
gulf.data$no_rings[which(gulf.data$no_rings>138)] <- 14

# Calculate the biological age (a.k.a., calendar age), derived to assign fish to the correct year class assuming a January 1st birthdate, used information on the relative stage of ring formation on the outer otolith margin and date of capture.  
gulf.data$biological_age <- ifelse(yday(gulf.data$date) < 91 & gulf.data$margin_code >= 3 ,gulf.data$no_rings + 1, gulf.data$no_rings)
gulf.data$fractional_age <- gulf.data$biological_age + (yday(gulf.data$date) - 91)/365


##########################
# Biological age asignment for the Atlantic
##########################
atlantic.data <- rbind.fill(sc.data, north.carolina.data, virginia.data)
atlantic.data$biological_age <-  atlantic.data$no_rings 
atlantic.data$fractional_age <-  atlantic.data$biological_age + (yday(atlantic.data$date) - 121)/365


#######################
# Compile all the data!
#######################
library(plyr)
gulf.data <- rbind.fill(gulf.data, atlantic.data)


########################
# SL, FL, TL Regressions
########################
library(lme4)
library(nlme)

fl.tl <- gulf.data[which(!is.na(gulf.data$FL_mm) & !is.na(gulf.data$TL_mm)),]
fl.tl.regression.hlm <- lme(FL_mm ~ TL_mm, random = ~TL_mm|state, data = fl.tl) # Don't need
fl.tl.regression <- lm(FL_mm ~ TL_mm, data = fl.tl)
fl.tl.intercept <- fl.tl.regression$coefficients[1]
fl.tl.slope <- fl.tl.regression$coefficients[2]
gulf.data$FL_mm[which(is.na(gulf.data$FL_mm))] <- fl.tl.intercept + fl.tl.slope * gulf.data$TL_mm[which(is.na(gulf.data$FL_mm))]

fl.sl <- gulf.data[which(!is.na(gulf.data$FL_mm) & !is.na(gulf.data$SL_mm)),]
fl.sl.regression.hlm <- lme(FL_mm ~ SL_mm, random = ~SL_mm|state, data = fl.sl)
fl.sl.regression <- lm(FL_mm ~ SL_mm, data = fl.sl)
fl.sl.intercept <- fl.sl.regression$coefficients[1]
fl.sl.slope <- fl.sl.regression$coefficients[2]
gulf.data$FL_mm[which(is.na(gulf.data$FL_mm))] <- fl.sl.intercept + fl.sl.slope * gulf.data$SL_mm[which(is.na(gulf.data$FL_mm))]



#########################
# Save
#########################
dat <- gulf.data

state.id <- data.frame(state = c("VA_ocean","VA_bay","NC","SC", "GA","FL_atlantic","FL_gulf","AL","MS","LA","TX"), state_no = c(1:11))
dat <- merge(dat, state.id , by = "state")

dat <- dat[which(!is.na(dat$state_no)),]

# Add spatial data
library(rgdal)
library(sp)
management_centroids <- rgdal::readOGR("Data/State boundaries", "coastal_management_centroids_w_ga")
management_centroids <- as.data.frame(management_centroids)[,c(16,15,7)]
management_centroids <- management_centroids[order(management_centroids$NAME),]
management_centroids$NAME <- sort(c("VA_ocean","VA_bay","NC","SC", "GA","FL_atlantic","FL_gulf","AL","MS","LA","TX"))
colnames(management_centroids) <- c("Lat","Lon","state")
dat <- merge(dat, management_centroids, by = "state")


save(dat, file = "Data/Sheepshead_data_GOM.RData" )
write.csv(dat, file = "Data/Sheepshead_VA_to_TX_biological_data.csv")
