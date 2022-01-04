library(tidyr)

source("Data/Data prep/VBGM_Data_prep.R")

# Sample size
summ_stat <- data.frame(table(dat$state_no, dat$sex))
summ_stat <- tidyr::spread(summ_stat, Var2, Freq)
colnames(summ_stat) <- c("state_no", "n_f", "n_m")

# Minimum size
min_size <- aggregate(FL_mm ~ state_no + sex, dat, function(x) min(x))
min_size <- tidyr::spread(min_size, sex, FL_mm)
colnames(min_size) <- c("state_no", "min_f", "min_m")

summ_vbgf <- merge(summ_stat, min_size, by = "state_no" )

# Add state names and reorder
state.id <- data.frame(state = c("VA_ocean","VA_bay","NC","SC", "GA","FL_atlantic","FL_gulf","AL","MS","LA","TX"), state_no = c(1:11))
summ_vbgf <- merge(summ_vbgf, state.id, by = "state_no" )
summ_vbgf <- summ_vbgf[order(summ_vbgf$state_no),]

# save
write.csv(summ_vbgf, file = "Table_3_VBGF_Summary.csv")
