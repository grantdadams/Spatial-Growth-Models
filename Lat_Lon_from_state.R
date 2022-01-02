# Load packages
library(rgeos)
library(fossil)
library(sp)
library(maptools)
library(rgdal)
library(spatstat)

# Load data
state_bounds <- rgdal::readOGR("State Boundaries/tl_2017_us_state", "tl_2017_us_state")
management_areas <- rgdal::readOGR("State Boundaries", "coastal_management_areas")


plot(state_bounds, xlim = c(-95,-80), ylim = c(20,40))
plot(management_areas, add = T, col = 2)
points(management_centroids, col = 3, pch = 2, cex =1)


