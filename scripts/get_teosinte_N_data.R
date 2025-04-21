library(dplyr)
donor_passport <-  read.csv("~/Desktop/Bzea_metadata.csv")
NIL_donor <-  read.csv("~/Desktop/NIL_donor.csv")

donor_data <- donor_passport %>% 
  select(donor_id,longitude,latitude,elevation) 

colnames(donor_passport)


gis_dir <- "/Volumes/rsstu/users/r/rrellan/sara/gisdata/soilgrids"

raster_file <- file.path(gis_dir,"nitrogen_0-5cm_mean_1000_lonlat.tif")
env_raster <- raster::raster(raster_file)
donor_data$`nitrogen_0-5cm_mean_1000` <- raster::extract(env_raster,donor_data[,c("longitude","latitude")])

raster_file <- "~/Desktop/NHx-deposition1993.tif"
env_raster <- raster::raster(raster_file)
donor_data$`NHx_deposition1993` <- raster::extract(env_raster,donor_data[,c("longitude","latitude")])

sum(is.na(donor_data$`nitrogen_0-5cm_mean_1000`))
write.csv(donor_data, file = "~/Desktop/donor_data.csv", quote=FALSE, na='',row.names = FALSE)

