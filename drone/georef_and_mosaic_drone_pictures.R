rm(list = ls())


## input parameters to set up
path_to_folder_source_pic<-"/home/ptaconet/Téléchargements/ANR_CoheSIoN_Bke_Photos_Drone"  # folder containing the drone pictures to mosaic
path_to_folder_output_pic<-"/home/ptaconet/Téléchargements/output_ANR_CoheSIoN_Bke_Photos_Drone" # output folder where the georeferenced pictures and the mosaic will be stored. The folder will be created if is does not exist. The mosaic will be created and named "mosaic.tif"
ground_sample_distance<-1.7 # (in cm / pixel)

############ Start #####################
if (!require(exifr)){
  install.packages("exifr")
}
if (!require(sf)){
  install.packages("sf")
}
if (!require(raster)){
  install.packages("raster")
}
if (!require(gdalUtils)){
  install.packages("gdalUtils")
}

require(exifr)
require(sf)
require(raster)
require(gdalUtils)

# source function to mosaic
source("https://raw.githubusercontent.com/ptaconet/r_react/master/outdated_scripts/mosaic_tif_images.R")
files <- list.files(path_to_folder_source_pic,pattern = "*.JPG",full.names = T)
dat <- read_exif(files)

# Create output folder if it does not exist
if (!dir.exists(path_to_folder_output_pic)){
  dir.create(path_to_folder_output_pic)
}

# formule pour avoir la ground_sample_distance en fonction des métadonnées présentes dans le fichier EXIF : (sensor height (mm) x flight height (m) x 100) / (focal lenght (mm) x image height (pixel))
# from https://support.pix4d.com/hc/en-us/articles/202560249-TOOLS-GSD-calculator
#sens_width<-2*unique(dat$FocalLength)*tan(0.5*unique(dat$FOV)/57.296)  # from  https://vfxcamdb.com/film-back-calculator/
#dat$resolution<-(sens_width*as.numeric(dat$RelativeAltitude)*100)/(dat$FocalLength*dat$ImageWidth) # en cm


dat$resolution<-ground_sample_distance

dat_sf<-st_as_sf(dat, coords = c('GPSLongitude', 'GPSLatitude'), crs = 4326)

utm_zone_number<-(floor((mean(dat$GPSLongitude) + 180)/6) %% 60) + 1
if(mean(dat$GPSLatitude)>0){ # if latitudes are North
  epsg<-as.numeric(paste0("326",utm_zone_number))
} else { # if latitude are South
  epsg<-as.numeric(paste0("325",utm_zone_number))
}

dat_sf_utm<-st_transform(dat_sf,epsg)

# Set bounding box the picturec
dat_sf_utm$x_min<- st_coordinates(dat_sf_utm)[,1]-floor(dat_sf_utm$ImageHeight/2)*dat_sf_utm$resolution/100
dat_sf_utm$x_max<- st_coordinates(dat_sf_utm)[,1]+floor(dat_sf_utm$ImageHeight/2)*dat_sf_utm$resolution/100
dat_sf_utm$y_min<- st_coordinates(dat_sf_utm)[,2]-floor(dat_sf_utm$ImageWidth/2)*dat_sf_utm$resolution/100
dat_sf_utm$y_max<- st_coordinates(dat_sf_utm)[,2]+floor(dat_sf_utm$ImageWidth/2)*dat_sf_utm$resolution/100

# georeference each picture
for (i in 1:nrow(dat_sf_utm)){
  cat(paste0("\nprocessing image",dat_sf_utm$SourceFile[i]))
  map <- stack(dat_sf_utm$SourceFile[i])

  map <- t(map)
  map <- flip(map,'x')

  xmin(map) <- dat_sf_utm$x_min[i]
  xmax(map) <- dat_sf_utm$x_max[i]
  ymin(map) <- dat_sf_utm$y_min[i]
  ymax(map) <- dat_sf_utm$y_max[i]
  crs(map) <- CRS(paste0("+init=epsg:",epsg))

  writeRaster(map, file.path(path_to_folder_output_pic,gsub('JPG','TIF',dat_sf_utm$FileName[i])), "GTiff",overwrite=TRUE,options="COMPRESS=LZW",datatype='INT2S')
}

# Mosaic the georef pictures
mosaic_tif_images(list.files(path_to_folder_output_pic,full.names = T),file.path(path_to_folder_output_pic,"mosaic.tif"))

cat('\nEnd')
