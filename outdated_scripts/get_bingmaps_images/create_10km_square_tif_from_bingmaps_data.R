######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = create_10km_square_tif_from_bingmaps, title = Create a 10km² georeferenced TIF + an OGC geopackage of a given AOI on Earth using Bing maps satellite imagery, abstract = Create a 10km² georeferenced TIF + a OGC Geopackage of a given AOI on Earth using Bing maps satellite imagery at best available zoom (19, which corresponds to a resolution of approx. 0.5 cm);
# wps.in: id = kml_path, type = string, title = Path to the kml covering the area of study., value = "/home/ptaconet/REACT_BF.kml";
# wps.in: id = apiKey, type = string, title = Bing API key. You need an account on the bing maps dev center to get a key. See https://www.bingmapsportal.com/ . To get the key: My account > my keys  and then copy and paste the key in a txt file , value = scan("/home/ptaconet/react/r_bingmaps/bingAPIkey.txt",what="");
# wps.in: id = imagerySet, type = string, title =  The type of imagery., value = " Aerial|AerialWithLabels|AerialWithLabelsOnDemand|CanvasDark|CanvasLight|CanvasGray|Road";
# wps.in: id = destFolder, type = string, title = Destination folder (to store the data downloaded and produced). Must be writable. , value = "/home/ptaconet/react/r_bingmaps";
# wps.out: id = , type = image, title = Georeferenced images in TIF format + OGC geopackage database. 1 image per 10km * 10km square., value = "";


# Local variable (to set by the user)
#kml_path="/home/ptaconet/Téléchargements/kml_bobo.kml"
lat_min<-11.06246
lon_min<--4.384058
lat_max<-11.237
lon_max<--4.249606

#apiKey = scan("/home/ptaconet/Documents/react/r_bingmaps/bingAPIkey.txt",what="")
apiKey="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
imagerySet="Aerial"
destFolder="/home/ptaconet/Documents/react/r_bingmaps"

# Global variable (for the time being, can not be changed. Could be set to local in future dev)
mapsize_width=1500  #max: 2000
mapsize_height=1500 # max: 1500
zoomLevel=19
cell_size=10000 # width/height of the big tifs to create (unit: meters)

# Install / call useful libraries
if(!require(sf)){
  install.packages("sf")
}
if(!require(raster)){
  install.packages("raster")
}
if(!require(rjson)){
  install.packages("rjson")
}
if(!require(rgdal)){
  install.packages("rgdal")
}
if(!require(gdalUtils)){
  install.packages("gdalUtils")
}
if(!require(raster)){
  install.packages("raster")
}

require(sf)
require(raster)
require(rjson)
require(rgdal)
require(gdalUtils)
require(raster)

# Source useful functions
source("https://raw.githubusercontent.com/ptaconet/r_react/master/outdated_scripts/get_bingmaps_images/download_and_georeference_bingmaps_data.R")
source("https://raw.githubusercontent.com/ptaconet/r_react/master/outdated_scripts/mosaic_tif_images.R")

## Get UTM WGS84 Zone number of the ROI. from https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
cat("Warning: ROIs overlapping more than 1 UTM zone are currently not adapted in this workflow\n")
utm_zone_number<-(floor((lon_min + 180)/6) %% 60) + 1
if(lat_max>0){ # if latitudes are North
  epsg_utm<-as.numeric(paste0("326",utm_zone_number))
} else { # if latitude are South
  epsg_utm<-as.numeric(paste0("325",utm_zone_number))
}

# Convert area of interest to sf object
df <- data.frame(lon=c(lon_min,lon_max,lon_max,lon_min,lon_min), 
                 lat=c(lat_min,lat_min,lat_max,lat_max,lat_min))

poly <- st_sf(st_sfc(st_polygon(list(as.matrix(df)))), crs = 4326)



# Convert to EPSG 32630 (to get the units in meters instead of degrees). EPSG 32630 is the one covering our area of study (Burkina Faso and north Cote d'Ivoire)
poly <- st_transform(poly,crs=epsg_utm)

# Make grid over the bounding box. We want 10 km square grids, so as to be able to use them on the tablet. 
grid_10km=st_make_grid(poly,what="polygons",cellsize = cell_size)

# We convert the grid to geopackage and we save it
st_write(st_transform(grid_10km,crs=4326),file.path(destFolder,"raster_10km.gpkg"))

# Loop on each 10km square tile
for (i in 2:length(grid_10km)){

  
cat(paste0("Creating 10 square kilometers grid n° ",i, " over ",length(grid_10km)))
  
# We can download 1500 * 1500 pixels image from bing maps servers. We know that at zoom 19, 1 pixel=0.30m, so 1500 px = 450m. So within each 10km grid we make 400m grids (to be sure that all the area is encompassed)
grid_400m=st_make_grid(grid_10km[i],what="centers",cellsize = 400)

grid_400m=st_transform(grid_400m,crs=4326)

dir.create(file.path(destFolder, i))
# Download all the tiles from bing maps servers and convert them as tif
for (j in 1:length(grid_400m)){
  cat(paste0("downloading the data n° ",j, " over ",length(grid_400m)))
  fileName=paste0(i,"_10km_",j)
  map_this_tile<-download_and_georeference_bingmaps_data(apiKey=apiKey,
                                                    center_x=st_coordinates(grid_400m[j])[1,1],
                                                    center_y=st_coordinates(grid_400m[j])[1,2],
                                                    mapsize_width=mapsize_width,
                                                    mapsize_height=mapsize_width,
                                                    imagerySet=imagerySet,
                                                    zoomLevel=zoomLevel,
                                                    destFolder=destFolder,
                                                    fileName=paste0(i,"_10km_",j)
  )
  file.copy(paste0(destFolder,"/",fileName,'.tif'),paste0(destFolder,"/",i,"/",fileName,'.tif'))
  file.remove(paste0(destFolder,"/",fileName,'.tif'))
}

# make a big tif out of all the small tifs
cat("Creating the 10km x 10km tif")
all_my_rasts <- list.files(path=file.path(destFolder,i), pattern =paste0(i,"_10km"), full.names=TRUE)
all_my_rasts=as.vector(all_my_rasts)

mosaic_tif_images(all_my_rasts,paste0(destFolder,"/",i,"_10km.tif"))

## Convert the tif created to a raster in a OGC Geopackage and build pyramids at various zoom levels .
cat("Converting the tif to a geopackage file")
system(paste0("gdal_translate -ot Byte -of GPKG ",destFolder,"/",i,"_10km.tif ",destFolder,"/",i,"_10km.gpkg -co APPEND_SUBDATASET=YES -co RASTER_TABLE=",i,"_10km"))
cat("Building pyramids for quick rendering on GIS")
system(paste0("gdaladdo --config OGR_SQLITE_SYNCHRONOUS OFF -r AVERAGE ",destFolder,"/",i,"_10km.gpkg 2 4 8 16 32 64 128 256"))

}


