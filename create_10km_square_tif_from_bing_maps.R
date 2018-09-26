# wps.des: id = create_10km_square_tif_from_bingmaps, title = Create a 10km square georeferenced .tif and an OGC geopackage using bing maps satellite imagery, abstract = Create a 10km square georeferenced tif and geopackage using bing maps satellite imagery at best available zoom (19, which corresponds to a resolution of approx. 0.5 cm);
# wps.in: id = kml_path, type = string, title = Path to the kml covering the area of study., value = "/home/ptaconet/REACT_BF.kml";
# wps.in: id = apiKey, type = string, title = Bing API key. You need an account on the bing maps dev center to get a key. See https://www.bingmapsportal.com/ . To get the key: My account > my keys  and then copy and paste the key in a txt file , value = scan("/home/ptaconet/react/r_bingmaps/bingAPIkey.txt",what="");
# wps.in: id = imagerySet, type = string, title =  The type of imagery., value = " Aerial|AerialWithLabels|AerialWithLabelsOnDemand|CanvasDark|CanvasLight|CanvasGray|Road";
# wps.in: id = destFolder, type = string, title = Destination folder (to store the data downloaded and produced). Must be writable. , value = "/home/ptaconet/react/r_bingmaps";
# wps.out: id = , type = image, title = Georeferenced images in TIF format + OGC geopackage database. 1 image per 10km * 10km square., value = "";


# Local variable (to set by the user)
kml_path="/home/ptaconet/REACT_BF.kml"
apiKey = scan("/home/ptaconet/react/r_bingmaps/bingAPIkey.txt",what="")
imagerySet="Aerial"
destFolder="/home/ptaconet/react/r_bingmaps"

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

require(sf)
require(raster)
require(rjson)
require(rgdal)
require(gdalUtils)


# Source useful function
source("https://raw.githubusercontent.com/ptaconet/r_react/master/download_and_georeference_bingmaps_data.R")

# read kml as sf object
kml_sf <- st_read(kml_path)

# Convert to EPSG 32630 (to get the units in meters instead of degrees). EPSG 32630 is the one covering our area of study (Burkina Faso and north Cote d'Ivoire)
kml_sf <- st_transform(kml_sf,crs=32630)

# Make grid over the bounding box. We want 10 km square grids, so as to be able to use them on the tablet. 
grid_10km=st_make_grid(kml_sf,what="polygons",cellsize = cell_size)

# We convert the grid to kml and we save it
# st_write(st_transform(grid_10km,crs=4326),"/home/ptaconet/grid_10km.kml")

# Loop on each 10km square tile
for (i in 1:length(grid_10km)){

cat(paste0("Creating 10 square kilometers grid n° ",i, " over ",length(grid_10km)))
  
# We can download 1500 * 1500 pixels image from bing maps servers. We know that at zoom 19, 1 pixel=0.30m, so 1500 px = 450m. So within each 10km grid we make 400m grids (to be sure that all the area is encompassed)
grid_400m=st_make_grid(grid_10km[i],what="centers",cellsize = 400)

grid_400m=st_transform(grid_400m,crs=4326)

# Download all the tiles from bing maps servers and convert them as tif
for (j in 1:length(grid_400m)){
  cat(paste0("downloading the data n° ",j, " over ",length(grid_400m)))
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
  
}

# make a big tif out of all the small tifs
cat("Creating the 10km x 10km tif")
all_my_rasts <- list.files(path=destFolder, pattern =paste0(i,"_10km"), full.names=TRUE)
all_my_rasts=as.vector(all_my_rasts)

e<-extent(as.numeric(st_bbox(st_transform(grid_10km,crs=4326))$xmin),as.numeric(st_bbox(st_transform(grid_10km,crs=4326))$xmax),as.numeric(st_bbox(st_transform(grid_10km,crs=4326))$ymin),as.numeric(st_bbox(st_transform(grid_10km,crs=4326))$ymax))
template <- raster(e)
projection(template) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
writeRaster(template, file=paste0(destFolder,"/",i,"_10km.tif"), format="GTiff",overwrite=TRUE)
mosaic_rasters(gdalfile=all_my_rasts,dst_dataset=paste0(destFolder,"/",i,"_10km.tif"),of="GTiff")

## Convert the tif created to a raster in a geopackage.
cat("Converting the tif to a geopackage file")
system(paste0("gdal_translate -ot Byte -of GPKG ",destFolder,"/",i,"_10km.tif ",destFolder,"/",i,"_10km.gpkg -co APPEND_SUBDATASET=YES -co RASTER_TABLE=",i,"_10km"))
system(paste0("gdaladdo --config OGR_SQLITE_SYNCHRONOUS OFF -r AVERAGE ",destFolder,"/",i,"_10km.gpkg 2 4 8 16 32 64 128 256"))

}


