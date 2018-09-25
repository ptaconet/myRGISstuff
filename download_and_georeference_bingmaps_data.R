# wps.des: id = download_and_georeference_bingmaps_data, title = Download and georeference Bing maps data, abstract = Download and georeference Bing maps data using the Bing API. More info here: https://msdn.microsoft.com/en-us/library/ff701724.aspx;
# wps.in: id = apiKey, type = string, title = Bing API key. You need an account on the bing maps dev center to get a key. See https://www.bingmapsportal.com/ . To get the key: My account > my keys  and then copy and paste the key in a txt file , value = scan("/home/ptaconet/react/r_bingmaps/bingAPIkey.txt",what="");
# wps.in: id = center_x, type = numeric, title = Coordinate of the longitude where the map is centered , value = -3.4534;
# wps.in: id = center_y, type = numeric, title = Coordinate of the latitude where the map is centered , value = 10.9203;
# wps.in: id = mapsize_width, type = numeric, title = The width in pixels of the static map output. Max is 2000 , value = 2000;
# wps.in: id = mapsize_height, type = numeric, title = The height in pixels of the static map output. Max is 1500 , value = 1500;
# wps.in: id = imagerySet, type = string, title =  The type of imagery., value = " Aerial|AerialWithLabels|AerialWithLabelsOnDemand|CanvasDark|CanvasLight|CanvasGray|Road";
# wps.in: id = zoomLevel, type = numeric, title = The level of zoom to display. An integer between 0 and 21. Note: Some imagery may not be available at all zoom levels for all locations, value = 19;
# wps.in: id = destFolder, type = string, title = Destination folder (to store the image)., value = "/home/ptaconet/react/r_bingmaps/";
# wps.out: id = , type = image, title = Georeferenced image in TIF format, value = "";

#apiKey = scan("/home/ptaconet/react/r_bingmaps/bingAPIkey.txt",what="")
#center_x=6.7098
#center_y=45.8739
#mapsize_width=2000  #max: 2000
#mapsize_height=1500 # max: 1500
#imagerySet="Aerial"
#zoomLevel=19
#destFolder="/home/ptaconet/react/r_bingmaps/"
#fileName="r_bingmaps"
  
download_and_georeference_bingmaps_data<-function(apiKey,center_x,center_y,mapsize_width,mapsize_height,imagerySet,zoomLevel,destFolder,fileName){
## Load useful libraries
library(raster)
library(rjson)
library(rgdal)

## First download image and metadata 
# Download image
y = paste0("http://dev.virtualearth.net/REST/v1/Imagery/Map/",imagerySet,"/",center_y,",",center_x,"/",zoomLevel,"?mapSize=",mapsize_width,",",mapsize_height,"&format=jpeg&key=",apiKey)
download.file(y,paste0(destFolder,fileName,'.jpg'), mode = 'wb')

# Download metadata of the image (including bouding box: add mmd=1  see https://msdn.microsoft.com/en-us/library/ff701724.aspx for more info)
y = paste0("http://dev.virtualearth.net/REST/v1/Imagery/Map/",imagerySet,"/",center_y,",",center_x,"/",zoomLevel,"?mapSize=",mapsize_width,",",mapsize_height,"&mmd=1&format=png&key=",apiKey)
download.file(y,paste0(destFolder,fileName,'.json'))

## Then georeference and save the image as tif
# open image
map <- stack(paste0(destFolder,fileName,'.jpg')  )

# open metadata
result <- fromJSON(file = paste0(destFolder,fileName,".json"))

# To get the bounding box: result$resourceSets[[1]]$resources[[1]]$bbox
xmin(map) <- result$resourceSets[[1]]$resources[[1]]$bbox[2]
xmax(map) <- result$resourceSets[[1]]$resources[[1]]$bbox[4]
ymin(map) <- result$resourceSets[[1]]$resources[[1]]$bbox[1]
ymax(map) <- result$resourceSets[[1]]$resources[[1]]$bbox[3]
crs(map) <- "+proj=longlat +datum=WGS84"

# Aggregate, as writeRaster disaggregates (why? do not know...)
map=aggregate(map,fact=2,fun=mean,expand=TRUE)
writeRaster(map, paste0(destFolder,fileName,".tiff"), "GTiff",overwrite=TRUE,options="COMPRESS=LZW",datatype='INT2S')
}