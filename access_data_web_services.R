## Given a couple of input parameters {date, position of catch (lat and lon)}, this script prepares the data to use as input of the predictive modelisation.

#### Workflow steps :
### Step x - MODIS LST
  ## x.1 - Download MODIS LST product
  ## x.2 - Prepare the data
  ## x.3 - Extract the data

# - download the datasets : MODIS LST, MODIS NDVI, GPM, Wind ERA
# - pre-process and prepare the datasets: 
  # - for MODIS LST
  # - for MODIS NDVI
# Extract the data (i.e. compute zonal statistics) on a buffer around the villages


########## Specific for ERA INTERIM (instructions to install the clients and use the python environment, more info here : https://dominicroye.github.io/en/2018/access-to-climate-reanalysis-data-from-r/#era-interim)
##Set virtual env and install the python ECMWF API
#reticulate::py_install("ecmwf-api-client") #(from https://community.rstudio.com/t/problem-installing-python-libraries-with-reticulate-py-install-error-pip-not-found/26561/2)
#system("virtualenv -p /usr/bin/python2 /home/ptaconet/.virtualenvs/py2-virtualenv")
##Also works with this: virtualenv_create("py3-virtualenv", python = "/usr/bin/python3")
#reticulate::use_virtualenv("py2-virtualenv")
##install the python ECMWF API
#reticulate::py_install("ecmwf-api-client", envname = "py2-virtualenv")
## For ERA-5 :
#system("pip install cdsapi")

#### Start Workflow
########################################################################################################################
############ Set Input parameters for the workflow ############
########################################################################################################################

path_to_processing_folder<-"/home/ptaconet/Documents/react/data_CIV"  #<Path to the processing folder (i.e. where all the data produced by the workflow will be stored)>
path_to_roi_vector="ROI.kml" #<Path to the Region of interest in KML format>
path_to_grassApplications_folder<-"/usr/lib/grass74" #<Can be retrieved with grass74 --config path . More info on the use of rgrass7 at https://grasswiki.osgeo.org/wiki/R_statistics/rgrass7
# If we extract the dataset of dates and locations of HLC from a database, provide path to DB and QSL query to execute :
path_to_database<-"/home/ptaconet/Bureau/react.db"  
path_to_sql_query_dates_loc_hlc<-"https://raw.githubusercontent.com/ptaconet/r_react/master/db_sql/hlc_dates_location.sql"
path_to_sql_query_households_loc_pop<-"https://raw.githubusercontent.com/ptaconet/r_react/master/db_sql/loc_pop_households.sql"
# Else if we extract the dataset of dates and locations of HLC from a csv file, provide the path to the csv file :
path_to_csv_dates_loc<-NULL
path_to_csv_population<-NULL
# Credential to the NASA servers (EarthData)
username_EarthData<-"ptaconet"  #<EarthData username>
password_EarthData<-"HHKcue51"  #<EarthData password>

# URLs or paths to static data
url_to_srtm_datasets<-c("http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/N09W006.SRTMGL1.hgt.zip","http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/N08W006.SRTMGL1.hgt.zip") # for BF : http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/N10W004.SRTMGL1.hgt.zip   # Check tiles at : https://dwtkns.com/srtm30m/
path_to_texture_inertia<-"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_4.TIF"  # this file was computed from the Spot 6 image using Orfeo Toolbox (application HaralickTextureExtraction). The inertia texture was computed on a 5 x 5 pixel moving window
threshold_built_areas<-3 # for CIV : 3 ; for BF : 2.2
path_to_landcover_dataset<-"Classification/classification_L3.gpkg"


# Buffer size, within which the raster statistics will be computed (radius in meters)
buffer_size_meters=2000 
## Number of lag days for each source (ie number of days separating the human catch landing date and the first date of interest for the given source)
lag_days_modis_lst<-40
lag_days_modis_veget<-40
lag_days_modis_evapo<-40
lag_days_gpm<-40



########################################################################################################################
############ Prepare workflow ############
########################################################################################################################

## Call useful libraries
library(raster)
library(sf)
library(sp)
library(gdalUtils)
library(rgdal)
require(rgeos)
require(httr)
require(RSQLite)
require(dplyr)
require(ncdf4)
require(reticulate)
require(lubridate)
library(stringr)
require(rgrass7)
require(landscapemetrics)
require(landscapetools)
#library(getSpatialData)
#library(MODIS)
#library(MODISTools)
#require(gapfill)
## Load packages for working on multi-core
#library(parallel)
#library(doParallel)
#library(foreach)

## Urls to the various OpenDAP servers
url_modis_opendap<-"https://opendap.cr.usgs.gov/opendap/hyrax"
url_gpm_opendap<-"https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.06"
url_imcce_webservice<-"http://vo.imcce.fr/webservices/miriade/ephemcc_query.php?"
url_noaa_nighttime_webservice<-"https://gis.ngdc.noaa.gov/arcgis/rest/services/NPP_VIIRS_DNB/Monthly_AvgRadiance/ImageServer/exportImage"
modis_lst_terra_collection<-"MOD11A1.006"
modis_lst_aqua_collection<-"MYD11A1.006"
modis_veget_terra_collection<-"MOD13Q1.006"
modis_veget_aqua_collection<-"MYD13Q1.006"
modis_evapo_terra_collection<-"MOD16A2.006"
modis_evapo_aqua_collection<-"MYD16A2.006"
modis_crs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

## Connection to the EarthData servers
httr::set_config(authenticate(user=username_EarthData, password=password_EarthData, type = "basic"))

## Parameters for download of ERA 5 data
#import python CDS-API
cdsapi <- reticulate::import('cdsapi')
#for this step there must exist the file .cdsapirc
server = cdsapi$Client() #start the connection

## User-defined functions
# To read the files of OpenDap indexes
fun_get_opendap_index<-function(path_to_opendap_index){
  opendap_indexes<-read.csv(path_to_opendap_index,skip = 1)
  opendap_indexes[1]<-NULL
  opendap_indexes<-colnames(opendap_indexes)
  opendap_indexes<-gsub("X\\.","-",opendap_indexes)
  opendap_indexes<-gsub("X","",opendap_indexes)
  opendap_indexes<-as.numeric(opendap_indexes)
  return(opendap_indexes)
}

# To open a MODIS dataset that was downloaded via OpenDap.
fun_preprocess_modis_product<-function(path_to_raw_modis,var_name){
  grid_nc<-raster(path_to_raw_modis,varname=var_name)
  projection(grid_nc)<-modis_crs
  extent(grid_nc)[1:2]<-extent(grid_nc)[1:2]+res(grid_nc)[1]/2
  extent(grid_nc)[3:4]<-extent(grid_nc)[3:4]-res(grid_nc)[1]/2
  #grid_nc <- projectRaster(grid_nc, crs = CRS(paste0("+init=epsg:",epsg)))
  return(grid_nc)
}

# To preprocess ERA-5 products
fun_preprocess_era_product<-function(path_to_output_erawind_data,variable_name,roi_sp_utm,epsg,size_output_grid){
  # Open netcdf
  nc <- nc_open(path_to_output_erawind_data)
  # Get lat and lon
  lat <- ncvar_get(nc,'latitude')
  lon <- ncvar_get(nc,'longitude')
  # Get time 
  #t <- ncvar_get(nc, "time")
  #time unit: hours since 1900-01-01
  #ncatt_get(nc,'time')
  #convert the hours into date + hour
  #as_datetime() function of the lubridate package needs seconds
  #timestamp <- as_datetime(c(t*60*60),origin="1900-01-01")
  
  # Get the variable
  variable <- ncvar_get(nc,variable_name)
  
  # Convert each layer to a raster
  brick_era<-NULL
  for (i in 1:dim(variable)[3]){
    variable_rast <- raster(t(variable[,,i]), xmn=min(lon)-0.125, xmx=max(lon)+0.125, ymn=min(lat)-0.125, ymx=max(lat)+0.125, crs=CRS("+init=epsg:4326"))
    variable_rast <- projectRaster(variable_rast, crs = CRS(paste0("+init=epsg:",epsg)))
    # Resample to size_output_grid size 
    #size_output_grid<-250
    r<-variable_rast
    res(r)<-c(size_output_grid,size_output_grid)
    variable_rast_resamp<-resample(variable_rast,r,method='bilinear')
    # Crop to ROI
    variable_rast_resamp<-crop(variable_rast_resamp,roi_sp_utm)
    # Add to brick
    brick_era<-c(brick_era,variable_rast_resamp)
  }
  
  return(brick_era)
}


## Set working directory
setwd(path_to_processing_folder)

## Set the paths of output folders / files and create them
path_to_dem_folder<-file.path(path_to_processing_folder,"DEM_SRTM") # Path to the folder where the DEM raw data will be stored
path_to_hrsl_folder<-file.path(path_to_processing_folder,"HRSL") 
path_to_modislst_folder<-file.path(path_to_processing_folder,"MODIS_LST")
path_to_modisveget_folder<-file.path(path_to_processing_folder,"MODIS_veget")
path_to_modisevapo_folder<-file.path(path_to_processing_folder,"MODIS_evapo")
path_to_gpm_folder<-file.path(path_to_processing_folder,"GPM")
path_to_nighttime_folder<-file.path(path_to_processing_folder,"VIIRS_nighttime")
path_to_erawind_folder<-file.path(path_to_processing_folder,"ERA_WIND")
path_to_imcce_folder<-file.path(path_to_processing_folder,"IMCCE_Moon")

directories<-list(path_to_dem_folder,path_to_hrsl_folder,path_to_modislst_folder,path_to_gpm_folder,path_to_nighttime_folder,path_to_erawind_folder,path_to_modisveget_folder,path_to_imcce_folder,path_to_modisevapo_folder)
lapply(directories, dir.create)

## Set ROI as sp SpatialPolygon object in epsg 4326
roi_sp_4326<-rgdal::readOGR(path_to_roi_vector)

## Get UTM WGS84 Zone number. from https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
cat("Warning: ROIs overlapping more than 1 UTM zone are currently not adapted in this workflow\n")
utm_zone_number<-(floor((bbox(roi_sp_4326)[1,1] + 180)/6) %% 60) + 1
if(bbox(roi_sp_4326)[2,1]>0){ # if latitudes are North
  epsg<-as.numeric(paste0("326",utm_zone_number))
} else { # if latitude are South
  epsg<-as.numeric(paste0("325",utm_zone_number))
}

## Get MODIS tile number(s) for the ROI
if(!file.exists(file.path(path_to_processing_folder,"modis_sin.kmz"))){
  download.file("https://modis.ornl.gov/files/modis_sin.kmz",destfile = file.path(path_to_processing_folder,"modis_sin.kmz"))
}

tiles = readOGR(file.path(path_to_processing_folder,"modis_sin.kmz"),"Features")
coordinates = SpatialPoints(t(bbox(roi_sp_4326)), CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
modis_tile = as.character(raster::extract(tiles,coordinates)$Name)

if(length(unique(modis_tile))>1){
  cat("Your ROI is overlapping more than 1 MODIS tile. This workflow is currently not adapted for this case\n")
} else {
  modis_tile<-unique(modis_tile)
  modis_tile=gsub(" ","",modis_tile)
  modis_tile=gsub(":","",modis_tile)
  for (i in 1:9){
    modis_tile<-gsub(paste0("h",i,"v"),paste0("h0",i,"v"),modis_tile)
  }
  if(nchar(modis_tile)!=6){
    modis_tile<-paste0(substr(modis_tile,1,4),"0",substr(modis_tile,5,5))
  }
}


## Set GRASS environment and database location 
loc <- rgrass7::initGRASS(path_to_grassApplications_folder, home=getwd(), gisDbase="GRASS_TEMP", override=TRUE,mapset = "PERMANENT" )
execGRASS("g.proj",flags="c",parameters = list(proj4=paste0("+proj=utm +zone=",utm_zone_number," +datum=WGS84 +units=m +no_defs")))

## Get buffer size in degrees
fun_convert_meters_to_degrees<-function(buffer_size_meters,mean_latitude){
  buffer_size_degrees <- buffer_size_meters / (111.32 * 1000 * cos(mean_latitude * ((pi / 180))))
return(buffer_size_degrees)
}

mean_latitude<-mean(bbox(roi_sp_4326)[2,])
buffer_size_degrees<-fun_convert_meters_to_degrees(buffer_size_meters,mean_latitude)

## Set ROI as sp SpatialPolygon object in UTM epsg and in MODIS projection
roi_sp_utm <- spTransform(roi_sp_4326, CRS(paste0("+init=epsg:",epsg)))
roi_sp_modis_project <- spTransform(roi_sp_4326, modis_crs)

## Retrieve the dataset with location and dates of Human Catch Landings
if (!is.null(path_to_database)){
# Connect to project database
react_db <- dbConnect(RSQLite::SQLite(),path_to_database)
# Get query to retrieve dates and locations of human landing catches and execute it. The query is stored on my github repository
sql_query_hlc_dates_location<-paste(readLines(path_to_sql_query_dates_loc_hlc), collapse="\n")
df_dates_locations_hlc<-dbGetQuery(react_db, sql_query_hlc_dates_location)
} else {
  df_dates_locations_hlc<-read.csv(path_to_csv_dates_loc)
}

## Filter data only for CIV
df_dates_locations_hlc <- df_dates_locations_hlc %>% filter (codepays_fk=="CI")
## Turn dataframe to HLC to sp SpatialPointsDataFrame
dates_locations_hlc_sp<-SpatialPointsDataFrame(coords=data.frame(df_dates_locations_hlc$longitude,df_dates_locations_hlc$latitude),data=df_dates_locations_hlc,proj4string=CRS("+init=epsg:4326"))
#dates_locations_hlc_sp_modis_proj<-spTransform(dates_locations_hlc_sp_4326_proj,modis_crs)

# Check that all the points are within the ROI 
bbox_hlc<-bbox(dates_locations_hlc_sp)
bbox_roi<-bbox(roi_sp_4326)
if (bbox_hlc[1,1]<bbox_roi[1,1]|bbox_hlc[1,2]>bbox_roi[2,1]|bbox_hlc[2,1]<bbox_roi[2,1]||bbox_hlc[2,2]>bbox_roi[2,2]){
  stop("Some points of your HLC data lie outside of the ROI")
}

# To delete when the coordinates of the points are validated
dates_locations_hlc_sp<-crop(dates_locations_hlc_sp,extent(roi_sp_4326))




###############################################################################################################
############################### A. Static data ################################
###############################################################################################################

##########################################################################
########### A.1 DEM and DEM-derivatives ###############
##########################################################################

# Set output paths
dem_output_path<-file.path(path_to_dem_folder,"DEM.tif")
slope_output_path<-file.path(path_to_dem_folder,"slope.tif")
aspect_output_path<-file.path(path_to_dem_folder,"aspect.tif")
accumulation_output_path<-file.path(path_to_dem_folder,"accumulation.tif")
tci_output_path<-file.path(path_to_dem_folder,"tci.tif") 
twi_output_path<-file.path(path_to_dem_folder,"twi.tif") 

#####################################
########### A.1.1 Download the data ###############
#####################################

for (i in 1:length(url_to_srtm_datasets)){
  srtm_tile_name<-gsub(".*(.*)/","\\1",url_to_srtm_datasets[i])
  if (!(file.exists(file.path(path_to_dem_folder,srtm_tile_name)))){
  httr::GET(url_to_srtm_datasets[i],write_disk(file.path(path_to_dem_folder,srtm_tile_name)))
  unzip(file.path(path_to_dem_folder,srtm_tile_name),exdir = path_to_dem_folder)
  }
}

# Merge, crop and reproject to UTM the rasters
if (!(file.exists(dem_output_path))){
  path_to_srtm_tiles<-list.files(path_to_dem_folder,pattern = "hgt",full.names = T)

  s <- lapply(path_to_srtm_tiles, stack)
  m <- do.call(merge, s)
  m <- crop(m,roi_sp_4326)

  writeRaster(m,file.path(path_to_dem_folder,"DEM.tif"))

  # Convert from EPSG 4326 (default SRTM EPSG) to UTM EPSG
  gdalwarp(srcfile=dem_output_path,dstfile=gsub("DEM.tif","DEM_temp.tif",dem_output_path),dstnodata=0,t_srs=paste0("+proj=utm +zone=",utm_zone_number," +datum=WGS84 +units=m +no_defs"))
  file.remove(dem_output_path)
  file.rename(gsub("DEM.tif","DEM_temp.tif",dem_output_path),dem_output_path)
  
}

#####################################
########### A.1.2 Prepare the data ###############
#####################################

# extract indices from the DEM : slope, aspect, flow accumulation, topographic convergence index ####

## We use GRASS, calling it in R using the "rgrass7" package. We use two GRASS applications: r.slope.aspect and r.terraflow . Grass must be installed on the computer.

if(!file.exists(slope_output_path)|!file.exists(aspect_output_path)|!file.exists(accumulation_output_path)|!file.exists(tci_output_path)|!file.exists(twi_output_path)){
  # Import DEM to GRASS and set region
  execGRASS("r.external", flags="o", parameters=list(input=file.path(path_to_dem_folder,"DEM.tif"), output="tmprast",band=1))
  execGRASS("g.region", parameters=list(raster="tmprast")) 

  # Compute slope and aspect and save to disk
  execGRASS("r.slope.aspect", flags="overwrite", parameters=list(elevation="tmprast", slope="slope",aspect="aspect",format="percent", precision="FCELL",zscale=1,min_slope=0))
  execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="slope", output=slope_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))
  execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="aspect", output=aspect_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))

  # Compute hydrograpy indices and save to disk
  execGRASS("r.terraflow", flags="overwrite", parameters=list(elevation="tmprast", accumulation="accumulation", tci="tci"))
  execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="accumulation", output=accumulation_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))
  execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="tci", output=tci_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))

  # Compute TWI indice
  execGRASS("r.topidx", flags="overwrite", parameters=list(input="tmprast", output="twi"))
  execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="twi", output=twi_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))
}

# open all rasters
dem_and_derivatives_rast <- brick(list(dem_output_path,slope_output_path,aspect_output_path,accumulation_output_path,tci_output_path,twi_output_path))
names(dem_and_derivatives_rast)[1]<-"elevation"

# calculate the stats within the buffer. First convert HLC dataset to UTM proj . EDIT : for the accumulation, might be interesting to also extract the maximum
dates_locations_hlc_sp <- spTransform(dates_locations_hlc_sp,CRS(paste0("+init=epsg:",epsg)))
dates_locations_hlc_sp <- raster::extract(dem_and_derivatives_rast, dates_locations_hlc_sp, buffer=buffer_size_meters,fun=mean, na.rm=TRUE, sp=TRUE,small=FALSE) 
##### WARNING : the function raster::extract can produce error that are not visible (ie no error is returned) if points are outside the raster extent. It is very important to check that at least one raster cell is intersected


##########################################################################
########### A.2 Settlements and population ###############
##########################################################################

## Variables : surface built, shape indices

#####################################
########### A.2.1 Open the data ###############
#####################################

## Surface built
inertia_texture<-raster(path_to_texture_inertia)
names(inertia_texture)<-"built_surf"


## Location and population of households :

# Output is a data.frame with columns:
# - population = population in the household (numeric)
# - latitude = latitude of the household (numeric)
# - longitude = longitude of the household (numeric)
# - village = village name (or code) (character)

# Query DB 
sql_query_households_loc_pop<-paste(readLines(path_to_sql_query_households_loc_pop), collapse="\n")
df_households_loc_pop<-dbGetQuery(react_db, sql_query_households_loc_pop)
# Tidy the dataset
df_households_loc_pop <- df_households_loc_pop %>% filter (!is.na(latitude)) %>% filter (latitude!="")
df_households_loc_pop$latitude <- as.numeric(gsub(",",".",df_households_loc_pop$latitude))
df_households_loc_pop$longitude <- as.numeric(gsub(",",".",df_households_loc_pop$longitude))
df_households_loc_pop$codemenage<-df_households_loc_pop$codepays_fk<-NULL


#####################################
########### A.2.1 Prepare the data ###############
#####################################

### Extract information : population in each village, population density, built surface, distance from each catch point to the edge of the village

## Turn location and population of households to SpatialPointsDataFrame
df_households_loc_pop_sp<-SpatialPointsDataFrame(coords=data.frame(df_households_loc_pop$longitude,df_households_loc_pop$latitude),data=df_households_loc_pop,proj4string=CRS("+init=epsg:4326"))
df_households_loc_pop_sp<-raster::crop(df_households_loc_pop_sp,extent(roi_sp_4326))
df_households_loc_pop<-as.data.frame(df_households_loc_pop_sp)
  
## Get population of each village (from the census dataset)
df_households_village_pop<-df_households_loc_pop %>% group_by(village) %>% summarize(population=sum(population)) 

## Get convex hull polygon for each village. We consider the convex hull as the edges of the village. The convex hull is the minimum polygon that encompasses all the locations of the households.  (code from https://stackoverflow.com/questions/25606512/create-convex-hull-polygon-from-points-and-save-as-shapefile)
villages<-unique(df_households_loc_pop$village)
ps_list<-NULL
for (i in 1:length(villages)){
  this_village<-villages[i]
  df_households_loc_pop_this_village<-df_households_loc_pop %>% filter(village==this_village)
  dat <- as.matrix(data.frame(df_households_loc_pop_this_village$longitude,df_households_loc_pop_this_village$latitude))
  ch<-chull(dat)
  coords <- dat[c(ch, ch[1]), ]  # closed polygon
  p = Polygon(coords)
  ps = Polygons(list(p),i)
  ps_list<-c(ps_list,ps)
}

spatialpoly<-as.SpatialPolygons.PolygonsList(ps_list,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
df_villages_loc_pop_sp<-SpatialPolygonsDataFrame(spatialpoly,df_households_village_pop)

## Get area (in m2) and population density (pers / m2) of each village
df_villages_loc_pop_sp$area<-raster::area(df_villages_loc_pop_sp)
df_villages_loc_pop_sp$pop_density<-df_villages_loc_pop_sp$population/df_villages_loc_pop_sp$area

## Get the surface built in each village
# We use the haralick correlation texture file calculated with the Spot 6 image. The pixels with value > threshold_built_areas (26000 for CIV) are considered as built surfaces (visual thresholding), other values are unbuilt surfaces.

df_villages_loc_pop_sp$built_surf<-0
for (i in 1:length(villages)){
  df_villages_loc_pop_sp_th_village<-df_villages_loc_pop_sp[which(df_villages_loc_pop_sp$village==villages[i]),]
  df_villages_loc_pop_sp_th_village<-spTransform(df_villages_loc_pop_sp_th_village,CRS(paste0("+init=epsg:",epsg)))
  inertia_this_village<-raster::crop(inertia_texture,df_villages_loc_pop_sp_th_village)
  inertia_this_village[inertia_this_village<=threshold_built_areas]<-NA
  inertia_this_village[inertia_this_village>threshold_built_areas]=1
  df_villages_loc_pop_sp_th_village <- raster::extract(inertia_this_village, df_villages_loc_pop_sp_th_village,fun=sum, na.rm=TRUE, df=TRUE,small=TRUE) 
  df_villages_loc_pop_sp$built_surf[which(df_villages_loc_pop_sp$village==villages[i])]<-df_villages_loc_pop_sp_th_village$built_surf
}

# Note: to do the same thing we could have done the way here-under, easier. However the raster is very big, hence the step inertia_texture[inertia_texture<=3]<-NA takes too long. Hence we do it the 'nasty' (ie loop) way
#inertia_texture[inertia_texture<=3]<-NA
#inertia_texture[inertia_texture>3]=1
#df_villages_loc_pop_sp <- raster::extract(inertia_texture, df_villages_loc_pop_sp,fun=sum, na.rm=TRUE, sp=TRUE,small=TRUE) 

# We put the built surface in m2. 1 pixel is 1.57 m * 1.57 m (res(inertia_texture))
df_villages_loc_pop_sp$built_surf<-df_villages_loc_pop_sp$built_surf*res(inertia_texture)[1]^2

# Compute built-up density
df_villages_loc_pop_sp$built_density<-df_villages_loc_pop_sp$built_surf/df_villages_loc_pop_sp$area

## Finally link to dataset of dates and locations of HLC
dates_locations_hlc_sp<-merge(dates_locations_hlc_sp,as.data.frame(df_villages_loc_pop_sp),by="village")



## Get distance from each catch point to the edge of the village (from https://stackoverflow.com/questions/28382949/finding-the-minimum-distance-between-all-points-and-the-polygon-boundary)
dates_locations_hlc_sp$dist_to_edge_vill<-0
for (i in 1:nrow(dates_locations_hlc_sp)){
  df_villages_loc_pop_sp_th_point<-df_villages_loc_pop_sp[which(df_villages_loc_pop_sp$village==dates_locations_hlc_sp$village[i]),]
  df_villages_loc_pop_sp_th_point<-spTransform(df_villages_loc_pop_sp_th_point,CRS(paste0("+init=epsg:",epsg)))
  dates_locations_hlc_sp_th_point<-spTransform(dates_locations_hlc_sp[i,],CRS(paste0("+init=epsg:",epsg)))
  dates_locations_hlc_sp$dist_to_edge_vill[i]<-round(as.numeric(rgeos::gDistance(dates_locations_hlc_sp_th_point, as(df_villages_loc_pop_sp_th_point, "SpatialLines"), byid = TRUE)))
}


## Get an indice of the shape of the village (is is dense ? spread ? ). For this we use the locations of the households



##########################################################################
########### A.4 Land use / land cover ###############
##########################################################################

#####################################
########### A.4.1 Prepare the data ###############
#####################################


###############################################################################################################
############################### B. Dynamic data ################################
###############################################################################################################

## Get all dates
all_dates_hlc<-unique(df_dates_locations_hlc$datecapture)


## For tests: retrieve all the lines for the first date
i=20
locations_hlc_sp_this_date<-dates_locations_hlc_sp[which(dates_locations_hlc_sp$datecapture==all_dates_hlc[i]),]
this_date_hlc<-as.Date(all_dates_hlc[i])


# Initiate cluster for paralell download 
#no_cores <- detectCores() - 1
#cl <- makeCluster(no_cores, type = "PSOCK")
#registerDoParallel(cl)

## Parameters for download of ERA INTERIM data
#import the python library ecmwfapi
#ecmwf <- reticulate::import('ecmwfapi')
#for this step there must exist the file .ecmwfapirc
#server = ecmwf$ECMWFDataServer() #start the connection


############ Download and prepare OpenDap time and spatial indexes to further download the MODIS and GPM datasets on opendap servers ############

## Set OpenDap URLs to each collection
url_opendap_modis_lst_terra<-paste0(url_modis_opendap,"/",modis_lst_terra_collection,"/",modis_tile)
url_opendap_modis_lst_aqua<-paste0(url_modis_opendap,"/",modis_lst_aqua_collection,"/",modis_tile)
url_opendap_modis_veget_terra<-paste0(url_modis_opendap,"/",modis_veget_terra_collection,"/",modis_tile)
url_opendap_modis_veget_aqua<-paste0(url_modis_opendap,"/",modis_veget_aqua_collection,"/",modis_tile)
url_opendap_modis_evapo_terra<-paste0(url_modis_opendap,"/",modis_evapo_terra_collection,"/",modis_tile)
url_opendap_modis_evapo_aqua<-paste0(url_modis_opendap,"/",modis_evapo_aqua_collection,"/",modis_tile)

## Download the time indexes on opendap for each collection
if (!(file.exists(file.path(path_to_modislst_folder,"opendap_modis_lst_terra_time_index.txt")))){
  httr::GET(paste0(url_opendap_modis_lst_terra,".ncml.ascii?time"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_terra_time_index.txt")))
}
if (!(file.exists(file.path(path_to_modislst_folder,"opendap_modis_lst_aqua_time_index.txt")))){
  httr::GET(paste0(url_opendap_modis_lst_aqua,".ncml.ascii?time"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_aqua_time_index.txt")))
}
if (!(file.exists(file.path(path_to_modisveget_folder,"opendap_modis_veget_terra_time_index.txt")))){
  httr::GET(paste0(url_opendap_modis_veget_terra,".ncml.ascii?time"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_terra_time_index.txt")))
}
if (!(file.exists(file.path(path_to_modisveget_folder,"opendap_modis_veget_aqua_time_index.txt")))){
  httr::GET(paste0(url_opendap_modis_veget_aqua,".ncml.ascii?time"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_aqua_time_index.txt")))
}
if (!(file.exists(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_terra_time_index.txt")))){
  httr::GET(paste0(url_opendap_modis_evapo_terra,".ncml.ascii?time"),write_disk(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_terra_time_index.txt")))
}
if (!(file.exists(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_aqua_time_index.txt")))){
  httr::GET(paste0(url_opendap_modis_evapo_aqua,".ncml.ascii?time"),write_disk(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_aqua_time_index.txt")))
}

## Get the time indexes as vectors
opendap_modis_lst_terra_time_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_terra_time_index.txt"))
opendap_modis_lst_aqua_time_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_aqua_time_index.txt"))
opendap_modis_veget_terra_time_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_terra_time_index.txt"))
opendap_modis_veget_aqua_time_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_aqua_time_index.txt"))
opendap_modis_evapo_terra_time_index<-fun_get_opendap_index(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_terra_time_index.txt"))
opendap_modis_evapo_aqua_time_index<-fun_get_opendap_index(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_aqua_time_index.txt"))

## Downlaad the spatial indexes on opendap for each collection
if (!(file.exists(file.path(path_to_modislst_folder,"opendap_modis_lst_XDim_index.txt")))){
  httr::GET(paste0(url_opendap_modis_lst_terra,".ncml.ascii?XDim"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_XDim_index.txt")))
}
if (!(file.exists(file.path(path_to_modislst_folder,"opendap_modis_lst_YDim_index.txt")))){
  httr::GET(paste0(url_opendap_modis_lst_terra,".ncml.ascii?YDim"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_YDim_index.txt")))
}    
if (!(file.exists(file.path(path_to_modisveget_folder,"opendap_modis_veget_XDim_index.txt")))){
  httr::GET(paste0(url_opendap_modis_veget_terra,".ncml.ascii?XDim"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_XDim_index.txt")))
}   
if (!(file.exists(file.path(path_to_modisveget_folder,"opendap_modis_veget_YDim_index.txt")))){
  httr::GET(paste0(url_opendap_modis_veget_terra,".ncml.ascii?YDim"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_YDim_index.txt")))
}   
if (!(file.exists(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_XDim_index.txt")))){
  httr::GET(paste0(url_opendap_modis_evapo_terra,".ncml.ascii?XDim"),write_disk(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_XDim_index.txt")))
}   
if (!(file.exists(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_YDim_index.txt")))){
  httr::GET(paste0(url_opendap_modis_evapo_terra,".ncml.ascii?YDim"),write_disk(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_YDim_index.txt")))
}   
if (!(file.exists(file.path(path_to_gpm_folder,"opendap_gpm_lon_index.txt")))){
  httr::GET(paste0(url_gpm_opendap,"/2018/02/3B-DAY.MS.MRG.3IMERG.20180201-S000000-E235959.V06.nc4.ascii?lon"),write_disk(file.path(path_to_gpm_folder,"opendap_gpm_lon_index.txt")))
}   
if (!(file.exists(file.path(path_to_gpm_folder,"opendap_gpm_lat_index.txt")))){
  httr::GET(paste0(url_gpm_opendap,"/2018/02/3B-DAY.MS.MRG.3IMERG.20180201-S000000-E235959.V06.nc4.ascii?lat"),write_disk(file.path(path_to_gpm_folder,"opendap_gpm_lat_index.txt")))
}

## Get the spatial indexes as vectors
opendap_modis_lst_XDim_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_XDim_index.txt"))
opendap_modis_lst_YDim_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_YDim_index.txt"))
opendap_modis_veget_XDim_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_XDim_index.txt"))
opendap_modis_veget_YDim_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_YDim_index.txt"))
opendap_modis_evapo_XDim_index<-fun_get_opendap_index(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_XDim_index.txt"))
opendap_modis_evapo_YDim_index<-fun_get_opendap_index(file.path(path_to_modisevapo_folder,"opendap_modis_evapo_YDim_index.txt"))
opendap_gpm_lon_index<-fun_get_opendap_index(file.path(path_to_gpm_folder,"opendap_gpm_lon_index.txt"))
opendap_gpm_lat_index<-fun_get_opendap_index(file.path(path_to_gpm_folder,"opendap_gpm_lat_index.txt"))


## Extract indexes for lat min, lat max, lon min and lon max for our bounding box
roi_bbox_4326<-bbox(roi_sp_4326) ## for MODIS 
roi_bbox_modisproject<-bbox(roi_sp_modis_project) # for GPM

index_opendap_gpm_lon_min<-which.min(abs(opendap_gpm_lon_index-roi_bbox_4326[1,1]))-4
index_opendap_gpm_lon_max<-which.min(abs(opendap_gpm_lon_index-roi_bbox_4326[1,2]))+4 
index_opendap_gpm_lat_min<-which.min(abs(opendap_gpm_lat_index-roi_bbox_4326[2,1]))-4
index_opendap_gpm_lat_max<-which.min(abs(opendap_gpm_lat_index-roi_bbox_4326[2,2]))+4 

index_opendap_modisveget_lon_min<-which.min(abs(opendap_modis_veget_XDim_index-roi_bbox_modisproject[1,1]))-1
index_opendap_modisveget_lon_max<-which.min(abs(opendap_modis_veget_XDim_index-roi_bbox_modisproject[1,2]))-1 
index_opendap_modisveget_lat_max<-which.min(abs(opendap_modis_veget_YDim_index-roi_bbox_modisproject[2,1]))-1
index_opendap_modisveget_lat_min<-which.min(abs(opendap_modis_veget_YDim_index-roi_bbox_modisproject[2,2]))-1 

index_opendap_modisevapo_lon_min<-which.min(abs(opendap_modis_evapo_XDim_index-roi_bbox_modisproject[1,1]))-1
index_opendap_modisevapo_lon_max<-which.min(abs(opendap_modis_evapo_XDim_index-roi_bbox_modisproject[1,2]))-1 
index_opendap_modisevapo_lat_max<-which.min(abs(opendap_modis_evapo_YDim_index-roi_bbox_modisproject[2,1]))-1
index_opendap_modisevapo_lat_min<-which.min(abs(opendap_modis_evapo_YDim_index-roi_bbox_modisproject[2,2]))-1 

index_opendap_modislst_lon_min<-which.min(abs(opendap_modis_lst_XDim_index-roi_bbox_modisproject[1,1]))-1
index_opendap_modislst_lon_max<-which.min(abs(opendap_modis_lst_XDim_index-roi_bbox_modisproject[1,2]))-1 
index_opendap_modislst_lat_max<-which.min(abs(opendap_modis_lst_YDim_index-roi_bbox_modisproject[2,1]))-1
index_opendap_modislst_lat_min<-which.min(abs(opendap_modis_lst_YDim_index-roi_bbox_modisproject[2,2]))-1 


## Indexes of the coordinates of the bounding box for OpenDap query to download GPM and modis vegetation data.
# OpenDap is a standard protocol that enables to downlaad the data only for a given ROI, instead of downloading the whole data which is provided by default on a given tile. 
# The coordinates of the bounding box are not provided in lat and lon. They are provided as indexes. 
#index_opendap_gpm_lon_min<-1741  # For CIV
#index_opendap_gpm_lon_max<-1745  # For CIV
#index_opendap_gpm_lat_min<-989  # For CIV
#index_opendap_gpm_lat_max<-995  # For CIV

#index_opendap_gpm_lon_min<-1763  # For BF
#index_opendap_gpm_lon_max<-1768  # For BF
#index_opendap_gpm_lat_min<-1005  # For BF
#index_opendap_gpm_lat_max<-1010  # For BF

# For CIV

#index_opendap_modisveget_lon_min<-2030
#index_opendap_modisveget_lon_max<-2240
#index_opendap_modisveget_lat_max<-510
#index_opendap_modisveget_lat_min<-200

#index_opendap_modislst_lon_min<-
#  index_opendap_modislst_lon_max<-
#  index_opendap_modislst_lat_max<-
#  index_opendap_modislst_lat_min<-

# For BF
#modis_tile<-"h17v07"
#index_opendap_modisveget_lon_min<-3080
#index_opendap_modisveget_lon_max<-3335
#index_opendap_modisveget_lat_max<-4530
#index_opendap_modisveget_lat_min<-4300


########################################################################################################################

########################################################################################################################
############ Start Workflow ############
########################################################################################################################
cat("Starting workflow")


#####################################
########### 1. Modis LST ###############
#####################################

# Il y a deux résolutions temporelles disponibles: 1 jour et 8 jours. Pour chaque résolution temporelle il y a 2 fichiers disponibles: le fichier issu de Terra et le fichier issu de Aqua (relevés distants d'à peu près 2/3 h entre Terra et Aqua). Enfin pour chaque fichier il y a 2 relevés: température diurne et température nocturne. La différence entre Terra et Aqua est l'heure de relevé, ainsi que l'étendue des données disponibles (ie couverture nuageuse au moment de l'acquisition)

### Questions à répondre :
# Prend-on le fichier 1 jour ou le fichier 8 jours ? Peut-etre prendre le fichier 1 jour pour la nuit de capture, et le fichier (ou 2 fichiers) 8 jours pour étudier le développement larvaire
# Prend-on le fichier Terra ou Acqua ou une composition des 2 (i.e. on rempli les trous de Terra par ceux de Aqua ou vice-versa) ? 

# Je prendrais: la résolution 1km/1jour et les variables suivantes: température diurne, température nocturne, différence températures diurne/nocturne, 
# Puis avec des cross corrélation maps on calcule la corrélation entre la variable à expliquer (eg incidence) et la moyenne (+stdev?) de chacune des variables 
# ensuite pour chaque village on stocke le résultat des CCM dans un tableau. On regarde ensuite quelles sont les différences entre les villages (i.e. quelles sont les variations selon les villages des meilleures corrélations entre la température et l'indicdence ? )

# Get modis names available with GetSpatialData. (See names here: https://modis.gsfc.nasa.gov/data/dataprod/ )
# product_names <- getSpatialData::getMODIS_names()

# MODIS LST products are listed here https://modis.gsfc.nasa.gov/data/dataprod/mod11.php) : 
# We are interested in daily L3 global 1km so we take : MOD11A1 , MYD11A1 
# Note that the 8-Day product (MOD11A2, MYD11A2) is simply an average of the 1 day products over the time period

# Apart from the GetSpatialData library, another option is to access MODIS data via OpenDAP. All the collections are available here:
# For instance for MOD11A1 v6 : https://opendap.cr.usgs.gov/opendap/hyrax/MYD11A1.006/contents.html
# Then we need to select the MODIS tile for our ROI. For instance for CIV the tile is h17v08, hence our data is: https://opendap.cr.usgs.gov/opendap/hyrax/MYD11A1.006/h17v08.ncml.html





## Thèse Nico :
# - 7 jours. Il prend les données de la semaine de capture + 2 semaines précédentes. 
# - stage Mader : moyenne sur les 7 jours aussi
# - Weiss et al : 1 month

## For the quality control of MODIS layers : https://www.r-bloggers.com/modis-qc-bits/

##############################################################
#### 1.1 - Download the data ####
##############################################################

dates_modis_lst <-seq(this_date_hlc,this_date_hlc-lag_days_modis_lst,-1)

# Function to build the OpenDap URL to dowload MODIS LST night and day 1km data. input parameters : date and satellite (terra or aqua)
fun_build_modis_lst_opendap_url<-function(date_catch,satellite){
  
  index_date<-as.integer(difftime(date_catch ,"2000-01-01" , units = c("days"))) # time is provided as number of days since 2000-01-01
  
  if (satellite=="terra"){
    vector_time_index<-opendap_modis_lst_terra_time_index
    url_opendap<-url_opendap_modis_lst_terra
  } else if (satellite=="aqua"){
    vector_time_index<-opendap_modis_lst_aqua_time_index
    url_opendap<-url_opendap_modis_lst_aqua
  }
  
  opendap_index_date<-which(vector_time_index==index_date)-1
  
  url_opendap_data<-paste0(url_opendap,".ncml.nc4?MODIS_Grid_Daily_1km_LST_eos_cf_projection,LST_Day_1km[",opendap_index_date,"][",index_opendap_modislst_lat_min,":",index_opendap_modislst_lat_max,"][",index_opendap_modislst_lon_min,":",index_opendap_modislst_lon_max,"],time[",opendap_index_date,"],YDim[",index_opendap_modislst_lat_min,":",index_opendap_modislst_lat_max,"],XDim[",index_opendap_modislst_lon_min,":",index_opendap_modislst_lon_max,"],LST_Night_1km[",opendap_index_date,"][",index_opendap_modislst_lat_min,":",index_opendap_modislst_lat_max,"][",index_opendap_modislst_lon_min,":",index_opendap_modislst_lon_max,"],time[",opendap_index_date,"],YDim[",index_opendap_modislst_lat_min,":",index_opendap_modislst_lat_max,"],XDim[",index_opendap_modislst_lon_min,":",index_opendap_modislst_lon_max,"]")
  
  return(url_opendap_data)
}


## Download the data 
for (i in 1:length(dates_modis_lst)){
  cat(paste0("Downloading MODIS LST Terra and Aqua (MOD11A1 and MYD11A1) for the ROI and for date ",dates_modis_lst[i],"\n"))
  for (j in c("terra","aqua")){
    opendap_modis<-fun_build_modis_lst_opendap_url(dates_modis_lst[i],j)
    path_to_dataset<-file.path(path_to_modislst_folder,paste0(gsub("-","",dates_modis_lst[i]),"_",j,".nc4"))
    if (!(file.exists(path_to_dataset))){
      httr::GET(opendap_modis,write_disk(path_to_dataset))
    } else {
      cat(paste0("Data already exist for date ",dates_modis_lst[i],"\n"))
    }
  }
}



##############################################################
#### 1.2 - Prepare the data ####
##############################################################

## For each date before the d date:
# - For each cell, extract the maximum day available temperature among Terra and Aqua. In case of NA on one of the two datasets, take the value of the other dataset.
# - For each cell, extract the minimum night available temperature among Terra and Aqua. In case of NA on one of the two datasets, take the value of the other dataset.
# - Extract the mean of the available LST cells on a buffer around each catch point
# - If no data are available for a given catch point (i.e. cells are all NA), set temperature of the previous night (recursively until a value is available)
# - If no data are available for the whole time series, ??? TODO
# - Create thermal amplitude value (i.e. difference between day and night)


# Get Terra and Aqua products to preprocess 
terra_names<-file.path(path_to_modislst_folder,paste0(gsub("-","",dates_modis_lst),"_terra.nc4"))
aqua_names<-file.path(path_to_modislst_folder,paste0(gsub("-","",dates_modis_lst),"_aqua.nc4"))


brick_lst_day<-NULL
brick_lst_night<-NULL

for (i in 1:length(dates_modis_lst)){

  cat(paste0("pre-processing LST for date ",dates_modis_lst[i],"\n"))

# Create raster of LST_day maximum value combining Terra and Aqua by taking the maximum available value for each pixel
rast_lst_day_terra<-fun_preprocess_modis_product(terra_names[i],"LST_Day_1km")
rast_lst_day_aqua<-fun_preprocess_modis_product(aqua_names[i],"LST_Day_1km")
rast_lst_day<-brick(rast_lst_day_terra,rast_lst_day_aqua)

rast_max_lst_day<-max(rast_lst_day,na.rm=TRUE)
brick_lst_day<-c(brick_lst_day,rast_max_lst_day)
names(brick_lst_day[[length(brick_lst_day)]])<-paste0("lst_max_",i-1)

# Create raster of LST_night minimum value combining Terra and Aqua by taking the minimum available value for each pixel
rast_lst_night_terra<-fun_preprocess_modis_product(terra_names[i],"LST_Night_1km")
rast_lst_night_aqua<-fun_preprocess_modis_product(aqua_names[i],"LST_Night_1km")
rast_lst_night<-brick(rast_lst_night_terra,rast_lst_night_aqua)

rast_min_lst_night<-min(rast_lst_night,na.rm=TRUE)
brick_lst_night<-c(brick_lst_night,rast_min_lst_night)
names(brick_lst_night[[length(brick_lst_night)]])<-paste0("lst_min_",i-1)

}

brick_lst_day<-brick(brick_lst_day)
brick_lst_night<-brick(brick_lst_night)

## Extract mean of LST in a buffer of 2 km for all the dates
# na.rm = TRUE means that NA values in the raster are not taken into account
# small=TRUE means that each time the buffer intersects a non NA cell it takes into account the cell value
# sp=TRUE means that the extracted values are added to the data.frame of the Spatial object
locations_hlc_sp_this_date_modis_proj <- raster::extract(brick_lst_day, locations_hlc_sp_this_date_modis_proj, buffer=buffer_size_meters,fun=count, na.rm=TRUE, sp=TRUE,small=TRUE) 
locations_hlc_sp_this_date_modis_proj <- raster::extract(brick_lst_night, locations_hlc_sp_this_date_modis_proj, buffer=buffer_size_meters,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 



#####################################
########### 2. Modis vegetation indices NDVI and EVI ###############
#####################################

## Modis NDVI and EVI are 250 m / 16 days resolutions. We take the data 

## MODIS Vegetations indices are quite huge files (approx. 110 MB for 1 single MODIS tile). We could download them through the GetSpatialData package but here we will prefer OpenDAP 
# The OpenDAP server where many MODIS data are stored is https://opendap.cr.usgs.gov/opendap/hyrax/ , the v6 collection for Modis Vegetation indices 250 m / 16 days Terra and Aqua are respectively MOD13Q1.006 and MYD13Q1.006
# The OpenDAP servers are hence https://opendap.cr.usgs.gov/opendap/hyrax/MOD13Q1.006 and https://opendap.cr.usgs.gov/opendap/hyrax/MYD13Q1.006
# On these servers the data are filtered by tile name. For us, for CIV the tile is h17v08 and for BF it's h17v07 . So the links to the ncml files are respectively : 
# Terra : https://opendap.cr.usgs.gov/opendap/hyrax/MOD13Q1.006/h17v08.ncml.html and Aqua : https://opendap.cr.usgs.gov/opendap/hyrax/MYD13Q1.006/h17v08.ncml.html

# What we do is :
# - first we query the metadata of the data available for our time frame using getSpatialData::getMODIS_query
# - from there we retrieve the date of aquisition (attribute AcquisitionDate) of the week including the catch + one week before + two weeks before
# - we convert this date to a integer "number of days since 2000-01-01" which is the time attribute for the openDap query and we divide by 16 to get the index since there is one acquisition every 16 days
# - we download the data using the indexes of the lat and lon coordinates provided as parameters


##############################################################
#### 2.1 - Download the data ####
##############################################################

# MODIS Veget are every 8 days
date_modis_veget<-seq(this_date_hlc,this_date_hlc-lag_days_modis_veget,-8)

## Retrieve date and satellite of the week including the catch (i.e. closest start aquisition date from the catch date) and build the link to download the data
# Function to build the OpenDap URL to dowload MODIS Vegetion indices NDVI and EVI 250m data. input parameters : date

fun_build_modis_veget_opendap_url<-function(date_catch){
  
  date_catch_julian<-as.integer(difftime(date_catch ,"2000-01-01" , units = c("days")))
  
  index_opendap_terra_closest_to_date<-which.min(abs(opendap_modis_veget_terra_time_index-date_catch_julian))
  index_opendap_aqua_closest_to_date<-which.min(abs(opendap_modis_veget_aqua_time_index-date_catch_julian))
  
  days_sep_terra_from_date<-opendap_modis_veget_terra_time_index[index_opendap_terra_closest_to_date]-date_catch_julian
  days_sep_aqua_from_date<-opendap_modis_veget_aqua_time_index[index_opendap_aqua_closest_to_date]-date_catch_julian
  
  if(days_sep_terra_from_date<0 & days_sep_aqua_from_date>0){
    opendap_collection_closest_date_to_catch<-"MOD13Q1.006"
    date_data<-as.Date("2000-01-01")+opendap_modis_veget_terra_time_index[index_opendap_terra_closest_to_date]
    opendap_index_closest_date_to_catch<-index_opendap_terra_closest_to_date-1 # We put minus 1 because on OpenDAP the first index is 0 and not 1
  } else if (days_sep_terra_from_date>0 & days_sep_aqua_from_date<0){
    opendap_collection_closest_date_to_catch<-"MYD13Q1.006"
    date_data<-as.Date("2000-01-01")+opendap_modis_veget_aqua_time_index[index_opendap_aqua_closest_to_date]
    opendap_index_closest_date_to_catch<-index_opendap_aqua_closest_to_date-1
  } else {
   if(abs(days_sep_terra_from_date)<abs(days_sep_aqua_from_date)){
      opendap_collection_closest_date_to_catch<-"MOD13Q1.006"
      date_data<-as.Date("2000-01-01")+opendap_modis_veget_terra_time_index[index_opendap_terra_closest_to_date]
      opendap_index_closest_date_to_catch<-index_opendap_terra_closest_to_date-1
   } else {
      opendap_collection_closest_date_to_catch<-"MYD13Q1.006"
      date_data<-as.Date("2000-01-01")+opendap_modis_veget_aqua_time_index[index_opendap_aqua_closest_to_date]
      opendap_index_closest_date_to_catch<-index_opendap_aqua_closest_to_date-1
   }
  }
  
  url_opendap<-paste0(url_modis_opendap,"/",opendap_collection_closest_date_to_catch,"/",modis_tile,".ncml.nc4?MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection,_250m_16_days_NDVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],_250m_16_days_EVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"]")
  
  return(list(url_opendap,opendap_collection_closest_date_to_catch,opendap_index_closest_date_to_catch,date_data))
  
}

## Download data

modisveget_list_paths<-NULL
modisveget_list_names<-NULL

for (i in 1:length(date_modis_veget)){
  cat(paste0("Downloading MODIS Vegetation indices data for the ROI and for available file closest to date ",date_modis_veget[i],"\n"))
  urls_infos<-fun_build_modis_veget_opendap_url(date_modis_veget[i])
  path<-file.path(path_to_modisveget_folder,paste0(urls_infos[[2]],"_",gsub("-","_",urls_infos[[4]]),".nc4"))
  modisveget_list_paths<-c(modisveget_list_paths,path)
  modisveget_list_names<-c(modisveget_list_names,as.character(abs(difftime(date_modis_veget[i] ,this_date_hlc , units = c("days")))))
  if (!(file.exists(path))){
    httr::GET(url = urls_infos[[1]], write_disk(path))
  } else {
    cat(paste0("Data already exist for date ",date_modis_veget[i],"\n"))
  }
}

#catch_week<-fun_build_modis_veget_opendap_url(this_date_hlc)
#one_week_before<-fun_build_modis_veget_opendap_url(this_date_hlc-8)
#two_week_before<-fun_build_modis_veget_opendap_url(this_date_hlc-16)

#path_catch_week<-file.path(path_to_modisveget_raw_folder,paste0(catch_week[[2]],"_",catch_week[[3]],".nc4"))
#path_one_week_before<-file.path(path_to_modisveget_raw_folder,paste0(one_week_before[[2]],"_",one_week_before[[3]],".nc4"))
#path_two_week_before<-file.path(path_to_modisveget_raw_folder,paste0(two_week_before[[2]],"_",two_week_before[[3]],".nc4"))

#if (!(file.exists(path_catch_week))){
#  httr::GET(url = catch_week[[1]], write_disk(path_catch_week))
#}
#if (!(file.exists(path_one_week_before))){
#  httr::GET(url = one_week_before[[1]], write_disk(path_one_week_before))
#}
#if (!(file.exists(path_two_week_before))){
#  httr::GET(url = two_week_before[[1]], write_disk(path_two_week_before))
#}


##############################################################
#### 2.2 - Prepare the data ####
##############################################################

brick_modisveget<-NULL
for (i in 1:length(modisveget_list_paths)){
  rast_catch_ndvi<-fun_preprocess_modis_product(modisveget_list_paths[i],"_250m_16_days_NDVI")
  rast_catch_evi<-fun_preprocess_modis_product(modisveget_list_paths[i],"_250m_16_days_EVI")
  brick_modisveget<-c(brick_modisveget,rast_catch_ndvi)
  names(brick_modisveget[[length(brick_modisveget)]])<-paste0("ndvi_",modisveget_list_names[i])
  brick_modisveget<-c(brick_modisveget,rast_catch_evi)
  names(brick_modisveget[[length(brick_modisveget)]])<-paste0("evi_",modisveget_list_names[i])
}

brick_modisveget<-brick(brick_modisveget)

# Extract mean of vegetation indices for each raster 
dates_locations_hlc_sp_modis_proj <- raster::extract(brick_modisveget, dates_locations_hlc_sp_modis_proj, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 



#####################################
########### 3. Modis evapotranspiration ###############
#####################################

##############################################################
#### 3.1 - Download the data ####
##############################################################

## Function to look for the closest date to the HLC and retrieve the corresponding OpenDAP time index
## MODIS evaporation are every 8 days, and both Terra and Aqua have the same dates (not as for Vegetation)

dates_modis_evapo<-seq(this_date_hlc,this_date_hlc-lag_days_modis_evapo,-8)

fun_build_modis_evapo_opendap_url<-function(date_catch,satellite){
  
  if (satellite=="terra"){
    opendap_time_index<-opendap_modis_evapo_terra_time_index
    url_opendap<-url_opendap_modis_evapo_terra
  } else if (satellite=="aqua"){
    opendap_time_index<-opendap_modis_evapo_aqua_time_index
    url_opendap<-url_opendap_modis_evapo_aqua
  }
  
  date_catch_julian<-as.integer(difftime(date_catch ,"2000-01-01" , units = c("days")))
  index_opendap_closest_to_date<-which.min(abs(opendap_time_index-date_catch_julian))
  days_sep_from_date<-opendap_time_index[index_opendap_closest_to_date]-date_catch_julian
  if(days_sep_from_date<=0){ 
    index_opendap_closest_to_date<-index_opendap_closest_to_date-1
  } else {
    index_opendap_closest_to_date<-index_opendap_closest_to_date-2
  }

  date<-as.Date("2000-01-01")+opendap_time_index[index_opendap_closest_to_date+1]

  url_opendap<-paste0(url_opendap,".ncml.nc4?MOD_Grid_MOD16A2_eos_cf_projection,ET_500m[",index_opendap_closest_to_date,"][",index_opendap_modisevapo_lat_min,":",index_opendap_modisevapo_lat_max,"][",index_opendap_modisevapo_lon_min,":",index_opendap_modisevapo_lon_max,"],time[",index_opendap_closest_to_date,"],YDim[",index_opendap_modisevapo_lat_min,":",index_opendap_modisevapo_lat_max,"],XDim[",index_opendap_modisevapo_lon_min,":",index_opendap_modisevapo_lon_max,"]")
  
  return(list(url_opendap,date))
}
  

dates_modis_evapo_real<-NULL

## Download the data 
for (i in 1:length(dates_modis_evapo)){
  cat(paste0("Downloading MODIS Evapotranspiration Terra and Aqua (MOD16A2 and MYD16A2) for the ROI and for the closest date to ",dates_modis_evapo[i],"\n"))
  for (j in c("terra","aqua")){
    opendap_modis<-fun_build_modis_evapo_opendap_url(dates_modis_evapo[i],j)
    path_to_dataset<-file.path(path_to_modisevapo_folder,paste0(gsub("-","",opendap_modis[[2]]),"_",j,".nc4"))
    if (!(file.exists(path_to_dataset))){
      httr::GET(opendap_modis[[1]],write_disk(path_to_dataset))
    } else {
      cat(paste0("Data already exist for date ",dates_modis_evapo[i],"\n"))
    }
  }
  dates_modis_evapo_real<-c(dates_modis_evapo_real,as.character(opendap_modis[[2]]))
}


##############################################################
#### 3.2 - Prepare the data ####
##############################################################

# Get Terra and Aqua products to preprocess 
terra_names<-file.path(path_to_modisevapo_folder,paste0(gsub("-","",dates_modis_evapo_real),"_terra.nc4"))
aqua_names<-file.path(path_to_modisevapo_folder,paste0(gsub("-","",dates_modis_evapo_real),"_aqua.nc4"))

brick_modisevapo<-NULL

for (i in 1:length(dates_modis_evapo_real)){
  
  cat(paste0("pre-processing LST for date ",dates_modis_evapo_real[i],"\n"))
  
  # Create raster of evapotranspiration by combining Terra and Aqua and taking the average between both values
  rast_evapo_terra<-fun_preprocess_modis_product(terra_names[i],"ET_500m")
  rast_evapo_aqua<-fun_preprocess_modis_product(aqua_names[i],"ET_500m")
  
  # Set pixel values >= 32760 (quality pixel values) to NA 
  rast_evapo_terra[rast_evapo_terra >= 32760] <- NA
  rast_evapo_aqua[rast_evapo_aqua >= 32760] <- NA
  
  rast_evapo<-brick(rast_evapo_terra,rast_evapo_aqua)
  rast_evapo<-mean(rast_evapo,na.rm=TRUE)

  brick_modisevapo<-c(brick_modisevapo,rast_evapo)
  names(brick_modisevapo[[length(brick_modisevapo)]])<-paste0("evapo_",i-1)
  
}

brick_modisevapo<-brick(brick_modisevapo)

locations_hlc_sp_this_date_modis_proj <- raster::extract(brick_modisevapo, locations_hlc_sp_this_date_modis_proj, buffer=buffer_size_meters,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 


#####################################
########### 4. GPM ###############
#####################################

## Websites
# TRMM/GPM doc: https://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B43/doc/README.TRMM_V7.pdf
# TRMM/GPM products :https://pmm.nasa.gov/data-access/downloads/trmm
# FAQ: https://pmm.nasa.gov/data-access#prodrt
# TRMM real time products: https://pmm.nasa.gov/sites/default/files/document_files/3B4XRT_doc_V7_180426.pdf
# TRMM offical FTPs: https://pmm.nasa.gov/data-access/data-sources#register
# TRMM opendap server: https://disc2.gesdisc.eosdis.nasa.gov/opendap/
# Doc to download data from opendap server: https://opendap.github.io/documentation/QuickStart.html#_an_easy_way_using_the_browser_based_opendap_server_dataset_access_form
# To access data in R via opendap protocol: https://publicwiki.deltares.nl/pages/viewpage.action?pageId=91848908
# R package to handle TRMM data: https://github.com/barryrowlingson/trmm
## Earthdata: 
# How to Download Data Files from HTTPS Service with wget: https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Download%20Data%20Files%20from%20HTTPS%20Service%20with%20wget

# from https://wiki.earthdata.nasa.gov/display/EL/How+to+access+data+with+R

# The GPM data at 1°/1 day resolution are available :
# through OpenDAP protocol here (ability to subset the dataset on our bounding box before downloading it): https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.05
# through standard http protocol here (impossible to subset the dataset, ie must download the whole dataset (approx. 20 mb / dataset)): https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDF.05

## Exhaustive info and links on the GPM product (including citation): https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDF_V06/summary?keywords=3imergdf
##############################################################
#### 4.1 - Download the data ####
##############################################################

# We retrieve the data 15 days before the date_epidemio_campain (source of the duration: thèse Nico)

time_range<-as.character(c(this_date_hlc-lag_days_gpm,this_date_hlc))
dates <-seq(as.Date(time_range[1]),as.Date(time_range[2]),1)

# Function to build the OpenDap URL to dowload GPM 1km data. input parameters : date

fun_build_gpm_opendap_url<-function(date_catch){
  
  year<-format(date_catch,'%Y')
  month<-format(date_catch,'%m')
  product_name<-paste0("3B-DAY.MS.MRG.3IMERG.",gsub("-","",date_catch),"-S000000-E235959.V06.nc4.nc4")
  
  # We use the precipitationCal variable, as explained in the technical doc (https://pps.gsfc.nasa.gov/Documents/IMERG_doc_190313.pdf , p.39) : Note  well  that  HQprecipitation only includes  microwave  data  (hence  “HQ”),  meaning  it  has significant gaps.  precipitationCal is the complete estimate that most users will want to access
  url_product<-paste0(url_gpm_opendap,"/",year,"/",month,"/",product_name,"?precipitationCal[0:0][",index_opendap_gpm_lon_min,":",index_opendap_gpm_lon_max,"][",index_opendap_gpm_lat_min,":",index_opendap_gpm_lat_max,"],lon[",index_opendap_gpm_lon_min,":",index_opendap_gpm_lon_max,"],lat[",index_opendap_gpm_lat_min,":",index_opendap_gpm_lat_max,"]")
  
  return(url_product)
}

## Download the data

gpm_list_paths<-NULL
#gpm_list_names<-NULL

dates<-rev(dates)

for (i in 1:length(dates)){
  cat(paste0("Downloading GPM data for the ROI and for date ",dates[i],"\n"))
  # build the url of the dataset
  url_opendap<-fun_build_gpm_opendap_url(dates[i])
  # Download the data
  path_to_output_gpm_data<-file.path(path_to_gpm_raw_folder,paste0(gsub("-","",dates[i]),".nc4"))
  gpm_list_paths<-c(gpm_list_paths,path_to_output_gpm_data)
  #gpm_list_names<-c(gpm_list_names,paste0("gpm_",as.character(abs(difftime(dates[i] ,this_date_hlc , units = c("days"))))))
  if (!(file.exists(path_to_output_gpm_data))){
    httr::GET(url = url_opendap, write_disk(path_to_output_gpm_data))
  } else {
    cat(paste0("Data is already existing for date ",dates[i],"\n"))
  }
}



# Another working solution to DL the data :
#system(paste0("wget --http-user=",username_EarthData," --http-password=",password_EarthData," --output-document=",file.path(path_to_gpm_raw_folder,product_name)," ",url_product))

## Only works with wget by a command system... we have tried many other options (via download.file, etc) but none work. Seems that R has problems with proxy

## Open the file as raster and reproject in right EPSG : 
# from https://rpubs.com/boyerag/297592
#
#grid_nc<-nc_open("/home/ptaconet/Documents/react/data_CIV/GPM/raw_data/3B-DAY.MS.MRG.3IMERG.20160119-S000000-E235959.V05.nc4.nc4")
#grid_nc <- ncvar_get(grid_nc, "HQprecipitation") 

#grid_nc <- raster(t(grid_nc), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+init=epsg:4326"))
#grid_nc <- projectRaster(grid_nc, crs = CRS(paste0("+init=epsg:",epsg)))

# Open GPM data as raster and pre-process


##############################################################
#### 4.2 - Prepare the data ####
##############################################################

brick_gpm<-NULL
brick_gpm_positive_precip<-NULL


for (i in 1:length(gpm_list_paths)){
  
  gpm_rast<-raster(gpm_list_paths[i])
  projection(gpm_rast)<-CRS("+init=epsg:4326")
  # The raster has to be flipped. Output was validated with the data from 2017-09-20 (see https://docserver.gesdisc.eosdis.nasa.gov/public/project/GPM/browse/GPM_3IMERGDF.png)
  gpm_rast <- t(gpm_rast)
  gpm_rast <- flip(gpm_rast,'y')
  gpm_rast <- flip(gpm_rast,'x')
  gpm_rast <- projectRaster(gpm_rast, crs = CRS(paste0("+init=epsg:",epsg)))

  ### Interpolate 
  size_output_grid<-250
  r<-gpm_rast
  res(r)<-c(size_output_grid,size_output_grid)
  
  ## Using Inverse distance weighting. More info : https://pro.arcgis.com/fr/pro-app/help/analysis/geostatistical-analyst/how-inverse-distance-weighted-interpolation-works.htm
  #grd <- as(r, 'SpatialGrid')
  #P<-rasterToPoints(gpm_rast, spatial = TRUE)
  # Interpolate the grid cells using a power value of 2 (idp=2.0). 
  #P.idw <- gstat::idw(P$layer ~ 1, P, newdata=grd, idp=2.0)
  # Convert to raster object
  #gpm_rast_resample_idw <- raster(P.idw)
  
  ## Using resample (bilinear resampling, cf. http://desktop.arcgis.com/fr/arcmap/latest/extensions/spatial-analyst/performing-analysis/cell-size-and-resampling-in-analysis.htm)
  gpm_rast<-resample(gpm_rast,r,method='bilinear')
  
  brick_gpm<-c(brick_gpm,gpm_rast)
  #names(brick_gpm[[length(brick_gpm)]])<-gpm_list_names[i]
  
  gpm_positive_precip_rast<-gpm_rast
  gpm_positive_precip_rast[gpm_positive_precip_rast > 0] <- 1
  gpm_positive_precip_rast[gpm_positive_precip_rast <= 0] <- 0
  brick_gpm_positive_precip<-c(brick_gpm_positive_precip,gpm_positive_precip_rast)
}

brick_gpm<-brick(brick_gpm)

## Sum of the precipitations for the n days (on a 3 days moving window ?)
precip_sum<-sum(brick_gpm,na.rm = T)
names(precip_sum)<-"precip_sum"
## Precipitation for the date of HLC
precip_0<-brick_gpm[[1]]
names(precip_0)<-"precip_0"
## Number of days with precipitations during the n days before the HLC
precip_ndays<-sum(brick(brick_gpm_positive_precip),na.rm = T)
names(precip_ndays)<-"precip_ndays"

# Extract mean of GPM for each raster
locations_hlc_sp_this_date <- raster::extract(precip_sum, locations_hlc_sp_this_date, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 
locations_hlc_sp_this_date <- raster::extract(precip_0, locations_hlc_sp_this_date, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 
locations_hlc_sp_this_date <- raster::extract(precip_ndays, locations_hlc_sp_this_date, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 



#####################################
########### 5. Night lights ###############
#####################################

# info on data : https://ngdc.noaa.gov/eog/viirs/download_dnb_composites.html and https://noaa.maps.arcgis.com/home/item.html?id=d7c95b2da6fd43cd9dec19b212f145db

# map viewer : https://maps.ngdc.noaa.gov/viewers/VIIRS_DNB_nighttime_imagery/index.html

# data service : https://gis.ngdc.noaa.gov/arcgis/rest/services/NPP_VIIRS_DNB/Monthly_AvgRadiance/ImageServer/

# NOAA VIIRS DNB Nighttime Lights Monthly Composites are monthly products. We download the data for the month of the HCL

##############################################################
#### 5.1 - Download the data ####
##############################################################

date_start<-as.Date(paste0(substr(this_date_hlc,1,7),"-01"))
date_end<-seq(date_start, by = "1 month", length = 2)[2]

time_start<-as.integer(difftime(date_start ,"1970-01-01" , units = c("secs")))*1000
time_end<-as.integer(difftime(date_end ,"1970-01-01" , units = c("secs")))*1000

url_product<-paste0(url_noaa_nighttime_webservice,"?bbox=",bbox_roi[1,1],",",bbox_roi[2,1],",",bbox_roi[1,2],",",bbox_roi[2,2],"&time=",time_start,",",time_end,"&format=tiff&f=image")
path_to_output_nighttime<-file.path(path_to_nighttime_folder,paste0(gsub("-","_",date_start),".tif"))

url_product_cloudcover<-gsub("Monthly_AvgRadiance","Monthly_CloudFreeCoverage",url_product)
path_to_output_nighttime_cloudcover<-file.path(path_to_nighttime_folder,paste0(gsub("-","_",date_start),"_cloudfreecov.tif"))

if (!file.exists(path_to_output_nighttime)){
  cat(paste0("Downloading Night lights data for the ROI and for month ",date_start,"\n"))
  GET(url_product,write_disk(path_to_output_nighttime))
  GET(url_product_cloudcover,write_disk(path_to_output_nighttime_cloudcover))
} else {
  cat(paste0("Data is already existing for month ",date_start,"\n"))
}


##############################################################
#### 5.2 - Prepare the data ####
##############################################################

# Open the data 
nighttime_rast<-raster(path_to_output_nighttime)
nighttime_cloudcover<-raster(path_to_output_nighttime_cloudcover)

# Quality control : exclude pixels with 0 cloud-free observations
nighttime_cloudcover[nighttime_cloudcover==0]<-NA
nighttime_rast <- mask(nighttime_rast, nighttime_cloudcover)

# Extract mean of DNB radiance. However change the buffer size (500 m) ? 
names(nighttime_rast)<-"nightligth_mean"
locations_hlc_sp_this_date <- raster::extract(nighttime_rast, locations_hlc_sp_this_date, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 


#####################################
########### 5. Wind (ERA) ###############
#####################################
# Check : https://dominicroye.github.io/en/2018/access-to-climate-reanalysis-data-from-r/

# ERA 5 : https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

# Description of the wind data: https://apps.ecmwf.int/codes/grib/param-db?id=165 and https://apps.ecmwf.int/codes/grib/param-db?id=166

##############################################################
#### 5.1 - Download the data ####
##############################################################

bbox_4326<-bbox(spTransform(dates_locations_hlc_sp,CRS("+init=epsg:4326")))
# extend a bit the size of the bbox
bbox_4326[,1]=bbox_4326[,1]-1
bbox_4326[,2]=bbox_4326[,2]+1

path_to_output_erawind_data_this_date<-file.path(path_to_erawind_folder,paste0(gsub("-","_",this_date_hlc),".nc"))
path_to_output_erawind_data_this_date_plus_one<-file.path(path_to_erawind_folder,paste0(gsub("-","_",this_date_hlc+1),".nc"))


#we create the query for the date of catch, hours 18h to 23h
query_this_date_hlc <- r_to_py(list(
  variable= c("10m_u_component_of_wind","10m_v_component_of_wind"),
  product_type= "reanalysis",
  year= format(this_date_hlc, "%Y"),
  month= format(this_date_hlc, "%m"), #formato: "01","01", etc.
  day= format(this_date_hlc, "%d"), #stringr::str_pad(1:31,2,"left","0"),   
  time= stringr::str_c(18:23,"00",sep=":")%>%str_pad(5,"left","0"),
  format= "netcdf",
  area = paste0(bbox_4326[2,2],"/",bbox_4326[1,1],"/",bbox_4326[2,1],"/",bbox_4326[1,2]) # North, West, South, East
))

#we create the query for the date of catch + 1, hours 00h to 08h
query_this_date_hlc_plus_one <- r_to_py(list(
  variable= c("10m_u_component_of_wind","10m_v_component_of_wind"),
  product_type= "reanalysis",
  year= format(this_date_hlc+1, "%Y"),
  month= format(this_date_hlc+1, "%m"), 
  day= format(this_date_hlc+1, "%d"), 
  time= stringr::str_c(00:08,"00",sep=":")%>%str_pad(5,"left","0"),
  format= "netcdf",
  area = paste0(bbox_4326[2,2],"/",bbox_4326[1,1],"/",bbox_4326[2,1],"/",bbox_4326[1,2])
))

#query the server to get the ncdf for date of catch and date of catch + 1
server$retrieve("reanalysis-era5-single-levels",
                query_this_date_hlc,
                path_to_output_erawind_data_this_date)

server$retrieve("reanalysis-era5-single-levels",
                query_this_date_hlc_plus_one,
                path_to_output_erawind_data_this_date_plus_one)

##############################################################
#### 5.2 - Prepare the data ####
##############################################################

# open path_to_output_erawind_data_this_date and preprocess each component (open each layer as raster (ie each hour), reproject to right epsg, resample, crop to ROI extent)
# open path_to_output_erawind_data_this_date and preprocess (open each layer as raster (ie each hour), reproject to right epsg, resample, crop to ROI extent)
# make a raster brick out of all the pre-processed rasters of path_to_output_erawind_data_this_date
# make a raster brick out of all the pre-processed rasters of path_to_output_erawind_data_this_date


brick_u10_this_date<-fun_preprocess_era_product(path_to_output_erawind_data_this_date,"u10",roi_sp_utm,epsg,size_output_grid=250)
brick_u10_this_date_plus_one<-fun_preprocess_era_product(path_to_output_erawind_data_this_date_plus_one,"u10",roi_sp_utm,epsg,size_output_grid=250)
brick_u10<-c(brick_u10_this_date,brick_u10_this_date_plus_one)
brick_u10<-brick(brick_u10)
wind_u10_mean<-mean(brick_u10,na.rm=TRUE)

brick_v10_this_date<-fun_preprocess_era_product(path_to_output_erawind_data_this_date,"v10",roi_sp_utm,epsg,size_output_grid=250)
brick_v10_this_date_plus_one<-fun_preprocess_era_product(path_to_output_erawind_data_this_date_plus_one,"v10",roi_sp_utm,epsg,size_output_grid=250)
brick_v10<-c(brick_v10_this_date,brick_v10_this_date_plus_one)
brick_v10<-brick(brick_v10)
wind_v10_mean<-mean(brick_v10,na.rm=TRUE)

# Wind speed is calculated by sqrt(u^2+v^2)
wind_speed<-sqrt(wind_u10_mean^2+wind_v10_mean^2)
names(wind_speed)="wind_speed"

# Extract wind speed
locations_hlc_sp_this_date <- raster::extract(wind_speed, locations_hlc_sp_this_date, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 


#################################################################
########### 6. Ephemeris of the Moon ###############
#################################################################

##############################################################
#### 6.1 - Download the data ####
##############################################################

## More info on IMCCE web services at http://vo.imcce.fr/webservices/miriade/?ephemcc

# Set output path
path_to_output_imcce_data<-file.path(path_to_imcce_folder,paste0(gsub("-","_",this_date_hlc),".csv"))

# Create URL to send as a web service
observer<-paste0(mean(locations_hlc_sp_this_date$longitude),"%20",mean(locations_hlc_sp_this_date$latitude),"%200.0") # %20 is to encode the space in the URL

url_imcce_moon<-paste0(url_imcce_webservice,"-name=s:Moon&-type=Satellite&-ep=",this_date_hlc,"T23:30:00&-nbd=1d&-step=1h&-tscale=UTC&-observer=",observer,"&-theory=INPOP&-teph=1&-tcoor=1&-mime=text/csv&-output=--jd&-extrap=0&-from=MiriadeDoc")

# Download the data
if (!(file.exists(path_to_output_imcce_data))){
  httr::GET(url = url_imcce_moon, write_disk(path_to_output_imcce_data))
} else {
  cat(paste0("Data is already existing for date ",this_date_hlc,"\n"))
}

##############################################################
#### 6.2 - Prepare the data ####
##############################################################

# Open the data
moon_magnitude<-read.csv(path_to_output_imcce_data,skip=10)
colnames(moon_magnitude)<-gsub("\\.","_",colnames(moon_magnitude))

# Add moon visual magnitude to data frame
locations_hlc_sp_this_date$moon_vmag <- moon_magnitude$V_Mag



########### Other dynamic data: 
# - MODIS yearly Land cover (https://cmr.earthdata.nasa.gov/search/concepts/C186286578-LPDAAC_ECS) (https://opendap.cr.usgs.gov/opendap/hyrax/MCD12Q1.006/contents.html)
# - wet zones (extracted from S1)



############ Static data: 
# - DEM and DEM-derivatives
# - land use / land cover (my data) . Check tutos at : https://www.r-craft.org/r-news/efficient-landscape-metrics-calculations-for-buffers-around-sampling-points/  and https://cran.rstudio.com/web/packages/landscapemetrics/vignettes/getstarted.html
# - built-up areas (facebook data) -> apply the same processings as for land use/land cover, 
# - population (either fb data or react data)
# - pedology
# - pistes -> comment on les choppe ? Uniquement les routes principales ou aussi les petits sentiers ? 

# - night ligths : either dynamic (VNP46 products) or static (NASA 2013 dataset)




###############################
####### Sentinel 1 ########
###############################

### Questions à répondre : Quels fichiers telecharger parmi tous ceux disponibles ? A voir en fonction de la méthodo pour calculer les données d'intérêt 

## Sentinel 1: approx. 6 days of revisit time. We progressively extend the time around the date until we find a valid record

i=0
records=NULL
while (is.null(records)){
  i=i+1
  timeframe_S1<-c(as.Date(date_epidemio_campain)-i,as.Date(date_epidemio_campain)+i)
  records <- getSentinel_query(time_range = as.character(timeframe_S1), platform = "Sentinel-1", verbose = F)
}

## Retrieve the dataset to download
# Cut by the AOI
# TODO

# Preview a single record
#getSentinel_preview(record = records[1,],show_aoi = TRUE)

## Download some datasets
#datasets <- getSentinel_data(records = records[c(1),],dir_out=path_to_sentinel1_raw_folder)




########################################################################################################################
############ Close the workflow ############
########################################################################################################################

# Remove grass temporary folder
system(paste0("rm -r ", file.path(getwd(),"GRASS_TEMP")))  
file.remove(file.path(getwd(),".grassrc7"))

# Close the connection to the database
dbDisconnect(react_db)

#### stop cluster ####
#stopCluster(cl)

########################################################################################################################
############ The end ############
########################################################################################################################
cat("End workflow")


