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




#### Start Workflow
########################################################################################################################
############ Set Input parameters for the workflow ############
########################################################################################################################

### Global variables used throughout the WF
path_to_processing_folder<-"/home/ptaconet/Documents/react/data_CIV"  #<Path to the processing folder (i.e. where all the data produced by the workflow will be stored)>
path_to_roi_vector="ROI.kml" #<Path to the Region of interest in KML format>
path_to_database<-"/home/ptaconet/Bureau/react.db"  #<Data frame containing 
epsg=32630
modis_tile<-"h17v08" # For CIV. For BF : modis_tile<-"h17v07" 
url_query_hlc_dates_location<-"https://raw.githubusercontent.com/ptaconet/r_react/master/db_sql/hlc_dates_location.sql"

## Credential to the various servers where the source data are stored
username_cophub<-"ptaconet" #<ESA Copernicus Scihub username>
password_cophub<-"HHKcue51" #<ESA Copernicus Scihub password>
username_USGS<-"ptaconet" #<USGS username>
password_USGS<-"HHKcue51" #<USGS password>
username_EarthData<-"ptaconet"  #<EarthData username>
password_EarthData<-"HHKcue51"  #<EarthData password>
  



## Buffer size, within which the raster statistics will be computed (radius in meters)
buffer_size=2000


########################################################################################################################
############ Prepare workflow ############
########################################################################################################################

### Call useful libraries
library(getSpatialData)
library(raster)
library(sf)
library(sp)
library(gdalUtils)
library(rgdal)
require(httr)
require(RSQLite)
#require(gapfill)
require(dplyr)
require(ncdf4)
require(reticulate)
require(lubridate)
## Load packages for working on multi-core
library(parallel)
library(doParallel)
library(foreach)

## Set working directory
setwd(path_to_processing_folder)

### Set the paths of output folders / files
path_to_sentinel1_products<-file.path(path_to_processing_folder,"Sentinel_1")
path_to_sentinel1_raw_folder<-file.path(path_to_sentinel1_products,"raw_data")
path_to_sentinel1_processed_folder<-file.path(path_to_sentinel1_products,"processed_data")
path_to_modislst_folder<-file.path(path_to_processing_folder,"MODIS_LST")
path_to_modislst_raw_folder<-file.path(path_to_modislst_folder,"raw_data")
path_to_modislst_processed_folder<-file.path(path_to_modislst_folder,"processed_data")
path_to_modisveget_folder<-file.path(path_to_processing_folder,"MODIS_veget")
path_to_modisveget_raw_folder<-file.path(path_to_modisveget_folder,"raw_data")
path_to_modisveget_processed_folder<-file.path(path_to_modisveget_folder,"processed_data")
path_to_gpm_folder<-file.path(path_to_processing_folder,"GPM")
path_to_gpm_raw_folder<-file.path(path_to_gpm_folder,"raw_data")
path_to_gpm_processed_folder<-file.path(path_to_gpm_folder,"processed_data")
path_to_erawind_folder<-file.path(path_to_processing_folder,"ERA_WIND")
path_to_erawind_raw_folder<-file.path(path_to_erawind_folder,"raw_data")
path_to_erawind_processed_folder<-file.path(path_to_erawind_folder,"processed_data")

## Create the output folders
directories<-list(path_to_sentinel1_products,path_to_sentinel1_raw_folder,path_to_sentinel1_processed_folder,path_to_modislst_folder,path_to_modislst_raw_folder,path_to_modislst_processed_folder,path_to_gpm_folder,path_to_gpm_raw_folder,path_to_gpm_processed_folder,path_to_erawind_folder,path_to_erawind_raw_folder,path_to_erawind_processed_folder,path_to_modisveget_folder,path_to_modisveget_raw_folder,path_to_modisveget_processed_folder)
lapply(directories, dir.create)

## Set connections for getSpatialData
getSpatialData::login_CopHub(username = username_cophub, password = password_cophub)
getSpatialData::login_USGS(username = username_USGS, password = password_USGS)

## Set urls to the various OpenDAP servers and connect
url_gpm_opendap<-"https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.06"
url_modis_opendap<-"https://opendap.cr.usgs.gov/opendap/hyrax"
modis_lst_terra_product<-"MOD11A1.006"
modis_lst_aqua_product<-"MYD11A1.006"
modis_veget_terra_product<-"MOD13Q1.006"
modis_veget_aqua_product<-"MYD13Q1.006"
modis_crs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

httr::set_config(authenticate(user=username_EarthData, password=password_EarthData, type = "basic"))

fun_get_opendap_index<-function(path_to_opendap_index){
  
  opendap_indexes<-read.csv(path_to_opendap_index,skip = 1)
  opendap_indexes[1]<-NULL
  opendap_indexes<-colnames(opendap_indexes)
  opendap_indexes<-gsub("X\\.","-",opendap_indexes)
  opendap_indexes<-gsub("X","",opendap_indexes)
  opendap_indexes<-as.numeric(opendap_indexes)
  
  return(opendap_indexes)
}


fun_preprocess_modis_product<-function(path_to_raw_modis,var_name){
  grid_nc<-raster(path_to_raw_modis,varname=var_name)
  projection(grid_nc)<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
  #extent(grid_nc)[1:2]<-extent(grid_nc)[1:2]+res(grid_nc)[1]/2
  #extent(grid_nc)[3:4]<-extent(grid_nc)[3:4]-res(grid_nc)[1]/2
  grid_nc <- projectRaster(grid_nc, crs = CRS(paste0("+init=epsg:",epsg)))
  return(grid_nc)
}

# initiate cluster for paralell download 
#no_cores <- detectCores() - 1
#cl <- makeCluster(no_cores, type = "PSOCK")
#registerDoParallel(cl)

## Set ROI as sf object
roi_sf <- read_sf(path_to_roi_vector)$geometry
set_aoi(roi_sf)

## Set ROI as sp SpatialPolygon object
roi_sp_4326<-rgdal::readOGR(path_to_roi_vector)
roi_sp_32630 <- spTransform(roi_sp, CRS(paste0("+init=epsg:",epsg)))
roi_sp_modis_project <- spTransform(roi_sp, modis_crs)

## Connect to DB
react_db <- dbConnect(RSQLite::SQLite(),path_to_database)

## Get query to retrieve dates and locations of human landing catches
sql_query_hlc_dates_location<-paste(readLines(url_query_hlc_dates_location), collapse="\n")

## Get data frame of dates and locations of human landing catches
df_dates_locations_hlc<-dbGetQuery(react_db, sql_query_hlc_dates_location)
## Turn to spatial dataset (sf and sp)
dates_locations_hlc_sf<-st_as_sf(df_dates_locations_hlc, coords = c("longitude", "latitude"), crs = 4326)
dates_locations_hlc_sf<-st_transform(dates_locations_hlc_sf, crs = epsg)
dates_locations_hlc_sp<-SpatialPointsDataFrame(coords=data.frame(df_dates_locations_hlc$longitude,df_dates_locations_hlc$latitude),data=df_dates_locations_hlc,proj4string=CRS("+init=epsg:4326"))
dates_locations_hlc_sp<-spTransform(dates_locations_hlc_sp,CRS(paste0("+init=epsg:",epsg)))


## Get all dates
all_dates_hlc<-unique(df_dates_locations_hlc$datecapture)

## For tests: retrieve all the lines for the first date
i=1
locations_hlc_sp_this_date<-dates_locations_hlc_sp[which(dates_locations_hlc_sp$datecapture==all_dates_hlc[i]),]
this_date_hlc<-as.Date(all_dates_hlc[i])

############# make menage file a geo dataset with positions #######

## Set date 
#date_epidemio_campain <- "2016-01-20"  # Date d'une capture au village DON au BF

#menages_csv<-read.csv(path_to_file_coord_menages)
#menages_csv<-menages_csv %>% filter (!(is.na(Part12Part2coordgpsLatitude))) %>% filter (!(is.na(Part12Part2coordgpsLongitude)))

#menages_sf<-st_as_sf(menages_csv, coords = c("Part12Part2coordgpsLongitude", "Part12Part2coordgpsLatitude"), crs = 4326)
#menages_sf<-st_transform(menages_sf, crs = 32630)

#menages_sp_4326<-SpatialPointsDataFrame(coords=data.frame(menages_csv$Part12Part2coordgpsLongitude,menages_csv$Part12Part2coordgpsLatitude),data=menages_csv,proj4string=CRS("+init=epsg:4326"))
#menages_sp<-spTransform(menages_sp,CRS("+init=epsg:32630"))

## Create 2 km buffer around the points
#buffer_2km <- st_buffer(menages_sf,2000) 
#buffer_2km<-st_union(buffer_2km,by_feature = TRUE)

## Get couple {date,coord}
#dates_loc_enquetes<-unique(menages_sf[c("dateenq", "geometry")])


########################################################################################################################
############ Start Workflow ############
########################################################################################################################
cat("Starting workflow")


########################################################################################################################
############ Download and prepare OpenDap time and spatial indexes to further download the MODIS and GPM datasets on opendap servers ############
########################################################################################################################

url_opendap_modis_lst_terra<-paste0(url_modis_opendap,"/",modis_lst_terra_product,"/",modis_tile)
url_opendap_modis_lst_aqua<-paste0(url_modis_opendap,"/",modis_lst_aqua_product,"/",modis_tile)
url_opendap_modis_veget_terra<-paste0(url_modis_opendap,"/",modis_veget_terra_product,"/",modis_tile)
url_opendap_modis_veget_aqua<-paste0(url_modis_opendap,"/",modis_veget_aqua_product,"/",modis_tile)

## Downlaad the time indexes on opendap for each opendap poduct
httr::GET(paste0(url_opendap_modis_lst_terra,".ncml.ascii?time"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_terra_time_index.txt")))
httr::GET(paste0(url_opendap_modis_lst_aqua,".ncml.ascii?time"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_aqua_time_index.txt")))
httr::GET(paste0(url_opendap_modis_veget_terra,".ncml.ascii?time"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_terra_time_index.txt")))
httr::GET(paste0(url_opendap_modis_veget_aqua,".ncml.ascii?time"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_aqua_time_index.txt")))

## Get the time indexes as vectors
opendap_modis_lst_terra_time_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_terra_time_index.txt"))
opendap_modis_lst_aqua_time_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_aqua_time_index.txt"))
opendap_modis_veget_terra_time_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_terra_time_index.txt"))
opendap_modis_veget_aqua_time_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_aqua_time_index.txt"))


## Downlaad the spatial indexes on opendap for each opendap poduct
httr::GET(paste0(url_opendap_modis_lst_terra,".ncml.ascii?XDim"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_XDim_index.txt")))
httr::GET(paste0(url_opendap_modis_lst_terra,".ncml.ascii?YDim"),write_disk(file.path(path_to_modislst_folder,"opendap_modis_lst_YDim_index.txt")))

httr::GET(paste0(url_opendap_modis_veget_terra,".ncml.ascii?XDim"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_XDim_index.txt")))
httr::GET(paste0(url_opendap_modis_veget_terra,".ncml.ascii?YDim"),write_disk(file.path(path_to_modisveget_folder,"opendap_modis_veget_YDim_index.txt")))

httr::GET(paste0(url_gpm_opendap,"/2018/02/3B-DAY.MS.MRG.3IMERG.20180201-S000000-E235959.V06.nc4.ascii?lon"),write_disk(file.path(path_to_gpm_folder,"opendap_gpm_lon_index.txt")))
httr::GET(paste0(url_gpm_opendap,"/2018/02/3B-DAY.MS.MRG.3IMERG.20180201-S000000-E235959.V06.nc4.ascii?lat"),write_disk(file.path(path_to_gpm_folder,"opendap_gpm_lat_index.txt")))


## Get the spatial indexes as vectors
opendap_modis_lst_XDim_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_XDim_index.txt"))
opendap_modis_lst_YDim_index<-fun_get_opendap_index(file.path(path_to_modislst_folder,"opendap_modis_lst_YDim_index.txt"))
opendap_modis_veget_XDim_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_XDim_index.txt"))
opendap_modis_veget_YDim_index<-fun_get_opendap_index(file.path(path_to_modisveget_folder,"opendap_modis_veget_YDim_index.txt"))
opendap_gpm_lon_index<-fun_get_opendap_index(file.path(path_to_gpm_folder,"opendap_gpm_lon_index.txt"))
opendap_gpm_lat_index<-fun_get_opendap_index(file.path(path_to_gpm_folder,"opendap_gpm_lat_index.txt"))


## Extract indexes for lat min, lat max, lon min and lon max for our bounding box

roi_bbox_4326<-bbox(roi_sp_4326) ## for MODIS 
roi_bbox_modisproject<-bbox(roi_sp_modis_project) # for GPM

index_opendap_gpm_lon_min<-which.min(abs(opendap_gpm_lon_index-roi_bbox_4326[1,1]))-1
index_opendap_gpm_lon_max<-which.min(abs(opendap_gpm_lon_index-roi_bbox_4326[1,2]))-1 
index_opendap_gpm_lat_min<-which.min(abs(opendap_gpm_lat_index-roi_bbox_4326[2,1]))-1
index_opendap_gpm_lat_max<-which.min(abs(opendap_gpm_lat_index-roi_bbox_4326[2,2]))-1 

index_opendap_modisveget_lon_min<-which.min(abs(opendap_modis_veget_XDim_index-roi_bbox_modisproject[1,1]))-1
index_opendap_modisveget_lon_max<-which.min(abs(opendap_modis_veget_XDim_index-roi_bbox_modisproject[1,2]))-1 
index_opendap_modisveget_lat_max<-which.min(abs(opendap_modis_veget_YDim_index-roi_bbox_modisproject[2,1]))-1
index_opendap_modisveget_lat_min<-which.min(abs(opendap_modis_veget_YDim_index-roi_bbox_modisproject[2,2]))-1 

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
  
#####################################
########### Modis LST ###############
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
##############################################################
#### 2.1 - Download data ####
##############################################################

time_range<-as.character(c(this_date_hlc-16,this_date_hlc))
dates <-seq(as.Date(time_range[1]),as.Date(time_range[2]),1)

fun_build_modis_lst_opendap_url<-function(date_hlc,satellite){
  
  index_date<-as.integer(difftime(date_hlc ,"2000-01-01" , units = c("days"))) # time is provided as number of days since 2000-01-01
  
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
for (i in 1:length(dates)){
  cat(paste0("Downloadind MODIS LST Terra and Aqua (MOD11A1 and MYD11A1) for the ROI and for date ",dates[i],"\n"))
  for (j in c("terra","aqua")){
    opendap_modis<-fun_build_modis_lst_opendap_url(dates[i],j)
    path_to_dataset<-file.path(path_to_modislst_raw_folder,paste0(gsub("-","",dates[i]),"_",j,".nc4"))
    if (!(file.exists(path_to_dataset))){
      httr::GET(opendap_modis,write_disk(path_to_dataset))
    }
  }
}



##############################################################
#### 2.2 - Prepare data ####
##############################################################

## For each date before the d date:
# - For each cell, extract the maximum day available temperature among Terra and Aqua. In case of NA on one of the two datasets, take the value of the other dataset.
# - For each cell, extract the minimum night available temperature among Terra and Aqua. In case of NA on one of the two datasets, take the value of the other dataset.
# - Extract the mean of the available LST cells on a buffer around each catch point
# - If no data are available for a given catch point (i.e. cells are all NA), set temperature of the previous night (recursively until a value is available)
# - If no data are available for the whole time series, ??? TODO
# - Create thermal amplitude value (i.e. difference between day and night)


# Get Terra and Aqua products to preprocess 
terra_names<-file.path(path_to_modislst_raw_folder,paste0(gsub("-","",dates),"_terra.nc4"))
aqua_names<-file.path(path_to_modislst_raw_folder,paste0(gsub("-","",dates),"_aqua.nc4"))

dates<-rev(dates)

list_lst_day<-list()
list_lst_nigth<-list()

for (i in 1:length(dates)){

  cat(paste0("pre-processing LST for date ",dates[i],"\n"))

# Create raster of LST_day maximum value combining Terra and Aqua by taking the maximum available value for each pixel
rast_lst_day_terra<-fun_preprocess_modis_product(terra_names[i],"LST_Day_1km")
rast_lst_day_aqua<-fun_preprocess_modis_product(aqua_names[i],"LST_Day_1km")
rast_lst_day<-brick(rast_lst_day_terra,rast_lst_day_aqua)

rast_max_lst_day<-max(rast_lst_day,na.rm=TRUE)
list_lst_day[[length(list_lst_day)+1]]<-rast_max_lst_day
names(list_lst_day[[length(list_lst_day)]])<-paste0("lst_max_d_",i)

# Create raster of LST_night minimum value combining Terra and Aqua by taking the minimum available value for each pixel
rast_lst_night_terra<-fun_preprocess_modis_product(terra_names[i],"LST_Night_1km")
rast_lst_night_aqua<-fun_preprocess_modis_product(aqua_names[i],"LST_Night_1km")
rast_lst_night<-brick(rast_lst_night_terra,rast_lst_night_aqua)

rast_min_lst_night<-min(rast_lst_night,na.rm=TRUE)
list_lst_nigth[[length(list_lst_nigth)+1]]<-rast_min_lst_night
names(list_lst_nigth[[length(list_lst_nigth)]])<-paste0("lst_min_d_",i)

}


lst_day<-brick(list_lst_day)
lst_night<-brick(list_lst_nigth)

## Extract mean of LST in a buffer of 2 km for all the dates
# na.rm = TRUE means that NA values in the raster are not taken into account
# small=TRUE means that each time the buffer intersects a non NA cell it takes into account the cell value
ex <- raster::extract(lst_day, locations_hlc_sp_this_date, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 
ex <- raster::extract(lst_night, ex, buffer=buffer_size,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 



#####################################
########### Modis NDVI and EVI ###############
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


time_range<-as.character(c(this_date_hlc-20,this_date_hlc+8))
dates <-seq(as.Date(time_range[1]),as.Date(time_range[2]),1)

#records_terra <- getMODIS_query(time_range = time_range, name = "MODIS_MOD13Q1_V6")
#records_aqua <- getMODIS_query(time_range = time_range, name = "MODIS_MYD13Q1_V6")
#records<-rbind(records_terra,records_aqua)

## Retrieve date and satellite of the week including the catch (i.e. closest start aquisition date from the catch date) and build the link to download the data

fun_build_modis_veget_opendap_url<-function(date_catch){
  
  date_catch_julian<-as.integer(difftime(date_catch ,"2000-01-01" , units = c("days")))
  
  index_opendap_terra_closest_to_date<-which.min(abs(opendap_modis_veget_terra_time_index-date_catch_julian))
  index_opendap_aqua_closest_to_date<-which.min(abs(opendap_modis_veget_aqua_time_index-date_catch_julian))
  
  days_sep_terra_from_date<-abs(opendap_modis_veget_terra_time_index[index_opendap_terra_closest_to_date]-date_catch_julian)
  days_sep_aqua_from_date<-abs(opendap_modis_veget_aqua_time_index[index_opendap_aqua_closest_to_date]-date_catch_julian)
  
  if(days_sep_terra_from_date<days_sep_aqua_from_date){
    opendap_collection_closest_date_to_catch<-"MOD13Q1.006"
    opendap_index_closest_date_to_catch<-index_opendap_terra_closest_to_date-1
  } else {
    opendap_collection_closest_date_to_catch<-"MYD13Q1.006"
    opendap_index_closest_date_to_catch<-index_opendap_aqua_closest_to_date-1
  }
  
  url_opendap<-paste0(url_modis_opendap,"/",opendap_collection_closest_date_to_catch,"/",modis_tile,".ncml.nc4?MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection,_250m_16_days_NDVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],_250m_16_days_EVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"]")
  
  return(list(url_opendap,opendap_collection_closest_date_to_catch,opendap_index_closest_date_to_catch))
  
}


#time_range<-as.character(c(this_date_hlc-20,this_date_hlc+8))
#records_terra <- getMODIS_query(time_range = time_range, name = "MODIS_MOD13Q1_V6")
#records_aqua <- getMODIS_query(time_range = time_range, name = "MODIS_MYD13Q1_V6")
#records<-rbind(records_terra,records_aqua)


#fun_build_modis_veget_opendap_url<-function(date_catch,records){

#index_closest_date_to_catch<-which(abs(as.Date(records$acquisitionDate)-this_date_hlc) == min(abs(as.Date(records$acquisitionDate) - this_date_hlc)))
#closest_date_to_catch<-min(records$acquisitionDate[index_closest_date_to_catch])
#sat_closest_date_to_catch<-records$product[which(records$acquisitionDate==closest_date_to_catch)]

#opendap_index_closest_date_to_catch<-as.integer(difftime(closest_date_to_catch ,"2000-01-01" , units = c("days")))
#opendap_index_closest_date_to_catch<-trunc(opendap_index_closest_date_to_catch/16)

#if (grepl("MOD13Q1",sat_closest_date_to_catch)){
#  opendap_collection_closest_date_to_catch<-"MOD13Q1.006"
#} else {
#  opendap_collection_closest_date_to_catch<-"MYD13Q1.006"
#}

#url_opendap<-paste0(url_modis_opendap,opendap_collection_closest_date_to_catch,"/",modis_tile,".ncml.nc4?MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection,_250m_16_days_NDVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],_250m_16_days_EVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"]")

#return(list(url_opendap,opendap_collection_closest_date_to_catch,opendap_index_closest_date_to_catch))

#}


catch_week<-fun_build_modis_veget_opendap_url(this_date_hlc,records)
one_week_before<-fun_build_modis_veget_opendap_url(as.character(this_date_hlc-8),records)
two_week_before<-fun_build_modis_veget_opendap_url(as.character(this_date_hlc-16),records)

path_catch_week<-file.path(path_to_modisveget_raw_folder,paste0(catch_week[[2]],"_",catch_week[[3]],".nc4"))
path_one_week_before<-file.path(path_to_modisveget_raw_folder,paste0(one_week_before[[2]],"_",one_week_before[[3]],".nc4"))
path_two_week_before<-file.path(path_to_modisveget_raw_folder,paste0(two_week_before[[2]],"_",two_week_before[[3]],".nc4"))

if (!(file.exists(path_catch_week))){
  httr::GET(url = catch_week[[1]], write_disk(path_catch_week))
}
if (!(file.exists(path_one_week_before))){
  httr::GET(url = one_week_before[[1]], write_disk(path_one_week_before))
}
if (!(file.exists(path_two_week_before))){
  httr::GET(url = two_week_before[[1]], write_disk(path_two_week_before))
}




rast_catch_week_ndvi<-fun_preprocess_modis_vegetation(path_catch_week,"_250m_16_days_NDVI")
rast_catch_week_evi<-fun_preprocess_modis_vegetation(path_catch_week,"_250m_16_days_EVI")

rast_catch_one_week_before_ndvi<-fun_preprocess_modis_vegetation(path_one_week_before,"_250m_16_days_NDVI")
rast_catch_one_week_before_evi<-fun_preprocess_modis_vegetation(path_one_week_before,"_250m_16_days_EVI")

rast_catch_two_week_before_ndvi<-fun_preprocess_modis_vegetation(path_two_week_before,"_250m_16_days_NDVI")
rast_catch_two_week_before_evi<-fun_preprocess_modis_vegetation(path_two_week_before,"_250m_16_days_EVI")


#grid_nc <- projectRaster(grid_nc, crs = CRS(paste0("+init=epsg:",epsg)))

#writeRaster(grid_nc,"/home/ptaconet/Bureau/test4.tif",overwrite=TRUE)




#datasets <- getMODIS_data(records = records[1,],dir_out=path_to_modisveget_raw_folder, force=FALSE) 



#####################################
########### GPM ###############
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


# The GPM data at 1°/1 day resolution are available :
# through OpenDAP protocol here (ability to subset the dataset on our bounding box before downloading it): https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.05
# through standard http protocol here (impossible to subset the dataset, ie must download the whole dataset (approx. 20 mb / dataset)): https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDF.05


# We retrieve the data 15 days before the date_epidemio_campain (source of the duration: thèse Nico)

time_range<-as.character(c(this_date_hlc-15,this_date_hlc))
dates <-seq(as.Date(time_range[1]),as.Date(time_range[2]),1)


for (i in 1:length(dates)){
  cat(paste0("Downloadind GPM data for the ROI and for date ",dates[i],"\n"))
  # build the url of the dataset
  year<-format(dates[i],'%Y')
  month<-format(dates[i],'%m')
  product_name<-paste0("3B-DAY.MS.MRG.3IMERG.",gsub("-","",dates[i]),"-S000000-E235959.V06.nc4.nc4")

  url_product<-paste0(url_gpm_opendap,"/",year,"/",month,"/",product_name,"?HQprecipitation[0:0][",index_opendap_gpm_lon_min,":",index_opendap_gpm_lon_max,"][",index_opendap_gpm_lat_min,":",index_opendap_gpm_lat_max,"],lon[",index_opendap_gpm_lon_min,":",index_opendap_gpm_lon_max,"],lat[",index_opendap_gpm_lat_min,":",index_opendap_gpm_lat_max,"]")

  # Download the data
  # from https://wiki.earthdata.nasa.gov/display/EL/How+to+access+data+with+R
  path_to_output_gpm_data<-file.path(path_to_gpm_raw_folder,paste0(gsub("-","",dates[i]),".nc4"))
  if (!(file.exists(path_to_output_gpm_data))){
    httr::GET(url = url_product, write_disk(path_to_output_gpm_data))
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

grid_nc<-raster(path_to_output_gpm_data)
projection(grid_nc)<-CRS("+init=epsg:4326")
grid_nc <- t(grid_nc)
grid_nc <- flip(grid_nc,'y')
grid_nc <- flip(grid_nc,'x')
grid_nc <- projectRaster(grid_nc, crs = CRS(paste0("+init=epsg:",epsg)))


#####################################
########### Wind (ERA) ###############
#####################################
# Check :
# https://dominicroye.github.io/en/2018/access-to-climate-reanalysis-data-from-r/

### Initialisation
# Set virtual env and install the python ECMWF API
reticulate::py_install("ecmwf-api-client") #(from https://community.rstudio.com/t/problem-installing-python-libraries-with-reticulate-py-install-error-pip-not-found/26561/2)
system("virtualenv -p /usr/bin/python2 /home/ptaconet/.virtualenvs/py2-virtualenv")
# Also works with this
#virtualenv_create("py3-virtualenv", python = "/usr/bin/python3")
reticulate::use_virtualenv("py2-virtualenv")
# install the python ECMWF API
reticulate::py_install("ecmwf-api-client", envname = "py2-virtualenv")

#import the python library ecmwfapi
ecmwf <- import::import('ecmwfapi')

#for this step there must exist the file .ecmwfapirc
server = ecmwf$ECMWFDataServer() #start the connection

bbox_4326<-bbox(menages_sp_4326)
# extend a bit the size of the bbox
bbox_4326[,1]=bbox_4326[,1]-1
bbox_4326[,2]=bbox_4326[,2]+1


#we create the query
query <-r_to_py(list(
  class='ei',
  dataset= "interim", #dataset
  date= "2016-01-20/to/2016-01-31", #time period
  expver= "1",
  grid= "0.3/0.3", #resolution
  levtype="sfc",
  param= "165.128/166.128", # wind vertical and horizontal components
  area=paste0(bbox_4326[2,2],"/",bbox_4326[1,1],"/",bbox_4326[2,1],"/",bbox_4326[1,2]), #N/W/S/E
  step= "0",
  stream="oper",
  time="00:00:00/06:00:00/12:00:00/18:00:00", #hours
  type="an",
  format= "netcdf", #format
  target='ERA_WIND/raw_data/ta2017.nc' #file name
))

#query to get the ncdf
server$retrieve(query)


nc <- nc_open("ta2017.nc")
lat <- ncvar_get(nc,'latitude')



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

#### stop cluster ####
stopCluster(cl)

########################################################################################################################
############ The end ############
########################################################################################################################
cat("End workflow")


