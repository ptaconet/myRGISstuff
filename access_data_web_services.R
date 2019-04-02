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
path_to_file_coord_menages<-"menage_coords.csv"  #<Data frame containing 
epsg="+init=epsg:32630"


## Credential to the various servers where the source data are stored
username_cophub<-"ptaconet" #<ESA Copernicus Scihub username>
password_cophub<-"HHKcue51" #<ESA Copernicus Scihub password>
username_USGS<-"ptaconet" #<USGS username>
password_USGS<-"HHKcue51" #<USGS password>
username_EarthData<-"ptaconet"  #<EarthData username>
password_EarthData<-"HHKcue51"  #<EarthData password>
  

## Indexes of the coordinates of the bounding box for OpenDap query to download GPM data.
# OpenDap is a standard protocol that enables to downlaad the data only for a given ROI, instead of downloading the whole GPM data which is provided by default on the whole globe. 
# The coordinates of the bounding box are not provided in lat and lon. They are provided as indexes. 
index_opendap_gpm_lon_min<-1741  # For CIV
index_opendap_gpm_lon_max<-1745  # For CIV
index_opendap_gpm_lat_min<-989  # For CIV
index_opendap_gpm_lat_max<-995  # For CIV

index_opendap_modisveget_lon_min<-2000  # For CIV
index_opendap_modisveget_lon_max<-2220  # For CIV
index_opendap_modisveget_lat_max<-510  # For CIV
index_opendap_modisveget_lat_min<-200  # For CIV

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

# initiate cluster for paralell download 
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = "PSOCK")
registerDoParallel(cl)

## Set ROI as sf object
roi_sf <- read_sf(path_to_roi_vector)$geometry
set_aoi(roi_sf)

## Set ROI as sp SpatialPolygon object
roi_sp<-rgdal::readOGR(path_to_roi_vector)
roi_sp <- spTransform(roi_sp, CRS(epsg))


############# make menage file a geo dataset with positions #######

## Set date 
date_epidemio_campain <- "2016-01-20"  # Date d'une capture au village DON au BF

menages_csv<-read.csv(path_to_file_coord_menages)
menages_csv<-menages_csv %>% filter (!(is.na(Part12Part2coordgpsLatitude))) %>% filter (!(is.na(Part12Part2coordgpsLongitude)))

menages_sf<-st_as_sf(menages_csv, coords = c("Part12Part2coordgpsLongitude", "Part12Part2coordgpsLatitude"), crs = 4326)
menages_sf<-st_transform(menages_sf, crs = 32630)

menages_sp_4326<-SpatialPointsDataFrame(coords=data.frame(menages_csv$Part12Part2coordgpsLongitude,menages_csv$Part12Part2coordgpsLatitude),data=menages_csv,proj4string=CRS("+init=epsg:4326"))
menages_sp<-spTransform(menages_sp,CRS("+init=epsg:32630"))

## Create 2 km buffer around the points
#buffer_2km <- st_buffer(menages_sf,2000) 
#buffer_2km<-st_union(buffer_2km,by_feature = TRUE)

## Get couple {date,coord}
dates_loc_enquetes<-unique(menages_sf[c("dateenq", "geometry")])


########################################################################################################################
############ Start Workflow ############
########################################################################################################################
cat("Starting workflow")

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

### Retrieve the metadata of the 1km/1-day resolution Terra and Aqua MODIS products 16 days before the campain date
time_range<-as.character(c(as.Date(date_epidemio_campain)-16,as.Date(date_epidemio_campain)))
records_terra <- getMODIS_query(time_range = time_range, name = "MODIS_MOD11A1_V6")
records_aqua <- getMODIS_query(time_range = time_range, name = "MODIS_MYD11A1_V6")
records<-rbind(records_terra,records_aqua)

# To preview the data :
# getMODIS_preview(records_terra[1,],show_aoi=TRUE)

### Download the data 

## Method 1 : linear dowload of all the files. Test: takes 1266 seconds for 18 records
# datasets <- getMODIS_data(records = records,dir_out=path_to_modislst_raw_folder) 

## Method 2 : parallel download. Test: much faster than linear download = takes 328 seconds for the same 18 records
# Code extracted from http://jxsw.de/getSpatialData/

# Download
files <- foreach(i = 1:nrow(records[]), 
                 .combine=c, 
                 .packages='getSpatialData') %dopar% {
                   getMODIS_data(records[i, ], dir_out = path_to_modislst_raw_folder, force=FALSE)
                 }



# Some data might not have been properly downloaded. Download linearly the files that have not been downloaded
# list available files on the processed folder
modis_lst_files<-list.files(path_to_modislst_raw_folder,recursive = TRUE,full.names = TRUE)
modis_lst_files_not_dl<-which(grepl("file",modis_lst_files))
if (length(modis_lst_files_not_dl)>0){
file.remove(modis_lst_files[modis_lst_files_not_dl])
datasets <- getMODIS_data(records = records[modis_lst_files_not_dl,],dir_out=path_to_modislst_raw_folder, force=FALSE) 
}


##############################################################
#### 2.2 - Prepare data ####
##############################################################

## For each date before the d date:
# - For each cell, extract the maximum day available temperature among Terra and Aqua. In case of NA on one of the two datasets, take the value of the other dataset.
# - For each cell, extract the minimum night available temperature among Terra and Aqua. In case of NA on one of the two datasets, take the value of the other dataset.
# - Extract the mean of the available LST cells on a buffer around each catch point
# - If no data are available for a given catch point (i.e. cells are all NA), set temperature of the previous night (recursively until a value is available)
# - If no data are available for the whole time series, ???
# - Create thermal amplitude value (i.e. difference between day and night)


# list available files on the output folder and retrieve path to the data to process
modis_lst_files<-list.files(path_to_modislst_raw_folder,recursive = TRUE,full.names = TRUE)

dates <-seq(as.Date(time_range[1]),as.Date(time_range[2]),1)
# convert dates to julian dates as they are provided this way in the MODIS products convention naming
dates_julian<-NULL
for (i in 1:length(dates)){
julian_day<-sprintf("%03d",julian(dates[i],origin = as.Date(x = paste0(format(dates[i],'%Y'),"-01-01"))-1))
dates_julian<-c(dates_julian,paste0(format(dates[i],'%Y'),julian_day))
}

terra_names<-paste0("MOD11A1.A",dates_julian)
aqua_names<-paste0("MYD11A1.A",dates_julian)


# Function to open the Aqua and Terra data for a given date and create a raster brick out of them
fun_open_project_crop_brick_modis_lst<-function(sds_index,sds_aqua,sds_terra,roi_sp,epsg){
  
  # Open Aqua (e.g. LST_day), reproject and crop following the ROI
  rast_aqua<-raster(sds_aqua[sds_index])
  rast_aqua <- projectRaster(rast_aqua, crs = CRS(epsg))
  rast_aqua<-crop(rast_aqua,roi_sp)
  
  # Same for Terra
  rast_terra<-raster(sds_terra[sds_index])
  rast_terra <- projectRaster(rast_terra, crs = CRS(epsg))
  rast_terra<-crop(rast_terra,roi_sp)
  
  # Create a raster brick out of the 2 rasters
  rast_brick<-brick(rast_aqua,rast_terra)
  
  return(rast_brick)
}


dates<-rev(dates)


list_lst_day<-list()
list_lst_nigth<-list()

for (i in 1:length(dates)){

  cat(paste0("pre-processing LST for date ",dates[i],"\n"))
  
### Open Aqua product
aqua_product_id<-modis_lst_files[which(grepl(aqua_names[i],modis_lst_files))]
# Get datasets available in the .hdf file (tells what subdatasets are within the hdf4 MODIS files and makes them into a list)
sds_aqua <- gdalUtils::get_subdatasets(aqua_product_id)

### Open Terra product
terra_product_id<-modis_lst_files[which(grepl(terra_names[i],modis_lst_files))]
# Get datasets available in the .hdf file
sds_terra <- gdalUtils::get_subdatasets(terra_product_id)


# Create raster of LST_day maximum value. sds_index for LST_day is 1
rast_lst_day<-fun_open_project_crop_brick_modis_lst(1,sds_aqua,sds_terra,roi_sp,epsg)
rast_max_lst_day<-max(rast_lst_day,na.rm=TRUE)
list_lst_day[[length(list_lst_day)+1]]<-rast_max_lst_day
names(list_lst_day[[length(list_lst_day)]])<-paste0("lst_max_d_",i)

# Create raster of LST_day miminum value. sds_index for LST_night is 5
rast_lst_night<-fun_open_project_crop_brick_modis_lst(5,sds_aqua,sds_terra,roi_sp,epsg)
rast_min_lst_night<-min(rast_lst_night,na.rm=TRUE)
list_lst_nigth[[length(list_lst_nigth)+1]]<-rast_min_lst_night
names(list_lst_nigth[[length(list_lst_nigth)]])<-paste0("lst_min_d_",i)

}


lst_day<-brick(list_lst_day)
lst_night<-brick(list_lst_nigth)

## Extract mean of LST in a buffer of 2 km for all the dates
# na.rm = TRUE means that NA values in the raster are not taken into account
# small=TRUE means that each time the buffer intersects a non NA cell it takes into account the cell value
ex <- raster::extract(lst_day, menages_sp, buffer=2000,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 
ex <- raster::extract(lst_night, ex, buffer=2000,fun=mean, na.rm=TRUE, sp=TRUE,small=TRUE) 


## 

#####  gapfill is a package to fill gaps on a raster with NA values. More about gapfill package: see https://gis.stackexchange.com/questions/214144/r-package-gapfill-how-to-convert-r-raster-stack-to-4-dimensional-array-and-then 
#r<-stack(brick(list_lst_day))

#tmp2 <- array(r, dim=c(dim(r[[1]])[2],dim(r[[1]])[1],2,16))
#input_array <- aperm(tmp2, c(2,1,3,4))


#gapfill::Image(input_array, col=rev(terrain.colors(100)),asRaster=TRUE)


#output <- gapfill::Gapfill(data=input_array, dopar = TRUE)
#output_array <- output$fill
#gapfill::Image(output_array)

## convert back
#output_stack <- stack(brick(array(output_array, c(dim(r[[1]])[1],dim(r[[1]])[2],16))))
#extent(output_stack)<-extent(r)
#crs(output_stack)<-crs(r)
#plot(output_stack)





# To convert to geotiff: https://stackoverflow.com/questions/36772341/reading-hdf-files-into-r-and-converting-them-to-geotiff-rasters

# Provides detailed data on hdf4 files but takes ages
# gdalUtils::gdalinfo(modis_lst_files[[1]])

# Tells what subdatasets are within the hdf4 MODIS files and makes them into a list

# Use gdal_translate to convert the first dataset to a .tif
# gdal_translate(sds[1], dst_dataset = "NPP2000.tif")

# Note: in Linux any sds can then be read directly using the raster function
#r <- raster(sds[1])

# For the MODIS 1-day and 8-day products (MOD11A1 , MYD11A1, MOD11A2 , MYD11A2) :
# sds[1] = LST_Day_1km
# sds[2] = QC_Day   
# sds[3] = Day_view_time
# sds[4] = Day_view_angl 
# sds[5] = LST_Night_1km  
# sds[6] = QC_Night      
# sds[7] = Night_view_time
# sds[8] = Night_view_angl
# sds[9] = Emis_31    
# sds[10] = Emis_32     
# sds[11] = Clear_day_cov 
# sds[12] = Clear_night_cov




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

httr::set_config(authenticate(user=username_EarthData, password=password_EarthData, type = "basic"))

url_modisveget_opendap<-"https://opendap.cr.usgs.gov/opendap/hyrax/"
modis_1km_tile<-"h17v08"

time_range<-as.character(c(as.Date(date_epidemio_campain)-20,as.Date(date_epidemio_campain)+8))
records_terra <- getMODIS_query(time_range = time_range, name = "MODIS_MOD13Q1_V6")
records_aqua <- getMODIS_query(time_range = time_range, name = "MODIS_MYD13Q1_V6")
records<-rbind(records_terra,records_aqua)

## Retrieve date and satellite of the week including the catch (i.e. closest start aquisition date from the catch date) and build the link to download the data

fun_build_modis_veget_opendap_url<-function(date_catch){

index_closest_date_to_catch<-which(abs(as.Date(records$acquisitionDate)-as.Date(date_catch)) == min(abs(as.Date(records$acquisitionDate) - as.Date(date_catch))))
closest_date_to_catch<-min(records$acquisitionDate[index_closest_date_to_catch])
sat_closest_date_to_catch<-records$product[which(records$acquisitionDate==closest_date_to_catch)]

opendap_index_closest_date_to_catch<-as.integer(difftime(closest_date_to_catch ,"2000-01-01" , units = c("days")))
opendap_index_closest_date_to_catch<-trunc(opendap_index_closest_date_to_catch/16)

if (grepl("MOD13Q1",sat_closest_date_to_catch)){
  opendap_collection_closest_date_to_catch<-"MOD13Q1.006"
} else {
  opendap_collection_closest_date_to_catch<-"MYD13Q1.006"
}

url_opendap<-paste0(url_modisveget_opendap,opendap_collection_closest_date_to_catch,"/",modis_1km_tile,".ncml.nc4?_250m_16_days_NDVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],_250m_16_days_EVI[",opendap_index_closest_date_to_catch,"][",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"][",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"],time[",opendap_index_closest_date_to_catch,"],YDim[",index_opendap_modisveget_lat_min,":",index_opendap_modisveget_lat_max,"],XDim[",index_opendap_modisveget_lon_min,":",index_opendap_modisveget_lon_max,"]")

return(list(url_opendap,opendap_collection_closest_date_to_catch,opendap_index_closest_date_to_catch))

}


catch_week<-fun_build_modis_veget_opendap_url(date_epidemio_campain)
one_week_before<-fun_build_modis_veget_opendap_url(as.character(as.Date(date_epidemio_campain)-8))
two_week_before<-fun_build_modis_veget_opendap_url(as.character(as.Date(date_epidemio_campain)-16))

path_catch_week<-file.path(path_to_modisveget_raw_folder,paste0(catch_week[[2]],"_",catch_week[[3]],".nc4"))
path_one_week_before<-file.path(path_to_modisveget_raw_folder,paste0(one_week_before[[2]],"_",one_week_before[[3]],".nc4"))
path_two_week_before<-file.path(path_to_modisveget_raw_folder,paste0(two_week_before[[2]],"_",two_week_before[[3]],".nc4"))

httr::GET(url = catch_week[[1]], write_disk(path_catch_week))
httr::GET(url = one_week_before[[1]], write_disk(path_one_week_before))
httr::GET(url = two_week_before[[1]], write_disk(path_two_week_before))


fun_preprocess_modis_vegetation<-function(path_to_raw_modis,var_name){
grid_nc<-raster(path_to_raw_modis,varname=var_name)
projection(grid_nc)<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
extent(grid_nc)[1:2]<-extent(grid_nc)[1:2]+res(grid_nc)[1]/2
extent(grid_nc)[3:4]<-extent(grid_nc)[3:4]-res(grid_nc)[1]/2
grid_nc <- projectRaster(grid_nc, crs = CRS("+init=epsg:32630"))
return(grid_nc)
}


rast_catch_week_ndvi<-fun_preprocess_modis_vegetation(path_catch_week,"_250m_16_days_NDVI")
rast_catch_week_evi<-fun_preprocess_modis_vegetation(path_catch_week,"_250m_16_days_EVI")

rast_catch_one_week_before_ndvi<-fun_preprocess_modis_vegetation(path_one_week_before,"_250m_16_days_NDVI")
rast_catch_one_week_before_evi<-fun_preprocess_modis_vegetation(path_one_week_before,"_250m_16_days_EVI")

rast_catch_two_week_before_ndvi<-fun_preprocess_modis_vegetation(path_two_week_before,"_250m_16_days_NDVI")
rast_catch_two_week_before_evi<-fun_preprocess_modis_vegetation(path_two_week_before,"_250m_16_days_EVI")


#grid_nc <- projectRaster(grid_nc, crs = CRS("+init=epsg:32630"))

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

url_gpm_data_server<-"https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.05"
httr::set_config(authenticate(user=username_EarthData, password=password_EarthData, type = "basic"))
  
# We retrieve the data 15 days before the date_epidemio_campain (source of the duration: thèse Nico)
GMP_data_to_retrieve<-as.Date(date_epidemio_campain)
for (i in 1:15){
  
# Get date to retrieve
GMP_data_to_retrieve<-as.Date(date_epidemio_campain)-i

# build the url of the dataset
year<-format(GMP_data_to_retrieve,'%Y')
month<-format(GMP_data_to_retrieve,'%m')
product_name<-paste0("3B-DAY.MS.MRG.3IMERG.",gsub("-","",GMP_data_to_retrieve),"-S000000-E235959.V05.nc4.ascii")

url_product<-paste0(url_gpm_data_server,"/",year,"/",month,"/",product_name,"?HQprecipitation[",index_opendap_gpm_lon_min,":",index_opendap_gpm_lon_max,"][",index_opendap_gpm_lat_min,":",index_opendap_gpm_lat_max,"],lon[",index_opendap_gpm_lon_min,":",index_opendap_gpm_lon_max,"],lat[",index_opendap_gpm_lat_min,":",index_opendap_gpm_lat_max,"]")

# from https://wiki.earthdata.nasa.gov/display/EL/How+to+access+data+with+R
if (!(file.exists(file.path(path_to_gpm_raw_folder,product_name)))){
httr::GET(url = url_product, write_disk(file.path(path_to_gpm_raw_folder,product_name)))
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
#grid_nc <- projectRaster(grid_nc, crs = CRS(epsg))

# Open GPM data as raster and pre-process
grid_nc<-raster("/home/ptaconet/Téléchargements/3B-DAY.MS.MRG.3IMERG.20170920-S000000-E235959.V05.nc4(1).nc4")
projection(grid_nc)<-CRS("+init=epsg:4326")
grid_nc <- t(grid_nc)
grid_nc <- flip(grid_nc,'y')
grid_nc <- flip(grid_nc,'x')
grid_nc <- projectRaster(grid_nc, crs = CRS(epsg))


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


