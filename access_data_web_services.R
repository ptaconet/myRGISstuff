


### Call useful libraries
library(getSpatialData)
library(raster)
library(sf)
library(sp)
library(gdalUtils)
library(rgdal)
require(httr)
## Load packages for working on multi-core
library(parallel)
library(doParallel)
library(foreach)

### Set Input parameters
## Working Directory
path_to_processing_folder<-"/home/ptaconet/Documents/react/data_CIV"  #<Path to the processing folder (i.e. where all the data produced by the workflow will be stored)>
path_to_roi_vector="ROI.kml" #<Path to the Region of interest in KML format>
#proj_srs="+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs" #<proj srs for the ROI>
epsg="+init=epsg:32630"

## Set date 
date_epidemio_campain <- "2017-01-20"  # Date d'une capture au village DON au BF
## Set credential to the various servers
username_cophub<-"ptaconet"
password_cophub<-"HHKcue51"
username_USGS<-"ptaconet"
password_USGS<-"HHKcue51"
username_EarthData<-"ptaconet"
password_EarthData<-"HHKcue51"  
  
## set login credentials for the Copernicus Scihub (Sentinel 1, 2, 3) and the EarthExplorer (MODIS)
login_CopHub(username = username_cophub, password = password_cophub)
login_USGS(username = username_USGS, password = password_USGS)

### Create directories for output files
path_to_sentinel1_products<-file.path(path_to_processing_folder,"Sentinel_1")
path_to_sentinel1_raw_folder<-file.path(path_to_sentinel1_products,"raw_data")
path_to_sentinel1_processed_folder<-file.path(path_to_sentinel1_products,"processed_data")
path_to_modislst_folder<-file.path(path_to_processing_folder,"MODIS_LST")
path_to_modislst_raw_folder<-file.path(path_to_modislst_folder,"raw_data")
path_to_modislst_processed_folder<-file.path(path_to_modislst_folder,"processed_data")
path_to_gpm_folder<-file.path(path_to_processing_folder,"GPM")
path_to_gpm_raw_folder<-file.path(path_to_gpm_folder,"raw_data")
path_to_gpm_processed_folder<-file.path(path_to_gpm_folder,"processed_data")

directories<-list(path_to_sentinel1_products,path_to_sentinel1_raw_folder,path_to_sentinel1_processed_folder,path_to_modislst_folder,path_to_modislst_raw_folder,path_to_modislst_processed_folder,path_to_gpm_folder,path_to_gpm_raw_folder,path_to_gpm_processed_folder)
lapply(directories, dir.create)


## Set working directory
setwd(path_to_processing_folder)

## set ROI 
roi_sf <- read_sf(path_to_roi_vector)$geometry
set_aoi(roi_sf)

## Get ROI as sp SpatialPolygon 
roi_sp<-rgdal::readOGR(path_to_roi_vector)
## Get ROI in the right EPSG
roi_sp <- spTransform(roi_sp, CRS(epsg))

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
# product_names <- GetSpatialData::getMODIS_names()

# MODIS LST products are listed here https://modis.gsfc.nasa.gov/data/dataprod/mod11.php) : 
# We are interested in daily L3 global 1km so we take : MOD11A1 , MYD11A1 
# Note that the 8-Day product (MOD11A2, MYD11A2) is simply an average of the 1 day products over the time period

# Downlaade 1km/1-day resolution Terra and Acqua MODIS records 8 days before the campain date (i.e. 8 )
time_range<-as.character(c(as.Date(date_epidemio_campain)-3,as.Date(date_epidemio_campain)))
records_terra <- getMODIS_query(time_range = time_range, name = "MODIS_MOD11A1_V6")
records_aqua <- getMODIS_query(time_range = time_range, name = "MODIS_MYD11A1_V6")

#getMODIS_preview(records_terra[1,],show_aoi=TRUE)

records<-rbind(records_terra,records_aqua)

### Download the data 
# Method 1 : linear dowload of all the files. Test: takes 1266 seconds for 18 records
# datasets <- getMODIS_data(records = records,dir_out=path_to_modislst_raw_folder) 

# Method 2 : parallel download. Test: much faster than linear download = takes 328 seconds for the same 18 records
# Code extracted from http://jxsw.de/getSpatialData/

# initiate cluster for paralell download #
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = "PSOCK")
registerDoParallel(cl)

# Download
files <- foreach(i = 1:nrow(records[]), 
                 .combine=c, 
                 .packages='getSpatialData') %dopar% {
                   getMODIS_data(records[i, ], dir_out = path_to_modislst_raw_folder)
                 }

#### stop cluster ####
stopCluster(cl)



# Convert to geotiff (from https://stackoverflow.com/questions/36772341/reading-hdf-files-into-r-and-converting-them-to-geotiff-rasters )

# list available files on the processed folder
modis_lst_files<-list.files(path_to_modislst_raw_folder,recursive = TRUE,full.names = TRUE)

# Provides detailed data on hdf4 files but takes ages
# gdalUtils::gdalinfo(modis_lst_files[[1]])

# Function to pre-process and prepare MODIS LST products 
# Algo : 
# - Open Aqua - day
# - Open Terra - day
# - Get location of NA values in Aqua
# - Replace by the values in Terra
# - Same for night
# - Create the day-night difference raster

dates <-seq(as.Date(time_range[1]),as.Date(time_range[2]),1)
# convert dates to julian dates as they are provided this way in the MODIS products convention naming
dates_julian<-NULL
for (i in 1:length(dates)){
julian_day<-sprintf("%03d",julian(dates[i],origin = as.Date(x = paste0(format(dates[i],'%Y'),"-01-01"))-1))
dates_julian<-c(dates_julian,paste0(format(dates[i],'%Y'),julian_day))
}

terra_names<-paste0("MOD11A1.A",dates_julian)
aqua_names<-paste0("MYD11A1.A",dates_julian)

# for (i in 1:length(dates)) # for each date
# fun_preprocess_modis_lst<-function(modis_lst_files,dates)
# i represents 1 date 
i=3



### Open Aqua product
aqua_product_id<-modis_lst_files[which(grepl(aqua_names[i],modis_lst_files))]
# Get datasets available in the .hdf file (tells what subdatasets are within the hdf4 MODIS files and makes them into a list)
sds_aqua <- gdalUtils::get_subdatasets(aqua_product_id)

### Open Terra product
terra_product_id<-modis_lst_files[which(grepl(terra_names[i],modis_lst_files))]
# Get datasets available in the .hdf file
sds_terra <- gdalUtils::get_subdatasets(terra_product_id)


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

# Create raster of LST_day maximum value. sds_index for LST_day is 1
rast_lst_day<-fun_open_project_crop_brick_modis_lst(1,sds_aqua,sds_terra,roi_sp,epsg)
rast_max_lst_day<-max(rast_lst_day,na.rm=TRUE)

# Create raster of LST_day miminum value. sds_index for LST_night is 5
rast_lst_night<-fun_open_project_crop_brick_modis_lst(5,sds_aqua,sds_terra,roi_sp,epsg)
rast_min_lst_night<-min(rast_lst_night,na.rm=TRUE)





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
########### GPM ###############
#####################################

## Websites
# TRMM doc: https://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B43/doc/README.TRMM_V7.pdf
# TRMM products :https://pmm.nasa.gov/data-access/downloads/trmm
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
# Through OpenDAP protocol here (possibility to subset the dataset on our bounding box before downloading it): https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.05
# Through standard http protocol here (impossible to subset the dataset): https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDF.05

url_gpm_data_server<"-https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDF.05"
set_config(authenticate(user=username_EarthData, password=password_EarthData, type = "basic"))
  
# We retrieve the data 15 days before the date_epidemio_campain (source: thèse Nico)
GMP_data_to_retrieve<-as.Date(date_epidemio_campain)
for (i in 1:15){
  
# Get date to retrieve
GMP_data_to_retrieve<-as.Date(date_epidemio_campain)-i

# build the link to the dataset
year<-format(GMP_data_to_retrieve,'%Y')
month<-format(GMP_data_to_retrieve,'%m')
product_name<-paste0("3B-DAY.MS.MRG.3IMERG.",gsub("-","",GMP_data_to_retrieve),"-S000000-E235959.V05.nc4")
url_product<-paste0(url_gpm_data_server,"/",year,"/",month,"/",product_name)

# from https://wiki.earthdata.nasa.gov/display/EL/How+to+access+data+with+R
httr::GET(url = url_product, write_disk(file.path(path_to_gpm_raw_folder,product_name), overwrite = TRUE))

}


# Another working solution to DL the data :
#system(paste0("wget --http-user=",username_EarthData," --http-password=",password_EarthData," --output-document=",file.path(path_to_gpm_raw_folder,product_name)," ",url_product))

## Only works with wget by a command system... we have tried many other options (via download.file, etc) but none work. Seems that R has problems with proxy



#####################################
########### Wind (ERA) ###############
#####################################
# Check :
# https://dominicroye.github.io/en/2018/access-to-climate-reanalysis-data-from-r/
# Install packages
if(!require("reticulate")) install.packages("reticulate")
require(reticulate)

### Initialisation
# Set virtual env and install the python ECMWF API
py_install("ecmwf-api-client") #(from https://community.rstudio.com/t/problem-installing-python-libraries-with-reticulate-py-install-error-pip-not-found/26561/2)
system("virtualenv -p /usr/bin/python2 /home/ptaconet/.virtualenvs/py2-virtualenv")
# Also works with this
#virtualenv_create("py3-virtualenv", python = "/usr/bin/python3")
use_virtualenv("py2-virtualenv")
# install the python ECMWF API
py_install("ecmwf-api-client", envname = "py2-virtualenv")
