#### Start Workflow
########################################################################################################################
############ Set Input parameters for the workflow ############
########################################################################################################################

rm(list = ls())

## Global parameters
path_to_processing_folder<-"/home/ptaconet/Documents/react/data_CIV"  #<Path to the processing folder (i.e. where all the data produced by the workflow will be stored)>
path_to_earthdata_credentials<-"credentials_earthdata.txt" # path to the file containing the credential to the NASA servers (EarthData)
path_to_grassApplications_folder<-"/usr/lib/grass74" #<Can be retrieved with grass74 --config path . More info on the use of rgrass7 at https://grasswiki.osgeo.org/wiki/R_statistics/rgrass7

## Path to the input dataset (date, lat, lon of HLC)
path_to_csv_hlc_dates_loc<-"df_hlc.csv"


### Set collections of data to use

## 1) Collections to use for the time series data
dataSources_timeSeries<-c("MOD11A1.006","MYD11A1.006","MOD13Q1.006","MYD13Q1.006","MOD16A2.006","MYD16A2.006","GPM_3IMERGDF","TAMSAT")
lagTime_timeSeries<-40  # lag days for the time series

## 2) Collections to use for the single date data (Data for the day of the catch)
dataSources_Dday<-c("GPM_3IMERGHH","ERA5","MIRIADE","VIIRS DNB")

## 3) Collections of data for no time series data (Data to calculate as static data)
dataSources_static<-c("SRTMGL1_v003","CGLS-LC100")  


buffer_sizes<-c(500,1000,2000)
gpmDay_resample_output_res=250
gpmHhour_resample_output_res=250
tamsatDay_resample_output_res=250
gpmHhour_hour_start<-"18"
gpmHhour_hour_end<-"08"

setwd(path_to_processing_folder)








library(raster)
library(sp)
require(sf)
require(tidyverse)
require(httr)
require(furrr)
require(rgeos)
require(lubridate)
source("/home/ptaconet/r_react/getData/getData_modis.R")
source("/home/ptaconet/r_react/getData/prepareData_modis.R")
source("/home/ptaconet/r_react/getData/getDataFromOpenDAP_ancillaryFunctions.R")
source("/home/ptaconet/r_react/getData/ancillaryFunctions.R")
source("/home/ptaconet/r_react/getData/getData_srtm.R")
source("/home/ptaconet/r_react/getData/getData_gpm.R")
source("/home/ptaconet/r_react/getData/prepareData_gpm.R")
source("/home/ptaconet/r_react/getData/getData_tamsat.R")
source("/home/ptaconet/r_react/getData/prepareData_tamsat.R")


## Connection to the EarthData servers
earthdata_credentials<-readLines(path_to_earthdata_credentials)
username_EarthData<-strsplit(earthdata_credentials,"=")[[1]][2]
password_EarthData<-strsplit(earthdata_credentials,"=")[[2]][2]
httr::set_config(authenticate(user=username_EarthData, password=password_EarthData, type = "basic"))


hlc_dates_loc<-read.csv(path_to_csv_hlc_dates_loc)

## Get ROI as sf object. We extend a bit the size of the bbox (of the max of the buffer size + 0.05°)
bbox_4326<-SpatialPointsDataFrame(coords=data.frame(hlc_dates_loc$longitude,hlc_dates_loc$latitude),data=hlc_dates_loc,proj4string=CRS("+init=epsg:4326")) %>% bbox
bbox_4326[,1]=bbox_4326[,1]-0.05-convertMetersToDegrees(max(buffer_sizes),mean(bbox_4326[2,])) # mean_latitude = mean(bbox_4326[2,])
bbox_4326[,2]=bbox_4326[,2]+0.05+convertMetersToDegrees(max(buffer_sizes),mean(bbox_4326[2,]))
roi<-rgeos::bbox2SP(bbox_4326[2,2],bbox_4326[2,1],bbox_4326[1,1],bbox_4326[1,2],proj4string=CRS("+init=epsg:4326")) %>% st_as_sf


##
dates_loc<- hlc_dates_loc %>% 
  mutate(date_capture=as.Date(date_capture)) %>%
  group_by(date_capture) %>%
  arrange(date_capture) %>%
  nest(latitude,longitude,idpointdecapture,village) %>%
  set_names(c("date_date","coords")) %>%
  mutate(sp_points=map(coords,~SpatialPointsDataFrame(coords=data.frame(.$longitude,.$latitude),data=.,proj4string=CRS("+init=epsg:4326")))) %>%
  dplyr::select(-coords) %>%
  mutate(date_numeric=as.integer(date_date))  %>%
  mutate(lagTime_timeSeries=lagTime_timeSeries)


dates<-dates_loc$date_date

dates_timesSeries<-dates_loc$date_date %>%
  map(~seq(.,.-lagTime_timeSeries,-1)) %>%
  unlist %>%
  unique %>%
  sort %>%
  as.Date(origin="1970-01-01")
  


plan(multiprocess) # use furrr for parallel computing
options(future.globals.maxSize= 10000*1024^2) # 10 GB for the max size to be exported for the furrr future expression (https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r)

extractVar_singleBuff<-function(rasts,names_rasts_to_use,spPoints,buffer_size){
  
  res<-rasts[names_rasts_to_use] %>%  # filter only the rasters between our dates of interest
    #rev %>%
    #future_   ## a priori ça ne sert à rien de paralleliser à ce niveau là, ca ne fait pas gagner de temps 
    map_dfr(~raster::extract(.,spTransform(spPoints,proj4string(.)),buffer=buffer_size,fun=mean, na.rm=TRUE, small=TRUE)) %>% # for each raster, calculate the stats
    #set_names(origin_date-as.numeric(names(.))) %>% ## to name by the number of days separating the date of interest from the the day of the data
    set_names(seq(0,ncol(.)-1,1)) %>%   ## to name by the lag index 
    mutate(idpointdecapture=as.character(spPoints$idpointdecapture)) %>%
    mutate(buffer=buffer_size) %>%
    gather(time_lag,val,-c(idpointdecapture,buffer))
  
  return(res) # to put in wide format : res <- res %>% unite(var,var,time_lag)) %>% spread(key=var,value=val)
}


extractVar<-function(buffer_sizes,dates_loc,names_rasts_to_use,rasters,var_name){
  
  res<-buffer_sizes %>% # for each buffer, calculate stats
    set_names %>%
    future_map_dfr(~pmap_dfr(list(names_rasts_to_use,dates_loc$sp_points,.),
                             ~extractVar_singleBuff(rasters,..1,..2,..3)))  %>%
    mutate(var=var_name) %>%
    dplyr::select(idpointdecapture,buffer,var,time_lag,val)
  # to put in wide format : lst_min <- lst_min %>% unite(var,var,time_lag) %>% spread(key=var,value=val)
  
  return(res)
}


# Get the names of the rasters to use for each date
getRasters_timeSeries<-function(dates_loc,xxxData_md,dataCollection){
  rastsNames<-map(dates_loc$date_numeric,~pluck(xxxData_md,as.character(.))) %>%
    map(.,pluck(dataCollection)) %>%
    map(.,pluck("name")) %>%
    map(.,as.character)
  return(rastsNames)
}

# Get the path of local datasets
getPaths<-function(xxxData_md,dataCollection){
  path<-modify(xxxData_md,dataCollection) %>%
    map(data.frame) %>%
    reduce(bind_rows) %>%
    dplyr::select(name,destfile) %>%
    unique
  return(path)
}

######################## 1/ MODIS #######################

## MODIS data are downloaded for all collections at the same time. 
## Then they are processed separately for each variable




modisCrs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
modisTile<-getMODIStileNames(roi)
roi_bbox<-st_bbox(st_transform(roi,modisCrs))

### Get Time and Space vectors of OpenDAP servers for each collection + spatial indices corresponding to the ROI boundaries
cat("Retrieving information to download MODIS data on the OpenDap servers...\n")
modisOpenDAP_md<- #data.frame(source=dataSources_timeSeries) %>%
  dataSources_timeSeries[which(dataSources_timeSeries %in% c("MOD11A1.006","MYD11A1.006","MOD11A2.006","MYD11A2.006","MOD13Q1.006","MYD13Q1.006","MOD16A2.006","MYD16A2.006"))] %>% 
  data.frame(stringsAsFactors = F) %>%
  set_names("source") %>%
  mutate(OpendapURL=paste0("https://opendap.cr.usgs.gov/opendap/hyrax","/",source,"/",modisTile,".ncml")) %>%
                  mutate(Opendap_timeVector=map(OpendapURL,~getOpenDAPvector(.,"time"))) %>%
                  mutate(Opendap_xVector=map(OpendapURL,~getOpenDAPvector(.,"XDim"))) %>%
                  mutate(Opendap_yVector=map(OpendapURL,~getOpenDAPvector(.,"YDim"))) %>%
                  mutate(Opendap_minLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox$xmin))-1)) %>%
                  mutate(Opendap_maxLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox$xmax))-1)) %>%
                  mutate(Opendap_minLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox$ymax))-1)) %>%
                  mutate(Opendap_maxLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox$ymin))-1)) %>%
                  mutate(roiSpatialIndexBound=pmap(list(Opendap_minLat,Opendap_maxLat,Opendap_minLon,Opendap_maxLon),.f=~c(..1,..2,..3,..4))) %>%
                  mutate(destFolder=file.path(path_to_processing_folder,source)) %>%
                  mutate(dimensionsToRetrieve=case_when(source %in% c("MOD11A1.006","MYD11A1.006","MOD11A2.006","MYD11A2.006") ~ list(c("LST_Day_1km","LST_Night_1km")),
                                                        source %in% c("MOD13Q1.006","MYD13Q1.006") ~ list(c("_250m_16_days_NDVI","_250m_16_days_EVI")),
                                                        source %in% c("MOD16A2.006","MYD16A2.006") ~ list(c("ET_500m"))
                                                        ))


## Build list of datasets to DL for all MODIS collection and for all dates
modisData_md<-dates %>%
  set_names(as.numeric(.)) %>% # names will be numeric format of the dates (days since 1970-01-01)
  map(~pmap(list(.,pluck(modisOpenDAP_md,"source"),pluck(modisOpenDAP_md,"destFolder"),pluck(modisOpenDAP_md,"dimensionsToRetrieve"),pluck(modisOpenDAP_md,"Opendap_timeVector"),pluck(modisOpenDAP_md,"roiSpatialIndexBound")),
            ~getData_modis(time_range = c(..1-lagTime_timeSeries,..1),
                           OpenDAPCollection=..2,
                           destFolder=..3,
                           modisTile=modisTile,
                           dimensionsToRetrieve=unlist(..4),
                           timeVector=unlist(..5),
                           roiSpatialIndexBound=unlist(..6))) %>%
        set_names(pluck(modisOpenDAP_md,"source")))


## Download datasets
# First create directories in the wd
directories<-modisOpenDAP_md$source %>%
  as.character() %>%
  as.list()  %>%
  lapply(dir.create,recursive = TRUE)

# Then download
df_DataToDL<-modisData_md %>%
  modify_depth(2, ~map2(.x=pluck(.,"url"),.y=pluck(.,"destfile"),cbind)) %>% 
  flatten %>% 
  flatten %>% 
  do.call(rbind,.) %>% 
  data.frame(stringsAsFactors = F) %>% 
  distinct %>% 
  set_names("url","destfile")

Dl_res<-downloadData(df_DataToDL$url,df_DataToDL$destfile,username_EarthData,password_EarthData,TRUE)

  
############ Process MODIS time series
  
  ## MODIS LST min and LST max

# Get the names of the rasters to use for each date
rastsNames_modLst<-getRasters_timeSeries(dates_loc,modisData_md,"MOD11A1.006")

  # Build paths to data
  path_to_mod<-getPaths(modisData_md,"MOD11A1.006")
  path_to_myd<-getPaths(modisData_md,"MYD11A1.006")
  path_to_mod_myd<-merge(path_to_mod,path_to_myd,by="name",suffixes = c("_mod","_myd"))
  
  # Pre-process
  rasts_modLstMin<-path_to_mod_myd %>%
    mutate(rast_mod=map(destfile_mod,~prepareData_modis(.,"LST_Night_1km"))) %>%
    mutate(rast_myd=map(destfile_myd,~prepareData_modis(.,"LST_Night_1km"))) %>%
    mutate(rast=map2(rast_mod,rast_myd,~min(.x,.y,na.rm = T))) %>%
    mutate(rast=map(rast,~.-273.15)) %>%
    pluck("rast") %>%
    set_names(path_to_mod_myd$name)
  
  rasts_modLstMax<-path_to_mod_myd %>%
    mutate(rast_mod=map(destfile_mod,~prepareData_modis(.,"LST_Day_1km"))) %>%
    mutate(rast_myd=map(destfile_myd,~prepareData_modis(.,"LST_Day_1km"))) %>%
    mutate(rast=map2(rast_mod,rast_myd,~max(.x,.y,na.rm = T))) %>%
    mutate(rast=map(rast,~.-273.15)) %>%
    pluck("rast") %>%
    set_names(path_to_mod_myd$name)
   
   # Extract
   cat("Extracting LstMin...\n")
   lstMin<-extractVar(buffer_sizes,dates_loc,rastsNames_modLst,rasts_modLstMin,"lstMin") # to put in wide format : lstMin <- lstMin %>% unite(var,var,time_lag) %>% spread(key=var,value=val)
   cat("Extracting LstMax...\n")
   lstMin<-extractVar(buffer_sizes,dates_loc,rastsNames_modLst,rasts_modLstMax,"lstMax")
   
   rm(rasts_modLstMin,rasts_modLstMax)
   
   ## MODIS vegetation indices
   
   # Get the names of the rasters to use for each date
   rastsNames_modVeget<-map2(getRasters_timeSeries(dates_loc,modisData_md,"MOD13Q1.006"),getRasters_timeSeries(dates_loc,modisData_md,"MYD13Q1.006"),c) %>%
     map(as.numeric) %>%
     map(~sort(.,decreasing=TRUE)) %>%
     map(as.character)
   
   # Build paths to data
   path_to_mod<-getPaths(modisData_md,"MOD13Q1.006")
   path_to_myd<-getPaths(modisData_md,"MYD13Q1.006")
   path_to_mod_myd<-rbind(path_to_mod,path_to_myd) %>% arrange(name)
   
   # Pre-process
   rasts_modNdvi<-path_to_mod_myd %>%
     mutate(rast=map(destfile,~prepareData_modis(.,"_250m_16_days_NDVI"))) %>%
     pluck("rast") %>%
     set_names(path_to_mod_myd$name) 
 
   rasts_modEvi<-path_to_mod_myd %>%
     mutate(rast=map(destfile,~prepareData_modis(.,"_250m_16_days_EVI"))) %>%
     pluck("rast") %>%
     set_names(path_to_mod_myd$name) 
     
   # Extract
   cat("Extracting ndvi...\n")
   ndvi<-extractVar(buffer_sizes,dates_loc,rastsNames_modVeget,rasts_modNdvi,"ndvi")
   cat("Extracting evi...\n")
   evi<-extractVar(buffer_sizes,dates_loc,rastsNames_modVeget,rasts_modEvi,"evi")
   
   rm(rasts_modNdvi,rasts_modEvi)
   
   ## MODIS evapotranspiration
   
   # Get the names of the rasters to use for each date
   rastsNames_modEt<-getRasters_timeSeries(dates_loc,modisData_md,"MOD16A2.006")
   
   # Build paths to data
   path_to_mod<-getPaths_modis(modisData_md,"MOD16A2.006")
   path_to_myd<-getPaths_modis(modisData_md,"MYD16A2.006")
   path_to_mod_myd<-merge(path_to_mod,path_to_myd,by="name",suffixes = c("_mod","_myd"))
   
   # Pre-process
   rasts_modEt<-path_to_mod_myd %>%
     mutate(rast_mod=map(destfile_mod,~prepareData_modis(.,"ET_500m"))) %>%
     mutate(rast_myd=map(destfile_myd,~prepareData_modis(.,"ET_500m"))) %>%
     mutate(rast_mod=map(rast_mod,~clamp(.x,upper=32760,useValues=FALSE))) %>% ## Set pixel values >= 32760 (quality pixel values) to NA 
     mutate(rast_myd=map(rast_myd,~clamp(.x,upper=32760,useValues=FALSE))) %>% ## Set pixel values >= 32760 (quality pixel values) to NA 
     mutate(rast=map2(rast_mod,rast_myd,~mean(.x,.y,na.rm = T))) %>%
     pluck("rast") %>%
     set_names(path_to_mod_myd$name)
     
   # Extract
   cat("Extracting et...\n")
   et<-extractVar(buffer_sizes,dates_loc,rastsNames_modEt,rasts_modEt,"et")
   rm(rasts_modEt)
   


   ######################## 2/ GPM #######################
   

   roi_bbox<-st_bbox(roi)
   
   ### Get Time and Space vectors of OpenDAP servers for GPM Daily
   cat("Retrieving information to download GPM data on the OpenDap servers...\n")
   gpmOpenDAP_md<- #data.frame(source=dataSources_timeSeries) %>%
     c("GPM_3IMERGDF.06","GPM_3IMERGHH.06") %>% 
     data.frame(stringsAsFactors = F) %>%
     set_names("source") %>%
     mutate(OpendapURL="https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGHH.06/2016/001/3B-HHR.MS.MRG.3IMERG.20160101-S000000-E002959.0000.V06B.HDF5") %>%
     mutate(Opendap_xVector=map(OpendapURL,~getOpenDAPvector(.,"lon"))) %>%
     mutate(Opendap_yVector=map(OpendapURL,~getOpenDAPvector(.,"lat"))) %>%
     mutate(Opendap_minLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox$xmin))-4)) %>%
     mutate(Opendap_maxLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox$xmax))+4)) %>%
     mutate(Opendap_minLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox$ymin))-4)) %>% ## careful, this line is not the same as for Modis. ymax has become ymin.
     mutate(Opendap_maxLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox$ymax))+4)) %>%
     mutate(roiSpatialIndexBound=pmap(list(Opendap_minLat,Opendap_maxLat,Opendap_minLon,Opendap_maxLon),.f=~c(..1,..2,..3,..4))) %>%
     mutate(destFolder=file.path(path_to_processing_folder,source)) %>%
     mutate(dimensionsToRetrieve=list(c("precipitationCal")))
                                        
   
   ## Build list of datasets to DL for GPM Daily and for all dates
   gpmOpenDAP_md_daily<-gpmOpenDAP_md %>% filter(source=="GPM_3IMERGDF.06")
   gpmData_md_daily<-dates %>%
     set_names(as.numeric(.)) %>% # names will be numeric format of the dates (days since 1970-01-01)
     map(~pmap(list(.,pluck(gpmOpenDAP_md_daily,"source"),pluck(gpmOpenDAP_md_daily,"destFolder"),pluck(gpmOpenDAP_md_daily,"dimensionsToRetrieve"),pluck(gpmOpenDAP_md_daily,"roiSpatialIndexBound")),
               ~getData_gpm(time_range = c(..1-lagTime_timeSeries,..1),
                              OpenDAPCollection=..2,
                              destFolder=..3,
                              dimensionsToRetrieve=unlist(..4),
                              roiSpatialIndexBound=unlist(..5))) %>%
           set_names(pluck(gpmOpenDAP_md_daily,"source")))
   
   ## Build list of datasets to DL for GPM half-hourly and for all dates
   gpmOpenDAP_md_hhourly<-gpmOpenDAP_md %>% filter(source=="GPM_3IMERGHH.06")
   gpmData_md_hhourly<-dates %>%
     set_names(as.numeric(.)) %>% # names will be numeric format of the dates (days since 1970-01-01)
     map(~pmap(list(.,pluck(gpmOpenDAP_md_hhourly,"source"),pluck(gpmOpenDAP_md_hhourly,"destFolder"),pluck(gpmOpenDAP_md_hhourly,"dimensionsToRetrieve"),pluck(gpmOpenDAP_md_hhourly,"roiSpatialIndexBound")),
               ~getData_gpm(time_range = c(paste0(as.Date(..1,origin="1970-01-01")," ",gpmHhour_hour_start,":00:00"),paste0(as.Date(..1,origin="1970-01-01")+1," ",gpmHhour_hour_end,":00:00")),
                            OpenDAPCollection=..2,
                            destFolder=..3,
                            dimensionsToRetrieve=unlist(..4),
                            roiSpatialIndexBound=unlist(..5))) %>%
           set_names(pluck(gpmOpenDAP_md_hhourly,"source")))
   
   ## Download datasets
   # First create directories in the wd
   directories<-gpmOpenDAP_md$source %>%
     as.character() %>%
     as.list()  %>%
     lapply(dir.create,recursive = TRUE)
   
   # Then download
   df_DataToDL<-gpmData_md_daily %>%
     append(gpmData_md_hhourly) %>%
     modify_depth(2, ~map2(.x=pluck(.,"url"),.y=pluck(.,"destfile"),cbind)) %>% 
     flatten %>% 
     flatten %>% 
     do.call(rbind,.) %>% 
     data.frame(stringsAsFactors = F) %>% 
     distinct %>% 
     set_names("url","destfile")
   
   Dl_res<-downloadData(df_DataToDL$url,df_DataToDL$destfile,username_EarthData,password_EarthData,TRUE)
   

   ############ Process GPM time series
   
   ## Daily rainfall
   # Get the names of the rasters to use for each date
   rastsNames_gpmDay<-getRasters_timeSeries(dates_loc,gpmData_md_daily,"GPM_3IMERGDF.06")
   
   # Build paths to data
   path_to_gpmDay<-getPaths(gpmData_md_daily,"GPM_3IMERGDF.06")
   
   # Pre-process TODO check quality
   rasts_gpmDay<-path_to_gpmDay %>%
     mutate(rast=future_map(destfile,~prepareData_gpm(.,"precipitationCal",TRUE,gpmDay_resample_output_res))) %>%  # TRUE means that we resample the dataset
     pluck("rast") %>%
     set_names(path_to_gpmDay$name)
   
  # Extract
   cat("Extracting daily precipitations (gpm)...\n")
   rain_gpmDay<-extractVar(buffer_sizes,dates_loc,rastsNames_gpmDay,rasts_gpmDay,"rain_gpmDay")
   rm(rasts_gpmDay)
   
   ## Rainfall on the HCL date 
   # Get the names of the rasters to use for each date
   rastsNames_gpmHhour<-getRasters_timeSeries(dates_loc,gpmData_md_hhourly,"GPM_3IMERGHH.06")
   
   # Build paths to data
   path_to_gpmHhour<-getPaths(gpmData_md_hhourly,"GPM_3IMERGHH.06")
   
   # Pre-process TODO check quality
   rasts_gpmHhour<-path_to_gpmHhour %>%
     mutate(rast=future_map(destfile,~prepareData_gpm(.,"precipitationCal",TRUE,gpmHhour_resample_output_res))) %>%  # TRUE means that we resample the dataset
     pluck("rast") %>%
     set_names(path_to_gpmHhour$name)
   
   # Extract
   cat("Extracting half hourly precipitations (gpm)...\n")
   rain_gpmHhour<-extractVar(buffer_sizes=10,dates_loc,rastsNames_gpmHhour,rasts_gpmHhour,"rain_gpmHhour") # we put buffer_sizes=10 to extract the half-hourly rainfall at the HLC position only (ie not in a buffer)
  # Output is the rainfall for each half hour. To summarize for the whole night : 
   rain_hhour<-rain_hhour %>%
     group_by(idpointdecapture,var,buffer) %>%
     summarise(val=sum(val))
   
   
   ######################## 3/ TAMSAT #######################
   
   tamsat_md<-data.frame(output_time_step=c("daily"), #,"monthly","monthly"),
                         output_product=c("rainfall_estimate"), #,"rainfall_estimate","anomaly"),
                         output_output=c("individual"), #,"individual","individual"),
                         stringsAsFactors = F) %>%
     mutate(source=paste(output_time_step,output_product,output_output,sep="_"))  %>%
    mutate(destFolder=file.path(path_to_processing_folder,"TAMSAT",source)) 
   
   tamsatData_md<-dates %>%
     set_names(as.numeric(.)) %>% # names will be numeric format of the dates (days since 1970-01-01)
     map(~pmap(list(.,pluck(tamsat_md,"destFolder"),pluck(tamsat_md,"output_time_step"),pluck(tamsat_md,"output_product"),pluck(tamsat_md,"output_output")),
               ~getData_tamsat(time_range = c(..1-lagTime_timeSeries,..1),
                            destFolder=..2,
                            output_time_step=..3,
                            output_product=..4,
                            output_output=..5)) %>%
                        set_names(pluck(tamsat_md,"source")))
   
   ## Download datasets
   # First create directories in the wd
   directories<-tamsat_md$destFolder %>%
     as.character() %>%
     as.list()  %>%
     lapply(dir.create,recursive = TRUE)
   
   # Then download
   df_DataToDL<-tamsatData_md %>%
     modify_depth(2, ~map2(.x=pluck(.,"url"),.y=pluck(.,"destfile"),cbind)) %>%
     flatten %>%
     flatten %>%
     do.call(rbind,.) %>%
     data.frame(stringsAsFactors = F) %>%
     distinct %>%
     set_names("url","destfile")
   
   Dl_res<-downloadData(df_DataToDL$url,df_DataToDL$destfile,parallelDL=TRUE)
   
   
   ############ Process TAMSAT time series
   
   ## Daily rainfall
   # Get the names of the rasters to use for each date
   rastsNames_tamsatRain<-getRasters_timeSeries(dates_loc,tamsatData_md,"daily_rainfall_estimate_individual")
   
   # Build paths to data
   path_to_tamsatRain<-getPaths(tamsatData_md,"daily_rainfall_estimate_individual")
   
   # Pre-process
   rasts_tamsatRain<-path_to_tamsatRain %>%
     mutate(rast=future_map(destfile,~prepareData_tamsat(.,roi,TRUE,tamsatDay_resample_output_res))) %>%  # TRUE means that we resample the dataset
     pluck("rast") %>%
     set_names(path_to_tamsatRain$name)
   
   # Extract
   cat("Extracting daily precipitations (tamsat)...\n")
   rain_tamsatDay<-extractVar(buffer_sizes,dates_loc,rastsNames_tamsatRain,rasts_tamsatRain,"rain_tamsatDay")
   
   
   ## compare gpm and tamsat
   #rain_gpm_tamsat<-merge(data.table(rain_gpmDay %>% filter(buffer==2000)),data.table(rain_tamsatDay %>% filter(buffer==2000)),by=c("idpointdecapture","buffer","time_lag"))
   