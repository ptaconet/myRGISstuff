roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T)
hlc_dates_loc<-read.csv("/home/ptaconet/Documents/react/data_CIV/df_hlc.csv")


require(dplyr)
require(sf)
require(stringr)
require(tidyverse)
lagTime_timeSeries<-20  # lag days for the time series
buffer_sizes<-c(500,1000,2000)
#buffer_sizes<-seq(200,2000,200)

dates_loc<- hlc_dates_loc %>% 
  mutate(date_capture=as.Date(date_capture)) %>%
  group_by(date_capture) %>%
  arrange(date_capture) %>%
  nest(latitude,longitude,idpointdecapture,village) %>%
  set_names(c("date_date","coords")) %>%
  mutate(sp_points=map(coords,~SpatialPointsDataFrame(coords=data.frame(.$longitude,.$latitude),data=.,proj4string=CRS("+init=epsg:4326")))) %>%
  select(-coords) %>%
  mutate(date_numeric=as.integer(date_date))  %>%
  mutate(lagTime_timeSeries=lagTime_timeSeries) #%>%
  #mutate(dates_timeSeries=map(date_numeric,~seq(.,.-lagTime_timeSeries,-1))) %>%
  #mutate(dates_timeSeries=map(date,~seq(.,.-lagTime_timeSeries,-1))) %>%
  #filter(date_date %in% c(as.Date("2016-09-21"),as.Date("2016-09-22")))


dates<-dates_loc$date_date

dates_timesSeries<-dates_loc$date_date %>%
  map(~seq(.,.-lagTime_timeSeries,-1)) %>%
  unlist %>%
  unique %>%
  sort %>%
  as.Date(origin="1970-01-01")
  


wd="/home/ptaconet/Documents/react/data_CIV"

## 1) time series data
dataSources_timeSeries<-c("MOD11A1.006","MYD11A1.006","MOD13Q1.006","MYD13Q1.006","MOD16A2.006","MYD16A2.006")#,"GPM_3IMERGDF","TAMSAT") # Data to calculate as time series
lagTime_timeSeries<-40  # lag days for the time series

## 2) single date data
dataSources_Dday<-c("MOD11A1.006","MYD11A1.006","MOD13Q1.006","MYD13Q1.006","MOD16A2.006","MYD16A2.006","TAMSAT","GPM_3IMERGHH","ERA5","MIRIADE","VIIRS DNB") # Data for the day of the catch

## 3) no time series data
dataSources_static<-c("SRTMGL1_v003","CGLS-LC100")  # Data to calculate as static data






if (nrow(dataSources_modis)>0){
  
modisCrs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
modisTile<-getMODIStileNames(roi)
roi_bbox_modisCRS<-st_bbox(st_transform(roi,modisCrs))

######################## 1/ MODIS #######################

## MODIS data are downloaded for all collections at the same time. 
## Then they are processed separately for each variable


### Get Time and Space vectors of OpenDAP servers for each collection + spatial indices corresponding to the ROI boundaries
modisOpenDAP_md<- #data.frame(source=dataSources_timeSeries) %>%
  dataSources_timeSeries[which(dataSources_timeSeries %in% c("MOD11A1.006","MYD11A1.006","MOD13Q1.006","MYD13Q1.006","MOD16A2.006","MYD16A2.006"))] %>% 
  data.frame(stringsAsFactors = F) %>%
  set_names("source") %>%
  mutate(OpendapURL=paste0("https://opendap.cr.usgs.gov/opendap/hyrax","/",source,"/",modisTile,".ncml")) %>%
                  mutate(Opendap_timeVector=map(OpendapURL,~getOpenDAPvector(.,"time"))) %>%
                  mutate(Opendap_xVector=map(OpendapURL,~getOpenDAPvector(.,"XDim"))) %>%
                  mutate(Opendap_yVector=map(OpendapURL,~getOpenDAPvector(.,"YDim"))) %>%
                  mutate(Opendap_minLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$xmin))-1)) %>%
                  mutate(Opendap_maxLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$xmax))-1)) %>%
                  mutate(Opendap_minLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$ymax))-1)) %>%
                  mutate(Opendap_maxLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$ymin))-1)) %>%
                  mutate(roiSpatialIndexBound=pmap(list(Opendap_minLat,Opendap_maxLat,Opendap_minLon,Opendap_maxLon),.f=~c(..1,..2,..3,..4))) %>%
                  mutate(destFolder=file.path(wd,source)) %>%
                  mutate(dimensionsToRetrieve=case_when(source %in% c("MOD11A1.006","MYD11A1.006") ~ list(c("LST_Day_1km","LST_Night_1km")),
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
directories<-timeSeries_dataSources$source %>%
  as.character() %>%
  as.list()  %>%
  lapply(dir.create)

# Then download

df_DataToDL<-modisData_md %>%
  modify_depth(2, ~map2(.x=pluck(.,"url"),.y=pluck(.,"destfile"),cbind)) %>% 
  flatten %>% 
  flatten %>% 
  do.call(rbind,.) %>% 
  data.frame(stringsAsFactors = F) %>% 
  distinct %>% 
  set_names("url","destfile")


Dl_res<-downloadData(df_DataToDL$url,df_DataToDL$destfile,"ptaconet","HHKcue51",TRUE)

  
############ Process time series
  
  ## LST min and LST max

  mod<-modify(modisData_md,c("MOD11A1.006"))
  myd<-modify(modisData_md,c("MYD11A1.006"))
    
  path_to_mod<-mod %>%
    map(data.frame) %>%
    reduce(bind_rows) %>%
    select(name,destfile) %>%
    unique

  path_to_myd<-myd %>%
    map(data.frame) %>%
    reduce(bind_rows) %>%
    select(name,destfile) %>%
    unique
  
  path_to_mod_myd<-merge(path_to_mod,path_to_myd,by="name",suffixes = c("_mod","_myd"))
  
  #require(furrr)
  #plan(multiprocess) # use furrr for parallel computing
  rasts_lst_max<-path_to_mod_myd %>%
    mutate(rast_mod=map(destfile_mod,~prepareData_modis_open(.,"LST_Day_1km"))) %>%
    mutate(rast_myd=map(destfile_myd,~prepareData_modis_open(.,"LST_Day_1km"))) %>%
    mutate(lst=map2(rast_mod,rast_myd,~max(.x,.y,na.rm = T))) %>%
    mutate(lst=map(lst,~.-273.15)) %>%
    pluck("lst") %>%
    set_names(path_to_mod_myd$name)

  rasts_lst_min<-path_to_mod_myd %>%
    mutate(rast_mod=map(destfile_mod,~prepareData_modis_open(.,"LST_Night_1km"))) %>%
    mutate(rast_myd=map(destfile_myd,~prepareData_modis_open(.,"LST_Night_1km"))) %>%
    mutate(lst=map2(rast_mod,rast_myd,~min(.x,.y,na.rm = T))) %>%
    mutate(lst=map(lst,~.-273.15)) %>%
    pluck("lst") %>%
    set_names(path_to_mod_myd$name)
    
   
   fun_extract<-function(list_rast,lagTime_timeSeries,spPoints,origin_date,covariateAbr,buffer_size){
     #list_rast=rasts_lst_min
     #spPoints=dates_loc$sp_points[1][[1]]
     #lagTime_timeSeries=40
     res<-list_rast[between(as.numeric(names(list_rast)),origin_date-lagTime_timeSeries,origin_date)] %>%  # filter only the rasters between our dates of interest
       rev %>%
       future_map_dfr(~raster::extract(.,spTransform(spPoints,proj4string(.)),buffer=buffer_size,fun=mean, na.rm=TRUE, small=FALSE)) %>% # for each raster, calculate the stats
       set_names(paste0(covariateAbr,"_",origin_date-as.numeric(names(.)))) %>%
       mutate(idpointdecapture=as.character(spPoints$idpointdecapture)) %>%
       mutate(date=origin_date) %>%
       mutate(buffer=buffer_size) #%>%
       #select(idpointdecapture,date,buffer,everything())
     
     return(res)
   }
   
   ## Add lstMin
   
   
   require(furrr)
   plan(multiprocess) # use furrr for parallel computing
   
   tic()
   lst_min<-buffer_sizes %>% # for each buffer, calculate stats
     set_names %>%
     future_map_dfr(~pmap_dfr(list(dates_loc$lagTime_timeSeries,dates_loc$sp_points,dates_loc$date_numeric,.),
                 ~fun_extract(rasts_lst_min,..1,..2,..3,"lstMin",..4))) %>%
     select(idpointdecapture,date,buffer,everything())
   toc()
   
   tic()
   lst_max<-buffer_sizes %>% # for each buffer, calculate stats
     set_names %>%
     future_map_dfr(~pmap_dfr(list(dates_loc$lagTime_timeSeries,dates_loc$sp_points,dates_loc$date_numeric,.),
                       ~fun_extract(rasts_lst_max,..1,..2,..3,"lstMax",..4))) %>%
     select(idpointdecapture,date,buffer,everything())
   toc()
   
   
   
   
   
}













fun_get_unique_dates<-function(dates,lag_days){

  dates_dl<-NULL
  for (i in 1:length(dates)){
    dates_dl<-c(dates_dl,seq(dates[i],dates[i]-lag_days,-1))
  }
  dates_dl<-unique(dates_dl)
  dates_dl<-as.character(as.Date("1970-01-01")+dates_dl)
    
  return(dates_dl)
  
}



















## GPM collections : Additional useful information for GPM products
roi_bbox_4326<-st_bbox(roi_sf)

# GPM collection : Get Space vector
# for this we pick up a random dataset (https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/gpm_daily_collection/2016/02/3B-DAY.MS.MRG.3IMERG.20160201-S000000-E235959.V06.nc4)
Opendap_xVector_gpmCollection<-getOpenDAPvector("https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.06/2016/02/3B-DAY.MS.MRG.3IMERG.20160201-S000000-E235959.V06.nc4",gpmCollection_x_vector)
Opendap_yVector_gpmCollection<-getOpenDAPvector("https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGDF.06/2016/02/3B-DAY.MS.MRG.3IMERG.20160201-S000000-E235959.V06.nc4",gpmCollection_y_vector)

# GPM collection : Get spatial indices of  corresponding to our boundaries
Opendap_minLon_gpmCollection<-which.min(abs(Opendap_xVector_gpmCollection-roi_bbox_4326$xmin)-4)
Opendap_maxLon_gpmCollection<-which.min(abs(Opendap_xVector_gpmCollection-roi_bbox_4326$xmax)+4)
Opendap_minLat_gpmCollection<-which.min(abs(Opendap_yVector_gpmCollection-roi_bbox_4326$ymin)-4)
Opendap_maxLat_gpmCollection<-which.min(abs(Opendap_yVector_gpmCollection-roi_bbox_4326$ymax)+4)



