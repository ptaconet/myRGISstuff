roi_sf<-st_as_sf(roi_sp_4326)
hlc_dates_loc<-read.csv("/home/ptaconet/Documents/react/data_CIV/df_hlc.csv")
dates<-as.Date(unique(hlc_dates_loc$date_capture))



modisCollection_crs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

openDap_metadata<-read.csv("/home/ptaconet/r_react/getData/opendap_info.csv")

openDap_metadata$GetData<-TRUE

openDap_metadata$lag_days<-40

## 


### MODIS : Get Time and Space vectors, and spatial indices corresponding to the ROI boundaries for each collection
modisCollection_tile<-getMODIStileNames(roi_sf)
roi_bbox_modisCRS<-st_bbox(st_transform(roi_sf,modis_crs))

modisCollections<-openDap_metadata %>% 
                  filter(Collection=="MODIS" & GetData==TRUE) %>%
                  mutate(OpendapURL=paste0(OpenDAPServerUrl,"/",CodeName,"/",modisCollection_tile,".ncml")) %>%
                  mutate(Opendap_timeVector=map2(.x=OpendapURL,.y=TimeVectorName,.f=getOpenDAPvector)) %>%
                  mutate(Opendap_xVector=map2(.x=OpendapURL,.y=SpatialXVectorName,.f=getOpenDAPvector)) %>%
                  mutate(Opendap_yVector=map2(.x=OpendapURL,.y=SpatialYVectorName,.f=getOpenDAPvector)) %>%
                  mutate(Opendap_minLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$xmin))-1)) %>%
                  mutate(Opendap_maxLon=map(Opendap_xVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$xmax))-1)) %>%
                  mutate(Opendap_maxLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$ymin))-1)) %>%
                  mutate(Opendap_minLat=map(Opendap_yVector,.f=~which.min(abs(.x-roi_bbox_modisCRS$ymax))-1))



modisCollections <- modisCollections %>%
  mutate(dates_to_dl=map(.x=lag_days,~fun_get_unique_dates(.x,dates))) %>%
  mutate(urls_opendap=pmap(list(dates_to_dl,), (date=)))
  
  modisCollections %>% select (Collection,,dates_to_dl) %>% unnest()                           
                 
fun_get_unique_dates<-function(lag_days,dates){

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



