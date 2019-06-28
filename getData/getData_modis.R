

# Build openDAP modis URL for a given collection, time range and ROI, to retrieve the data in the ROI and for the closest available dates to the one provided
getData_modis<-function(time_range=c("2010-01-01","2010-01-30"), # mandatory. either a time range (e.g. c(date_start,date_end) ) or a single date e.g. ( date_start )
                        roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T), # either provide roi (sf point or polygon) or provide roiSpatialIndexBound. if roiSpatialIndexBound is not provided, it will be calculated from roi
                        username=NULL, # EarthData user name
                        password=NULL, # EarthData password
                        OpenDAPCollection="MOD11A1.006", # mandatory
                        #download=FALSE, # TRUE will download the file and return a list with : the URL, the path to the output file, a boolean wether the dataset was properly downloaded or not. FALSE will return a list with the URL only
                        destFolder=NULL,
                        modisTile=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                        dimensionsToRetrieve=c("LST_Day_1km","LST_Night_1km"), # mandatory
                        timeVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                        XVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                        YVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                        roiSpatialIndexBound=NULL,# optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                        ...
){
  
  OpenDAPServerUrl="https://opendap.cr.usgs.gov/opendap/hyrax"
  TimeVectorName="time"
  SpatialXVectorName="XDim"
  SpatialYVectorName="YDim"
  modisCollection_crs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
  
  if (OpenDAPCollection %in% c("MOD11A1.006","MYD11A1.006")){
    gridDimensionName<-"MODIS_Grid_Daily_1km_LST_eos_cf_projection"
  } else if (OpenDAPCollection %in% c("MOD13Q1.006","MYD13Q1.006")){
    gridDimensionName<-"MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection"
  } else if (OpenDAPCollection %in% c("MOD16A2.006","MYD16A2.006")){
    gridDimensionName<-"MOD_Grid_MOD16A2_eos_cf_projection"
  }
  
  # Calculate modisTile if not provided
  if(is.null(modisTile)){
    modisTile<-getMODIStileNames(roi)
  }
  
  OpenDAPModisURL<-paste0(OpenDAPServerUrl,"/",OpenDAPCollection,"/",modisTile,".ncml")
  
  # Calculate TimeVector if not provided
  if(is.null(timeVector)){
    timeVector<-getOpenDAPvector(OpenDAPModisURL,TimeVectorName)
  }
  # Calculate XVector if not provided
  if(is.null(XVector) & is.null(roiSpatialIndexBound)){
    XVector<-getOpenDAPvector(OpenDAPModisURL,SpatialXVectorName)
  }
  # Calculate YVector if not provided
  if(is.null(YVector) & is.null(roiSpatialIndexBound)){
    YVector<-getOpenDAPvector(OpenDAPModisURL,SpatialYVectorName)
  }
  # Calculate roiSpatialIndexBound if not provided
  if(is.null(roiSpatialIndexBound)){
    roi_bbox_modisCRS<-sf::st_bbox(st_transform(roi,modisCollection_crs))
    Opendap_minLon<-which.min(abs(XVector-roi_bbox_modisCRS$xmin))-1
    Opendap_maxLon<-which.min(abs(XVector-roi_bbox_modisCRS$xmax))-1
    Opendap_maxLat<-which.min(abs(YVector-roi_bbox_modisCRS$ymin))-1
    Opendap_minLat<-which.min(abs(YVector-roi_bbox_modisCRS$ymax))-1
    roiSpatialIndexBound<-c(Opendap_minLat,Opendap_maxLat,Opendap_minLon,Opendap_maxLon)
  }
  

  # Get openDAP time indices for the time frame of interest
  time_range<-as.Date(time_range,origin="1970-01-01")
  if (length(time_range)==1){
    time_range=c(time_range,time_range)
  }
  
  revisit_time<-timeVector[2]-timeVector[1]
  
  timeIndices_of_interest<-seq(time_range[2],time_range[1],-revisit_time) %>% 
    map(~getOpenDAPtimeIndex_modis(.,timeVector)) %>% 
    do.call(rbind.data.frame,.) %>%
    set_names("ideal_date","date_closest_to_ideal_date","days_sep_from_ideal_date","index_opendap_closest_to_date") %>%
    mutate(ideal_date=as.Date(ideal_date,origin="1970-01-01")) %>%
    mutate(date_closest_to_ideal_date=as.Date(date_closest_to_ideal_date,origin="1970-01-01"))
  
  # Build URL to download data in NetCDF format
  
  table_urls<-timeIndices_of_interest %>%
    mutate(dimensions_url=map(.x=index_opendap_closest_to_date,.f=~getOpenDapURL_dimensions(dimensionsToRetrieve,.x,roiSpatialIndexBound,TimeVectorName,SpatialXVectorName,SpatialYVectorName))) %>%
    mutate(url=paste0(OpenDAPModisURL,".nc4?",gridDimensionName,",",dimensions_url))

  
  urls<-table_urls$url

  destfiles<-paste0(file.path(destFolder,paste0(OpenDAPCollection,"_",gsub("-","",table_urls$date_closest_to_ideal_date))),".nc4")
  
  names<-as.numeric(timeIndices_of_interest$date_closest_to_ideal_date)
  
  #if (download){
  #  res<-downloadData(urls,destfiles,username,password,parallel)
  #} else {
  #  res<-1:length(urls);NA
  #}
  
  
  return(list(name=names,url=urls,destfile=destfiles))
  
}

