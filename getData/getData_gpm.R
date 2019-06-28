getData_gpm<-function(time_range=c("2010-01-01","2010-01-30"), # mandatory. either a time range (e.g. c(date_start,date_end) ) or a single date e.g. ( date_start )
                      roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T), # either provide roi (sf point or polygon) or provide roiSpatialIndexBound. if roiSpatialIndexBound is not provided, it will be calculated from roi
                      username=NULL, # EarthData user name
                      password=NULL, # EarthData password
                      OpenDAPCollection="GPM_3IMERGHH.06", # mandatory
                      download=FALSE, # TRUE will download the file and return a list with : the URL, the path to the output file, a boolean wether the dataset was properly downloaded or not. FALSE will return a list with the URL only
                      destFolder=NULL,
                      modisTile=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      dimensionsToRetrieve=c("LST_Day_1km","LST_Night_1km"), # mandatory
                      timeVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      XVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      YVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      roiSpatialIndexBound=NULL,# optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      ...
){
  
  
  OpenDAPServerUrl="https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3"
  #TimeVectorName="time"
  SpatialXVectorName="lon"
  SpatialYVectorName="lat"
  #modisCollection_crs="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
  
  #if (OpenDAPCollection %in% c("MOD11A1.006","MYD11A1.006")){
  #  gridDimensionName<-"MODIS_Grid_Daily_1km_LST_eos_cf_projection"
  #} else if (OpenDAPCollection %in% c("MOD13Q1.006","MYD13Q1.006")){
  #  gridDimensionName<-"MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection"
  #} else if (OpenDAPCollection %in% c("MOD16A2.006","MYD16A2.006")){
  #  gridDimensionName<-"MOD_Grid_MOD16A2_eos_cf_projection"
  #}
  
  # Calculate modisTile if not provided
  #if(is.null(modisTile)){
  #  modisTile<-getMODIStileNames(roi)
  #}
  
  if(OpenDAPCollection=="GPM_3IMERGHH.06"){
    
  } else if(OpenDAPCollection=="GPM_3IMERGDF.06"){
  modisTile="2016/02/3B-DAY.MS.MRG.3IMERG.20160201-S000000-E235959.V06B.HDF5"
  
  #OpenDAPModisURL<-paste0(OpenDAPServerUrl,"/",OpenDAPCollection,"/",modisTile,".ncml")
  OpenDAPModisURL<-paste0(OpenDAPServerUrl,"/",OpenDAPCollection,"/",modisTile,".html")
  
  # Calculate TimeVector if not provided
  #if(is.null(timeVector)){
  #  timeVector<-getOpenDAPvector(OpenDAPModisURL,TimeVectorName)
  #}
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
  
  
  ## Build URLs
  
  
  
  
  
  
  
  
  
  
  
}