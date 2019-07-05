getData_gpm<-function(time_range=as.Date(c("2010-01-01","2010-01-30")), # mandatory. either a time range (e.g. c(date_start,date_end) ) or a single date e.g. ( date_start ) / or a as.POSIXlt single date or time range (e.g. "2010-01-01 18:00:00")
                      roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T), # either provide roi (sf point or polygon) or provide roiSpatialIndexBound. if roiSpatialIndexBound is not provided, it will be calculated from roi
                      username=NULL, # EarthData user name
                      password=NULL, # EarthData password
                      OpenDAPCollection="GPM_3IMERGHH.06", # mandatory
                      #download=FALSE, # TRUE will download the file and return a list with : the URL, the path to the output file, a boolean wether the dataset was properly downloaded or not. FALSE will return a list with the URL only
                      destFolder=NULL,
                      dimensionsToRetrieve=c("precipitationCal"), # mandatory
                      XVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      YVector=NULL, # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      roiSpatialIndexBound=NULL # optional. providing it will fasten the processing time. if not provided it will be calculated automatically
                      ){
  
  #require(lubridate)
  OpenDAPServerUrl="https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3"
  SpatialXVectorName="lon"
  SpatialYVectorName="lat"
  
  # Retrieve info to build url
  if(OpenDAPCollection=="GPM_3IMERGHH.06"){

    #times_gpm_hhourly<-seq(from=as.POSIXlt(paste0(this_date_hlc," ",hh_rainfall_hour_begin,":00:00")),to=as.POSIXlt(as.POSIXlt(paste0(this_date_hlc+1," ",hh_rainfall_hour_end,":00:00"))),by="30 min")
    time_range=as.POSIXlt(time_range)
    
    datesToRetrieve<-seq(from=time_range[2],to=time_range[1],by="-30 min") %>%
      data.frame(stringsAsFactors = F) %>%
      set_names("date") %>%
      mutate(date_character=as.character(as.Date(date))) %>%
      mutate(year=format(date,'%Y')) %>%
      mutate(month=format(date,'%m')) %>%
      mutate(day=sprintf("%03d",lubridate::yday(date))) %>%
      mutate(hour_start=paste0(sprintf("%02d",hour(date)),sprintf("%02d",minute(date)),sprintf("%02d",second(date)))) %>%
      mutate(hour_end=date+minutes(29)+seconds(59)) %>%
      mutate(hour_end=paste0(sprintf("%02d",hour(hour_end)),sprintf("%02d",minute(hour_end)),sprintf("%02d",second(hour_end)))) %>%
      mutate(number_minutes_from_start_day=sprintf("%04d",difftime(date,as.POSIXlt(paste0(as.Date(date)," 00:00:00")),units="mins")))
      
    urls<-datesToRetrieve %>% 
      mutate(product_name=paste0("3B-HHR.MS.MRG.3IMERG.",gsub("-","",date_character),"-S",hour_start,"-E",hour_end,".",number_minutes_from_start_day,".V06B.HDF5")) %>%
      mutate(url_product=paste(OpenDAPServerUrl,OpenDAPCollection,year,day,product_name,sep="/"))
    
  } else if(OpenDAPCollection=="GPM_3IMERGDF.06"){

    time_range=as.Date(time_range,origin="1970-01-01")
    
    datesToRetrieve<-seq(time_range[2],time_range[1],-1) %>%
      data.frame(stringsAsFactors = F) %>%
      set_names("date") %>%
      mutate(date_character=as.character(as.Date(date))) %>%
      mutate(year=format(date,'%Y')) %>%
      mutate(month=format(date,'%m'))
      
    urls<-datesToRetrieve %>% 
      mutate(product_name=paste0("3B-DAY.MS.MRG.3IMERG.",gsub("-","",date_character),"-S000000-E235959.V06.nc4")) %>%
      mutate(url_product=paste(OpenDAPServerUrl,OpenDAPCollection,year,month,product_name,sep="/"))
    
  }
  
  # To retrieve spatial indices
  OpenDAPURL<-urls$product_name[1]
    # Calculate XVector if not provided
  if(is.null(XVector) & is.null(roiSpatialIndexBound)){
    XVector<-getOpenDAPvector(OpenDAPURL,SpatialXVectorName)
  }
  # Calculate YVector if not provided
  if(is.null(YVector) & is.null(roiSpatialIndexBound)){
    YVector<-getOpenDAPvector(OpenDAPURL,SpatialYVectorName)
  }
  # Calculate roiSpatialIndexBound if not provided
  if(is.null(roiSpatialIndexBound)){
    roi_bbox<-sf::st_bbox(st_transform(roi,4326))
    Opendap_minLon<-which.min(abs(XVector-roi_bbox$xmin))-4
    Opendap_maxLon<-which.min(abs(XVector-roi_bbox$xmax))+4
    Opendap_minLat<-which.min(abs(YVector-roi_bbox$ymin))-4
    Opendap_maxLat<-which.min(abs(YVector-roi_bbox$ymax))+4
    roiSpatialIndexBound<-c(Opendap_minLat,Opendap_maxLat,Opendap_minLon,Opendap_maxLon)
  }
  
  # Build URL to download data in NetCDF format
  
  dim<-dimensionsToRetrieve %>%
    map(~paste0(.x,"[0:0][",roiSpatialIndexBound[3],":",roiSpatialIndexBound[4],"][",roiSpatialIndexBound[1],":",roiSpatialIndexBound[2],"],",SpatialYVectorName,"[",roiSpatialIndexBound[1],":",roiSpatialIndexBound[2],"],",SpatialXVectorName,"[",roiSpatialIndexBound[3],":",roiSpatialIndexBound[4],"]")) %>%
    unlist() %>%
    paste(collapse=",")
  
  table_urls<-urls %>%
    mutate(url=paste0(url_product,".nc4","?",dim)) %>%
    mutate(destfile=file.path(destFolder,paste0(OpenDAPCollection,product_name,".nc4"))) #%>%
    #mutate(names=)
  
  urls<-table_urls$url
  
  destfiles<-table_urls$destfile
  
  names<-table_urls$product_name
  
  #return(list(name=names,url=urls,destfile=destfiles))
  return(list(name=names,url=urls,destfile=destfiles))
  }
  
  
  
  
  
  
