getData_viirsDnb<-function(time_range=as.Date(c("2017-01-01","2017-01-30")), # mandatory. either a time range (e.g. c(date_start,date_end) ) or a single date e.g. ( date_start )
                           roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T), # either provide roi (sf point or polygon) or provide roiSpatialIndexBound. if roiSpatialIndexBound is not provided, it will be calculated from roi
                           dimensionsToRetrieve=c("Monthly_AvgRadiance","Monthly_CloudFreeCoverage"), # mandatory
                           destFolder=NULL
){
  
  url_noaa_nighttime_webservice<-"https://gis.ngdc.noaa.gov/arcgis/rest/services/NPP_VIIRS_DNB/"
  
  roi_bbox<-sf::st_bbox(st_transform(roi,4326))
  
  time_range=as.Date(time_range,origin="1970-01-01")
  
  if(length(time_range)==1){
    time_range=c(time_range,time_range %m+% days(1)) 
  }
  
  datesToRetrieve<-seq(from=time_range[2],to=time_range[1],by="-1 month") %>%
    data.frame(stringsAsFactors = F) %>%
    set_names("date") %>%
    mutate(date_character=as.character(as.Date(date))) %>%
    mutate(year=format(date,'%Y')) %>%
    mutate(month=format(date,'%m')) %>%
    mutate(date_start=as.Date(paste(year,month,"01",sep="-"))) %>%
    mutate(date_end=date_start %m+% months(1))  %>%
    mutate(time_start=as.integer(difftime(date_start ,"1970-01-01" , units = c("secs")))*1000) %>%
    mutate(time_end=as.integer(difftime(date_end ,"1970-01-01" , units = c("secs")))*1000)
    
  table_urls<-datesToRetrieve %>%
    select(year,month,time_start,time_end) %>%
    slice(rep(1:n(), each = length(dimensionsToRetrieve))) %>%
    mutate(dimensionsToRetrieve=rep(dimensionsToRetrieve,n()/2)) %>%
    mutate(url=paste0(url_noaa_nighttime_webservice,dimensionsToRetrieve,"/ImageServer/exportImage?bbox=",roi_bbox$xmin,",",roi_bbox$ymin,",",roi_bbox$xmax,",",roi_bbox$ymax,"&time=",format(time_start,scientific=FALSE),",",format(time_end,scientific=FALSE),"&format=tiff&f=image")) %>%
    mutate(product_name=paste0(dimensionsToRetrieve,"_",year,month)) %>%
    mutate(destfile=file.path(destFolder,paste0(product_name,".tif")))
  
    
  urls<-table_urls$url
  
  destfiles<-table_urls$destfile
  
  names<-table_urls$product_name
  
  return(list(name=names,url=urls,destfile=destfiles))
  
  
}
