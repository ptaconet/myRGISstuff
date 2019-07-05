

getData_era5<-function(time_range=as.Date(c("2010-01-01","2010-01-30")), # mandatory. either a time range (e.g. c(date_start,date_end) ) or a single date e.g. ( date_start )
                        roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T), # either provide roi (sf point or polygon) or provide roiSpatialIndexBound. if roiSpatialIndexBound is not provided, it will be calculated from roi
                        username=NULL, # EarthData user name
                        password=NULL, # EarthData password
                        destFolder=NULL,
                        dimensionsToRetrieve=c("10m_u_component_of_wind","10m_v_component_of_wind") # mandatory
){
  
  
  # Check : https://dominicroye.github.io/en/2018/access-to-climate-reanalysis-data-from-r/
  
  # ERA 5 : https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
  
  # Description of the wind data: https://apps.ecmwf.int/codes/grib/param-db?id=165 and https://apps.ecmwf.int/codes/grib/param-db?id=166
  
  
  roi_bbox<-sf::st_bbox(st_transform(roi,4326))

  time_range=as.POSIXlt(time_range)
  
  
  datesToRetrieve<-seq(from=time_range[2],to=time_range[1],by="-1 hour") %>%
    data.frame(stringsAsFactors = F) %>%
    set_names("date") %>%
    mutate(date_character=as.character(as.Date(date))) %>%
    mutate(year=format(date,'%Y')) %>%
    mutate(month=format(date,'%m')) %>%
    mutate(day=format(date,"%d")) %>%
    mutate(hour=format(date,"%H"))
    #mutate(hour_era5=as_datetime(c(t*60*60),origin="1900-01-01")) %>%

  #we create the query for the date of catch, hours 18h to 23h
  query <- r_to_py(list(
    variable= dimensionsToRetrieve,
    product_type= "reanalysis",
    year= sort(unique(datesToRetrieve$year)),
    month=sort(unique(datesToRetrieve$month)), #formato: "01","01", etc.
    day= sort(unique(datesToRetrieve$day)), #stringr::str_pad(1:31,2,"left","0"),   
    time= sort(unique(datesToRetrieve$hour)),
    format= "netcdf",
    area = paste0(roi_bbox$ymax+1,"/",roi_bbox$xmin-1,"/",roi_bbox$ymin-1,"/",roi_bbox$xmax+1) # North, West, South, East
  ))
  
  name=gsub("-","",min(as.character(as.Date(datesToRetrieve$date))))
  destfile=file.path(destFolder,paste0(name,".nc"))
    
  return(list(name=name,url=query,destfile=destfile))
  
  #query the server to get the ncdf for date of catch and date of catch + 1
  #server$retrieve("reanalysis-era5-single-levels",
   #               query_this_date_hlc,
    #              "/home/ptaconet/Documents/react/data_CIV/ERA_WIND/test.nc")
  
  
}