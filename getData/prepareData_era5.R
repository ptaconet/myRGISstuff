
prepareData_era5<-function(path_to_raw_era5,var_name,resample=FALSE,resample_output_res=NULL){
  
  # Open netcdf
  nc <- nc_open(path_to_raw_era5)
  # Get lat and lon
  lat <- ncvar_get(nc,'latitude')
  lon <- ncvar_get(nc,'longitude')
  # Get time 
  t <- ncvar_get(nc, "time")
  #time unit: hours since 1900-01-01
  #ncatt_get(nc,'time')
  #convert the hours into date + hour
  #as_datetime() function of the lubridate package needs seconds
  timestamp <- as_datetime(c(t*60*60),origin="1900-01-01")
  
  # Get the variable
  variable <- ncvar_get(nc,var_name)
  
  brick_era<-NULL
  for (i in 1:dim(variable)[3]){
    rast <- raster(t(variable[,,i]), xmn=min(lon)-0.125, xmx=max(lon)+0.125, ymn=min(lat)-0.125, ymx=max(lat)+0.125, crs=CRS("+init=epsg:4326"))
  # Resample to size_output_grid size 
  if(resample){
    resample_output_res<-convertMetersToDegrees(resample_output_res,latitude_4326=mean(c(extent(rast)[3],extent(rast)[4])))
    r<-rast
    res(r)<-resample_output_res
    rast<-raster::resample(rast,r,method='bilinear')
  }
    brick_era<-c(brick_era,variable_rast)
}
  
  
}

