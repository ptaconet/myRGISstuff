# To open a GPM dataset that was downloaded via OpenDap
prepareData_gpm<-function(path_to_raw_gpm,var_name,resample=FALSE,resample_output_res=NULL){
  
  gpm_rast<-raster(path_to_raw_gpm,varname=var_name)
  projection(gpm_rast)<-CRS("+init=epsg:4326")
  # The raster has to be flipped. Output was validated with the data from 2017-09-20 (see https://docserver.gesdisc.eosdis.nasa.gov/public/project/GPM/browse/GPM_3IMERGDF.png)
  gpm_rast <- t(gpm_rast)
  gpm_rast <- flip(gpm_rast,'y')
  gpm_rast <- flip(gpm_rast,'x')
  
  if(resample){
    resample_output_res<-convertMetersToDegrees(resample_output_res,latitude_4326=mean(c(extent(gpm_rast)[3],extent(gpm_rast)[4])))
    r<-gpm_rast
    res(r)<-resample_output_res
    gpm_rast<-resample(gpm_rast,r,method='bilinear')
  }
  
  return(gpm_rast)
  
}