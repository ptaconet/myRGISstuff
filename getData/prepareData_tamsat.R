# To prepare a TAMSAT dataset 
prepareData_tamsat<-function(path_to_raw_tamsat,roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T),resample=FALSE,resample_output_res=NULL){
  
  rast<-raster(path_to_raw_tamsat)
  
  # extend a bit the size of the bbox
  bbox_tamsat<-extent(roi)
  bbox_tamsat[1]=bbox_tamsat[1]-0.5
  bbox_tamsat[2]=bbox_tamsat[2]+0.5
  bbox_tamsat[3]=bbox_tamsat[3]-0.5
  bbox_tamsat[4]=bbox_tamsat[4]+0.5
  
  # Crop to the bbox
  rast<-raster(path_to_raw_tamsat) %>% 
    crop(bbox_tamsat)

  if(resample){
    resample_output_res<-convertMetersToDegrees(resample_output_res,latitude_4326=mean(c(extent(rast)[3],extent(rast)[4])))
    r<-rast
    res(r)<-resample_output_res
    rast<-raster::resample(rast,r,method='bilinear')
  }
  
  return(rast)
  
}