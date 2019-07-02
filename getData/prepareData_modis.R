
# To open a MODIS dataset that was downloaded via OpenDap
prepareData_modis<-function(path_to_raw_modis,var_name){
  grid_nc<-raster(path_to_raw_modis,varname=var_name)
  projection(grid_nc)<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
  extent(grid_nc)[1:2]<-extent(grid_nc)[1:2]+res(grid_nc)[1]/2
  extent(grid_nc)[3:4]<-extent(grid_nc)[3:4]-res(grid_nc)[1]/2
  return(grid_nc)
}