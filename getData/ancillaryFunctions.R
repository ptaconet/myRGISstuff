# Set of functions to calculate various information given a ROI 


getUTMepsg<-function(roi){
  
  bbox<-st_bbox(roi)
  #  cat("Warning: ROIs overlapping more than 1 UTM zone are currently not adapted in this workflow\n")
  utm_zone_number<-(floor((bbox$xmin + 180)/6) %% 60) + 1
  if(bbox$ymin>0){ # if latitudes are North
    epsg<-as.numeric(paste0("326",utm_zone_number))
  } else { # if latitude are South
    epsg<-as.numeric(paste0("325",utm_zone_number))
  }
  
  return(epsg)
}


getMODIStileNames<-function(roi){
  
  modis_tile = read_sf("https://modis.ornl.gov/files/modis_sin.kmz") %>% 
    st_intersection(roi) %>% 
    as.data.frame() %>%
    dplyr::select(Name) %>%
    as.character()
  
  if(length(unique(modis_tile))>1){
    stop("Your ROI is overlapping more than 1 MODIS tile. This workflow is currently not adapted for this case\n")
  } else {
    modis_tile<-modis_tile %>%
      unique() %>%
      stringr::str_replace_all(c(" "="",":"=""))
    for (i in 1:9){
      modis_tile<-gsub(paste0("h",i,"v"),paste0("h0",i,"v"),modis_tile)
    }
    if(nchar(modis_tile)!=6){
      modis_tile<-paste0(substr(modis_tile,1,4),"0",substr(modis_tile,5,5))
    }
  }
  
  return(modis_tile)
}


getSRTMtileNames<-function(roi){
  
  srtm_tiles <- geojsonsf::geojson_sf("http://dwtkns.com/srtm30m/srtm30m_bounding_boxes.json")  %>%
    sf::st_intersection(roi) %>%
    as.data.frame()

  SRTMtileNames<-substr(srtm_tiles$dataFile,1,7)
  
  return(SRTMtileNames)
}


convertMetersToDegrees<-function(length_meters,
                                 latitude_4326){

  length_degrees <- length_meters / (111.32 * 1000 * cos(latitude_4326 * ((pi / 180))))
  
  return(length_degrees)
}


downloadData<-function(urls,destfiles,username=NULL,password=NULL,parallelDL=FALSE){
  
  # check which data is already downloaded
  data_dl<-data.frame(url=urls,destfile=destfiles,stringsAsFactors = F) %>%
      mutate(fileExist=map_lgl(destfile,file.exists)) %>%
      mutate(status=ifelse(fileExist==TRUE,3,NA))
  
  # data already downloaded
  data_already_exist<-data_dl %>%
    filter(fileExist==TRUE)
    
  # data to download
  data_to_download<-data_dl %>%
    filter(fileExist==FALSE)

    # download data
  #for (i in 1:nrow(data_to_download)){
  #    httr::GET(data_to_download$url[i],httr::authenticate(username,password),write_disk(data_to_download$destfile[i]))
  # }
  
  if(is.null(username)){
    username<-password<-"no_auth"
  }
   dl_func<-function(url,output,username,password) {httr::GET(url,httr::authenticate(username,password),httr::write_disk(output),httr::progress())}
  
  if (parallelDL){
  require(parallel)
  cl <- makeCluster(detectCores())
  clusterMap(cl, dl_func, url=data_to_download$url,output=data_to_download$destfile,username=username,password=password,
             .scheduling = 'dynamic')
  stopCluster(cl)
  } else {
    for (i in 1:nrow(data_to_download)){
      dl_func(url=data_to_download$url[i],output=data_to_download$destfile[i],username=username,password=password)
    }
  }
  
    data_dl<-data_to_download %>%
    mutate(fileExist=map_lgl(destfile,file.exists)) %>%
    mutate(status=ifelse(fileExist==TRUE,1,2))  %>%
    rbind(data_already_exist)
  
  # 1 : download ok
  # 2 : download error
  # 3 : data already existing in output folder

    return(data_dl)
}
