getData_srtm<-function(roi,
                       download=FALSE, # TRUE will download the file and return a list with : the URL, the path to the output file, a boolean wether the dataset was properly downloaded or not. FALSE will return a list with the URL only
                       destFolder=NULL,
                       ){
  
  url_srtm_server<-"http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/"
  srtm_tiles<-getSRTMtileNames(roi)
  urls<-paste0(url_srtm_server,srtm_tiles,".SRTMGL1.hgt.zip")
    
  if (download){
    destfiles<-file.path(destFolder,paste0(srtm_tiles,".SRTMGL1.hgt.zip"))
    res<-downloadData(urls,destfiles)
    urls<-res$urls
    destfiles<-res$destfiles
    status<-res$status
    # unzip  
    destfiles %>% map(~unzip(.,exdir = destFolder))
  } else {
    destfiles<-res<-NULL
  }

return(list(urls,destfiles,res))
}
