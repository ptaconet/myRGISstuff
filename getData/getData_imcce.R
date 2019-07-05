getData_imcce<-function(time_range=as.Date(c("2010-01-01","2010-01-30")),
                        roi=st_read("/home/ptaconet/r_react/getData/ROI_test.kml",quiet=T),
                        destFolder=NULL
){
  
  url_imcce_webservice<-"http://vo.imcce.fr/webservices/miriade/ephemcc_query.php?"
  
  if(length(time_range)==1){
    time_range=c(time_range,time_range) 
  }
  
  roi_bbox<-sf::st_bbox(st_transform(roi,4326))
  
  datesToRetrieve<-seq(from=time_range[2],to=time_range[1],by=-1) %>%
    as.data.frame() %>%
    set_names("date")
    
  table_urls<-datesToRetrieve %>%
    mutate(url=paste0(url_imcce_webservice,"-name=s:Moon&-type=Satellite&-ep=",date,"T23:30:00&-nbd=1d&-step=1h&-tscale=UTC&-observer=",mean(c(roi_bbox$xmin,roi_bbox$xmax)),"%20",mean(c(roi_bbox$ymin,roi_bbox$ymax)),"%200.0&-theory=INPOP&-teph=1&-tcoor=1&-mime=text/csv&-output=--jd&-extrap=0&-from=MiriadeDoc")) %>%
    mutate(product_name=gsub("-","",date)) %>%
    mutate(destfile=file.path(destFolder,paste0(product_name,".csv")))
  
  urls<-table_urls$url
  
  destfiles<-table_urls$destfile
  
  names<-table_urls$product_name
  
  return(list(name=names,url=urls,destfile=destfiles))
  
  
}