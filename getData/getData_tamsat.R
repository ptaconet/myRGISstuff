getData_tamsat<-function(time_range=as.Date(c("2017-01-01","2017-01-30")), # mandatory. either a time range (e.g. c(date_start,date_end) ) or a single date e.g. ( date_start ) 
                         destFolder=NULL,
                         output_time_step="daily", # {daily,monthly}
                         output_product="rainfall_estimate", # {rainfall_estimate,anomaly (only if output_time_step==monthly) ,climatology (only if output_time_step==monthly)}
                         output_output="individual" # {individual,yearly} 
                         
){
  # see https://www.tamsat.org.uk/data/archive
  
  url_tamsat_data<-"https://www.tamsat.org.uk/public_data/TAMSAT3"
  
  time_range<-as.Date(time_range, origin="1970-01-01")
  
  datesToRetrieve<-seq(time_range[2],time_range[1],-1) %>%
    data.frame(stringsAsFactors = F) %>%
    set_names("date") %>%
    mutate(date_character=as.character(as.Date(date))) %>%
    mutate(year=format(date,'%Y')) %>%
    mutate(month=format(date,'%m')) %>%
    mutate(day=format(date,'%d'))
    
  urls<-datesToRetrieve %>%
    mutate(product_name_daily_rain_individual=paste0("rfe",year,"_",month,"_",day,".v3.nc")) %>%
    mutate(url_product_daily_rain_individual=paste0(url_tamsat_data,"/",year,"/",month,"/",product_name_daily_rain_individual)) %>%
    mutate(product_name_daily_rain_yearly=paste0("TAMSATv3.0_rfe_daily_",year,".zip")) %>%
    mutate(url_product_daily_rain_yearly=paste0(url_tamsat_data,"/zip/",product_name_daily_rain_yearly)) %>%
    mutate(product_name_monthly_rain_individual=paste0("rfe",year,"_",month,".v3.nc")) %>%
    mutate(url_product_monthly_rain_individual=paste0(url_tamsat_data,"/",year,"/",month,"/",product_name_monthly_rain_individual)) %>%
    mutate(product_name_monthly_rain_yearly=paste0("TAMSATv3.0_rfe_monthly_",year,".zip")) %>%
    mutate(url_product_monthly_rain_yearly=paste0(url_tamsat_data,"/zip/",product_name_monthly_rain_yearly)) %>%
    mutate(product_name_monthly_anomaly_individual=paste0("rfe",year,"_",month,"_anom.v3.nc")) %>%
    mutate(url_product_monthly_anomaly_individual=paste0(url_tamsat_data,"/",year,"/",month,"/",product_name_monthly_anomaly_individual)) #%>%
    #mutate(product_name_monthly_climatology_individual=paste0("rfe",year,"-",year,"_",month,"_clim.v3.nc")) %>%
    #mutate(url_product_monthly_climatology_individual=paste0(url_tamsat_data,"/clim","/",month,"/",product_name_monthly_climatology_individual))


  if (output_time_step=="daily" & output_product=="rainfall_estimate" & output_output=="individual"){
    urls <- urls %>% select(product_name_daily_rain_individual,url_product_daily_rain_individual)
  } else if (output_time_step=="daily" & output_product=="rainfall_estimate" & output_output=="yearly"){
    urls <- urls %>% select(product_name_daily_rain_yearly,url_product_daily_rain_yearly)
  } else if (output_time_step=="monthly" & output_product=="rainfall_estimate" & output_output=="individual"){
    urls <- urls %>% select(product_name_monthly_rain_individual,url_product_monthly_rain_individual)
  } else if (output_time_step=="monthly" & output_product=="rainfall_estimate" & output_output=="yearly"){
    urls <- urls %>% select(product_name_monthly_rain_yearly,url_product_monthly_rain_yearly)
  } else if (output_time_step=="monthly" & output_product=="anomaly" & output_output=="individual"){
    urls <- urls %>% select(product_name_monthly_anomaly_individual,url_product_monthly_anomaly_individual)
  } #else if (output_time_step=="monthly" & output_product=="climatology" & output_output=="individual"){
    #urls <- urls %>% select(product_name_monthly_climatology_individual,url_product_monthly_climatology_individual)
  #}
  
  urls$destfile <- file.path(destFolder,urls[,1])
  
  return(list(name=urls[,1],url=urls[,2],destfile=urls[,3]))
  
}