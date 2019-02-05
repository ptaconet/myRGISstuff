extract_zonal_statistics<-function(path_to_input_vector_data,indices_to_compute,path_to_input_rasters,methods_to_compute){
  
  loc <- rgrass7::initGRASS(path_to_grassApplications_folder, home=getwd(), gisDbase="GRASS_TEMP", override=TRUE,mapset = "PERMANENT" )
  execGRASS("g.proj",flags="c",parameters = list(proj4="+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs"))
  #execGRASS("v.in.ogr", flags="o", parameters=list(input=path_to_output_vector_group_by_landcover_type, output="tmpvect",min_area=0.0001, snap=-1.0))
  execGRASS("v.in.ogr", flags=c("o","overwrite"), parameters=list(input=path_to_input_vector_data, output="tmpvect",min_area=0.0001, snap=-1.0))
  
  for (i in 1:length(indices_to_compute)){
    cat(paste0("Computing zonal statistics for indice ",indices_to_compute[i],"\n"))
    execGRASS("r.external", flags="overwrite", parameters=list(input=path_to_input_rasters[i], output="tmprast",band=1))
    execGRASS("g.region", parameters=list(raster="tmprast")) 
    execGRASS("v.rast.stats", flags=c("c","verbose"), parameters=list(map="tmpvect", raster="tmprast",column_prefix=indices_to_compute[i],method=methods_to_compute,percentile=90))
    #stats_df<-as.data.frame(readVECT("tmpvect"))
  }
  
  writeOGR(readVECT("tmpvect"),path_to_ground_truth_training_stats,driver = "GPKG",layer="ground_truth_training_stats")
  
  # Remove temorary folder
  system(paste0("rm -r ", file.path(getwd(),"GRASS_TEMP")))  
  file.remove(file.path(getwd(),".grassrc7"))

}