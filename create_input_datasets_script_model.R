require(RSQLite)
path_to_gpkg_database<-"/home/ptaconet/Documents/react/react_db.gpkg"
react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)
query_hlc_dates_loc<-paste0("SELECT idpointdecapture,codevillage as village,date_capture,X as longitude,Y as latitude FROM hlc_dates_loc_times WHERE codepays='",codepays,"'")
query_households_population<-paste0("SELECT codevillage as village,population,X as longitude,Y as latitude FROM households_loc_pop WHERE codepays='",codepays,"'")
if (codepays=="CI"){
  table_landcover_dataset<-""
  table_pedology_dataset<-"pedology_civ"
  table_road_network<-""
} else if (codepays=="BF"){
  table_landcover_dataset<-""
  table_pedology_dataset<-"pedology_bf"
  table_road_network<-""
}
hlc_dates_loc<-dbGetQuery(react_gpkg, query_hlc_dates_loc)
households_population<-dbGetQuery(react_gpkg, query_households_population)

write.csv(hlc_dates_loc,file.path(path_to_processing_folder,path_to_csv_hlc_dates_loc),row.names = F)
write.csv(households_population,file.path(path_to_processing_folder,path_to_csv_households_population),row.names = F)

gdal_translate(path_to_gpkg_database,path_to_pedology_dataset,ot="Float32",of="GTiff",oo=paste0("TABLE=",table_pedology_dataset))


dbDisconnect(react_gpkg)
