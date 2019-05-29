rm(list = ls())
require(RSQLite)
require(dplyr)
require(sf)
require(rgdal)
require(raster)
require(gdalUtils)

path_to_amal_database<-"/home/ptaconet/Bureau/React_dbase_V7.db"  
path_to_gpkg_database<-"/home/ptaconet/Bureau/react_db.gpkg"  # Empty gpkg template is available here : http://www.geopackage.org/data/empty.gpkg
path_to_metadata_table<-"/home/ptaconet/r_react/database/metadata.csv"
path_to_metadata_mapping_table<-"/home/ptaconet/r_react/database/metadata_mapping.csv"

download.file("http://www.geopackage.org/data/empty.gpkg",path_to_gpkg_database)

## Connect to DBs
amal_db <- dbConnect(RSQLite::SQLite(),path_to_amal_database)
react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)

## Open and load tables metadata and metadata_mapping
metadata<-read.csv(path_to_metadata_table)
metadata_mapping<-read.csv(path_to_metadata_mapping_table)

dbWriteTable(react_gpkg,"metadata",metadata,append=TRUE)
dbWriteTable(react_gpkg,"metadata_mapping",metadata_mapping,append=TRUE)

## Create and load tables :

# raw_supervcapture
source("/home/ptaconet/r_react/database/raw_supervcapture.R") # créé le data.frame nommé "all_data" contenant les données brutes de supervision capture
supervcapture_sf<-st_as_sf(all_data,coords =  c("longitude", "latitude"), crs = 4326 )
supervcapture_sf<-cbind(supervcapture_sf,st_coordinates(supervcapture_sf))
st_write(supervcapture_sf, path_to_gpkg_database, "raw_supervcapture", update = TRUE)

# raw_villages
query<-"SELECT * FROM village"
villages<-dbGetQuery(amal_db, query)
villages$Latitude <- as.numeric(gsub(",",".",villages$Latitude))
villages$Longitude <- as.numeric(gsub(",",".",villages$Longitude))
villages$Latitude[which(is.na(villages$Latitude))]<-9
villages$Longitude[which(is.na(villages$Longitude))]<--5
villages_sf<-st_as_sf(villages,coords =  c("Longitude", "Latitude"), crs = 4326 )
villages_sf<-cbind(villages_sf,st_coordinates(villages_sf))
st_write(villages_sf, path_to_gpkg_database, "raw_villages", update = TRUE)

# raw_menages
query<-"SELECT * FROM menage"
menages<-dbGetQuery(amal_db, query)
menages$Latitude <- as.numeric(gsub(",",".",menages$Latitude))
menages$Longitude <- as.numeric(gsub(",",".",menages$Longitude))
menages$Latitude[which(is.na(menages$Latitude))]<-9
menages$Longitude[which(is.na(menages$Longitude))]<--5
menages_sf<-st_as_sf(menages,coords =  c("Longitude", "Latitude"), crs = 4326 )
menages_sf<-cbind(menages_sf,st_coordinates(menages_sf))
st_write(menages_sf, path_to_gpkg_database, "raw_menages", update = TRUE)

# raw_individus
query<-"SELECT * FROM individu"
individus<-dbGetQuery(amal_db, query)
dbWriteTable(react_gpkg,"raw_individus",individus)

# raw_capturedeterm
query<-"SELECT * FROM capturedeterm"
df_capturedeterm<-dbGetQuery(amal_db, query)
query<-"SELECT * FROM capturedeterm_ci_niv1"
df_capturedeterm_ci_niv1<-dbGetQuery(amal_db, query)
# on aligne les noms de colonne
df_capturedeterm_ci_niv1$identifiant<-NA
df_capturedeterm_ci_niv1$baro_id<-NA
df_capturedeterm_ci_niv1$row_id_pk<-NULL
colnames(df_capturedeterm)[which(colnames(df_capturedeterm)=="idmoustique_pk")]="idmoustique"
df_capturedeterm<-rbind(df_capturedeterm,df_capturedeterm_ci_niv1)
df_capturedeterm$postedecapture<-gsub("int","i",df_capturedeterm$postedecapture)
df_capturedeterm$postedecapture<-gsub("ext","e",df_capturedeterm$postedecapture)
df_capturedeterm$idpostedecapture<-paste0(df_capturedeterm$enquete,df_capturedeterm$codevillage_fk,df_capturedeterm$pointdecapture,df_capturedeterm$postedecapture)
dbWriteTable(react_gpkg,"raw_capturedeterm",df_capturedeterm)

# households_loc_pop
source("/home/ptaconet/r_react/database/households_loc_pop.R")
df_households_loc_pop_sf<-st_as_sf(df_households_loc_pop,coords =  c("longitude", "latitude"), crs = 4326 )
df_households_loc_pop_sf<-cbind(df_households_loc_pop_sf,st_coordinates(df_households_loc_pop_sf))
st_write(df_households_loc_pop_sf, path_to_gpkg_database, "households_loc_pop", update = TRUE)

# villages_loc_pop
source("/home/ptaconet/r_react/database/villages_loc_pop.R")
df_villages_loc_pop_sf<-st_as_sf(df_villages_loc_pop,coords =  c("X", "Y"), crs = 4326 )
df_villages_loc_pop_sf<-cbind(df_villages_loc_pop_sf,st_coordinates(df_villages_loc_pop_sf))
st_write(df_villages_loc_pop_sf, path_to_gpkg_database, "villages_loc_pop", update = TRUE)

# raw_bf_dates_hlc
query<-"SELECT * FROM paul_dates_captures_par_village_bf"
raw_bf_dates_hlc<-dbGetQuery(amal_db, query)
raw_bf_dates_hlc$date_de_captures<-as.character(as.Date(raw_bf_dates_hlc$date_de_captures, format="%d/%m/%Y"))
dbWriteTable(react_gpkg,"raw_bf_dates_hlc",raw_bf_dates_hlc)

# LU/LC training and validation parcels (raw and segmented)
ground_truth_data_civ_raw<-st_read("/home/ptaconet/Documents/react/data_CIV/Ground_truth/civ_groundtruth_vector_32630.gpkg")
st_write(ground_truth_data_civ_raw, path_to_gpkg_database, "lu_lc_ground_truth_civ_raw", update = TRUE)
ground_truth_data_civ_revised<-st_read("/home/ptaconet/Documents/react/data_CIV/Ground_truth/civ_groundtruth_objects_segmentation.gpkg")
st_write(ground_truth_data_civ_revised, path_to_gpkg_database, "lu_lc_ground_truth_civ_rev", update = TRUE)

ground_truth_data_bf_raw<-st_read("/home/ptaconet/Documents/react/data_BF/Ground_truth/bf_groundtruth_vector_32630.gpkg")
st_write(ground_truth_data_bf_raw, path_to_gpkg_database, "lu_lc_ground_truth_bf_raw", update = TRUE)
ground_truth_data_bf_revised<-st_read("/home/ptaconet/Documents/react/data_BF/Ground_truth/groundtruth_bf.gpkg")
st_write(ground_truth_data_bf_revised, path_to_gpkg_database, "lu_lc_ground_truth_bf_rev", update = TRUE)

# Pedology (raster)
path_to_pedology_civ<-"/home/ptaconet/Documents/react/data_CIV/pedology/pedo_final_32630.tif"
path_to_pedology_bf<-"/home/ptaconet/Documents/react/data_BF/pedology/pedo_final_32630.tif"
gdal_translate(path_to_pedology_civ,path_to_gpkg_database,ot="Float32",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=pedology_civ"))
gdal_translate(path_to_pedology_bf,path_to_gpkg_database,ot="Float32",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=pedology_bf"))



# hlc_locations_dates
source("/home/ptaconet/r_react/database/hlc_locations_dates.R")




#dbDisconnect(react_gpkg)
