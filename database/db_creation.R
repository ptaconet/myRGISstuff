rm(list = ls())
require(RSQLite)
require(dplyr)
require(sf)
require(rgdal)
require(raster)
require(gdalUtils)
require(readxl)
require(lubridate)
require(stringr)

path_to_amal_database<-"/home/ptaconet/Documents/react/miscellaneous_data/React_dbase_V7.db"  
path_to_gpkg_database<-"/home/ptaconet/Documents/react/react_db.gpkg"  # Empty gpkg template is available here : http://www.geopackage.org/data/empty.gpkg
path_to_metadata_table<-"/home/ptaconet/r_react/database/metadata.csv"
path_to_metadata_mapping_table<-"/home/ptaconet/r_react/database/metadata_mapping.csv"
upload_lulc_rasters<-TRUE

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
source("/home/ptaconet/r_react/database/raw_villages.R") 
villages_sf<-st_as_sf(villages,coords =  c("Longitude", "Latitude"), crs = 4326 )
villages_sf<-cbind(villages_sf,st_coordinates(villages_sf))
st_write(villages_sf, path_to_gpkg_database, "raw_villages", update = TRUE)

# raw_menages
source("/home/ptaconet/r_react/database/raw_menages.R") 
menages_sf<-st_as_sf(menages,coords =  c("Longitude", "Latitude"), crs = 4326 )
menages_sf<-cbind(menages_sf,st_coordinates(menages_sf))
st_write(menages_sf, path_to_gpkg_database, "raw_menages", update = TRUE)

# raw_individus
source("/home/ptaconet/r_react/database/raw_individus.R") 
individus <- cbind(fid = 1:nrow(individus), individus)
dbWriteTable(react_gpkg,"raw_individus",individus)

# raw_capturedeterm
source("/home/ptaconet/r_react/database/raw_capturedeterm.R") 
df_capturedeterm <- cbind(fid = 1:nrow(df_capturedeterm), df_capturedeterm)
dbWriteTable(react_gpkg,"raw_capturedeterm",df_capturedeterm)

# villages_households_loc_pop
source("/home/ptaconet/r_react/database/villages_households_loc_pop.R")
df_households_loc_pop_sf<-st_as_sf(df_households_loc_pop,coords =  c("longitude", "latitude"), crs = 4326 )
df_households_loc_pop_sf<-cbind(df_households_loc_pop_sf,st_coordinates(df_households_loc_pop_sf))
st_write(df_households_loc_pop_sf, path_to_gpkg_database, "villages_households_loc_pop", update = TRUE)

# villages_loc_pop
source("/home/ptaconet/r_react/database/villages_loc_pop.R")
df_villages_loc_pop_sf<-st_as_sf(df_villages_loc_pop,coords =  c("X", "Y"), crs = 4326 )
df_villages_loc_pop_sf<-cbind(df_villages_loc_pop_sf,st_coordinates(df_villages_loc_pop_sf))
st_write(df_villages_loc_pop_sf, path_to_gpkg_database, "villages_loc_pop", update = TRUE)

# raw_dates_hlc
source("/home/ptaconet/r_react/database/raw_dates_hlc.R")
raw_dates_hlc <- cbind(fid = 1:nrow(raw_dates_hlc), raw_dates_hlc)
dbWriteTable(react_gpkg,"raw_dates_hlc",raw_dates_hlc)

# hlc_dates_loc_times
source("/home/ptaconet/r_react/database/hlc_dates_loc_times.R")
hlc_dates_loc_times_sf<-st_as_sf(hlc_dates_loc_times,coords =  c("longitude", "latitude"), crs = 4326 )
hlc_dates_loc_times_sf<-cbind(hlc_dates_loc_times_sf,st_coordinates(hlc_dates_loc_times_sf))
st_write(hlc_dates_loc_times_sf, path_to_gpkg_database, "hlc_dates_loc_times", update = TRUE)

# ancillary_africa_countries  downloaded here : https://data.humdata.org/dataset/west-and-central-africa-administrative-boundaries-levels
adm_bound_sf<-read_sf("/home/ptaconet/Documents/react/miscellaneous_data/wca_adm0/wca_adm0.shp")
adm_bound_sf<-st_transform(adm_bound_sf,crs=4326)
st_write(adm_bound_sf, path_to_gpkg_database, "ancillary_africa_countries", update = TRUE)

# ancillary_bound_project
roi_civ_sf<-read_sf("/home/ptaconet/Documents/react/data_CIV/ROI.kml")
roi_civ_sf$codepays<-"CI"
roi_civ_sf<-roi_civ_sf[,"codepays"]
roi_civ_sf<-st_cast(roi_civ_sf,"POLYGON")
roi_bf_sf<-read_sf("/home/ptaconet/Documents/react/data_BF/ROI.kml")
roi_bf_sf$codepays<-"BF"
roi_bf_sf<-roi_bf_sf[,"codepays"]
roi_bf_sf<-st_zm(roi_bf_sf,drop = TRUE, what = "ZM")
roi<-rbind(roi_civ_sf,roi_bf_sf)
roi<-st_transform(roi,crs=32630)
st_write(roi, path_to_gpkg_database, "ancillary_bound_project", update = TRUE)

# ancillary_africa_cities
ancillary_africa_cities<-read_sf("/home/ptaconet/Documents/react/miscellaneous_data/africa_places/places.shp")
ancillary_africa_cities<-st_intersection(ancillary_africa_cities,adm_bound_sf)
st_write(ancillary_africa_cities, path_to_gpkg_database, "ancillary_africa_cities", update = TRUE)

# LU/LC training and validation parcels (raw and segmented)
ground_truth_data_civ_raw<-st_read("/home/ptaconet/Documents/react/data_CIV/Ground_truth/civ_groundtruth_vector_32630.gpkg")
ground_truth_data_civ_raw <- cbind(pk = 1:nrow(ground_truth_data_civ_raw), ground_truth_data_civ_raw)
st_write(ground_truth_data_civ_raw, path_to_gpkg_database, "landcover_civ_groundtruth_raw", update = TRUE)
ground_truth_data_civ_revised<-st_read("/home/ptaconet/Documents/react/data_CIV/Ground_truth/civ_groundtruth_objects_segmentation_v_classes_update.gpkg")
st_write(ground_truth_data_civ_revised, path_to_gpkg_database, "landcover_civ_groundtruth_processed", update = TRUE)

ground_truth_data_bf_raw<-st_read("/home/ptaconet/Documents/react/data_BF/Ground_truth/bf_groundtruth_vector_32630.gpkg")
st_write(ground_truth_data_bf_raw, path_to_gpkg_database, "landcover_bf_groundtruth_raw", update = TRUE)
ground_truth_data_bf_revised<-st_read("/home/ptaconet/Documents/react/data_BF/Ground_truth/groundtruth_bf_v_classes_update.gpkg")
st_write(ground_truth_data_bf_revised, path_to_gpkg_database, "landcover_bf_groundtruth_processed", update = TRUE)

## LU/LC maps 
if (upload_lulc_rasters){
# BF
  cat("loading BF LU/LC rasters...\n")
path_to_LU_L1_bf<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L1.tif"
path_to_LU_L2_bf<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L2.tif"
path_to_LU_L3_bf<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L3.tif"
path_to_LU_L4_bf<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L4.tif"
path_to_LU_L5_bf<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L5.tif"
gdal_translate(path_to_LU_L1_bf,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_bf_L1"))
gdal_translate(path_to_LU_L2_bf,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_bf_L2"))
gdal_translate(path_to_LU_L3_bf,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_bf_L3"))
gdal_translate(path_to_LU_L4_bf,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_bf_L4"))
gdal_translate(path_to_LU_L5_bf,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_bf_L5"))

path_to_LU_L1_classes<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L1.csv"
path_to_LU_L2_classes<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L2.csv"
path_to_LU_L3_classes<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L3.csv"
path_to_LU_L4_classes<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L4.csv"
path_to_LU_L5_classes<-"/home/ptaconet/Documents/react/data_BF/Classification/classification_L5.csv"

LU_L1_classes<-read.csv(path_to_LU_L1_classes)
LU_L2_classes<-read.csv(path_to_LU_L2_classes)
LU_L3_classes<-read.csv(path_to_LU_L3_classes)
LU_L4_classes<-read.csv(path_to_LU_L4_classes)
LU_L5_classes<-read.csv(path_to_LU_L5_classes)

LU_L1_classes$classif_level<-"classification_L1"
LU_L2_classes$classif_level<-"classification_L2"
LU_L3_classes$classif_level<-"classification_L3"
LU_L4_classes$classif_level<-"classification_L4"
LU_L5_classes$classif_level<-"classification_L5"

LU_classes<-rbind(LU_L1_classes,LU_L2_classes,LU_L3_classes,LU_L4_classes,LU_L5_classes)
LU_classes <- LU_classes %>% arrange(classif_level,pixval)
LU_classes <- cbind(fid = 1:nrow(LU_classes), LU_classes)
dbWriteTable(react_gpkg,"landcover_bf_pixval2class",LU_classes)

# CIV
cat("loading CIV LU/LC rasters...\n")
path_to_LU_L1_civ<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L1.tif"
path_to_LU_L2_civ<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L2.tif"
path_to_LU_L3_civ<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L3.tif"
path_to_LU_L4_civ<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L4.tif"
path_to_LU_L5_civ<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L5.tif"
gdal_translate(path_to_LU_L1_civ,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_civ_L1"))
gdal_translate(path_to_LU_L2_civ,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_civ_L2"))
gdal_translate(path_to_LU_L3_civ,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_civ_L3"))
gdal_translate(path_to_LU_L4_civ,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_civ_L4"))
gdal_translate(path_to_LU_L5_civ,path_to_gpkg_database,ot="UInt16",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=landcover_civ_L5"))

path_to_LU_L1_classes<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L1.csv"
path_to_LU_L2_classes<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L2.csv"
path_to_LU_L3_classes<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L3.csv"
path_to_LU_L4_classes<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L4.csv"
path_to_LU_L5_classes<-"/home/ptaconet/Documents/react/data_CIV/Classification/classification_L5.csv"

LU_L1_classes<-read.csv(path_to_LU_L1_classes)
LU_L2_classes<-read.csv(path_to_LU_L2_classes)
LU_L3_classes<-read.csv(path_to_LU_L3_classes)
LU_L4_classes<-read.csv(path_to_LU_L4_classes)
LU_L5_classes<-read.csv(path_to_LU_L5_classes)

LU_L1_classes$classif_level<-"classification_L1"
LU_L2_classes$classif_level<-"classification_L2"
LU_L3_classes$classif_level<-"classification_L3"
LU_L4_classes$classif_level<-"classification_L4"
LU_L5_classes$classif_level<-"classification_L5"

LU_classes<-rbind(LU_L1_classes,LU_L2_classes,LU_L3_classes,LU_L4_classes,LU_L5_classes)
LU_classes <- LU_classes %>% arrange(classif_level,pixval)
LU_classes <- cbind(fid = 1:nrow(LU_classes), LU_classes)
dbWriteTable(react_gpkg,"landcover_civ_pixval2class",LU_classes)
}

# Pedology (raster)
path_to_pedology_civ<-"/home/ptaconet/Documents/react/data_CIV/pedology/pedo_final_32630.tif"
path_to_pedology_bf<-"/home/ptaconet/Documents/react/data_BF/pedology/pedo_final_32630.tif"
gdal_translate(path_to_pedology_civ,path_to_gpkg_database,ot="Float32",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=pedology_civ"))
gdal_translate(path_to_pedology_bf,path_to_gpkg_database,ot="Float32",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=pedology_bf"))


dbSendQuery(react_gpkg,"VACUUM") # It is very important to Vacuum. Not vacuuming may prevent the DB to be opened. 
dbDisconnect(react_gpkg)
