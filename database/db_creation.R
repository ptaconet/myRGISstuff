rm(list = ls())
require(RSQLite)
require(dplyr)
require(sf)
require(rgdal)
require(raster)
require(gdalUtils)
require(readxl)

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
all_data <- cbind(pk = 1:nrow(all_data), all_data)
supervcapture_sf<-st_as_sf(all_data,coords =  c("longitude", "latitude"), crs = 4326 )
supervcapture_sf<-cbind(supervcapture_sf,st_coordinates(supervcapture_sf))
st_write(supervcapture_sf, path_to_gpkg_database, "raw_supervcapture", update = TRUE)

# raw_villages
query<-"SELECT * FROM village"
villages<-dbGetQuery(amal_db, query)
villages$Latitude <- as.numeric(gsub(",",".",villages$Latitude))
villages$Longitude <- as.numeric(gsub(",",".",villages$Longitude))
villages$Latitude[which(villages$codevillage_pk=="NAM")]<-8.8845
villages$Longitude[which(villages$codevillage_pk=="NAM")]<--5.75
villages$Latitude[which(villages$codevillage_pk=="KOL")]<-9.288
villages$Longitude[which(villages$codevillage_pk=="KOL")]<--5.524
villages$Latitude[which(villages$codevillage_pk=="BLA")]<-8.948
villages$Longitude[which(villages$codevillage_pk=="BLA")]<--5.652
villages$Latitude[which(is.na(villages$Latitude))]<-9
villages$Longitude[which(is.na(villages$Longitude))]<--5
villages$nomvillage[which(villages$codevillage_pk=="BLA")]<-"Blawara"
villages$nomvillage[which(villages$codevillage_pk=="NAM")]<-"Namasselikaha"
colnames(villages)<-gsub("_fk","",colnames(villages))
colnames(villages)<-gsub("_pk","",colnames(villages))
villages_sf<-st_as_sf(villages,coords =  c("Longitude", "Latitude"), crs = 4326 )
villages_sf<-cbind(villages_sf,st_coordinates(villages_sf))
st_write(villages_sf, path_to_gpkg_database, "raw_villages", update = TRUE)

# raw_menages
query<-"SELECT * FROM menage"
menages<-dbGetQuery(amal_db, query)
menages$Latitude <- as.numeric(gsub(",",".",menages$Latitude))
menages$Longitude <- as.numeric(gsub(",",".",menages$Longitude))
menages$codevillage_fk[which(menages$codevillage_fk=="KOL" & menages$codemenage_pk %in% c("LOK062","LOK012"))]<-"LOK"
index_menages_nam<-which(is.na(menages$Latitude) & menages$codevillage_fk=="NAA")
codemenage_menages_nam<-menages$codemenage_pk[index_menages_nam]
menages$codemenage_pk[index_menages_nam]<-gsub("NAA","NAM",menages$codemenage_pk[index_menages_nam])
menages$codevillage_fk[index_menages_nam]="NAM"
index_menages_bla<-which(is.na(menages$Latitude) & menages$codevillage_fk=="KOL")
codemenage_menages_bla<-menages$codemenage_pk[index_menages_bla]
menages$codemenage_pk[index_menages_bla]<-gsub("KOL","BLA",menages$codemenage_pk[index_menages_bla])
menages$codevillage_fk[index_menages_bla]="BLA"
menages$Latitude[which(is.na(menages$Latitude))]<-9
menages$Longitude[which(is.na(menages$Longitude))]<--5
colnames(menages)<-gsub("_fk","",colnames(menages))
colnames(menages)<-gsub("_pk","",colnames(menages))
menages_sf<-st_as_sf(menages,coords =  c("Longitude", "Latitude"), crs = 4326 )
menages_sf<-cbind(menages_sf,st_coordinates(menages_sf))
st_write(menages_sf, path_to_gpkg_database, "raw_menages", update = TRUE)

# raw_individus
query<-"SELECT * FROM individu"
individus<-dbGetQuery(amal_db, query)
individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_nam)]<-gsub("NAA","NAM",individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_nam)])
individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_nam)]<-gsub("NAA","NAM",individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_nam)])
individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_bla)]<-gsub("KOL","BLA",individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_bla)])
individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_bla)]<-gsub("KOL","BLA",individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_bla)])
individus <- cbind(fid = 1:nrow(individus), individus)
colnames(individus)<-gsub("_fk","",colnames(individus))
colnames(individus)<-gsub("_pk","",colnames(individus))
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
df_capturedeterm <- cbind(fid = 1:nrow(df_capturedeterm), df_capturedeterm)
colnames(df_capturedeterm)<-gsub("_fk","",colnames(df_capturedeterm))
colnames(df_capturedeterm)<-gsub("_pk","",colnames(df_capturedeterm))
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

# raw_dates_hlc
raw_bf_dates_hlc<-read_excel("/home/ptaconet/Bureau/reprise_data_amal/Heure de captures REACT-BF_M1_M7.xls")
colnames(raw_bf_dates_hlc)<-c("n_mission","nomvillage","codevillage","CSPS","date_de_captures","heure_de_debut","heure_de_fin","n_sac","observations")
raw_bf_dates_hlc$date_de_captures<-as.character(as.Date(raw_bf_dates_hlc$date_de_captures, format="%d/%m/%Y"))
raw_bf_dates_hlc$codepays<-"BF"
raw_bf_dates_hlc$heure_de_debut<-gsub("H",":",raw_bf_dates_hlc$heure_de_debut)
raw_bf_dates_hlc$heure_de_debut<-gsub("h",":",raw_bf_dates_hlc$heure_de_debut)
raw_bf_dates_hlc$heure_de_fin<-gsub("H",":",raw_bf_dates_hlc$heure_de_fin)
raw_bf_dates_hlc$heure_de_fin<-gsub("h",":",raw_bf_dates_hlc$heure_de_fin)
raw_bf_dates_hlc$heure_de_debut<-paste0(raw_bf_dates_hlc$heure_de_debut,":00")
raw_bf_dates_hlc$heure_de_fin<-paste0(raw_bf_dates_hlc$heure_de_fin,":00")
raw_bf_dates_hlc$date_heure_debut<-paste(raw_bf_dates_hlc$date_de_captures,raw_bf_dates_hlc$heure_de_debut,sep=" ")
raw_bf_dates_hlc$date_heure_fin<-paste(as.Date(raw_bf_dates_hlc$date_de_captures)+1,raw_bf_dates_hlc$heure_de_fin,sep=" ")
raw_civ_dates_hlc<-read_excel("/home/ptaconet/Bureau/reprise_data_amal/Dates_capture entomo cote  d'Ivoire.xlsx")
colnames(raw_civ_dates_hlc)<-c("date_de_captures","nomvillage","n_mission")
raw_civ_dates_hlc$date_de_captures<-as.character(raw_civ_dates_hlc$date_de_captures)
raw_civ_dates_hlc$nomvillage <- raw_civ_dates_hlc$nomvillage %>% str_replace_all(c("Penatiguikaha" = "Penatiguikaha_Gopko", "Logaha" = "Lokaha","Kolékaha"= "Kolekaha","Lagomokaha"= "Lagomounkaha","Yenessonkaha"= "Yenessonkaha_Gofionkaha","Narlougokaha"= "Nalourgokala","Katiorpokaha"= "Katiorkpo","Kogninguekaha"= "Koguin","Tagbarakaha"= "Tagbara","Nongotakaha"= "Nongotanakaha","Karafiné"= "Karafine","Nongowélékaha"= "Nangowelekaha","Kougniguékaha"= "Koungniguékaha","Félékaha"= "Felekaha","Tahouélékaha"= "Tahouelekaha","Blaouara"= "Blawara" ))
query<-"SELECT codepays, nomvillage ,codevillage FROM raw_villages WHERE codepays='CI'"
villages<-dbGetQuery(react_gpkg, query)
raw_civ_dates_hlc<-left_join(raw_civ_dates_hlc,villages)
raw_civ_dates_hlc$CSPS<-raw_civ_dates_hlc$heure_de_debut<-raw_civ_dates_hlc$heure_de_fin<-raw_civ_dates_hlc$n_sac<-raw_civ_dates_hlc$observations<-raw_civ_dates_hlc$date_heure_debut<-raw_civ_dates_hlc$date_heure_fin<-NA
raw_dates_hlc<-rbind(raw_bf_dates_hlc,raw_civ_dates_hlc)
raw_dates_hlc <- cbind(fid = 1:nrow(raw_dates_hlc), raw_dates_hlc)
dbWriteTable(react_gpkg,"raw_dates_hlc",raw_dates_hlc)

# LU/LC training and validation parcels (raw and segmented)
ground_truth_data_civ_raw<-st_read("/home/ptaconet/Documents/react/data_CIV/Ground_truth/civ_groundtruth_vector_32630.gpkg")
ground_truth_data_civ_raw <- cbind(pk = 1:nrow(ground_truth_data_civ_raw), ground_truth_data_civ_raw)
st_write(ground_truth_data_civ_raw, path_to_gpkg_database, "lu_lc_ground_truth_civ_raw", update = TRUE)
ground_truth_data_civ_revised<-st_read("/home/ptaconet/Documents/react/data_CIV/Ground_truth/civ_groundtruth_objects_segmentation.gpkg")
st_write(ground_truth_data_civ_revised, path_to_gpkg_database, "lu_lc_ground_truth_civ_rev", update = TRUE)

ground_truth_data_bf_raw<-st_read("/home/ptaconet/Documents/react/data_BF/Ground_truth/bf_groundtruth_vector_32630.gpkg")
st_write(ground_truth_data_bf_raw, path_to_gpkg_database, "lu_lc_ground_truth_bf_raw", update = TRUE)
ground_truth_data_bf_revised<-st_read("/home/ptaconet/Documents/react/data_CIV/Ground_truth/civ_groundtruth_objects_segmentation.gpkg")
st_write(ground_truth_data_bf_revised, path_to_gpkg_database, "lu_lc_ground_truth_bf_rev", update = TRUE)

# Pedology (raster)
path_to_pedology_civ<-"/home/ptaconet/Documents/react/data_CIV/pedology/pedo_final_32630.tif"
path_to_pedology_bf<-"/home/ptaconet/Documents/react/data_BF/pedology/pedo_final_32630.tif"
gdal_translate(path_to_pedology_civ,path_to_gpkg_database,ot="Float32",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=pedology_civ"))
gdal_translate(path_to_pedology_bf,path_to_gpkg_database,ot="Float32",of="GPKG",b=1,co=c("APPEND_SUBDATASET=YES","RASTER_TABLE=pedology_bf"))



# hlc_locations_dates
source("/home/ptaconet/r_react/database/hlc_locations_dates.R")




#dbDisconnect(react_gpkg)
