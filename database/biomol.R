require(readxl)
require(dplyr)
require(RSQLite)
require(lubridate)

path_to_xls_civ<-"/home/ptaconet/Documents/react/miscellaneous_data/Rapportdactivité_LMB_IRD-IRP_10au16Décembre2018.xls"
path_to_xls_bf<-"/home/ptaconet/Documents/react/miscellaneous_data/Version  06_Trans_Entomo_REACT-BF.xls"
path_to_amal_database<-"/home/ptaconet/Documents/react/miscellaneous_data/React_dbase_V7.db"  
path_to_gpkg_database<-"/home/ptaconet/Documents/react/react_db.gpkg"  # Empty gpkg template is available here : http://www.geopackage.org/data/empty.gpkg

amal_db <- dbConnect(RSQLite::SQLite(),path_to_amal_database)
react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)

# rapatrier les fichiers excel
biomol_civ<-read_excel(path_to_xls_civ, sheet = "Rapport_Biomol_20181220")
biomol_bf<-read_excel(path_to_xls_bf, sheet = "Anopheles sp_Lab Analyse_OK_upd")

# rapatrier la table biomol de la base d'amal
query<-"SELECT * FROM biomol"
df_biomol_db<-dbGetQuery(amal_db, query)

# rapatrier les villages de la nouvelle base de données
query<-"SELECT codepays,codevillage FROM raw_villages"
df_villages<-dbGetQuery(react_gpkg, query)

# nb de lignes de biomol par pays
df_biomol_db$codevillage<-substr(df_biomol_db$idmoustique_pk_fk,2,4)
df_biomol_db<-left_join(df_biomol_db,df_villages)
df_biomol_db %>% group_by(codepays) %>% summarise(count=n())

# on harmonise les noms de colonne des 2 fichiers excel
colnames(biomol_civ)[which(colnames(biomol_civ)=="code")]<-"nummoustique"
colnames(biomol_civ)[which(colnames(biomol_civ)=="enquete")]<-"nummission"
colnames(biomol_civ)[which(colnames(biomol_civ)=="date")]<-"date_determination"
colnames(biomol_civ)[which(colnames(biomol_civ)=="village")]<-"codevillage"
colnames(biomol_civ)[which(colnames(biomol_civ)=="point capture")]<-"pointdecapture"
colnames(biomol_civ)[which(colnames(biomol_civ)=="int/ext")]<-"postedecapture"
colnames(biomol_civ)[which(colnames(biomol_civ)=="heure")]<-"heuredecapture"
colnames(biomol_civ)[which(colnames(biomol_civ)=="Traitement")]<-"traitement"
colnames(biomol_civ)[which(colnames(biomol_civ)=="Parturité")]<-"parturite"
colnames(biomol_civ)[which(colnames(biomol_civ)=="Obs.")]<-"observations"
colnames(biomol_civ)[which(colnames(biomol_civ)=="anophele vecteur")]<-"espece"
colnames(biomol_civ)[which(colnames(biomol_civ)=="Espèces")]<-"pcr_espece"
colnames(biomol_civ)[which(colnames(biomol_civ)=="Plasmodium sp.")]<-"pcr_pf"
colnames(biomol_civ)[which(colnames(biomol_civ)=="kdr")]<-"kdrw"
colnames(biomol_civ)[which(colnames(biomol_civ)=="Formes moléc.")]<-"forme_moleculaire"
colnames(biomol_civ)[which(colnames(biomol_civ)=="Ace1")]<-"ace1"

colnames(biomol_bf)[which(colnames(biomol_bf)=="NumMoustique")]<-"nummoustique"
colnames(biomol_bf)[which(colnames(biomol_bf)=="NumMission")]<-"nummission"
colnames(biomol_bf)[which(colnames(biomol_bf)=="DateDétermination")]<-"date_determination"
colnames(biomol_bf)[which(colnames(biomol_bf)=="CodeVillage")]<-"codevillage"
colnames(biomol_bf)[which(colnames(biomol_bf)=="PointDeCapture")]<-"pointdecapture"
colnames(biomol_bf)[which(colnames(biomol_bf)=="PosteDeCapture")]<-"postedecapture"
colnames(biomol_bf)[which(colnames(biomol_bf)=="HeureDeCapture")]<-"heuredecapture"
colnames(biomol_bf)[which(colnames(biomol_bf)=="Traitement")]<-"traitement"
colnames(biomol_bf)[which(colnames(biomol_bf)=="TauxParturite")]<-"parturite"
colnames(biomol_bf)[which(colnames(biomol_bf)=="Remarques_SDD")]<-"observations"
colnames(biomol_bf)[which(colnames(biomol_bf)=="EspeceAnophele")]<-"espece"
colnames(biomol_bf)[which(colnames(biomol_bf)=="PcrEspece")]<-"pcr_espece"
colnames(biomol_bf)[which(colnames(biomol_bf)=="PcrPf")]<-"pcr_pf"
colnames(biomol_bf)[which(colnames(biomol_bf)=="kdrw")]<-"kdrw"
colnames(biomol_bf)[which(colnames(biomol_bf)=="kdre")]<-"kdre"
colnames(biomol_bf)[which(colnames(biomol_bf)=="IdentifiantMoustiq")]<-"idmoustique"

# on reconstitue l'idmoustique depuis les données des fichiers excel
#biomol_bf$idmoustique2<-paste0(biomol_bf$nummission,biomol_bf$codevillage,biomol_bf$pointdecapture,biomol_bf$postedecapture,sprintf("%02d",biomol_bf$heuredecapture),substr(biomol_bf$Genre,1,2),substr(biomol_bf$espece,1,1),substr(biomol_bf$EtatAbdomen,1,1),biomol_bf$nummoustique)
#a<-which(biomol_bf$idmoustique!=biomol_bf$idmoustique2)

biomol_bf$postedecapture<-substr(biomol_bf$postedecapture,1,1)
biomol_civ$postedecapture<-substr(biomol_civ$postedecapture,1,1)
biomol_civ$heuredecapture<-substr(biomol_civ$heuredecapture,1,2)

## on remplace les occurences dans les colonnes de biomol_civ, pour aligner les noms dans les colonnes de biomol_bf
biomol_civ$postedecapture <- biomol_bf$datefin %>% str_replace_all(c("" = "", ))


