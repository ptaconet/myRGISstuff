#require(RSQLite)
#require(dplyr)
#amal_db <- dbConnect(RSQLite::SQLite(),path_to_amal_database)

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
colnames(df_capturedeterm)<-gsub("_fk","",colnames(df_capturedeterm))
colnames(df_capturedeterm)<-gsub("_pk","",colnames(df_capturedeterm))
