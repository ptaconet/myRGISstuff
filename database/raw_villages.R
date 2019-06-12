#require(RSQLite)
#require(dplyr)
#amal_db <- dbConnect(RSQLite::SQLite(),path_to_amal_database)

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
