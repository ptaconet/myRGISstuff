#require(RSQLite)
#require(dplyr)
#amal_db <- dbConnect(RSQLite::SQLite(),path_to_amal_database)

query<-"SELECT * FROM menage"
menages<-dbGetQuery(amal_db, query)
menages$Latitude <- as.numeric(gsub(",",".",menages$Latitude))
menages$Longitude <- as.numeric(gsub(",",".",menages$Longitude))
menages$codevillage_fk[which(menages$codevillage_fk=="KOL" & menages$codemenage_pk %in% c("LOK062","LOK012"))]<-"LOK"
menages$codevillage_fk[which(menages$codemenage_pk=="PEN035")]<-"PEN"
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
