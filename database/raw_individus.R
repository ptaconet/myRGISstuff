#require(RSQLite)
#require(dplyr)
#amal_db <- dbConnect(RSQLite::SQLite(),path_to_amal_database)

query<-"SELECT * FROM individu"
individus<-dbGetQuery(amal_db, query)
individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_nam)]<-gsub("NAA","NAM",individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_nam)])
individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_nam)]<-gsub("NAA","NAM",individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_nam)])
individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_bla)]<-gsub("KOL","BLA",individus$codeindividu_pk[which(individus$codemenage_fk %in% codemenage_menages_bla)])
individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_bla)]<-gsub("KOL","BLA",individus$codemenage_fk[which(individus$codemenage_fk %in% codemenage_menages_bla)])
colnames(individus)<-gsub("_fk","",colnames(individus))
colnames(individus)<-gsub("_pk","",colnames(individus))
