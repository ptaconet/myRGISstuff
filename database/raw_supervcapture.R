require(readxl)
require(dplyr)
require(stringr)


##############################################################################################
# CIV  (BF en dessous)
##############################################################################################
path_to_xls_civ<-"/home/ptaconet/Bureau/reprise_data_amal/DataSupervCaptureCI_jointes_orig_all.xlsx"

lot1<-read_excel(path_to_xls_civ, sheet = "Lot1")
lot2<-read_excel(path_to_xls_civ, sheet = "Lot2")
lot3<-read_excel(path_to_xls_civ, sheet = "Lot3")
lot4<-read_excel(path_to_xls_civ, sheet = "Lot4")
lot5<-read_excel(path_to_xls_civ, sheet = "Lot5")
lot6<-read_excel(path_to_xls_civ, sheet = "Lot6")
lot7<-read_excel(path_to_xls_civ, sheet = "Lot7")
lot8<-read_excel(path_to_xls_civ, sheet = "Lot8",col_types = "text")
lot9<-read_excel(path_to_xls_civ, sheet = "lot9")
lot10<-read_excel(path_to_xls_civ, sheet = "lot10")

# Align colum names with lot1
### Lot 1, 3, 4, 5, 6, 8  are treated the same way
# lot3 and lot4 have exactly the same column names and numbers
# lot1 and lot5  have exactly the same column names and numbers

all_column_names<-sort(unique(c(colnames(lot1),colnames(lot2),colnames(lot3),colnames(lot4),colnames(lot5),colnames(lot6),colnames(lot7),colnames(lot8),colnames(lot9),colnames(lot10))))


# Keep good columns only
lot1<-lot1 %>% filter(metainstanceid==parentuid)
colnames(lot3)[which(colnames(lot3)=="meta-instanceID")]<-"instanceID"
lot3<-lot3 %>% filter(instanceID==PARENT_KEY)
colnames(lot4)[which(colnames(lot4)=="meta-instanceID")]<-"instanceID"
lot4<-lot4 %>% filter(instanceID==PARENT_KEY)
lot5<-lot5 %>% filter(metainstanceid==parentuid)

lot6$postedecapture[which(lot6$postedecapture=="Ext")]<-"e"
lot6$postedecapture[which(lot6$postedecapture=="Int")]<-"i"
lot6$idpostedecapture<-paste0(lot6$NumMission,lot6$codevillage,lot6$pointdecapture,lot6$postedecapture)

## MANQUE LE NUMERO DE MISSION DANS LOT8 ET LOT9 ET LOT10
lot8$postedecapture[which(lot8$postedecapture=="Ext")]<-"e"
lot8$postedecapture[which(lot8$postedecapture=="Int")]<-"i"
#lot8$idpostedecapture<-paste0(lot8$NumMission,lot8$codevillage,lot8$pointdecapture,lot8$postedecapture)

lot9$postedecapture[which(lot9$postedecapture=="Ext")]<-"e"
lot9$postedecapture[which(lot9$postedecapture=="Int")]<-"i"
#lot9$idpostedecapture<-paste0(lot9$NumMission,lot9$codevillage,lot9$pointdecapture,lot9$postedecapture)

lot10$postedecapture[which(lot10$postedecapture=="Ext")]<-"e"
lot10$postedecapture[which(lot10$postedecapture=="Int")]<-"i"
#lot9$idpostedecapture<-paste0(lot9$NumMission,lot9$codevillage,lot9$pointdecapture,lot9$postedecapture)


# lot2
colnames(lot2)[which(colnames(lot2)=="meta-instanceName")]<-"instanceName"
colnames(lot2)[which(colnames(lot2)=="coordgpspointdecapture-Latitude")]<-"coordgpspointdecaptureLatitude"
colnames(lot2)[which(colnames(lot2)=="coordgpspointdecapture-Longitude")]<-"coordgpspointdecaptureLongitude"
colnames(lot2)[which(colnames(lot2)=="coordgpspointdecapture-Altitude")]<-"coordgpspointdecaptureAltitude"
colnames(lot2)[which(colnames(lot2)=="coordgpspointdecapture-Accuracy")]<-"coordgpspointdecaptureAccuracy"
lot2$F25<-NA
lot2$parentuid<-NA
lot2$metainstanceid<-NA
lot2$metamodelversion<-NA
lot2$metauiversion<-NA
lot2$metasubmissiondate<-NA
lot2$metaiscomplete<-NA
lot2$metadatemarkedascomplete<-NA
lot2$instanceID<-NA
lot2$deviceid<-NA


# lot3
colnames(lot3)[which(colnames(lot3)=="meta-instanceName")]<-"instanceName"
colnames(lot3)[which(colnames(lot3)=="SubmissionDate")]<-"metasubmissiondate"
colnames(lot3)[which(colnames(lot3)=="meta-instanceID")]<-"instanceID"
colnames(lot3)[which(colnames(lot3)=="device_id")]<-"deviceid"
colnames(lot3)[which(colnames(lot3)=="coordgpspointdecapture-Latitude")]<-"coordgpspointdecaptureLatitude"
colnames(lot3)[which(colnames(lot3)=="coordgpspointdecapture-Longitude")]<-"coordgpspointdecaptureLongitude"
colnames(lot3)[which(colnames(lot3)=="coordgpspointdecapture-Altitude")]<-"coordgpspointdecaptureAltitude"
colnames(lot3)[which(colnames(lot3)=="coordgpspointdecapture-Accuracy")]<-"coordgpspointdecaptureAccuracy"
lot3$'SupervCaptureLot3part1_SET-OF-Part2'<-NULL
lot3$'SupervCaptureLot3part1_KEY'<-NULL
lot3$PARENT_KEY<-NULL
lot3$SupervCaptureLot3part2_KEY<-NULL
lot3$'SupervCaptureLot3part2_SET-OF-Part2'<-NULL
lot3$F25<-NA
lot3$parentuid<-NA
lot3$metainstanceid<-NA
lot3$metamodelversion<-NA
lot3$metauiversion<-NA
lot3$metaiscomplete<-NA
lot3$metadatemarkedascomplete<-NA


# lot4
colnames(lot4)[which(colnames(lot4)=="meta-instanceName")]<-"instanceName"
colnames(lot4)[which(colnames(lot4)=="SubmissionDate")]<-"metasubmissiondate"
colnames(lot4)[which(colnames(lot4)=="meta-instanceID")]<-"instanceID"
colnames(lot4)[which(colnames(lot4)=="device_id")]<-"deviceid"
colnames(lot4)[which(colnames(lot4)=="coordgpspointdecapture-Latitude")]<-"coordgpspointdecaptureLatitude"
colnames(lot4)[which(colnames(lot4)=="coordgpspointdecapture-Longitude")]<-"coordgpspointdecaptureLongitude"
colnames(lot4)[which(colnames(lot4)=="coordgpspointdecapture-Altitude")]<-"coordgpspointdecaptureAltitude"
colnames(lot4)[which(colnames(lot4)=="coordgpspointdecapture-Accuracy")]<-"coordgpspointdecaptureAccuracy"
lot4$'SupervCaptureLot7part1_KEY'<-NULL
lot4$'SupervCaptureLot7part2_KEY'<-NULL
lot4$PARENT_KEY<-NULL
lot4$SupervCaptureLot3part2_KEY<-NULL
lot4$'SupervCaptureLot7part1_SET-OF-Part2'<-NULL
lot4$'SupervCaptureLot7part2_SET-OF-Part2'<-NULL
lot4$F25<-NA
lot4$parentuid<-NA
lot4$metainstanceid<-NA
lot4$metamodelversion<-NA
lot4$metauiversion<-NA
lot4$metaiscomplete<-NA
lot4$metadatemarkedascomplete<-NA

# lot5
# nothing to do, same columns

# lot6
colnames(lot6)[which(colnames(lot6)=="NumMission")]<-"nummission"
lot6$F25<-NA
lot6$parentuid<-NA
lot6$rcpidpostedecapture<-NA
lot6$rcpcoordgpsptdecapture<-NA
lot6$utilisationbonsac<-NA

# lot7
colnames(lot7)[which(colnames(lot7)=="meta-instanceName")]<-"instanceName"
colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Latitude")]<-"coordgpspointdecaptureLatitude"
colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Longitude")]<-"coordgpspointdecaptureLongitude"
colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Altitude")]<-"coordgpspointdecaptureAltitude"
colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Accuracy")]<-"coordgpspointdecaptureAccuracy"
colnames(lot7)<-gsub("Part1-","",colnames(lot7))
colnames(lot7)<-gsub("Part2-","",colnames(lot7))
lot7$metainstanceid<-NA
lot7$metamodelversion<-NA
lot7$metauiversion<-NA
lot7$metasubmissiondate<-NA
lot7$metaiscomplete<-NA
lot7$metadatemarkedascomplete<-NA
lot7$instanceID<-NA
lot7$deviceid<-NA
lot7$nummission<-NA
lot7$F25<-NA
lot7$parentuid<-NA
lot7$rcpidpostedecapture<-NA
lot7$rcpcoordgpsptdecapture<-NA
lot7$utilisationbonsac<-NA


#lot8
colnames(lot8)[which(colnames(lot8)=="coordgpspointdecapture-Latitude")]<-"coordgpspointdecaptureLatitude"
colnames(lot8)[which(colnames(lot8)=="coordgpspointdecapture-Longitude")]<-"coordgpspointdecaptureLongitude"
colnames(lot8)[which(colnames(lot8)=="coordgpspointdecapture-Altitude")]<-"coordgpspointdecaptureAltitude"
colnames(lot8)[which(colnames(lot8)=="coordgpspointdecapture-Accuracy")]<-"coordgpspointdecaptureAccuracy"
lot8$KEY<-NULL
colnames(lot8)[which(colnames(lot8)=="SubmissionDate")]<-"metasubmissiondate"
colnames(lot8)[which(colnames(lot8)=="meta-instanceID")]<-"instanceID"
colnames(lot8)[which(colnames(lot8)=="meta-instanceName")]<-"instanceName"
colnames(lot8)[which(colnames(lot8)=="device_id")]<-"deviceid"
lot8$metainstanceid<-NA
lot8$metamodelversion<-NA
lot8$metauiversion<-NA
lot8$metaiscomplete<-NA
lot8$metadatemarkedascomplete<-NA
lot8$nummission<-NA
lot8$F25<-NA
lot8$parentuid<-NA
lot8$rcpidpostedecapture<-NA
lot8$rcpcoordgpsptdecapture<-NA
lot8$utilisationbonsac<-NA


#lot9
colnames(lot9)[which(colnames(lot9)=="coordgpspointdecapture-Latitude")]<-"coordgpspointdecaptureLatitude"
colnames(lot9)[which(colnames(lot9)=="coordgpspointdecapture-Longitude")]<-"coordgpspointdecaptureLongitude"
colnames(lot9)[which(colnames(lot9)=="coordgpspointdecapture-Altitude")]<-"coordgpspointdecaptureAltitude"
colnames(lot9)[which(colnames(lot9)=="coordgpspointdecapture-Accuracy")]<-"coordgpspointdecaptureAccuracy"
lot9$metainstanceid<-NA
lot9$metamodelversion<-NA
lot9$metauiversion<-NA
lot9$metasubmissiondate<-NA
lot9$metaiscomplete<-NA
lot9$metadatemarkedascomplete<-NA
lot9$instanceID<-NA
lot9$instanceName<-NA
lot9$deviceid<-NA
lot9$nummission<-NA
lot9$F25<-NA
lot9$parentuid<-NA
lot9$rcpidpostedecapture<-NA
lot9$rcpcoordgpsptdecapture<-NA
lot9$utilisationbonsac<-NA

#lot10
lot10$heurecapture2<-NULL  # heurecapture2 n'a que des NULL
lot10$metainstanceid<-NA
lot10$metamodelversion<-NA
lot10$metauiversion<-NA
lot10$metasubmissiondate<-NA
lot10$metaiscomplete<-NA
lot10$metadatemarkedascomplete<-NA
lot10$instanceID<-NA
lot10$instanceName<-NA
lot10$deviceid<-NA
lot10$nummission<-NA
lot10$F25<-NA
lot10$parentuid<-NA
lot10$rcpidpostedecapture<-NA
lot10$rcpcoordgpsptdecapture<-NA
lot10$utilisationbonsac<-NA
lot10$datedebut<-lot10$datecapture
lot10$datefin<-lot10$datecapture



## Set date formats
lot1$datedebut <- as.POSIXct(lot1$datedebut, format="%m/%d/%Y %H:%M:%OS")
lot1$datefin <- as.POSIXct(lot1$datefin, format="%m/%d/%Y %H:%M:%OS")

lot2$datedebut <- lot2$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot2$datedebut <- as.POSIXct(lot2$datedebut, format="%d/%m/%Y %H:%M:%OS")
lot2$datefin <- lot2$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot2$datefin <- as.POSIXct(lot2$datefin, format="%d/%m/%Y %H:%M:%OS")

lot3$datedebut <- lot3$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot3$datedebut <- as.POSIXct(lot3$datedebut, format="%d/%m/%Y %H:%M:%OS")
lot3$datefin <- lot3$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot3$datefin <- as.POSIXct(lot3$datefin, format="%d/%m/%Y %H:%M:%OS")

lot4$datedebut <- lot4$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot4$datedebut <- as.POSIXct(lot4$datedebut, format="%d/%m/%Y %H:%M:%OS")
lot4$datefin <- lot4$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot4$datefin <- as.POSIXct(lot4$datefin, format="%d/%m/%Y %H:%M:%OS")

lot5$datedebut <- as.POSIXct(lot5$datedebut, format="%m/%d/%Y %H:%M:%OS")
lot5$datefin <- as.POSIXct(lot5$datefin, format="%m/%d/%Y %H:%M:%OS")

lot6$datedebut <- as.POSIXct(lot6$datedebut, format="%m/%d/%Y %H:%M:%OS")
lot6$datefin <- as.POSIXct(lot6$datefin, format="%m/%d/%Y %H:%M:%OS")

lot7$datedebut <- lot7$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot7$datedebut <- as.POSIXct(lot7$datedebut, format="%d/%m/%Y %H:%M:%OS")
lot7$datefin <- lot7$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" ,"avr." = "/04/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot7$datefin <- as.POSIXct(lot7$datefin, format="%d/%m/%Y %H:%M:%OS")


# lot8
# isolate bad dates
lot8_bad_dates<- lot8[which(substr(lot8$datedebut,1,2)=="42"),]
lot8_good_dates<- lot8[which(substr(lot8$datedebut,1,2)!="42"),]

lot8_good_dates$datedebut <- lot8_good_dates$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "fÃ©vr." = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "dÃ©c." = "/12/" ,"avr." = "/04/", "fÃ©vr." = "/02/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot8_good_dates$datedebut <- as.POSIXct(lot8_good_dates$datedebut, format="%d/%m/%Y %H:%M:%OS")
lot8_good_dates$datefin <- lot8_good_dates$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "fÃ©vr." = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "dÃ©c." = "/12/" ,"avr." = "/04/", "fÃ©vr." = "/02/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot8_good_dates$datefin <- as.POSIXct(lot8_good_dates$datefin, format="%d/%m/%Y %H:%M:%OS")

lot8_bad_dates$datedebut<- as.POSIXlt("1900-01-01 00:00:00", tz="GMT") + as.numeric(lot8_bad_dates$datedebut)*3600*24
lot8_bad_dates$datefin<- as.POSIXlt("1900-01-01 00:00:00", tz="GMT") + as.numeric(lot8_bad_dates$datefin)*3600*24

lot8<-rbind(lot8_good_dates,lot8_bad_dates)

# lot9
# isolate bad dates
lot9_bad_dates<- lot9[which(substr(lot9$datedebut,1,2)=="42"),]
lot9_good_dates<- lot9[which(substr(lot9$datedebut,1,2)!="42"),]

lot9_good_dates$datedebut <- lot9_good_dates$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "févr." = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "déc." = "/12/" ,"avr." = "/04/", "fÃ©vr." = "/02/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot9_good_dates$datedebut <- as.POSIXct(lot9_good_dates$datedebut, format="%d/%m/%Y %H:%M:%OS")
lot9_good_dates$datefin <- lot9_good_dates$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "févr." = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "déc." = "/12/" ,"avr." = "/04/", "fÃ©vr." = "/02/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
lot9_good_dates$datefin <- as.POSIXct(lot9_good_dates$datefin, format="%d/%m/%Y %H:%M:%OS")

lot9_bad_dates$datedebut<- as.POSIXlt("1900-01-01 00:00:00", tz="GMT") + as.numeric(lot9_bad_dates$datedebut)*3600*24
lot9_bad_dates$datefin<- as.POSIXlt("1900-01-01 00:00:00", tz="GMT") + as.numeric(lot9_bad_dates$datefin)*3600*24

lot9<-rbind(lot9_good_dates,lot9_bad_dates)


# lot10 ok

colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Latitude")]<-"coordgpspointdecaptureLatitude"
colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Longitude")]<-"coordgpspointdecaptureLongitude"
colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Altitude")]<-"coordgpspointdecaptureAltitude"
colnames(lot7)[which(colnames(lot7)=="coordgpspointdecapture-Accuracy")]<-"coordgpspointdecaptureAccuracy"
all_data_civ<-rbind(lot1,lot2,lot3,lot4,lot5,lot6,lot7,lot8,lot9,lot10)

all_data_civ$latitude<-as.numeric(all_data_civ$'coordgpspointdecaptureLatitude')
all_data_civ$longitude<-as.numeric(all_data_civ$'coordgpspointdecaptureLongitude')
all_data_civ$altitude<-as.numeric(all_data_civ$'coordgpspointdecaptureAltitude')
all_data_civ$accuracy<-as.numeric(all_data_civ$'coordgpspointdecaptureAccuracy')
all_data_civ$pointdecapture<-as.numeric(all_data_civ$pointdecapture)

all_data_civ$datecapture<-as.Date(all_data_civ$datedebut)

# Subset columns useful for us
all_data_civ_subset<-all_data_civ  %>% dplyr::select(datedebut,datefin,datecapture,codesuperviseur,codevillage,nummission,pointdecapture,idpointdecapture,postedecapture,idpostedecapture,idappdemes,latitude,longitude,altitude,accuracy,idcaptureur,heurecapture,respectpointdecapture,captureurauposte,captureureveille,posicaptureuradequate,habitcaptureuradequat,materielcapturedispo,sacsheureprecedtecollectes,utilisationbonsac,presenceappdemes,bnedispoappdemes,observations)

all_data_civ_subset<-unique(all_data_civ_subset)

# Get the data where we do not have the nummission and set the mission number for these data
all_data_civ_no_nummission<-all_data_civ_subset %>% filter(is.na(nummission) | nummission>8)
all_data_civ_no_nummission$nummission[which(all_data_civ_no_nummission$datecapture<="2016-10-23")]<-1
all_data_civ_no_nummission$nummission[which(all_data_civ_no_nummission$datecapture>="2016-11-16" & all_data_civ_no_nummission$datecapture<="2016-12-07")]<-2
all_data_civ_no_nummission$nummission[which(all_data_civ_no_nummission$datecapture>="2017-02-14" & all_data_civ_no_nummission$datecapture<="2017-03-07")]<-3
all_data_civ_no_nummission$nummission[which(all_data_civ_no_nummission$datecapture>="2017-03-21" & all_data_civ_no_nummission$datecapture<="2017-04-29")]<-4
all_data_civ_no_nummission$nummission[which(all_data_civ_no_nummission$datecapture>="2017-11-05" & all_data_civ_no_nummission$datecapture<="2017-11-07")]<-5
all_data_civ_no_nummission$nummission[which(all_data_civ_no_nummission$datecapture>="2018-03-17" & all_data_civ_no_nummission$datecapture<="2018-03-17")]<-8

all_data_civ_with_nummission<-all_data_civ_subset %>% filter(!(is.na(nummission) | nummission>8))
all_data_civ_subset<-rbind(all_data_civ_with_nummission,all_data_civ_no_nummission)

all_data_civ_subset$idpostedecapture<-paste0(all_data_civ_subset$nummission,all_data_civ_subset$codevillage,all_data_civ_subset$pointdecapture,all_data_civ_subset$postedecapture)

# nombre de lignes par poste de capture (ie nombre de passages par poste de capture)
number_lines_by_postecapture<-all_data_civ_subset %>% group_by(idpostedecapture) %>% summarise(n=n())
hist(number_lines_by_postecapture$n)


## On veut vérifier les coordonnées des points de capture
# Pour cela 1) on calcule la mediane des coordonéees de chaque idpostedecapture puis 2) on calcule la médiane de la distance à ces coordonnées médianes

#coords_median_postecapture <- all_data_civ_subset %>% group_by(idpostedecapture) %>% summarize(median_lat=median(latitude,na.rm = T),median_lon=median(longitude,na.rm = T)) #summarise(stdev_dist=sd(sqrt(latitude-longitude)^2)*(111.32 * 1000 * cos(mean(all_data_bf_subset$latitude) * ((pi / 180)))),count=n())  #  
#all_data_civ_subset<-left_join(all_data_civ_subset,coords_median_postecapture,by="idpostedecapture")
#all_data_civ_subset$dist_to_median<-sqrt((all_data_civ_subset$latitude-all_data_civ_subset$median_lat)^2+(all_data_civ_subset$longitude-all_data_civ_subset$median_lon)^2)
#median_dist_to_median_coords<-all_data_civ_subset %>% group_by(idpostedecapture) %>% summarise(median_dist_to_median_coords=median(dist_to_median)*(111.32 * 1000 * cos(mean(all_data_civ_subset$latitude,na.rm=TRUE)) * ((pi / 180))))

#coords_postecapture$lat[which(coords_postecapture$idpostedecapture == "3LAT2i")]<-coords_postecapture$lat[which(coords_postecapture$idpostedecapture == "3LAT2e")]
#coords_postecapture$lon[which(coords_postecapture$idpostedecapture == "3LAT2i")]<-coords_postecapture$lon[which(coords_postecapture$v == "3LAT2e")]


##############################################################################################
# BF
##############################################################################################

path_to_xls_bf<-"/home/ptaconet/Bureau/reprise_data_amal/All_DAta_SupervCaptureBF.xlsx"
path_to_data_mission1<-"/home/ptaconet/Bureau/reprise_data_amal/SupervCapture1BF.csv"
path_to_data_data_h_capt<-"/home/ptaconet/Bureau/reprise_data_amal/Heure de captures REACT-BF_M1_M7.xls"

Lot1part1 <- read_excel(path_to_xls_bf, sheet = "Lot1part1")
Lot1part2 <- read_excel(path_to_xls_bf, sheet = "Lot1part2")
Lot1_merge<-merge(Lot1part1,Lot1part2,by="SET-OF-Part2")

Lot2part1 <- read_excel(path_to_xls_bf, sheet = "Lot2part1")
Lot2part1$datedebut<-as.character(Lot2part1$datedebut)
Lot2part1$datefin<-as.character(Lot2part1$datefin)
Lot2part1$datecapture<-as.character(Lot2part1$datecapture)
Lot2part2 <- read_excel(path_to_xls_bf, sheet = "Lot2part2")
Lot2part1$rcpcoordgpsptdecapture <- paste0(Lot2part1$coordgpspointdecaptureLatitude," ",Lot2part1$coordgpspointdecaptureLongitude," ",Lot2part1$coordgpspointdecaptureAltitude,".0 ",Lot2part1$coordgpspointdecaptureAccuracy,".0")
Lot2_merge<-merge(Lot2part1,Lot2part2,by="rcpcoordgpsptdecapture")

Lot3part1 <- read_excel(path_to_xls_bf, sheet = "Lot3part1")
Lot3part2 <- read_excel(path_to_xls_bf, sheet = "Lot3part2")
Lot3_merge<-merge(Lot3part1,Lot3part2,by="SET-OF-Part2")

Lot4part1 <- read_excel(path_to_xls_bf, sheet = "Lot4part1")
Lot4part2 <- read_excel(path_to_xls_bf, sheet = "Lot4part2")
Lot4_merge<-merge(Lot4part1,Lot4part2,by.x="metainstanceid",by.y="parentuid")

# Lot1_merge et Lot3_merge ont les memes noms de colonne. On aligne Lot2_merge et Lot4_merge sur ces noms de colonnes
colnames(Lot2_merge)[which(colnames(Lot2_merge)=="instanceName")]<-"meta-instanceName"
colnames(Lot2_merge)[which(colnames(Lot2_merge)=="coordgpspointdecaptureLatitude")]<-"coordgpspointdecapture-Latitude"
colnames(Lot2_merge)[which(colnames(Lot2_merge)=="coordgpspointdecaptureLongitude")]<-"coordgpspointdecapture-Longitude"
colnames(Lot2_merge)[which(colnames(Lot2_merge)=="coordgpspointdecaptureAltitude")]<-"coordgpspointdecapture-Altitude"
colnames(Lot2_merge)[which(colnames(Lot2_merge)=="coordgpspointdecaptureAccuracy")]<-"coordgpspointdecapture-Accuracy"
Lot2_merge$'SET-OF-Part2'<-NA
Lot2_merge$SubmissionDate<-NA
Lot2_merge$'meta-instanceID'<-NA
Lot2_merge$KEY.x<-NA
Lot2_merge$KEY.y<-NA
Lot2_merge$PARENT_KEY<-NA
Lot2_merge$device_id<-NA


colnames(Lot4_merge)[which(colnames(Lot4_merge)=="metainstanceid")]<-"meta-instanceID"
colnames(Lot4_merge)[which(colnames(Lot4_merge)=="coordgpspointdecaptureLatitude")]<-"coordgpspointdecapture-Latitude"
colnames(Lot4_merge)[which(colnames(Lot4_merge)=="coordgpspointdecaptureLongitude")]<-"coordgpspointdecapture-Longitude"
colnames(Lot4_merge)[which(colnames(Lot4_merge)=="coordgpspointdecaptureAltitude")]<-"coordgpspointdecapture-Altitude"
colnames(Lot4_merge)[which(colnames(Lot4_merge)=="coordgpspointdecaptureAccuracy")]<-"coordgpspointdecapture-Accuracy"
Lot4_merge$'SET-OF-Part2'<-NA
Lot4_merge$SubmissionDate<-NA
Lot4_merge$'meta-instanceID'<-NA
Lot4_merge$KEY.x<-NA
Lot4_merge$KEY.y<-NA
Lot4_merge$PARENT_KEY<-NA
Lot4_merge$device_id<-NA
Lot4_merge$'meta-instanceName'<-NA
Lot4_merge$..25<-NULL
Lot4_merge$metamodelversion<-NULL
Lot4_merge$metauiversion<-NULL
Lot4_merge$metasubmissiondate<-NULL
Lot4_merge$metaiscomplete<-NULL
Lot4_merge$metadatemarkedascomplete<-NULL
Lot4_merge$instanceID<-NULL
Lot4_merge$instanceName<-NULL
Lot4_merge$deviceid<-NULL

## Correction des format de dates. On passe tout au format Y-M-D H:M.S
#lot1merge et lot3merge
Lot1_merge$datedebut <- Lot1_merge$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
Lot1_merge$datedebut <- as.POSIXct(Lot1_merge$datedebut, format="%d/%m/%Y %H:%M:%OS")
Lot1_merge$datefin <- Lot1_merge$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
Lot1_merge$datefin <- as.POSIXct(Lot1_merge$datefin, format="%d/%m/%Y %H:%M:%OS")

Lot3_merge$datedebut <- Lot3_merge$datedebut %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
Lot3_merge$datedebut <- as.POSIXct(Lot3_merge$datedebut, format="%d/%m/%Y %H:%M:%OS")
Lot3_merge$datefin <- Lot3_merge$datefin %>% str_replace_all(c("janv." = "/01/", "nov." = "/11/","sept."= "/09/", "f\\?®vr" = "/02/", "juil." = "/07/" , "ao\\?\\?t" = "/08/","oct." = "/10/", "d\\?®c." = "/12/" )) %>% str_replace_all(c(" /" = "/", "/ " = "/", "/. " = "/"))
Lot3_merge$datefin <- as.POSIXct(Lot3_merge$datefin, format="%d/%m/%Y %H:%M:%OS")

# lot4merge
Lot4_merge$datedebut <- as.POSIXct(Lot4_merge$datedebut, format="%m/%d/%Y %H:%M:%OS")
Lot4_merge$datefin <- as.POSIXct(Lot4_merge$datefin, format="%m/%d/%Y %H:%M:%OS")

Lot2_merge$datedebut <- as.POSIXct(Lot2_merge$datedebut)
Lot2_merge$datefin <- as.POSIXct(Lot2_merge$datefin)


# open mission 1 file
df_mission1<-read.csv(path_to_data_mission1,header=T)
colnames(df_mission1)<-gsub("Part1.","",colnames(df_mission1))
colnames(df_mission1)<-gsub("Part2.","",colnames(df_mission1))
df_mission1$KEY.x<-NA
df_mission1$KEY.y<-NA
df_mission1$PARENT_KEY<-NA
df_mission1$datedebut<-df_mission1$datecapture
df_mission1$datefin<-df_mission1$datecapture
df_mission1$Num<-NULL
df_mission1$'SET-OF-Part2'<-NA
df_mission1$SubmissionDate<-NA
df_mission1$'meta-instanceID'<-NA
df_mission1$'meta-instanceName'<-NA
df_mission1$rcpcoordgpsptdecapture<-NA
df_mission1$utilisationbonsac<-NA
df_mission1$device_id<-NA
df_mission1$rcpidpostedecapture<-NA

df_mission1$idpostedecapture<-paste0(df_mission1$idpointdecapture,df_mission1$postedecapture)
df_mission1$idpostecapture<-NULL

colnames(df_mission1)[which(colnames(df_mission1)=="coordgpspointdecapture.Latitude")]<-"coordgpspointdecapture-Latitude"
colnames(df_mission1)[which(colnames(df_mission1)=="coordgpspointdecapture.Longitude")]<-"coordgpspointdecapture-Longitude"
colnames(df_mission1)[which(colnames(df_mission1)=="coordgpspointdecapture.Altitude")]<-"coordgpspointdecapture-Altitude"
colnames(df_mission1)[which(colnames(df_mission1)=="coordgpspointdecapture.Accuracy")]<-"coordgpspointdecapture-Accuracy"


df_mission1$datedebut <- df_mission1$datedebut %>% str_replace_all(c("Wed Jan 18 00:00:00 UTC 2017" = "2017-01-18", 
                                                                     "Thu Jan 19 00:00:00 UTC 2017" = "2017-01-19", 
                                                                     "Fri Jan 20 00:00:00 UTC 2017" = "2017-01-20", 
                                                                     "Sat Jan 21 00:00:00 UTC 2017" = "2017-01-21", 
                                                                     "Sun Jan 22 00:00:00 UTC 2017" = "2017-01-22", 
                                                                     "Wed Feb 22 00:00:00 UTC 2017" = "2017-01-22", 
                                                                     "Mon Jan 23 00:00:00 UTC 2017" = "2017-01-23", 
                                                                     "Thu Nov 23 00:00:00 UTC 2017" = "2017-01-23", 
                                                                     "Tue Jan 24 00:00:00 UTC 2017" = "2017-01-24", 
                                                                     "Sun Jan 24 00:00:00 UTC 2017" = "2017-01-24", 
                                                                     "Wed Jan 25 00:00:00 UTC 2017" = "2017-01-25", 
                                                                     "Fri Feb 24 00:00:00 UTC 2017" = "2017-01-24", 
                                                                     "Thu Jan 26 00:00:00 UTC 2017" = "2017-01-26", 
                                                                     "Fri Jan 27 00:00:00 UTC 2017" = "2017-01-27", 
                                                                     "Mon Feb 27 00:00:00 UTC 2017" = "2017-01-27", 
                                                                     "Sat Jan 28 00:00:00 UTC 2017" = "2017-01-28", 
                                                                     "Sun Jan 29 00:00:00 UTC 2017" = "2017-01-29"))

df_mission1$datedebut<-as.POSIXct(df_mission1$datedebut)                                                                  
df_mission1$datefin<-df_mission1$datedebut


#bind everything
Lot1_merge$lot<-1
Lot2_merge$lot<-2
Lot3_merge$lot<-3
Lot4_merge$lot<-4
df_mission1$lot<-NA

all_data_bf<-rbind(Lot1_merge,Lot2_merge,Lot3_merge,Lot4_merge,df_mission1)

all_data_bf$lot<-NULL
all_data_bf<-unique(all_data_bf)

all_data_bf$latitude<-as.numeric(all_data_bf$'coordgpspointdecapture-Latitude')
all_data_bf$longitude<-as.numeric(all_data_bf$'coordgpspointdecapture-Longitude')
all_data_bf$altitude<-as.numeric(all_data_bf$'coordgpspointdecapture-Altitude')
all_data_bf$accuracy<-as.numeric(all_data_bf$'coordgpspointdecapture-Accuracy')
all_data_bf$pointdecapture<-as.numeric(all_data_bf$pointdecapture)
all_data_bf$datecapture<-as.Date(all_data_bf$datedebut)

all_data_bf$codevillage[which(all_data_bf$codevillage=="NPI")]<-"NIP"
all_data_bf$idpointdecapture<-gsub("NPI","NIP",all_data_bf$idpointdecapture)
all_data_bf$idpostedecapture<-gsub("NPI","NIP",all_data_bf$idpostedecapture)


## On corrige les numéros de mission
df_dates_mission<-read_excel(path_to_data_data_h_capt, sheet = "Heure de capture")
df_dates_mission<-df_dates_mission[,c(1,3,5)]
colnames(df_dates_mission)<-c("num_mission","codevillage","datecapture")
df_dates_mission$datecapture<-as.Date(df_dates_mission$datecapture)       
df_dates_mission2<-df_dates_mission
df_dates_mission2$datecapture<-df_dates_mission2$datecapture+1
df_dates_mission<-rbind(df_dates_mission,df_dates_mission2)

all_data_bf<-left_join(all_data_bf,df_dates_mission)

#histogram of number of lines / month
hist(all_data_bf$datedebut,"month")

# nombre de lignes par poste de capture (ie nombre de passages par poste de capture)
number_lines_by_postecapture<-all_data_bf %>% group_by(idpostedecapture) %>% summarise(n=n())
hist(number_lines_by_postecapture$n)

## on ne conserve que les colonneset lignes qui nous intéressent
all_data_bf_subset<-all_data_bf %>% dplyr::select(datedebut,datefin,datecapture,codesuperviseur,codevillage,nummission,pointdecapture,idpointdecapture,postedecapture,idpostedecapture,idappdemes,latitude,longitude,altitude,accuracy,idcaptureur,heurecapture,respectpointdecapture,captureurauposte,captureureveille,posicaptureuradequate,habitcaptureuradequat,materielcapturedispo,sacsheureprecedtecollectes,utilisationbonsac,presenceappdemes,bnedispoappdemes,observations)


## On veut vérifier les coordonnées des points de capture
# Pour cela 1) on calcule la mediane des coordonéees de chaque idpostedecapture puis 2) on calcule la médiane de la distance à ces coordonnées médianes

#coords_median_postecapture <- all_data_bf_subset %>% group_by(idpostedecapture) %>% summarize(median_lat=median(latitude,na.rm = T),median_lon=median(longitude,na.rm = T)) #summarise(stdev_dist=sd(sqrt(latitude-longitude)^2)*(111.32 * 1000 * cos(mean(all_data_bf_subset$latitude) * ((pi / 180)))),count=n())  #  
#all_data_bf_subset<-left_join(all_data_bf_subset,coords_median_postecapture,by="idpostedecapture")
#all_data_bf_subset$dist_to_median<-sqrt((all_data_bf_subset$latitude-all_data_bf_subset$median_lat)^2+(all_data_bf_subset$longitude-all_data_bf_subset$median_lon)^2)
#median_dist_to_median_coords<-all_data_bf_subset %>% group_by(idpostedecapture) %>% summarise(median_dist_to_median_coords=median(dist_to_median)*(111.32 * 1000 * cos(mean(all_data_bf_subset$latitude) * ((pi / 180)))))

#on corrige les coords de 4PER1, 4PER2, 4PER3, 4PER4 qui sont mauvaises

#coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "4PER1")]<-coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "1PER1")]
#coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "4PER1")]<-coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "1PER1")]
#coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "4PER2")]<-coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "1PER2")]
#coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "4PER2")]<-coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "1PER2")]
#coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "4PER3")]<-coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "1PER3")]
#coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "4PER3")]<-coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "1PER3")]
#coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "4PER4")]<-coords_postecapture$lat[which(coords_postecapture$idpointdecapture == "1PER4")]
#coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "4PER4")]<-coords_postecapture$lon[which(coords_postecapture$idpointdecapture == "1PER4")]
#coords_postecapture$lat[which(coords_postecapture$idpostedecapture == "6PAL3e")]<-coords_postecapture$lat[which(coords_postecapture$idpostedecapture == "6PAL3i")]
#coords_postecapture$lon[which(coords_postecapture$idpostedecapture == "6PAL3e")]<-coords_postecapture$lon[which(coords_postecapture$idpostedecapture == "6PAL3i")]



############################
# Bind both datasets and load in the DB
############################
all_data_bf_subset$codepays<-"BF"
all_data_civ_subset$codepays<-"CI"

all_data<-rbind(all_data_bf_subset,all_data_civ_subset)

all_data$datedebut<-as.character(all_data$datedebut)
all_data$datefin<-as.character(all_data$datefin)
all_data$datecapture<-as.character(all_data$datecapture)
all_data$nummission<-as.integer(all_data$nummission)
all_data$pointdecapture<-as.integer(all_data$pointdecapture)

all_data <- all_data %>% filter (!(is.na(latitude)))

# on corrige les données erronnés
all_data$codevillage[which(all_data$datecapture %in% c("2017-04-15","2017-04-16") & all_data$codevillage=="NAL")]="NAK"
all_data$codevillage[which(all_data$datecapture %in% c("2017-11-09","2017-11-10","2018-01-17","2018-01-18") & all_data$codevillage=="PER")]="PAL"
all_data$nummission[which(all_data$nummission==1 & all_data$datecapture=="2018-01-17" & all_data$codevillage=="LAG")]=7
all_data$codevillage[which(all_data$codevillage=="NAV" & all_data$codepays=="CI")]="NAA" # a cause de homonimie avec BF
all_data$codevillage[which(all_data$latitude<9 & all_data$codevillage=="NAA" & all_data$codepays=="CI")]="NAM"
all_data$codevillage[which(all_data$latitude<9 & all_data$codevillage=="KOL" & all_data$codepays=="CI")]="BLA"
all_data$codevillage[which(all_data$codevillage=="NOK" & all_data$idpostedecapture=="6NOK2e")]="NAK"
all_data$codevillage[which(all_data$codevillage=="LAG" & all_data$idpostedecapture=="3LAG2i")]="KAT"
all_data$codevillage[which(all_data$codevillage=="NON" & all_data$idpostedecapture=="6NON1e")]="LAG"
all_data$postedecapture[which(all_data$idpostedecapture=="2TAK2NA")]<-"e"
all_data$codevillage[which(all_data$codevillage=="KOU" & all_data$codepays=="CI")]<-"KON" # a cause de homonimie avec BF
all_data$codevillage[which(all_data$codevillage=="NON" & all_data$codepays=="CI" & all_data$latitude>=9.48)]<-"NOW"
all_data$codevillage[which(all_data$codevillage=="KOL" & all_data$codepays=="CI" & all_data$nummission==7)]<-"NOW"
all_data$codevillage[which(all_data$codevillage=="KAM")]<-"PAL"
all_data$codevillage[which(all_data$codevillage=="BAR")]<-"PAL"
all_data<-all_data %>% filter(!(codevillage %in% c("YOU","WUL")))


#num_rows_by_week_year_civ<-all_data %>% filter(codepays=="CI") %>% group_by(week = lubridate::week(datedebut),year=year(datedebut),nummission) %>% summarise(value = n())

all_data$idpointdecapture<-paste0(all_data$nummission,all_data$codevillage,all_data$pointdecapture)
all_data$idpostedecapture<-paste0(all_data$nummission,all_data$codevillage,all_data$pointdecapture,all_data$postedecapture)

