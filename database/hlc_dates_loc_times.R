#require(RSQLite)
#require(dplyr)
#require(lubridate)
#react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)

query<-"SELECT * FROM raw_supervcapture"

## Execute query
df<-dbGetQuery(react_gpkg, query)

# erase data from mission 3 from BF and other bad rows
df <- df %>% filter(!(is.na(nummission))) %>% filter(nummission<=8) %>% filter(!(nummission==3 & codepays=="BF"))


coords_median_postecapture <- df %>% group_by(idpointdecapture,codepays) %>% dplyr::summarize(median_lat=median(Y,na.rm = T),median_lon=median(X,na.rm = T)) 

## On veut vérifier les coordonnées des points de capture
#Pour cela 1) on calcule la mediane des coordonéees de chaque idpointdecapture puis 2) on calcule la médiane de la distance à ces coordonnées médianes
df<-left_join(df,coords_median_postecapture,by=c("idpointdecapture","codepays"))
df$dist_to_median<-sqrt((df$Y-df$median_lat)^2+(df$X-df$median_lon)^2)
median_dist_to_median_coords<-df %>% group_by(idpointdecapture) %>% summarise(median_dist_to_median_coords=median(dist_to_median)*(111.32 * 1000 * cos(mean(df$Y) * ((pi / 180)))))

# grâce à ça on a identifié les points de capture qui posent problème. On les corrige à la main les coords qui sont mauvaises
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpointdecapture == "2KAT1")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpointdecapture == "1KAT1")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpointdecapture == "2KAT1")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpointdecapture == "1KAT1")]

# D'abord on identifie les points de capture manquants dans les données de tablettes
coords_median_postecapture<-as.data.frame(lapply(coords_median_postecapture,rep,2))
coords_median_postecapture$postedecapture<-rep(c("i","e"),nrow(coords_median_postecapture)/2)
coords_median_postecapture$idpostedecapture<-paste0(coords_median_postecapture$idpointdecapture,coords_median_postecapture$postedecapture)
coords_median_postecapture<- coords_median_postecapture %>% arrange(idpostedecapture)
coords_median_postecapture$quality_flag_position<-1

## Pour la mission 3 au BF on a perdu les données. On prend donc pour coordonnées des postes de captures la moyenne des coordonnées de chaque poste de capture des autres missions
coords_median_postecapture$postedecapture<-substr(coords_median_postecapture$idpostedecapture,2,6)
coords_median_postecapture_mission3_bf<- coords_median_postecapture %>% filter(codepays=="BF") %>% group_by(postedecapture) %>% dplyr::summarize(median_lat=mean(median_lat),median_lon=mean(median_lon))
coords_median_postecapture_mission3_bf$idpostedecapture<-paste0("3",coords_median_postecapture_mission3_bf$postedecapture)
coords_median_postecapture_mission3_bf$postedecapture<-NULL
coords_median_postecapture_mission3_bf$codepays<-"BF"
coords_median_postecapture_mission3_bf$quality_flag_position<-2
coords_median_postecapture$postedecapture<-coords_median_postecapture$idpointdecapture<-NULL
coords_median_postecapture<-rbind(as.data.frame(coords_median_postecapture),as.data.frame(coords_median_postecapture_mission3_bf))


## On vérifie et éventuellement on corrige les dates pour le BF
# On récupère pour chaque idpostedecapture le minimum de la colonne 'datedecapture'
dates_captures_from_supervcapture_bf<-df %>% filter(codepays=="BF") %>% group_by(nummission,codevillage,datecapture) %>% summarise(n=n()) 
dates_captures_from_supervcapture_bf<-dates_captures_from_supervcapture_bf %>% filter(n>1) %>% group_by(nummission,codevillage) %>% summarise(date_from_supervcapture=min(datecapture)) 
#dates_captures_from_supervcapture_bf$idpostedecapture<-substr(dates_captures_from_supervcapture_bf$idpostedecapture,1,4)
dates_captures_from_supervcapture_bf$date_from_supervcapture<-as.Date(dates_captures_from_supervcapture_bf$date_from_supervcapture)

# on compare pour le BF avec le fichier de dieudonné
query<-"SELECT * FROM raw_dates_hlc"
df_dates_captures_par_villages<-dbGetQuery(react_gpkg, query)
df_dates_captures_par_villages_bf<-df_dates_captures_par_villages %>% filter(codepays=="BF")
#df_dates_captures_par_villages_bf$idpostedecapture<-paste0(df_dates_captures_par_villages_bf$n_mission,df_dates_captures_par_villages_bf$code_village)
df_dates_captures_par_villages_bf$date_de_captures<-as.Date(df_dates_captures_par_villages_bf$date_de_captures)
df_dates_captures_par_villages_bf<-df_dates_captures_par_villages_bf[,c("n_mission","codevillage","date_de_captures","date_heure_debut","date_heure_fin")]
colnames(df_dates_captures_par_villages_bf)<-c("nummission","codevillage","date_from_dieudo","date_heure_debut","date_heure_fin")
dates_captures_from_supervcapture_bf<-right_join(dates_captures_from_supervcapture_bf,df_dates_captures_par_villages_bf)
dates_captures_from_supervcapture_bf$corresp<-dates_captures_from_supervcapture_bf$date_from_supervcapture-dates_captures_from_supervcapture_bf$date_from_dieudo

# Pour la grande majorité des données on a les mêmes dates entre le fichier de Dieudo et le fichier de supervcapture. Pour les 2 dates qui ne correspondent pas on a vérifié que les bonnes dates sont celles du fichier de Dieudo
coords_median_postecapture_bf<-coords_median_postecapture %>% filter(codepays=="BF")
coords_median_postecapture_bf$idpostedecapture2<-substr(coords_median_postecapture_bf$idpostedecapture,1,4)
dates_captures_from_supervcapture_bf$idpostedecapture2<-paste0(dates_captures_from_supervcapture_bf$nummission,dates_captures_from_supervcapture_bf$codevillage)
coords_median_postecapture_bf<-left_join(coords_median_postecapture_bf,dates_captures_from_supervcapture_bf)
coords_median_postecapture_bf$idpostedecapture2<-coords_median_postecapture_bf$date_from_supervcapture<-coords_median_postecapture_bf$corresp<-NULL

# On ajoute les lignes qui manquent (5NAV1i, 5NAV2e, 5NAV2i, 5NAV3e, 5NAV3i, 5NAV4e, 5NAV4i)
coords_median_postecapture_bf <-coords_median_postecapture_bf %>% dplyr::select(idpostedecapture,codepays,median_lat,median_lon,nummission,codevillage,date_from_dieudo,date_heure_debut,date_heure_fin,quality_flag_position)
coords_median_postecapture_bf <- coords_median_postecapture_bf %>% filter (!(codevillage=="NAV" & nummission==5))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV1e","BF",10.89252,-3.329525,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV1i","BF",10.89252,-3.329525,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV2e","BF",10.89294,-3.329959,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV2i","BF",10.89294,-3.329959,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV3e","BF",10.89194,-3.329896,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV3i","BF",10.89194,-3.329896,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV4e","BF",10.8925,-3.328943,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV4i","BF",10.8925,-3.328943,5,"NAV","2017-12-14","2017-12-14 17:00:00","2017-12-15 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB1e","BF",10.74350,-3.399661,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB1i","BF",10.74350,-3.399661,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB2e","BF",10.74495,-3.399188,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB2i","BF",10.74495,-3.399188,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB3e","BF",10.74457,-3.398204,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB3i","BF",10.74457,-3.398204,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB4e","BF",10.74427,-3.397680,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("1DMB4i","BF",10.74427,-3.397680,1,"DMB","2017-01-20","2017-01-20 18:00:00","2017-01-21 09:00:00",2))

coords_median_postecapture_bf$quality_flag_horaires<-1
coords_median_postecapture_bf$poindecapture<-substr(coords_median_postecapture_bf$idpostedecapture,5,5)
coords_median_postecapture_bf$idpointdecapture<-substr(coords_median_postecapture_bf$idpostedecapture,1,5)
coords_median_postecapture_bf$idpostedecapture<-NULL
coords_median_postecapture_bf<-unique(coords_median_postecapture_bf)

colnames(coords_median_postecapture_bf)<-c("codepays","latitude","longitude","nummission","codevillage","date_capture","date_heure_debut","date_heure_fin","quality_flag_position","quality_flag_horaires","pointdecapture","idpointdecapture")
coords_median_postecapture_bf<-coords_median_postecapture_bf %>% mutate_all(as.character)

## On vérifie et éventuellement on corrige les dates pour la CIV
dates_captures_from_supervcapture_civ<-df %>% filter(codepays=="CI") %>% group_by(nummission,codevillage,pointdecapture,datecapture) %>% summarise(n=n()) %>% arrange(nummission,codevillage,pointdecapture,datecapture)
dates_captures_from_supervcapture_civ$datecapture<-as.Date(dates_captures_from_supervcapture_civ$datecapture)
#dates_captures_from_supervcapture_civ$idpointdecapture<-paste0(dates_captures_from_supervcapture_civ$nummission,dates_captures_from_supervcapture_civ$codevillage,dates_captures_from_supervcapture_civ$pointdecapture)
#dates_captures_from_supervcapture_civ2<-dates_captures_from_supervcapture_civ %>% group_by(idpointdecapture) %>% summarise(n_dates=n())
#dates_captures_from_supervcapture_civ<-left_join(dates_captures_from_supervcapture_civ,dates_captures_from_supervcapture_civ2)

# n donne le nombre de passages pour un point de capture donné (exterieur et intérieur compris) à une date donnée. On conserve à partir de 4 passages. Pour les données sans date, on verra au cas par cas
dates_captures_from_supervcapture_civ<-dates_captures_from_supervcapture_civ %>% filter(n>=5) %>% group_by(nummission,codevillage,pointdecapture) %>% summarise(date_from_supervcapture=min(datecapture)) 
dates_captures_from_supervcapture_civ$idpostedecapture2<-paste0(dates_captures_from_supervcapture_civ$nummission,dates_captures_from_supervcapture_civ$codevillage,dates_captures_from_supervcapture_civ$pointdecapture)

coords_median_postecapture_civ<-coords_median_postecapture %>% filter(codepays=="CI")
coords_median_postecapture_civ$idpostedecapture2<-substr(coords_median_postecapture_civ$idpostedecapture,1,5)
coords_median_postecapture_civ<-left_join(coords_median_postecapture_civ,dates_captures_from_supervcapture_civ,by="idpostedecapture2")
coords_median_postecapture_civ$idpostedecapture2<-NULL
coords_median_postecapture_civ$nummission<-as.numeric(substr(coords_median_postecapture_civ$idpostedecapture,1,1))
coords_median_postecapture_civ$codevillage<-substr(coords_median_postecapture_civ$idpostedecapture,2,4)
coords_median_postecapture_civ$pointdecapture<-as.numeric(substr(coords_median_postecapture_civ$idpostedecapture,5,5))

## Pour ces données, on croise avec les données de Barnabas pour vérifier les dates
df_dates_captures_par_villages_civ<-df_dates_captures_par_villages %>% filter(codepays=="CI")
df_dates_captures_par_villages_civ$date_de_captures<-as.Date(df_dates_captures_par_villages_civ$date_de_captures)
df_dates_captures_par_villages_civ<-df_dates_captures_par_villages_civ[,c("n_mission","codevillage","date_de_captures")]
colnames(df_dates_captures_par_villages_civ)<-c("nummission","codevillage","date_from_barnabas")
dates_captures_from_supervcapture_civ<-right_join(dates_captures_from_supervcapture_civ,df_dates_captures_par_villages_civ)
dates_captures_from_supervcapture_civ$corresp<-as.Date(dates_captures_from_supervcapture_civ$date_from_supervcapture)-dates_captures_from_supervcapture_civ$date_from_barnabas
# Pour la grande majorité des données on a les mêmes dates entre le fichier de Barbanabas et le fichier de supervcapture. Pour les qq dates qui ne correspondent pas on a vérifié que les bonnes dates sont celles du fichier de Barbabas
dates_captures_from_supervcapture_civ$date_from_supervcapture<-dates_captures_from_supervcapture_civ$idpostedecapture2<-dates_captures_from_supervcapture_civ$corresp<-NULL
coords_median_postecapture_civ$date_from_supervcapture<-NULL
coords_median_postecapture_civ<-right_join(coords_median_postecapture_civ,dates_captures_from_supervcapture_civ)
# En plus de tous les NA, il manque les points de captures n°1 et n°2 de LAT à la mission 3   et   les points de capture n°3 et n°4 de FEL à la mission 1. On les complete à la main
# Pour les coordonnées on prend la moyenne des coordonnées des points de capture lors des autres missions
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.210729,-5.756261,"3LAT1e",2,3,"LAT",1,"2017-03-04"))
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.210729,-5.756261,"3LAT1i",2,3,"LAT",1,"2017-03-04"))
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.210018,-5.756201,"3LAT2e",2,3,"LAT",2,"2017-03-04"))
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.210018,-5.756201,"3LAT2i",2,3,"LAT",2,"2017-03-04"))
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.448971,-5.470124,"1FEL3e",2,1,"FEL",3,"2016-10-08"))
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.448971,-5.470124,"1FEL3i",2,1,"FEL",3,"2016-10-08"))
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.449943,-5.470300,"1FEL4e",2,1,"FEL",4,"2016-10-08"))
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,c("CI",9.449943,-5.470300,"1FEL4i",2,1,"FEL",4,"2016-10-08"))

## Pour les points de capture dont on n'a aucune info dans les données de tablette : 
# - pour les dates, on utilise le tableau de Barnabas
# - pour les coordonnées geographiques, on prend la moyenne des coordonnées issues des autres missions
coords_median_postecapture_civ_miss_coord<-coords_median_postecapture_civ %>% filter(is.na(idpostedecapture))
coords_median_postecapture_civ<-coords_median_postecapture_civ %>% filter(!(is.na(idpostedecapture)))
coords_median_postecapture_civ_miss_coord<-as.data.frame(lapply(coords_median_postecapture_civ_miss_coord,rep,4))
coords_median_postecapture_civ_miss_coord$pointdecapture<-rep(c(1,2,3,4),nrow(coords_median_postecapture_civ_miss_coord)/4)
coords_median_postecapture_civ_miss_coord<-as.data.frame(lapply(coords_median_postecapture_civ_miss_coord,rep,2)) %>% arrange(nummission,codevillage,pointdecapture)
coords_median_postecapture_civ_miss_coord$postedecapture<-rep(c("i","e"),nrow(coords_median_postecapture_civ_miss_coord)/2)
coords_median_postecapture_civ_miss_coord$idpostedecapture<-paste0(coords_median_postecapture_civ_miss_coord$nummission,coords_median_postecapture_civ_miss_coord$codevillage,coords_median_postecapture_civ_miss_coord$pointdecapture,coords_median_postecapture_civ_miss_coord$postedecapture)
coords_median_postecapture_civ_miss_coord$codepays<-"CI"
coords_median_postecapture_civ_miss_coord$quality_flag_position<-2

coords_mean<-coords_median_postecapture_civ %>% group_by(codevillage,pointdecapture) %>% summarise(mean_lat=mean(as.numeric(median_lat)),mean_lon=mean(as.numeric(median_lon)))
coords_mean$pointdecapture<-as.numeric(coords_mean$pointdecapture)
coords_median_postecapture_civ_miss_coord<-left_join(coords_median_postecapture_civ_miss_coord,coords_mean)
coords_median_postecapture_civ_miss_coord$median_lat<-coords_median_postecapture_civ_miss_coord$mean_lat
coords_median_postecapture_civ_miss_coord$median_lon<-coords_median_postecapture_civ_miss_coord$mean_lon
coords_median_postecapture_civ_miss_coord$mean_lat<-coords_median_postecapture_civ_miss_coord$mean_lon<-NULL

coords_median_postecapture_civ$postedecapture<-substr(coords_median_postecapture_civ$idpostedecapture,6,6)
coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,coords_median_postecapture_civ_miss_coord)


# On ajoute les horaires de début de fin
horaires_civ<-df %>% filter(codepays=="CI")
horaires_civ$datedebut=as.POSIXct(horaires_civ$datedebut)
horaires_civ$datefin=as.POSIXct(horaires_civ$datefin)
horaires_civ <- horaires_civ %>% group_by(nummission,codevillage,pointdecapture) %>% summarise(date_min=min(datedebut),date_max=max(datefin),count=n())
horaires_civ$date_min=as.POSIXlt(horaires_civ$date_min)
horaires_civ$date_max=as.POSIXlt(horaires_civ$date_max)
horaires_civ$diff_time<-difftime(horaires_civ$date_max,horaires_civ$date_min,units="hours")
horaires_civ$heure_debut <- hour(horaires_civ$date_min)
horaires_civ$heure_fin <- hour(horaires_civ$date_max)
horaires_civ$minute_debut <- minute(horaires_civ$date_min)
horaires_civ$minute_fin <- minute(horaires_civ$date_max)
horaires_civ$date_min=as.character(horaires_civ$date_min)
horaires_civ$date_max=as.character(horaires_civ$date_max)
# on supprime toutes les données qui ne sont pas plausibles
horaires_civ_filt<-horaires_civ %>% filter(between(count,15,30)) %>% filter(between(diff_time,9,18))
horaires_civ_filt$count<-horaires_civ_filt$diff_time<-horaires_civ_filt$heure_debut<-horaires_civ_filt$heure_fin<-horaires_civ_filt$minute_debut<-horaires_civ_filt$minute_fin<-NULL
horaires_civ_filt$nummission<-as.character(horaires_civ_filt$nummission)
horaires_civ_filt$pointdecapture<-as.character(horaires_civ_filt$pointdecapture)
# on joint avec les données de coordonnées geographiques
coords_median_postecapture_civ<-left_join(coords_median_postecapture_civ,horaires_civ_filt)
coords_median_postecapture_civ$date_min[which(coords_median_postecapture_civ$codevillage=="KOG" & coords_median_postecapture_civ$nummission==3)] <- gsub("2017-03-06","2017-03-04",coords_median_postecapture_civ$date_min[which(coords_median_postecapture_civ$codevillage=="KOG" & coords_median_postecapture_civ$nummission==3)])
coords_median_postecapture_civ$date_min[which(coords_median_postecapture_civ$codevillage=="TAH" & coords_median_postecapture_civ$nummission==3)] <- gsub("2017-03-04","2017-03-02",coords_median_postecapture_civ$date_min[which(coords_median_postecapture_civ$codevillage=="KOG" & coords_median_postecapture_civ$nummission==3)])
# pour toutes les données on l'on n'a pas les horaires de captures, on prend la moyenne des horaires pour le village et l'enquete
coords_postecapture_civ_no_horaire<-coords_median_postecapture_civ %>% filter (is.na(date_min))
coords_median_postecapture_civ<-coords_median_postecapture_civ %>% filter (!(is.na(date_min)))
coords_median_postecapture_civ$quality_flag_horaires<-1
coords_median_postecapture_civ$date_min=as.POSIXlt(coords_median_postecapture_civ$date_min)
coords_median_postecapture_civ$date_max=as.POSIXlt(coords_median_postecapture_civ$date_max)
coords_median_postecapture_civ$heure_debut <- hour(coords_median_postecapture_civ$date_min) + minute(coords_median_postecapture_civ$date_min)/60
coords_median_postecapture_civ$heure_fin <- hour(coords_median_postecapture_civ$date_max) + minute(coords_median_postecapture_civ$date_max)/60
mean_horaires_by_vill_miss<-coords_median_postecapture_civ %>% group_by(codevillage,nummission) %>% summarise(mean_heure_min_vill_mission=mean(heure_debut),mean_heure_max_vill_mission=mean(heure_fin))
# To see the standard deviation of horaires (get an idea on if the starting and ending times are more or less equal for a given village and enquete) : sd_horaires_by_vill_miss<-coords_median_postecapture_civ %>% group_by(codevillage,nummission) %>% summarise(sd_heure_min=sd(heure_debut),sd_heure_max=sd(heure_fin))
coords_postecapture_civ_no_horaire<-left_join(coords_postecapture_civ_no_horaire,mean_horaires_by_vill_miss)
# Il nous manque encore des horaires. pour ces horaires manquants on prend les horaires moyens toutes missions confondues pour ces villages (en vérifiant avec qu'ils ne différent pas tropavec l'écart type) avec: sd_horaires_by_vill<-coords_median_postecapture_civ %>% group_by(codevillage) %>% summarise(sd_heure_min=sd(heure_debut),sd_heure_max=sd(heure_fin))
mean_horaires_by_vill<-coords_median_postecapture_civ %>% group_by(codevillage) %>% summarise(mean_heure_min_vill=mean(heure_debut),mean_heure_max_vill=mean(heure_fin))
coords_postecapture_civ_no_horaire<-left_join(coords_postecapture_civ_no_horaire,mean_horaires_by_vill)
# Il nous manque encore des horaires. pour ces horaires manquants on prend les horaires moyens tous villages et toutes missions confondues
coords_postecapture_civ_no_horaire$mean_heure_min<-mean(coords_median_postecapture_civ$heure_debut)
coords_postecapture_civ_no_horaire$mean_heure_max<-mean(coords_median_postecapture_civ$heure_fin)
# Horaires quality flag
coords_postecapture_civ_no_horaire$quality_flag_horaires<-2
coords_postecapture_civ_no_horaire$quality_flag_horaires[which(is.na(coords_postecapture_civ_no_horaire$mean_heure_min_vill_mission))]<-3
coords_postecapture_civ_no_horaire$quality_flag_horaires[which(is.na(coords_postecapture_civ_no_horaire$mean_heure_min_vill))]<-4

coords_postecapture_civ_no_horaire$heure_debut<-coords_postecapture_civ_no_horaire$mean_heure_min_vill_mission
coords_postecapture_civ_no_horaire$heure_fin<-coords_postecapture_civ_no_horaire$mean_heure_max_vill_mission
coords_postecapture_civ_no_horaire$heure_debut[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==3)]<-coords_postecapture_civ_no_horaire$mean_heure_min_vill[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==3)]
coords_postecapture_civ_no_horaire$heure_fin[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==3)]<-coords_postecapture_civ_no_horaire$mean_heure_max_vill[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==3)]
coords_postecapture_civ_no_horaire$heure_debut[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==4)]<-coords_postecapture_civ_no_horaire$mean_heure_min[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==4)]
coords_postecapture_civ_no_horaire$heure_fin[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==4)]<-coords_postecapture_civ_no_horaire$mean_heure_max[which(coords_postecapture_civ_no_horaire$quality_flag_horaires==4)]


coords_postecapture_civ_no_horaire$heure_debut<-paste0(floor(coords_postecapture_civ_no_horaire$heure_debut),":",sprintf("%02d",round(coords_postecapture_civ_no_horaire$heure_debut%%1*60)),":00")
coords_postecapture_civ_no_horaire$heure_fin<-paste0(floor(coords_postecapture_civ_no_horaire$heure_fin),":",sprintf("%02d",round(coords_postecapture_civ_no_horaire$heure_fin%%1*60)),":00")
coords_postecapture_civ_no_horaire$date_min<-paste0(coords_postecapture_civ_no_horaire$date_from_barnabas," ",coords_postecapture_civ_no_horaire$heure_debut)
coords_postecapture_civ_no_horaire$date_max<-paste0(coords_postecapture_civ_no_horaire$date_from_barnabas+1," ",coords_postecapture_civ_no_horaire$heure_fin)

coords_median_postecapture_civ$heure_debut<-coords_median_postecapture_civ$heure_fin<-NULL
coords_postecapture_civ_no_horaire$mean_heure_min_vill_mission<-coords_postecapture_civ_no_horaire$mean_heure_max_vill_mission<-coords_postecapture_civ_no_horaire$mean_heure_min_vill<-coords_postecapture_civ_no_horaire$mean_heure_max_vill<-coords_postecapture_civ_no_horaire$mean_heure_min<-coords_postecapture_civ_no_horaire$mean_heure_max<-coords_postecapture_civ_no_horaire$heure_debut<-coords_postecapture_civ_no_horaire$heure_fin<-NULL

coords_median_postecapture_civ<-rbind(coords_median_postecapture_civ,coords_postecapture_civ_no_horaire)

coords_median_postecapture_civ$date_min<-as.character(coords_median_postecapture_civ$date_min)
coords_median_postecapture_civ$date_max<-as.character(coords_median_postecapture_civ$date_max)


coords_median_postecapture_civ$idpointdecapture<-substr(coords_median_postecapture_civ$idpostedecapture,1,5)
coords_median_postecapture_civ$idpostedecapture<-coords_median_postecapture_civ$postedecapture<-NULL
coords_median_postecapture_civ<-unique(coords_median_postecapture_civ)

colnames(coords_median_postecapture_civ)<-c("codepays","latitude","longitude","quality_flag_position","nummission","codevillage","pointdecapture","date_capture","date_heure_debut","date_heure_fin","quality_flag_horaires","idpointdecapture")
coords_median_postecapture_civ<-coords_median_postecapture_civ %>% mutate_all(as.character)

hlc_dates_loc_times<-rbind(coords_median_postecapture_bf,coords_median_postecapture_civ)

hlc_dates_loc_times$quality_flag_dates<-"1"
hlc_dates_loc_times <-hlc_dates_loc_times %>% dplyr::select(idpointdecapture,nummission,codevillage,pointdecapture,codepays,longitude,latitude,date_capture,date_heure_debut,date_heure_fin,quality_flag_position,quality_flag_dates,quality_flag_horaires) %>% arrange(codepays,idpointdecapture)

#### Vérification via qq graphiques
#hlc_dates_loc_times$nummission<-as.numeric(hlc_dates_loc_times$nummission)
#hlc_dates_loc_times$pointdecapture<-as.numeric(hlc_dates_loc_times$pointdecapture)
#hlc_dates_loc_times$date_capture<-as.Date(hlc_dates_loc_times$date_capture)
#hlc_dates_loc_times$date_heure_debut<-as.POSIXct(hlc_dates_loc_times$date_heure_debut)
#hlc_dates_loc_times$date_heure_fin<-as.POSIXct(hlc_dates_loc_times$date_heure_fin)

#dat<-hlc_dates_loc_times %>% filter(codepays=="BF")

#hist(dat$date_heure_debut,"week") # répartition dans le temps des jours de capture
#hist(as.numeric(dat$date_heure_fin-dat$date_heure_debut)) # durée des captures
#nb_pts_captures_villages_mission<-dat %>% group_by(codepays,codevillage,nummission) %>% summarise(count=n())



##### Verification avec les données de capturedeterm
#query<-"SELECT raw_capturedeterm.*,raw_villages.codepays FROM raw_capturedeterm JOIN raw_villages ON raw_villages.codevillage=raw_capturedeterm.codevillage"
#df_capturedeterm<-dbGetQuery(react_gpkg, query)
#df_capturedeterm$date<-as.Date(df_capturedeterm$date)-1
## Pour le BF
#df_capturedeterm_bf<-df_capturedeterm %>% filter(codepays=="BF")
#df_capturedeterm_bf<-full_join(df_capturedeterm_bf,coords_median_postecapture_bf,by="idpostedecapture")
#df_capturedeterm_bf$diff_date<-df_capturedeterm_bf$date-df_capturedeterm_bf$date_from_dieudo
## Pour le BF on a l'air OK !!
## Pour la CIV
#df_capturedeterm_civ<-df_capturedeterm %>% filter(codepays=="CI")
#df_capturedeterm_civ<-full_join(df_capturedeterm_civ,coords_median_postecapture_civ,by="idpostedecapture")
#df_capturedeterm_civ$diff_date<-df_capturedeterm_civ$date-df_capturedeterm_civ$date_from_barnabas
#idpostedecapture_missing<-df_capturedeterm_civ %>% filter(is.na(nummission))
#idpostedecapture_missing<-sort(unique(idpostedecapture_missing$idpostedecapture))