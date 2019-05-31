require(RSQLite)
require(dplyr)
react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)

query<-"SELECT * FROM raw_supervcapture"

## Execute query
df<-dbGetQuery(react_gpkg, query)

# erase data from mission 3 from BF and other bad rows
df <- df %>% filter(!(is.na(nummission))) %>% filter(nummission<=8) %>% filter(!(nummission==3 & codepays=="BF"))


coords_median_postecapture <- df %>% group_by(idpostedecapture,codepays) %>% dplyr::summarize(median_lat=median(Y,na.rm = T),median_lon=median(X,na.rm = T)) 

## On veut vérifier les coordonnées des points de capture
#Pour cela 1) on calcule la mediane des coordonéees de chaque idpostedecapture puis 2) on calcule la médiane de la distance à ces coordonnées médianes
df<-left_join(df,coords_median_postecapture,by=c("idpostedecapture","codepays"))
df$dist_to_median<-sqrt((df$Y-df$median_lat)^2+(df$X-df$median_lon)^2)
median_dist_to_median_coords<-df %>% group_by(idpostedecapture) %>% summarise(median_dist_to_median_coords=median(dist_to_median)*(111.32 * 1000 * cos(mean(df$Y) * ((pi / 180)))))

# grâce à ça on a identifié les points de capture qui posent problème. On les corrige à la main les coords qui sont mauvaises
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "3LAT2i")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "3LAT2e")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "3LAT2i")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "3LAT2e")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "2KAT1e")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1KAT1e")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "2KAT1e")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1KAT1e")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "2YL2e")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "2YL2i")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "2YL2e")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "2YL2i")]

## Pour la mission 3 au BF on a perdu les données. On prend donc pour coordonnées des postes de captures la moyenne des coordonnées de chaque poste de capture des autres missions
coords_median_postecapture$postedecapture<-substr(coords_median_postecapture$idpostedecapture,2,6)
coords_median_postecapture_mission3_bf<- coords_median_postecapture %>% filter(codepays=="BF") %>% group_by(postedecapture) %>% dplyr::summarize(median_lat=mean(median_lat),median_lon=mean(median_lon))
coords_median_postecapture_mission3_bf$idpostedecapture<-paste0("3",coords_median_postecapture_mission3_bf$postedecapture)
coords_median_postecapture_mission3_bf$postedecapture<-NULL
coords_median_postecapture_mission3_bf$codepays<-"BF"
coords_median_postecapture$postedecapture<-NULL
coords_median_postecapture<-rbind(as.data.frame(coords_median_postecapture),as.data.frame(coords_median_postecapture_mission3_bf))

## On vérifie et éventuellement on corrige les dates pour le BF
# On récupère pour chaque idpostedecapture le minimum de la colonne 'datedecapture'
dates_captures_from_supervcapture_bf<-df %>% filter(codepays=="BF") %>% group_by(nummission,codevillage,datecapture) %>% summarise(n=n()) 
dates_captures_from_supervcapture_bf<-dates_captures_from_supervcapture_bf %>% filter(n>1) %>% group_by(nummission,codevillage) %>% summarise(date_from_supervcapture=min(datecapture)) 
#dates_captures_from_supervcapture_bf$idpostedecapture<-substr(dates_captures_from_supervcapture_bf$idpostedecapture,1,4)
dates_captures_from_supervcapture_bf$date_from_supervcapture<-as.Date(dates_captures_from_supervcapture_bf$date_from_supervcapture)

# on compare pour le BF avec le fichier de dieudonné
query<-"SELECT * FROM raw_bf_dates_hlc"
df_dates_captures_par_villages_bf<-dbGetQuery(react_gpkg, query)
#df_dates_captures_par_villages_bf$idpostedecapture<-paste0(df_dates_captures_par_villages_bf$n_mission,df_dates_captures_par_villages_bf$code_village)
df_dates_captures_par_villages_bf$date_de_captures<-as.Date(df_dates_captures_par_villages_bf$date_de_captures)
df_dates_captures_par_villages_bf<-df_dates_captures_par_villages_bf[,c("n_mission","code_village","date_de_captures")]
colnames(df_dates_captures_par_villages_bf)<-c("nummission","codevillage","date_from_dieudo")
dates_captures_from_supervcapture_bf<-right_join(dates_captures_from_supervcapture_bf,df_dates_captures_par_villages_bf)
dates_captures_from_supervcapture_bf$corresp<-dates_captures_from_supervcapture_bf$date_from_supervcapture-dates_captures_from_supervcapture_bf$date_from_dieudo

# Pour la grande majorité des données on a les mêmes dates entre le fichier de Dieudo et le fichier de supervcapture. Pour les 2 dates qui ne correspondent pas on a vérifié que les bonnes dates sont celles du fichier de Dieudo
coords_median_postecapture_bf<-coords_median_postecapture %>% filter(codepays=="BF")
coords_median_postecapture_bf$idpostedecapture2<-substr(coords_median_postecapture_bf$idpostedecapture,1,4)
dates_captures_from_supervcapture_bf$idpostedecapture2<-paste0(dates_captures_from_supervcapture_bf$nummission,dates_captures_from_supervcapture_bf$codevillage)
coords_median_postecapture_bf<-left_join(coords_median_postecapture_bf,dates_captures_from_supervcapture_bf)
coords_median_postecapture_bf$idpostedecapture2<-coords_median_postecapture_bf$date_from_supervcapture<-coords_median_postecapture_bf$corresp<-NULL

# On ajoute les lignes qui manquent (5NAV1i, 5NAV2e, 5NAV2i, 5NAV3e, 5NAV3i, 5NAV4e, 5NAV4i)
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV1e","BF",10.89252,-3.329525,5,"NAV","2017-12-14"))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV1i","BF",10.89252,-3.329534,5,"NAV","2017-12-14"))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV2e","BF",10.89294,-3.329959,5,"NAV","2017-12-14"))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV2i","BF",10.89292,-3.329994,5,"NAV","2017-12-14"))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV3e","BF",10.89194,-3.329896,5,"NAV","2017-12-14"))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV3i","BF",10.89199,-3.329933,5,"NAV","2017-12-14"))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV4e","BF",10.8925,-3.328943,5,"NAV","2017-12-14"))
coords_median_postecapture_bf<-rbind(coords_median_postecapture_bf,c("5NAV4i","BF",10.89249,-3.328945,5,"NAV","2017-12-14"))


## On vérifie et éventuellement on corrige les dates pour la CIV
dates_captures_from_supervcapture_civ<-df %>% filter(codepays=="CI") %>% group_by(nummission,codevillage,pointdecapture,datecapture) %>% summarise(n=n()) %>% arrange(nummission,codevillage,pointdecapture,datecapture)
#dates_captures_from_supervcapture_civ$idpointdecapture<-paste0(dates_captures_from_supervcapture_civ$nummission,dates_captures_from_supervcapture_civ$codevillage,dates_captures_from_supervcapture_civ$pointdecapture)
#dates_captures_from_supervcapture_civ2<-dates_captures_from_supervcapture_civ %>% group_by(idpointdecapture) %>% summarise(n_dates=n())
#dates_captures_from_supervcapture_civ<-left_join(dates_captures_from_supervcapture_civ,dates_captures_from_supervcapture_civ2)

# n donne le nombre de passages pour un point de capture donné (exterieur et intérieur compris) à une date donnée sur toute la nuit. On conserve à partir de 4 passages. Pour les données sans date, on verra au cas par cas
dates_captures_from_supervcapture_civ<-dates_captures_from_supervcapture_civ %>% filter(n>=5) %>% group_by(nummission,codevillage,pointdecapture) %>% summarise(date_hlc=min(datecapture)) 
dates_captures_from_supervcapture_civ$idpostedecapture2<-paste0(dates_captures_from_supervcapture_civ$nummission,dates_captures_from_supervcapture_civ$codevillage,dates_captures_from_supervcapture_civ$pointdecapture)

coords_median_postecapture_civ<-coords_median_postecapture %>% filter(codepays=="CI")
coords_median_postecapture_civ$idpostedecapture2<-substr(coords_median_postecapture_civ$idpostedecapture,1,5)
coords_median_postecapture_civ<-left_join(coords_median_postecapture_civ,dates_captures_from_supervcapture_civ,by="idpostedecapture2")
coords_median_postecapture_civ$idpostedecapture2<-NULL
coords_median_postecapture_civ$nummission<-substr(coords_median_postecapture_civ$idpostedecapture,1,1)
coords_median_postecapture_civ$codevillage<-substr(coords_median_postecapture_civ$idpostedecapture,2,4)
coords_median_postecapture_civ$pointdecapture<-substr(coords_median_postecapture_civ$idpostedecapture,5,5)

# On corrige à la main les dates manquantes
#coords_median_postecapture_civ$date_hlc[which(coords_median_postecapture_civ$codevillage=="GUE" & coords_median_postecapture_civ$nummission==1)]<-"2016-09-21"
#coords_median_postecapture_civ$date_hlc[which(coords_median_postecapture_civ$codevillage=="KOL" & coords_median_postecapture_civ$nummission==1)]<-"2016-09-23"
#coords_median_postecapture_civ$date_hlc[which(coords_median_postecapture_civ$codevillage=="KOU" & coords_median_postecapture_civ$nummission==3)]<-"2017-05-25"



##### Verification avec les données de capturedeterm
query<-"SELECT raw_capturedeterm.*,raw_villages.codepays_fk FROM raw_capturedeterm JOIN raw_villages ON raw_villages.codevillage_pk=raw_capturedeterm.codevillage_fk"
df_capturedeterm<-dbGetQuery(react_gpkg, query)
df_capturedeterm$date<-as.Date(df_capturedeterm$date)-1
## Pour le BF
df_capturedeterm_bf<-df_capturedeterm %>% filter(codepays_fk=="BF")
df_capturedeterm_bf<-full_join(df_capturedeterm_bf,coords_median_postecapture_bf,by="idpostedecapture")
df_capturedeterm_bf$diff_date<-df_capturedeterm_bf$date-df_capturedeterm_bf$date_from_dieudo
## Pour le BF on a l'air OK !!
## Pour la CIV
df_capturedeterm_civ<-df_capturedeterm %>% filter(codepays_fk=="CI")
df_capturedeterm_civ<-left_join(df_capturedeterm_civ,coords_median_postecapture_civ,by="idpostedecapture")
# idpostedecapture pour lesquels on a des données d'identification de moustique mais pas de données de supervision capture
idpostedecapture_missing<-df_capturedeterm_civ %>% filter(is.na(nummission))
idpostedecapture_missing<-sort(unique(idpostedecapture_missing$idpostedecapture))