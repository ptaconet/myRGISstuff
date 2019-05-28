require(RSQLite)
require(dplyr)
#react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)

query<-"SELECT * FROM raw_supervcapture"

## Execute query
df<-dbGetQuery(react_gpkg, query)

# erase data from mission 3 from BF and other bad rows
df <- df %>% filter(!(is.na(nummission))) %>% filter(nummission<=8) %>% filter(!(nummission==3 & codepays=="BF"))

coords_median_postecapture <- df %>% group_by(idpostedecapture,codepays) %>% summarize(median_lat=median(Y,na.rm = T),median_lon=median(X,na.rm = T)) 

## On veut vérifier les coordonnées des points de capture
#Pour cela 1) on calcule la mediane des coordonéees de chaque idpostedecapture puis 2) on calcule la médiane de la distance à ces coordonnées médianes
#all_data_bf_subset<-left_join(all_data_bf_subset,coords_median_postecapture,by="idpostedecapture")
#all_data_bf_subset$dist_to_median<-sqrt((all_data_bf_subset$latitude-all_data_bf_subset$median_lat)^2+(all_data_bf_subset$longitude-all_data_bf_subset$median_lon)^2)
#median_dist_to_median_coords<-all_data_bf_subset %>% group_by(idpostedecapture) %>% summarise(median_dist_to_median_coords=median(dist_to_median)*(111.32 * 1000 * cos(mean(all_data_bf_subset$latitude) * ((pi / 180)))))

# grâce à ça on corrige à la main les coords de 4PER1, 4PER2, 4PER3, 4PER4, 6PAL3e (BF) et 3LAT2i (CIV) qui sont mauvaises

coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER1i")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER1i")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER1i")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER1i")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER1e")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER1e")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER1e")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER1e")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER2i")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER2i")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER2i")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER2i")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER2e")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER2e")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER2e")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER2e")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER3i")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER3i")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER3i")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER3i")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER3e")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER3e")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER3e")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER3e")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER4i")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER4i")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER4i")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER4i")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "4PER4e")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "1PER4e")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "4PER4e")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "1PER4e")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "6PAL3e")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "6PAL3i")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "6PAL3e")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "6PAL3i")]
coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "3LAT2i")]<-coords_median_postecapture$median_lat[which(coords_median_postecapture$idpostedecapture == "3LAT2e")]
coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "3LAT2i")]<-coords_median_postecapture$median_lon[which(coords_median_postecapture$idpostedecapture == "3LAT2e")]

## Pour la mission 3 au BF on a perdu les données. On prend donc pour coordonnées des postes de captures la moyenne des coordonnées de chaque poste de capture des autres missions
coords_median_postecapture$postedecapture<-substr(coords_median_postecapture$idpostedecapture,2,6)
coords_median_postecapture_mission3_bf<- coords_median_postecapture %>% filter(codepays=="BF") %>% group_by(postedecapture) %>% summarize(median_lat=mean(median_lat),median_lon=mean(median_lon))
coords_median_postecapture_mission3_bf$idpostedecapture<-paste0("3",coords_median_postecapture_mission3_bf$postedecapture)
coords_median_postecapture_mission3_bf$postedecapture<-NULL
coords_median_postecapture_mission3_bf$codepays<-"BF"
coords_median_postecapture$postedecapture<-NULL
coords_median_postecapture<-rbind(as.data.frame(coords_median_postecapture),as.data.frame(coords_median_postecapture_mission3_bf))

## On vérifie les dates
# On récupère pour chaque idpostedecapture le minimum de la colonne 'datedecapture'
dates_captures_from_supervcapture_bf<-df %>% filter(codepays=="BF") %>% group_by(nummission,codevillage,datecapture) %>% summarise(n=n()) 
dates_captures_from_supervcapture_bf<-dates_captures_from_supervcapture_bf %>% filter(n>1) %>% group_by(nummission,codevillage) %>% summarise(date_from_supervcapture=min(datecapture)) 
#dates_captures_from_supervcapture_bf$idpostedecapture<-substr(dates_captures_from_supervcapture_bf$idpostedecapture,1,4)
dates_captures_from_supervcapture_bf$date_from_supervcapture<-as.Date(dates_captures_from_supervcapture_bf$date_from_supervcapture)

# on compare pour le BF avec le fichier de dieudonné
query<-"SELECT * FROM raw_bf_dates_hlc"
df_dates_captures_par_villages_bf<-dbGetQuery(react_gpkg, query)
#df_dates_captures_par_villages_bf$idpostedecapture<-paste0(df_dates_captures_par_villages_bf$n_mission,df_dates_captures_par_villages_bf$code_village)
df_dates_captures_par_villages_bf$date_de_captures<-as.Date(df_dates_captures_par_villages_bf$date_de_captures, format="%d/%m/%Y")
df_dates_captures_par_villages_bf<-df_dates_captures_par_villages_bf[,c("n_mission","code_village","date_de_captures")]
colnames(df_dates_captures_par_villages_bf)<-c("nummission","codevillage","date_from_dieudo")
dates_captures_from_supervcapture_bf<-left_join(dates_captures_from_supervcapture_bf,df_dates_captures_par_villages_bf)
dates_captures_from_supervcapture_bf$corresp<-dates_captures_from_supervcapture_bf$date_from_supervcapture-dates_captures_from_supervcapture_bf$date_from_dieudo

# Pour la grande majorité des données on a les mêmes dates entre le fichier de Dieudo et le fichier de supervcapture. Pour les dates qui ne correspondent pas on prend les dates issues du fichier de Dieudo.
coords_median_postecapture_bf<-coords_median_postecapture %>% filter(codepays=="BF")
coords_median_postecapture_bf$idpostedecapture2<-substr(coords_median_postecapture_bf$idpostedecapture,1,4)
dates_captures_from_supervcapture_bf$idpostedecapture2<-paste0(dates_captures_from_supervcapture_bf$nummission,dates_captures_from_supervcapture_bf$codevillage)
coords_median_postecapture_bf<-left_join(coords_median_postecapture_bf,dates_captures_from_supervcapture_bf)
coords_median_postecapture_bf$idpostedecapture2<-coords_median_postecapture_bf$date_from_supervcapture<-coords_median_postecapture_bf$corresp<-NULL

# pour la CIV on a uniquement le fichier de supervcapture donc c'est celui ci qu"on utilise pour déterminer les dates
dates_captures_from_supervcapture_civ<-df %>% filter(codepays=="CI") %>% group_by(nummission,codevillage,pointdecapture,datecapture) %>% summarise(n=n()) %>% arrange(nummission,codevillage,pointdecapture,datecapture)
#dates_captures_from_supervcapture_civ$idpointdecapture<-paste0(dates_captures_from_supervcapture_civ$nummission,dates_captures_from_supervcapture_civ$codevillage,dates_captures_from_supervcapture_civ$pointdecapture)
#dates_captures_from_supervcapture_civ2<-dates_captures_from_supervcapture_civ %>% group_by(idpointdecapture) %>% summarise(n_dates=n())
#dates_captures_from_supervcapture_civ<-left_join(dates_captures_from_supervcapture_civ,dates_captures_from_supervcapture_civ2)

# n donne le nombre de passages pour un point de capture donné (exterieur et intérieur compris) à une date donnée sur toute la nuit. C'est trop peu. On conserve à partir de 4 points de passage. Pour les données sans date, on verra au cas par cas
dates_captures_from_supervcapture_civ<-dates_captures_from_supervcapture_civ %>% filter(n>=4) %>% group_by(nummission,codevillage,pointdecapture) %>% summarise(date_hlc=min(datecapture)) 
dates_captures_from_supervcapture_civ$idpostedecapture2<-paste0(dates_captures_from_supervcapture_civ$nummission,dates_captures_from_supervcapture_civ$codevillage,dates_captures_from_supervcapture_civ$pointdecapture)

coords_median_postecapture_civ<-coords_median_postecapture %>% filter(codepays=="CI")
coords_median_postecapture_civ$idpostedecapture2<-substr(coords_median_postecapture_civ$idpostedecapture,1,5)
coords_median_postecapture_civ<-left_join(coords_median_postecapture_civ,dates_captures_from_supervcapture_civ,by="idpostedecapture2")
coords_median_postecapture_civ$idpostedecapture2<-NULL


# On vérifie que globalement, pour chaque village et chaque mission, on a 1 seule date
coords_median_postecapture_civ_n_dates<-coords_median_postecapture_civ %>% group_by(nummission,codevillage,date_hlc) %>% summarise(n=n())  %>% arrange(nummission,codevillage,date_hlc)