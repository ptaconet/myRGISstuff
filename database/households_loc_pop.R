# Script to correct the households data from the orignial data stored in the DB v7
require(RSQLite)
require(dplyr)
react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)
sql_query_households_loc_pop<-"SELECT
raw_menages.codemenage as codemenage, count(codeindividu) as population, raw_menages.Y as latitude,raw_menages.X as longitude, raw_villages.codevillage as village,raw_villages.codepays
FROM
raw_individus
JOIN raw_menages on raw_individus.codemenage=raw_menages.codemenage
JOIN raw_villages on raw_villages.codevillage=raw_menages.codevillage
GROUP BY raw_individus.codemenage,raw_menages.codemenage,raw_menages.X,raw_menages.Y,raw_villages.codevillage,raw_villages.codepays
ORDER BY raw_villages.codevillage"

df_households_loc_pop<-dbGetQuery(react_gpkg, sql_query_households_loc_pop)

# Tidy the dataset

# Attention les données de localisation des ménages pour BLA et NAM sont FAUSSES !! on ne les intégre ici que pour avoir l'étendue du village.
# Cependant les données de population sont justes.
df_households_loc_pop$latitude[which(df_households_loc_pop$village=="BLA")]<-8.948
df_households_loc_pop$longitude[which(df_households_loc_pop$village=="BLA")]<--5.652
df_households_loc_pop$latitude[which(df_households_loc_pop$village=="NAM")]<-8.8845
df_households_loc_pop$longitude[which(df_households_loc_pop$village=="NAM")]<--5.75

df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="BLA028")]<-8.94788
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="BLA028")]<--5.65341
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="BLA030")]<-8.94994
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="BLA030")]<--5.65308
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="BLA032")]<-8.94965
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="BLA032")]<--5.65098
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="BLA033")]<-8.94786
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="BLA033")]<--5.65146

df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM044")]<-8.8838
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM044")]<--5.7512
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM045")]<-8.8879
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM045")]<--5.7504
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM046")]<-8.8880
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM046")]<--5.7485
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM047")]<-8.8863
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM047")]<--5.7476
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM048")]<-8.8842
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM048")]<--5.7490
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM049")]<-8.8844
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM049")]<--5.7497
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM050")]<-8.8823
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM050")]<--5.7501
df_households_loc_pop$latitude[which(df_households_loc_pop$codemenage=="NAM051")]<-8.8846
df_households_loc_pop$longitude[which(df_households_loc_pop$codemenage=="NAM051")]<--5.7513


df_households_loc_pop$latitude[which(df_households_loc_pop$latitude==9)]<-NA
df_households_loc_pop$longitude[which(df_households_loc_pop$longitude==-5)]<-NA

df_households_loc_pop_mean_pos_by_vill<-df_households_loc_pop %>% filter(!is.na(df_households_loc_pop$latitude)) %>% group_by(village) %>%  summarise(latitude=mean(latitude),longitude=mean(longitude))

df_households_loc_pop_without_coords<-df_households_loc_pop %>% filter(is.na(latitude)) %>% dplyr::select(-c(latitude,longitude))
df_households_loc_pop<-df_households_loc_pop %>% filter(!is.na(df_households_loc_pop$latitude))

df_households_loc_pop_without_coords<-merge(df_households_loc_pop_without_coords,df_households_loc_pop_mean_pos_by_vill,by="village")

df_households_loc_pop<-rbind(df_households_loc_pop,df_households_loc_pop_without_coords)

colnames(df_households_loc_pop)[which(colnames(df_households_loc_pop)=="village")]<-"codevillage"

# Correction des points abérrants

df_households_loc_pop[571,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="KAR"),c("latitude","longitude")]
df_households_loc_pop[581,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="KAR"),c("latitude","longitude")]
df_households_loc_pop[1521,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="NBO"),c("latitude","longitude")]
df_households_loc_pop[1734,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="NOT"),c("latitude","longitude")]
df_households_loc_pop[309,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="FEL"),c("latitude","longitude")]

# Ajout attention
df_households_loc_pop$note_importante<-NA
df_households_loc_pop$note_importante[which(df_households_loc_pop$codevillage %in% c("BLA","NAM"))]<-"Attention les données de localisation (X et Y) des ménages pour BLA et NAM sont FAUSSES !! on ne les a intégré ici que pour avoir l'étendue spatiale maximale du village. Cependant les données de population (nombre d’habitants) sont justes."
  
