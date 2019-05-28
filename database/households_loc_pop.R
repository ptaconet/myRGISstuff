# Script to correct the households data from the orignial data stored in the DB v7
require(RSQLite)
require(dplyr)
#react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_empty_database)
sql_query_households_loc_pop<-"SELECT
raw_menages.codemenage_pk as codemenage, count(codeindividu_pk) as population, raw_menages.Y as latitude,raw_menages.X as longitude, raw_villages.codevillage_pk as village,raw_villages.codepays_fk
FROM
raw_individus
JOIN raw_menages on raw_individus.codemenage_fk=raw_menages.codemenage_pk
JOIN raw_villages on raw_villages.codevillage_pk=raw_menages.codevillage_fk
GROUP BY raw_individus.codemenage_fk,raw_menages.codemenage_pk,raw_menages.X,raw_menages.Y,raw_villages.codevillage_pk,raw_villages.codepays_fk
ORDER BY raw_villages.codevillage_pk"

df_households_loc_pop<-dbGetQuery(react_gpkg, sql_query_households_loc_pop)

# Tidy the dataset
df_households_loc_pop$latitude[which(df_households_loc_pop$latitude==9)]<-NA
df_households_loc_pop$longitude[which(df_households_loc_pop$longitude==-5)]<-NA

df_households_loc_pop_mean_pos_by_vill<-df_households_loc_pop %>% filter(!is.na(df_households_loc_pop$latitude)) %>% group_by(village) %>%  summarise(latitude=mean(latitude),longitude=mean(longitude))

df_households_loc_pop_without_coords<-df_households_loc_pop %>% filter(is.na(latitude)) %>% dplyr::select(-c(latitude,longitude))
df_households_loc_pop<-df_households_loc_pop %>% filter(!is.na(df_households_loc_pop$latitude))

df_households_loc_pop_without_coords<-merge(df_households_loc_pop_without_coords,df_households_loc_pop_mean_pos_by_vill,by="village")

df_households_loc_pop<-rbind(df_households_loc_pop,df_households_loc_pop_without_coords)

# Correction des points abérrants
df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="KOL"),c("latitude","longitude")]<-c(9.28812,-5.52354)

df_households_loc_pop[571,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="KAR"),c("latitude","longitude")]
df_households_loc_pop[581,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="KAR"),c("latitude","longitude")]
df_households_loc_pop[1521,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="NBO"),c("latitude","longitude")]
df_households_loc_pop[1734,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="NOT"),c("latitude","longitude")]
df_households_loc_pop[2769:2778,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="KOL"),c("latitude","longitude")]
df_households_loc_pop[686:687,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="KOL"),c("latitude","longitude")]
df_households_loc_pop[309,c("latitude","longitude")]<-df_households_loc_pop_mean_pos_by_vill[which(df_households_loc_pop_mean_pos_by_vill$village=="FEL"),c("latitude","longitude")]

