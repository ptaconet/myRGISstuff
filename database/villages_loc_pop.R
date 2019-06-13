#require(RSQLite)
#require(dplyr)
#react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)

query<-"select sum(villages_households_loc_pop.population) as population, villages_households_loc_pop.codevillage, raw_villages.codepays, raw_villages.X, raw_villages.Y
from villages_households_loc_pop
left join raw_villages on  villages_households_loc_pop.codevillage=raw_villages.codevillage
group by villages_households_loc_pop.codevillage
order by raw_villages.codepays, population desc
"
df_villages_loc_pop<-dbGetQuery(react_gpkg, query)
