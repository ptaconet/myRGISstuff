require(RSQLite)
require(dplyr)
#react_gpkg <- dbConnect(RSQLite::SQLite(),path_to_gpkg_database)

query<-"select sum(households_loc_pop.population) as population, households_loc_pop.village, raw_villages.codepays, raw_villages.X, raw_villages.Y
from households_loc_pop
left join raw_villages on  households_loc_pop.village=raw_villages.codevillage
group by households_loc_pop.village
order by raw_villages.codepays, population desc
"
df_villages_loc_pop<-dbGetQuery(react_gpkg, query)
