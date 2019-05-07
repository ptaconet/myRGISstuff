-- Location, population and village of each household. Data extracted from the census made in the frame of the project.
-- output columns : 
-- codemenage : code of the household
-- population : population of the household
-- latitude: latitude the household
-- longitude: longitude of the household
-- village : code of the village
-- codepays_fk : code of the country
-- author : Paul Taconet
-- date : 2019-05-06
-- db version : 5

SELECT
menage.codemenage_pk as codemenage, count(codeindividu_pk) as population, menage.Latitude as latitude,menage.Longitude as longitude, village.codevillage_pk as village,village.codepays_fk
FROM
individu
JOIN menage on individu.codemenage_fk=menage.codemenage_pk
JOIN village on village.codevillage_pk=menage.codevillage_fk
GROUP BY individu.codemenage_fk,menage.codemenage_pk,menage.Latitude,menage.Longitude,village.codevillage_pk,village.codepays_fk
ORDER BY village.codevillage_pk