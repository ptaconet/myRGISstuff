-- Dates and locations of each human landing catch (HLC)
-- output columns : 
-- datecapture : date of human landing catch
-- idpostedecapture : catch place identifier
-- latitude : average of collected latitude
-- longitude: average of collected longitude
-- nummission : campain number 
-- codepays_fk : country identifier 

-- author : Paul Taconet
-- date : 2019-04-11
-- db version : 5

SELECT DISTINCT 
datecapture,
idpostedecapture,
avg(supervcapture.Latitude) as latitude,
avg(supervcapture.longitude) as longitude,
nummission,
codepays_fk 
FROM supervcapture
JOIN village ON village.codevillage_pk=supervcapture.codevillage_fk
group by datecapture,nummission,idpostedecapture,codepays_fk
ORDER BY codepays_fk,nummission,datecapture

