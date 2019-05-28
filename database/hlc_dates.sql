-- Dates of each human landing catch (HLC) and the associated mission number for each country
-- output columns : 
-- datecapture : date of human landing catch
-- nummission : campain number 
-- codepays_fk : country identifier 

-- author : Paul Taconet
-- date : 2019-04-11
-- db version : 5

SELECT DISTINCT
datecapture,
nummission,
codepays_fk 
FROM supervcapture
JOIN village ON village.codevillage_pk=supervcapture.codevillage_fk
ORDER BY codepays_fk,nummission,datecapture
