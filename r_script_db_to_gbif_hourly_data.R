library(DBI)
library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)

### connect to the database
path_to_db <- "/home/ptaconet/phd/data/react_db/react_db.gpkg" 
react_gpkg <- DBI::dbConnect(RSQLite::SQLite(),dbname = path_to_db) 

villages <- dbReadTable(react_gpkg, 'recensement_villages_l1') %>% dplyr::select(codevillage,nomvillage)

### create 'events' table
# final column names : EventID	parentEventID	samplingProtocol	samplingEffort	sampleSizeValue	sampleSizeUnit	eventDate	EventTime	LocationID	decimalLatitude	decimalLongitude	geodeticDatum	country	countryCode	locality stateProvince	continent	coordinatePrecision	institutionCode
events <- dbReadTable(react_gpkg, 'entomo_csh_metadata_l1') %>% dplyr::select(-geom) %>% filter(!(nummission %in% c("11","12","13","15"))) %>%
  mutate(year = year(date_capture), month = month(date_capture), day = day(date_capture)) %>%
  mutate(date_heure_fin=ifelse(date_heure_fin=="2016-12-05 08:60:00","2016-12-05 09:00:00",date_heure_fin)) %>%
  mutate(date_heure_debut=as_datetime(date_heure_debut), date_heure_fin = as_datetime(date_heure_fin)) %>%
  mutate(date_heure_debut=round_date(date_heure_debut, unit = "hour"), date_heure_fin = round_date(date_heure_fin, unit = "hour")) %>%
  mutate(sampleSizeValue = 2*as.numeric(date_heure_fin - date_heure_debut)) %>%
  mutate(date_heure_debut = gsub(" ","T",date_heure_debut),date_heure_fin = gsub(" ","T",date_heure_fin)) %>%
  mutate(date_heure_debut = paste0(date_heure_debut,"Z/"),date_heure_fin = paste0(date_heure_fin,"Z")) %>%
  mutate(eventDate = paste0(date_heure_debut,date_heure_fin)) %>%
  mutate(eventTime = paste0(substr(eventDate,12,21),substr(eventDate,33,41))) %>%
  mutate(eventID = idpointdecapture) %>%
  mutate( verbatimEventDate = paste0("survey n°",nummission)) %>%
  left_join(villages) %>%
  rename(countryCode=codepays, LocationID = idpointdecapture, locality = nomvillage, decimalLatitude = Y, decimalLongitude = X) %>%
  mutate(samplingProtocol="Human Landing Catch (both indoors and outdoors)", samplingEffort=paste0("1 Night of Human Landing Catch indoors and outdoors (2 x ",sampleSizeValue/2," observer-hours)"), sampleSizeUnit="hour",geodeticDatum="WGS84", country=ifelse(countryCode=="BF","Burkina Faso","Côte d'Ivoire"),continent="Africa",coordinatePrecision="0.00001",habitat="rural",institutionCode=ifelse(countryCode=="BF","IRD | IRSS","IRD | IPR"), stateProvince=ifelse(countryCode=="BF","Diébougou","Korhogo")) %>%
  dplyr::select(eventID, samplingProtocol,	samplingEffort,	sampleSizeValue,	sampleSizeUnit,	eventDate,eventTime,	year, month, day,  verbatimEventDate, LocationID,	decimalLatitude,	decimalLongitude,	geodeticDatum,	country,	countryCode,	locality,	continent,	coordinatePrecision,	institutionCode,	stateProvince,	habitat) %>%
  arrange(country, LocationID, eventDate )
  

### create 'occurrence' table    eventID	occurrenceID	basisOfRecord	individualCount	organismQuantity	organismQuantityType	occurrenceStatus	scientificName	kingdom	phylum	class	order	family	taxonRank	dateIdentified
occurrence <- dbReadTable(react_gpkg, 'entomo_idmoustiques_l0') %>%
  filter(genre == "Anopheles", !(nummission %in% c("11","12","13","15"))) %>%
  rename(occurrenceID = idmoustique, eventID = idpointdecapture) %>%
  left_join(events %>% dplyr::select(eventID,eventDate,decimalLatitude,	decimalLongitude,	geodeticDatum,	countryCode)) %>% 
  mutate(basisOfRecord = "HumanObservation",
       organismQuantityType = "individuals", 
       occurrenceStatus = "present",
       individualCount = 1,
       organismQuantity = 1,
       kingdom = "Animalia",
       phylum	= "Arthropoda",
       class = "Insecta",
       order = "Diptera",
       family = "Culicidae",
       genus = "Anopheles",
       lifeStage = "adult",
       sex = "female",
       behavior = "host seeking",
       scientificName = ifelse(!is.na(pcr_espece),pcr_espece,especeanoph),
       taxonRank = ifelse(!is.na(pcr_espece),"species","speciesAggregate")) %>%
  filter(!is.na(scientificName)) %>%
  mutate(genericName = case_when(scientificName=="An.funestus_ss" ~ "Anopheles funestus",
         scientificName=="An.gambiae s.l." ~ "Anopheles gambiae sensu lato",
         scientificName=="An.coluzzii" ~ "Anopheles coluzzii",
         scientificName=="An.pharoensis" ~ "Anopheles pharoensi",
         scientificName=="An.gambiae_ss" ~ "Anopheles gambiae sensu stricto",
         scientificName=="An.nili" ~ "Anopheles nili",
         scientificName=="An.ziemanni" ~ "Anopheles ziemanni",
         scientificName=="An.arabiensis" ~ "Anopheles arabiensis",
         scientificName=="An.coustani" ~ "Anopheles coustani",
         scientificName=="An.rufipes" ~ "Anopheles rufipes",
         scientificName=="An.squamosus" ~ "Anopheles squamosus",
         scientificName=="An.flavicosta" ~ "Anopheles flavicosta")) %>%
  mutate(identificationRemarks = ifelse(!is.na(pcr_espece),"morphological identification followed by polymerase chain reaction (PCR)", "morphological identification using morphological keys")) %>% 
  mutate(scientificName = case_when(scientificName=="An.funestus_ss" ~ "Anopheles (Cellia) funestus Giles, 1900",
                                    scientificName=="An.gambiae s.l." ~ "Anopheles gambiae complex (White 1985)",
                                    scientificName=="An.coluzzii" ~ "Anopheles (Cellia) coluzzii Coetzee & Wilkerson, 2013",
                                    scientificName=="An.pharoensis" ~ "Anopheles (Cellia) pharoensis Theobald, 1901",
                                    scientificName=="An.gambiae_ss" ~ "Anopheles (Cellia) gambiae Giles, 1902",
                                    scientificName=="An.nili" ~ "Anopheles (Cellia) nili (Theobald, 1904)",
                                    scientificName=="An.ziemanni" ~ "Anopheles (Cellia) ziemanni Grünberg, 1902",
                                    scientificName=="An.arabiensis" ~ "Anopheles (Cellia) arabiensis Patton, 1905",
                                    scientificName=="An.coustani" ~ "Anopheles (Cellia) coustani Laveran, 1900",
                                    scientificName=="An.rufipes" ~ "Anopheles (Cellia) rufipes (Gough, 1910)",
                                    scientificName=="An.squamosus" ~ "Anopheles (Cellia) squamosus Theobald, 1901",
                                    scientificName=="An.flavicosta" ~ "Anopheles (Cellia) flavicosta Edwards, 1911")) %>%
  mutate(identificationReferences = case_when(identificationRemarks == "morphological identification using morphological keys" ~ "Gillies M, Meillon D. 1968. The Anophelinae of Africa south of the Sahara. Publication of the South African Institute for Medical Research, Johannesburg, 54, 1–343. | Gillies MT. 1987. A supplement to the Anophelinae of Africa south of the Sahara (Afrotropical Region). Publications of the South African Institute for Medical Research, 55, 1–143",
                                              identificationRemarks == "morphological identification followed by polymerase chain reaction (PCR)" ~ "https://doi.org/10.1046/j.1365-2583.2001.00236.x | https://doi.org/10.4269/ajtmh.1993.49.520 | https://doi.org/10.4269/ajtmh.2002.66.804")) %>%
  mutate(eventDate = substr(eventDate,1,10)) %>%
  mutate(nameAccordingTo = "Harbach, R.E. 2013. Mosquito Taxonomic Inventory, https://mosquito-taxonomic-inventory.myspecies.info/, accessed on 25 april 2023") %>%
  mutate(identifiedBy = ifelse(countryCode=="BF","Institut de Recherche en Sciences de la Santé (IRSS)","Institut Pierre Richet (IPR)")) %>%
  mutate(occurrenceID = paste0(idpostedecapture,substr(occurrenceID,9,nchar(occurrenceID)))) %>%
  mutate(occurrenceID = ifelse(codevillage=="KOG" & nummission==8, paste0(occurrenceID,fid) , occurrenceID)) %>%
  arrange(occurrenceID, scientificName)
  
  

n_occur <- data.frame(table(occurrence$occurrenceID), stringsAsFactors = F)
doublons <- n_occur[n_occur$Freq > 1,]$Var1

occurrence <- occurrence %>% 
  filter(!(occurrenceID %in% doublons))


#########
### create 'extended_measurement_fact' table to store attributes for mosquito data at individual level (see https://rs.gbif.org/extension/obis/extended_measurement_or_fact.xml )
#########

extendedMeasurement <- occurrence %>%
  mutate(measurementDeterminedBy = ifelse(countryCode=="BF","Institut de Recherche en Sciences de la Santé (IRSS)","Institut Pierre Richet (IPR)")) %>%
  dplyr::select(eventID, occurrenceID, measurementDeterminedBy, postedecapture, parturite, pcr_pf, kdrw, kdre, ace1, countryCode) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(-c("eventID","occurrenceID","measurementDeterminedBy","countryCode")) %>%
  filter(!is.na(value)) %>%
  rename(measurementValue = value, measurementType = name) %>%
  mutate(measurementType = case_when(measurementType == "postedecapture" ~ "place of collection",
                                     measurementType == "parturite" ~ "parity state",
                                     measurementType == "pcr_pf" ~ "plasmodium falciparum infection",
                                     measurementType == "kdrw" ~ "L1014F (kdr-w) mutation",
                                     measurementType == "kdre" ~ "L1014S (kdr-e) mutation",
                                     measurementType == "ace1" ~ "G119S (ace-1) mutation")) %>%
  mutate(measurementValue = case_when(measurementValue == "e" ~ "exterior of the household",
                                      measurementValue == "i" ~ "interior of the household",
                                      measurementValue == "P" ~ "parous",
                                      measurementValue == "N" ~ "nulliparous",
                                      measurementValue == "0" ~ "FALSE",
                                      measurementValue == "1" ~ "TRUE",
                                      measurementValue == "RR" ~ "homozygous resistant (RR)",
                                      measurementValue == "RS" ~ "heterozygote (RS)",
                                      measurementValue == "SS" ~ "homozygous susceptible (SS)",
                                      TRUE ~ measurementValue
  )) %>%
  mutate(measurementMethod = case_when(measurementType == "parity state" ~ "dissection (Detinova TS; 1962; 47: 13–191. PMID: 13885800)",
                                       measurementType == "plasmodium falciparum infection" ~ "quantitative polymerase chain reaction (qPCR) | https://doi.org/10.1371/journal.pone.0054820",
                                       measurementType == "L1014F (kdr-w) mutation" & countryCode == "BF" ~ "polymerase chain reaction (PCR) | https://doi.org/10.1046/j.1365-2583.1998.72062.x",
                                       measurementType == "L1014F (kdr-w) mutation" & countryCode == "CI" ~ "quantitative polymerase chain reaction (qPCR) | https://doi.org/10.1186/1475-2875-6-111",
                                       measurementType == "G119S (ace-1) mutation" ~ "quantitative polymerase chain reaction (qPCR) | https://doi.org/10.1186/1475-2875-12-404" ,
                                       measurementType == "L1014S (kdr-e) mutation" & countryCode == "BF" ~ "polymerase chain reaction (PCR) | https://doi.org/10.1046/j.1365-2583.2000.00209.x" )) %>%
    filter(!(measurementValue %in% c("ND","IC"))) %>%
  dplyr::select(-countryCode) %>%
  mutate(measurementDeterminedBy = ifelse(measurementType == "place of collection",NA,measurementDeterminedBy)) %>%
  group_by(occurrenceID) %>%
  mutate(index2 = row_number()) %>%
  ungroup() %>%
  mutate(measurementID = paste0(occurrenceID,"_",index2)) %>%
  dplyr::select(-index2) %>%
  relocate(measurementDeterminedBy, .after = measurementMethod ) %>%
  relocate(measurementID, 1) %>%
  arrange(measurementID)


occurrence <- occurrence %>%
  dplyr::select(occurrenceID, eventID, basisOfRecord, eventDate, 	scientificName, genericName, nameAccordingTo,  identificationRemarks, identifiedBy, identificationReferences,	kingdom,	phylum,	class,	order,	family,genus, taxonRank, individualCount,organismQuantity,organismQuantityType,occurrenceStatus, lifeStage, sex, behavior, decimalLatitude,	decimalLongitude,	geodeticDatum,	countryCode, identificationReferences)

#########
### create 'measurement_fact' table for micro-climatic data + quality control data (see https://rs.gbif.org/extension/dwc/measurements_or_facts_2022-02-02.xml )
#########

# spatial-only explanatory variables
env_spatial <- dbReadTable(react_gpkg,'env_spatial') %>% dplyr::select(-fid) %>% dplyr::rename(idpointdecapture = id) %>% filter( var!="POH") %>% filter(buffer == 2000)
env_spatial <- env_spatial %>%
  dplyr::select(-buffer) %>%
  rename(eventID = idpointdecapture) %>%
  mutate(val = as.character(val))


BCH_POP <- dbReadTable(react_gpkg,'env_spatial') %>% dplyr::select(-fid) %>% dplyr::rename(idpointdecapture = id) %>% filter(var %in% c("BCH","POP","POH"), buffer==500) %>% mutate(buffer=2000) %>% mutate(val = as.character(val)) %>%
  dplyr::select(-buffer) %>%
  rename(eventID = idpointdecapture)

env_spatial <- rbind(env_spatial,BCH_POP)


# non-spatial explanatory variables
env_static <- dbReadTable(react_gpkg, 'env_static') %>% dplyr::select(-fid) %>% dplyr::rename(eventID = id) %>%
  mutate(val = case_when(val %in% c("Ctrle","IEC") ~ "LLIN",
                         #VCM == "IEC" ~ "LLIN + IEC",
                         val == "IRS" ~ "LLIN + IRS",
                         val == "IVM" ~ "LLIN + IVM",
                         val == 'Larvicide' ~ 'LLIN + Larv.'))


# variables for the night of catch
env_nightcatch <- dbReadTable(react_gpkg, 'env_nightcatch') %>% dplyr::select(-fid) %>% dplyr::rename(eventID = id) %>% dplyr::filter(var != "WDR") %>%
  mutate(val = as.character(val))


# data landcover
th_env_nightcatch_postedecapture <- read.csv("/home/ptaconet/Bureau/data_anglique.csv") %>%
  dplyr::select(idpostedecapture,date_time,temperature_hygro,humidity_hygro,pressure_baro,luminosite_hobo) %>%
  dplyr::rename(NMT = temperature_hygro, NMH = humidity_hygro, NMA = pressure_baro, NML = luminosite_hobo) %>%
  mutate(date = as.Date(date_time), heuredecapture = hour(date_time)) %>%
  mutate(idpointdecapture = substr(idpostedecapture,1,5)) %>%
  dplyr::group_by(idpointdecapture) %>%
  dplyr::summarise(NMT = mean(NMT,na.rm = T), NMH = mean(NMH,na.rm = T), NMA=mean(NMA,na.rm = T), NML=mean(NML,na.rm = T)) %>%
  as_tibble() %>%
  mutate(NMT = round(NMT,1), NMH = round(NMH,2), NMA = round(NMA,2), NML = round(NML,2)) %>%
  mutate(idpointdecapture=as.character(idpointdecapture)) %>%
  mutate(NMT = ifelse(NMT<10,NA,NMT),NMA = ifelse(NMA<800,NA,NMA),NML = ifelse(NML>2000,2000,NML))

th_env_nightcatch_postedecapture <- th_env_nightcatch_postedecapture %>% pivot_longer(-idpointdecapture) %>% rename(val = value, var = name, eventID = idpointdecapture) %>%
  mutate(val = as.character(val))

env_landcover <-  dbReadTable(react_gpkg, 'env_landcover') %>% dplyr::select(-fid) %>% dplyr::rename(idpointdecapture = id) %>% filter(layer_id %in% c(3,8)) %>% filter(buffer == 2000, metric == "pland")

env_landcover_classes <- dbReadTable(react_gpkg, 'lco_metadata') %>% filter(layer_id %in% c(3,8)) %>% dplyr::select(pixval, layer_id, pixlabel_english) 

env_landcover <- env_landcover %>% 
  left_join(env_landcover_classes) %>% 
  mutate(var = paste0("% of landscape occupied by the class \'",tolower(pixlabel_english),"\' in a 2-km radius buffer zone around the sampling site")) %>%
  dplyr::select(idpointdecapture, var, val) %>%
  rename(eventID = idpointdecapture) %>%
  mutate(val = as.character(val))


# data meteo
env_spatiotemporal <- dbReadTable(react_gpkg, 'env_spatiotemporal') %>% dplyr::select(-fid) %>% mutate(date = as.Date(date)) %>% dplyr::rename(idpointdecapture = id) %>% filter(var %in% c("TMAX1","TMIN1","RFD1F"), buffer == 2000, lag_time >= 0 , lag_time <= 42)

fun_summarize_week <- function(var_to_summarize){
  
  if(var_to_summarize=="RFD1F"){
    env_spatiotemporal_summarize <- env_spatiotemporal %>%
      filter(var==var_to_summarize) %>%
      group_by(idpointdecapture,buffer,lag_n = lubridate::week(date)) %>%
      summarise(val=sum(val, na.rm = T),date = min(date)) %>%
      group_by(idpointdecapture,buffer) %>%
      mutate(lag_n=seq(n()-1,0,-1)) %>%
      mutate(var = gsub("1","7",var_to_summarize)) %>%
      as_tibble()
  } else {
    env_spatiotemporal_summarize <- env_spatiotemporal %>%
      filter(var==var_to_summarize) %>%
      group_by(idpointdecapture,buffer,lag_n = lubridate::week(date)) %>%
      summarise(val=mean(val, na.rm = T),date = min(date)) %>%
      group_by(idpointdecapture,buffer) %>%
      mutate(lag_n=seq(n()-1,0,-1)) %>%
      mutate(var = gsub("1","7",var_to_summarize)) %>%
      as_tibble()
  }
  return(env_spatiotemporal_summarize)
  
}

env_meteo <- fun_summarize_week("TMAX1") %>%
  bind_rows(fun_summarize_week("TMIN1")) %>%
  bind_rows(fun_summarize_week("RFD1F"))

env_meteo <- env_meteo %>% 
  mutate(var = case_when(var=="TMAX7" ~ "Daytime land surface temperature",
                         var=="TMIN7" ~ "Nighttime land surface temperature",
                         var=="RFD7F" ~ "Rainfall")) %>%
  mutate(var = ifelse(lag_n == 0, paste0(var," the week preceding the event"), paste0(var," ",lag_n," weeks preceding the event"))) %>%
  filter(lag_n <= 6) %>%
  dplyr::select(-c("buffer","lag_n","date")) %>%
  rename(eventID = idpointdecapture) %>%
  mutate(val = as.character(val))





source("/home/ptaconet/phd/r_scripts/data_analysis_tests/functions_script_data_analysis.R")
entomo_csh_metadata_l1 <- dbReadTable(react_gpkg, 'entomo_csh_metadata_l1') %>% dplyr::select(-geom) %>% filter(!(nummission %in% c("11","12","13","15")))

hum_behav_bf <- load_hmnbehav_data("BF", entomo_csh_metadata_l1)
LUS_bf = hum_behav_bf[[1]]
hum_behav_4_earlylatebiting_bf = hum_behav_bf[[3]]

hum_behav_ci <- load_hmnbehav_data("CI", entomo_csh_metadata_l1)
LUS_ci = hum_behav_ci[[1]]
hum_behav_4_earlylatebiting_ci = hum_behav_ci[[3]]

hum_behav <- rbind(hum_behav_4_earlylatebiting_bf,hum_behav_4_earlylatebiting_ci) %>% 
  rename(eventID = idpointdecapture) %>%
  pivot_longer(-eventID) %>%
  rename(var = name, val = value) %>%
  mutate(val = as.character(val))

LUS <- rbind(LUS_bf,LUS_ci) %>% 
  rename(eventID = idpointdecapture) %>%
  rename(val = LUS) %>%
  mutate(var = "LUS") %>%
  mutate(val = as.character(val))


# load time since vector control measure
time_since_vc_bf <- load_time_since_vc("BF", entomo_csh_metadata_l1)
time_since_vc_ci <- load_time_since_vc("CI", entomo_csh_metadata_l1)

time_since_vc <- rbind(time_since_vc_bf,time_since_vc_ci) %>%
  rename(eventID = idpointdecapture) %>%
  rename(val = VCT) %>%
  mutate(var = "VCT") %>%
  dplyr::select(-VCT3) %>%
  mutate(val = as.character(val))


googlesheets4::sheets_deauth()
prediction_vars <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1dIeSOa2WinXvOQGLmIjA0gFdsHnb6zMMsDME-G5pyMc/edit?usp=sharing", sheet = "var_explication", col_types="c")
prediction_vars <- prediction_vars %>%
  dplyr::select(code, long_name, unit) %>%
  rename(var = code)

measurements <- env_spatial %>%
  bind_rows(env_static) %>%
  bind_rows(env_nightcatch) %>%
  bind_rows(th_env_nightcatch_postedecapture) %>%
  bind_rows(env_landcover) %>%
  bind_rows(env_meteo) %>%
  bind_rows(hum_behav) %>%
  bind_rows(LUS) %>%
  bind_rows(time_since_vc)


measurements <- measurements %>%
  filter(eventID %in% events$eventID, !is.na(val)) %>%
  mutate(c = nchar(var)) %>%
  #filter(c > 5 | var %in% c("RFH","WLS","WMD","BDE","BCH","NMT","NMH","NML","WSP","LUS","VCT","HBI2","HBB2","VCM","POP","LMN")) %>%
  filter(c > 5, !(var %in% c("% of landscape occupied by the class 'main tracks' in a 2-km radius buffer zone around the sampling site","% of landscape occupied by the class 'main road' in a 2-km radius buffer zone around the sampling site"))) %>%
  left_join(prediction_vars) %>%
  mutate(measurementType = ifelse(is.na(long_name),var,long_name)) %>%
  rename(measurementValue = val, measurementUnit = unit) %>%
  mutate(measurementType = case_when(measurementType == "Total length of the hydrographic network" ~ "Total length of the hydrographic network in a 2-km radius buffer zone around the event site", 
                                     measurementType == "Total human population (census/ground data)" ~ "Human population in the sampling village", 
                                     measurementType == "Degree of clustering or ordering of the households" ~ "Degree of clustering or ordering of the households in the sampling village", 
                                     measurementType == "Vector control tool type (LLIN / IRS / anti-larval fighting / communication strategies)" ~ "Vector control tool in use in the village at the time of sampling", 
                                     measurementType == "Proportion of half-hours with positive precipitation for the whole duration of the sampling" ~ "Proportion of half-hours with positive precipitation during the night of sampling",
                                     measurementType == "Mean wind speed direction during the night of catch" ~ "Mean wind speed direction during the night of sampling", 
                                     measurementType == "Mean temperature during the night of catch" ~ "Mean temperature during the night of sampling", 
                                     measurementType == "Mean humidity during the night of catch" ~ "Mean humidity during the night of sampling", 
                                     measurementType == "Mean luminosity during the night of catch" ~ "Mean luminosity during the night of sampling", 
                                     measurementType == "HBI2" ~ "Average time spent indoors by the population of the village during the night of sampling", 
                                     measurementType == "HBB2" ~ "Average time spent under LLIN by the population of the village during the night of sampling", 
                                     measurementType == "LUS" ~ "Average LLIN use rate among the population of the village on the night of sampling", 
                                     measurementType == "Numer of days since introduction of Vector control tools (to study rebounding residual transmission, cf kileen)" ~ "Time since last LLIN distribution", 
                                     measurementType == "Apparent magnitude of the Moon" ~ "Average apparent magnitude of the Moon during the night of sampling", 
                                     TRUE ~ measurementType
  )) %>%
  mutate(measurementRemarks = case_when(grepl("landscape|Degree of clustering|hydrographic network", measurementType) ~ "Landscape",
                                        grepl("land surface temperature|Rainfall", measurementType) ~ "Meteorology preceding the event",
                                        grepl("Mean wind speed|Mean temperature|Mean humidity|Mean luminosity|half-hours with positive precipitation|Moon", measurementType) ~ "Micro-climate during the event",
                                        grepl("Time since last LLIN distribution|Vector control tool in use", measurementType) ~ "Vector control",
                                        grepl("Average time spent|Average LLIN use rate|Human population", measurementType) ~ "Host availability",
  )) %>%
  mutate(measurementUnit = case_when(grepl("landscape", measurementType) ~ "% of landscape", 
                                     grepl("hydrographic network", measurementType) ~ "meters", 
                                     grepl("wind speed", measurementType) ~ "meters/second", 
                                     grepl("Time since last LLIN distribution", measurementType) ~ "number of months", 
                                     grepl("Average time spent under LLIN|Average time spent indoors", measurementType) ~ "number of hours", 
                                     grepl("Human population in the sampling village", measurementType) ~ "number of persons", 
                                     grepl("land surface temperature", measurementType) ~ "degree Celsius", 
                                     grepl("Rainfall", measurementType) ~ "cumulated millimeters", 
                                     TRUE ~ measurementUnit)) %>%
  mutate(measurementMethod = case_when(grepl("landscape|Degree of clustering|hydrographic network", measurementType) ~ "Derived from satellite products (detailed method : https://doi.org/10.1186/s13071-021-04851-x)",
                                       grepl("land surface temperature|Rainfall", measurementType) ~ "Derived from satellite products (detailed method : https://doi.org/10.1186/s13071-021-04851-x)",
                                       grepl("Average time spent|Average LLIN use rate|Human population", measurementType) ~ "Derived from population census and human behaivour surveys (detailed method : https://doi.org/10.1186/s12889-021-10304-y)",
                                       measurementRemarks %in% c("Micro-climate during the event") ~ "Measured with data loggers during the entomological collections (detailed method : https://doi.org/10.1101/2022.08.20.504631)",
                                       measurementRemarks %in% c("Vector control") ~ "https://doi.org/10.1101/2022.08.20.504631"
  )) %>%
  # mutate(measurementDeterminedBy = case_when(measurementRemarks %in% c("Landscape","Meteorology preceding the event") ~ "Paul Taconet",
  #                                            grepl("Mean wind speed|half-hours with positive precipitation|Moon", measurementType) ~ "Paul Taconet",
  #                                            grepl("Mean wind speed|Mean temperature|Mean humidity|Mean luminosity|half-hours with positive precipitation|Moon", measurementType) ~ "Micro-climate during the event",
  #                                            
  #   
  # )) %>%
  arrange(eventID,measurementRemarks) %>% 
  group_by(substr(eventID, 1, 5)) %>%
  mutate(index2 = row_number()) %>%
  ungroup() %>%
  mutate(measurementID = paste0(substr(eventID, 1, 5),"_",index2)) %>%
  dplyr::select(measurementID, eventID,measurementType, measurementRemarks, measurementValue, measurementUnit, measurementMethod)




write.csv(events,"events.csv", quote = F, row.names = F)
write.table(occurrence,"occurrence.csv", quote = F, row.names = F, sep = ";", dec = ".")
write.csv(measurements,"measurements.csv", quote = F, row.names = F)
write.csv(extendedMeasurement,"extendedMeasurement.csv", quote = F, row.names = F)

