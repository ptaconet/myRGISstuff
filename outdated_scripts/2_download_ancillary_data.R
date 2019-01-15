######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = download_ancillary_data, title = Download from the internet a set ancillary data that will be used for the classification or the modeling, abstract = This script downloads a set of ancillary data from the internet that will afterwards be used for the work (classification, modeling) on both ROIs (BF and CIV). 

Sentinel2_products_id<-c("bc6bafd7-d44f-4d62-8754-3cc4ba4e8cc0","18895056-852f-4a4f-a3aa-ca7882fe79de")
Sentinel2_output_folder_path<-"/home/ptaconet/Documents/react/data_BF/HR_Sentinel2/raw_data"

source("functions_data_download/download_copernicus_scihub_products.R")

#### Sentinel optical and radar images
# First we went on the Copernicus scihub (https://scihub.copernicus.eu/dhus/#/home) to identify the images of interest. Then we download them through the Scihub API using their UUID. The following link [https://scihub.copernicus.eu/userguide/BatchScripting] presents how to use the Scihub API in bash to query and download data. Example : wget --content-disposition --continue --user=**** --password=**** "https://scihub.copernicus.eu/dhus/odata/v1/Products('19521a0a-996b-4579-90eb-c0a64c46e4ba')/\$value"

## Sentinel 2 
# S2 data for the BF area are available through 2 tiles (T30PVS and T30PVT). We chose the images following the criteria: no cloud + dates of the ground throuth data collection
# T30PVS uuid : bc6bafd7-d44f-4d62-8754-3cc4ba4e8cc0 and T30PVT uuid : 18895056-852f-4a4f-a3aa-ca7882fe79de
for (i in 1:length(Sentinel2_products_id)){
cat("Downloading Sentinel 2 products...")
res<-download_copernicus_scihub_products("ptaconet","***",Sentinel2_products_id[i],Sentinel2_output_folder_path)
cat(res)
}


# S2 data for the CIV area are available through 1 tile (T29PRL). We chose the images following the criteria: no cloud + dates of the ground throuth data collection
# T29PRL uuid : 6236fb46-41c1-4950-b6c8-602c48b90049
system("wget --content-disposition --continue --user=**** --password=**** \"https://scihub.copernicus.eu/dhus/odata/v1/Products('6236fb46-41c1-4950-b6c8-602c48b90049')/\\$value\\")  ## NOT RUN need to specify the username and the password 

## Sentinel 1



#### SRTM Digital Elevation Model. 
# First we went on the USGS EarthExplorer () to identify the tiles corresponding to our ROIs and related products URLs. Then we download them 

# SRTM is available through 2 tiles for the BF area (N10W004 and N11W004) and through 2 tiles for the CIV area (N08W006 and N09W006).
# The bash script enables to download the data from EarthExplorer using the API. The script 
system("/home/ptaconet/r_react/object_based_image_analysis/2_download_SRTM_from_earthexplorer.sh")