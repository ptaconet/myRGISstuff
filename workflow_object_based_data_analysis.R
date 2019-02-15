######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = workflow_object_based_data_analysi, title = Workflow for the land cover classification using object based image analysis, abstract = This script is a workflow for the object based image analysis work used to create the land cover maps. It performs all the pre-processing, data preparation, classification and post-processing steps. The script uses applications of the Orfeo Toolbox v6.6.1 (https://www.orfeo-toolbox.org/), the GRASS library v7.4 (https://grass.osgeo.org)  and the GDAL library v2.2.1 ()
# wps.in: id = path_to_otbApplications_folder, type = string, title = Path to the folder containing the OTB applications. , value = "/home/ptaconet/OTB-6.6.1-Linux64/bin";
# wps.in: id = path_to_grassApplications_folder, type = string, title = Path to the folder containing the GRASS applications. , value = "/usr/lib/grass74";
# wps.in: id = path_to_processing_folder, type = string, title = Path to working folder , value = "/home/ptaconet/Documents/react/data_CIV";
# wps.in: id = path_to_roi_vector, type = string, title = Path to the ROI in kml or shp format, value = "ROI.kml";
# wps.in: id = path_to_spot67_raw_folder, type = string, title = Path to the folder where input data (i.e. Spot6/7 products as tar.gz files) are stored , value = "VHR_SPOT6/raw_data";
# wps.in: id = path_to_dem_raw_folder, type = string, title = Path to the folder containing the SRTM Digital Elevation Model files. value = "DEM_SRTM/raw_data";
# wps.out: id = output_zip, type = text/zip, title = Folder path_to_spot67_preprocessed_folder containing the data pre-processed: for each input product, the orthorectified multispectral image (extracted on the ROI), the orthorectified panchromatic image (extracted on the ROI), the orthorectified pansharpened image (extracted on the ROI). In addition, if there are multiple input products, the previously mentioned products will be mosaiced 

#### Workflow :
### Step 1 - Download the SRTM tiles for the ROI
### Step 2 - Pre-process the Spot6/7 image(s) :
  ## 2.1 - fusion the tiles of the panchromatic image
  ## 2.2 - convert the multispectral and panchromatic images from digital numbers to TOA reflectance
  ## 2.3 - orthorectifying the multispectral and panchromatic images
  ## 2.4 - extract the ROI
  ## 2.5 - pansharpen
  ## 2.6 - mosaic the various tiles covering the ROI (if relevant)
### Step 3 - Download the ancillary data :
  ## 3.1 - Download Sentinel 2 image(s)
### Step 4 - Preprocess the ancillary data :
  ## 4.1 - preprocess the DEM : mosaic the various tiles covering the ROI (if relevant), and then extract the ROI
  ## 4.2 - preprocess the Sentinel 2 image(s) : mosaic the various images covering the ROI (if relevant), and then extract the ROI
### Step 5 - Prepare the data for the classification : 
  ## 5.1 - extract indices from the DEM : slope, aspect, flow accumulation, flow direction, topographic convergence index
  ## 5.2 - extract textural indices from the Spot6/7 panchromatic image : 
  ## 5.3 - extract radiometric indices from the Spot6/7 pansharpened image : 
  ## 5.4 - extract radiometric indices from the S2 image : 
  ## 5.5 - Split the bands of the Spot6/7 image
### Step 6 - Segment the Spot6/7 image 
### Step 7 - Extract zonal statistics (radiometric, textural, ancillary, shape, contextual)
  ## 7.1 - Extract zonal statistics for the ground truth dataset
  ## 7.2 - Extract zonal statistics for the segmented objects dataset
### Step 8 - Train the classifier : generate the classification model using random forest
  ## 8.1 - Parameterize the model: features to keep, optimum number of trees, optimum number of variables to test at each division
### Step 9 - Classify
  ## 9.1 Classify the segmented objects using the model prepared at step 9
  ## 9.2 Rasterize the output
### Step 10 - Post - process
  ## 10.1 Get clean raster dataset
  ## 10.2 Get clear vector dataset

########################################################################################################################
############ Set Input parameters for the workflow ############
########################################################################################################################

### Global variables used throughout the WF
path_to_otbApplications_folder<-"/home/ptaconet/OTB-6.6.1-Linux64/bin"
path_to_grassApplications_folder<-"/usr/lib/grass74" #<Can be retrieved with grass74 --config path . More info on the use of rgrass7 at https://grasswiki.osgeo.org/wiki/R_statistics/rgrass7
path_to_processing_folder<-"/home/ptaconet/Documents/react/data_BF"  #<Path to the processing folder (i.e. where all the data produced by the workflow will be stored)>
path_to_roi_vector="ROI.kml" #<Path to the Region of interest in KML format>


### Parameters for step 3 : Downloading ancillary data
copernicus_scihub_username="ptaconet" #<Copenicus scihub username>
copernicus_scihub_password="****"  #<Copenicus scihub password>
Sentinel2_products_id<-c("bc6bafd7-d44f-4d62-8754-3cc4ba4e8cc0","18895056-852f-4a4f-a3aa-ca7882fe79de") #<Ids of the products to download in the Copernicus Scihub>
# BF: Sentinel2_products_id<-c("bc6bafd7-d44f-4d62-8754-3cc4ba4e8cc0","18895056-852f-4a4f-a3aa-ca7882fe79de")
# CIV: Sentinel2_products_id<-c("6236fb46-41c1-4950-b6c8-602c48b90049")

### Parameters for step 4
proj_srs="+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs" #<proj srs for the ROI>
threshold_accumulation_raster<-800   #<Threshold for the water accumulation raster file. All cells above this threshold are considered as the hydrographic network.>
xrad=c(5,9,17)  #<size of the moving windows in the x direction for the computation of the textures>
yrad=c(5,9,17)  #<size of the moving windows in the y direction for the computation of the textures>
nbbin=64  #<number of bins in the output image for the computation of the textures>
# BF: threshold_accumulation_raster<-1000
# CIV: threshold_accumulation_raster<-800

### Parameters for step 5
radiometric_indices_list_spot67<-c("Vegetation:NDVI","Water:NDWI","Soil:BI") #<Radiometric indices to compute for the Spot6/7 image>

### Parameters for step 6
segmentation_threshold=100       #<segmentation scale parameter> 
segmentation_cw=0.1              #<segmentation shape parameter> 
segmentation_sw=0.9             #<segmentation spectral parameter>
# BF: segmentation_threshold<-100
# CIV: threshold_accumulation_raster<-160
# BF: segmentation_cw<-0.1
# CIV: segmentation_cw<-0.1
# BF: segmentation_sw<-0.9
# CIV: segmentation_sw<-0.8

### Parameters for step 7
path_to_groundtruth_folder<-file.path(path_to_processing_folder,"Ground_truth")
path_to_ground_truth_data<-file.path(path_to_groundtruth_folder,"groundtruth_bf.gpkg") #<Path to the ground truth dataset. The geometry column must be named "geom">
methods_to_compute<-"average,stddev" #<methods_to_compute for the primitives. Available methods are: "minimum,maximum,range,average,stddev,variance,coeff_var,first_quartile,median,third_quartile,percentile">
indices_for_classif_labels<-c("DEM",
                              "slope",
                              "accumulation",
                              "NDVI_SPOT6",
                              "NDWI_SPOT6",
                              "BI_SPOT6",
                              "B02_S2",
                              "B03_S2",
                              "B04_S2",
                              "B05_S2",
                              "B06_S2",
                              "B07_S2",
                              "B08_S2",
                              "B8A_S2",
                              "B11_S2",
                              "B12_S2",
                              "B0_SPOT6",
                              "B1_SPOT6",
                              "B2_SPOT6",
                              "B3_SPOT6",
                              "PAN_SPOT6",
                              "text_energy_5",
                              "text_entropy_5",
                              "text_correlation_5",
                              "text_inertia_5",
                              "text_haralickcorellation_5",
                              "text_mean_5",
                              "text_energy_9",
                              "text_entropy_9",
                              "text_correlation_9",
                              "text_inertia_9",
                              "text_haralickcorellation_9",
                              "text_mean_9",
                              "text_energy_17",
                              "text_entropy_17",
                              "text_correlation_17",
                              "text_inertia_17",
                              "text_haralickcorellation_17",
                              "text_mean_17",
                              "NDVI_S2",
                              "NDWI_S2",
                              "BRI_S2",
                              "MNDWI_S2",
                              "MNDVI_S2",
                              "RNDVI_S2"
                              ) # parameter indices_for_classif_paths sets the path to each band
indices_for_classif_paths<-c(file.path(path_to_processing_folder,"DEM_SRTM/processed_data/DEM.tif"),
                             file.path(path_to_processing_folder,"DEM_SRTM/processed_data/slope.tif"),
                             file.path(path_to_processing_folder,"DEM_SRTM/processed_data/accumulation.tif"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/NDVI.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/NDWI.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/BI.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B02.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B03.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B04.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B05.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B06.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B07.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B08.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B8A.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B11.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B12.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_3.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PAN.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_4.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_7.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_advanced_5_5_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_4.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_7.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_advanced_9_9_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_4.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_7.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_advanced_17_17_0.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/NDVI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/NDWI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/BRI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/MNDWI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/MNDVI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/RNDVI.tif")
                             
)
                             
### Parameters for step 8
column_names_lc_classes<-c("level_1_fr","level_2_fr","level_3_fr","level_4_fr") #<Names of the columns of land cover classes in the ground truth dataset. eg : c("type_1","type_2"). Typically type_1 is the most aggregated land cover, type_2 is a less aggregated classification, etc.>
column_name_lc_classification<-"level_1_fr" #<Name of the column of land cover class that we want to classify in the ground truth dataset. Column type must be character string>

########################################################################################################################
############ Prepare workflow ############
########################################################################################################################

### Provide useful functions

#all_functions<-list.files("/home/ptaconet/r_react/functions",full.names = T)
#lapply(all_functions, source)

## Function to convert digital numbers to Top of Atmosphere reflectance. Details: https://www.orfeo-toolbox.org/Applications/OpticalCalibration.html
convert_dn_to_toa_reflectance<-function(input_path,output_path,otbapplications_folder_path){
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OpticalCalibration")," -in ",input_path," -out ",gsub(".TIF","_temp.TIF",output_path)," float -level toa")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BandMathX")," -il ",gsub(".TIF","_temp.TIF",output_path)," -exp \"10000*im1\" -out ",output_path," uint16")
  system(otb_appli)
  file.remove(gsub(".TIF","_temp.TIF",output_path))
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}

## Function to extract a region of interest defined by a user. Details: https://www.orfeo-toolbox.org/Applications/ExtractROI.html 
extract_roi<-function(input_path,output_path,roi_vector_path,otbapplications_folder_path,dem_folder_path){
  otb_appli<-paste0(file.path(otbapplications_folder_path,"otbcli_ExtractROI")," -in ",input_path," -out ",output_path," uint16 -mode fit -mode.fit.vect ",roi_vector_path," -elev.dem ",dem_folder_path)
  system(otb_appli)
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}

## Function to mosaic multiple remote sensing images. Details: https://github.com/remicres/otb-mosaic
mosaic<-function(input_path,output_path,otbapplications_folder_path){
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(input_path, collapse = ' ')," -out ",output_path," uint16")
  system(otb_appli)
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}

## Function to mosaic and extract ROI
mosaic_and_extract_roi<-function(input_path_files,output_path_file,roi_vector_path,otbapplications_folder_path,dem_folder_path){
  
  cat("Mosaicing the tiles...")
  output_path_temp_file=gsub("\\.","_temporary.",output_path_file)
  res<-mosaic(input_path_files,output_path_temp_file,otbapplications_folder_path)
  if (res=="Done"){
    cat ("Mosaic OK")
  }
  
  cat("Extracting the ROI...")
  res<-extract_roi(output_path_temp_file,output_path_file,roi_vector_path,otbapplications_folder_path,dem_folder_path)
  if (res=="Done"){
    cat ("ROI Extraction OK")
  }
  
  file.remove(output_path_temp_file)
  
  if (file.exists(output_path_file)){
    return("Done")
  } else {
    return("Error")
  }
}

## Function to orthorectify images. Details: https://www.orfeo-toolbox.org/Applications/OrthoRectification.html
orthorectify_spot67<-function(input_path,output_path,otbapplications_folder_path,dem_folder_path){
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OrthoRectification")," -io.in \"",input_path,"?&skipcarto=true\""," -io.out ",output_path," uint16 -elev.dem ",dem_folder_path)
  system(otb_appli)
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}


### Call useful packages
require(raster)
require(rgrass7)
require(RSAGA)
require(rgdal)
require(sf)
require(fasterize)
require(dplyr)
require(reshape)
require(randomForest)
require(caret)

### Set working directory
setwd(path_to_processing_folder)

## Set paths of output folders / files
path_to_spot67_raw_folder<-file.path(path_to_processing_folder,"VHR_SPOT6/raw_data") # Path to the folder where the Spot6/7 products are stored. Within that folder, there must be 1 folder / product. Each of these folder contains 2 files: the Panchromatic and mutlispectral .tar.gz files.
# Step 1
path_to_dem_raw_folder=file.path(path_to_processing_folder,"DEM_SRTM/raw_data") # Path to the folder where the DEM raw data will be stored
# Step 2
path_to_spot67_preprocessed_folder=gsub("raw","processed",path_to_spot67_raw_folder) # Path to the folder where the datasets extracted from the Spot6/7 will be stored
path_to_spot67_preprocessed_ms=file.path(path_to_spot67_preprocessed_folder,"MS.TIF") # Path to the output pre-processed Spot6/7 multispectral image
path_to_spot67_preprocessed_pan=file.path(path_to_spot67_preprocessed_folder,"PAN.TIF") # Path to the output pre-processed Spot6/7 panchromatic image
path_to_spot67_preprocessed_pansharpen=file.path(path_to_spot67_preprocessed_folder,"PANSHARPEN.TIF") # Path to the output pre-processed Spot6/7 pansharpened image
# Step 3
path_to_sentinel2_raw_folder<-file.path(path_to_processing_folder,"HR_Sentinel2/raw_data") # Path to the folder where the Sentinel 2 raw data will be stored
# Step 4
path_to_dem_preprocessed_folder=file.path(path_to_processing_folder,"DEM_SRTM/processed_data") # Path to the folder where the DEM processed data and related data (slope, accumulation, etc.) will be stored
path_to_sentinel2_preprocessed_folder<-file.path(path_to_processing_folder,"HR_Sentinel2/processed_data") # Path to the folder where the Sentinel 2 processed data will be stored
dem_output_path_file<-file.path(path_to_dem_preprocessed_folder,"DEM.tif") # output DEM file path
# Step 5
path_to_simple_textural_indices<-file.path(path_to_spot67_preprocessed_folder,"HaralickTextures_simple.TIF") # Path to textural indice
path_to_advanced_textural_indices<-file.path(path_to_spot67_preprocessed_folder,"HaralickTextures_advanced.TIF")
path_to_output_accumulation_threshold<-file.path(path_to_dem_preprocessed_folder,"accumulation_treshold_vector.gpkg")
# Step 6
path_to_segmentation_folder<-file.path(path_to_processing_folder,"Segmentation") # Path to the folder where the outputs of the segmentation process will be stored
path_to_segmented_dataset<-file.path(path_to_segmentation_folder,"segmentation_vector.gpkg")
# Step 7
path_to_ground_truth_stats<-file.path(path_to_groundtruth_folder,"ground_truth_stats.gpkg") # Path to the ground truth datasets with zonal statistics
path_to_zonalstatistics_folder<-file.path(path_to_processing_folder,"Zonal_statistics") # Path to the folder where the output of the zonal statistics process will be stored
# Step 8
path_to_segmented_dataset_stats<-file.path(path_to_segmentation_folder,"segmented_dataset_stats.gpkg") # Path to the object segmented datasets with zonal statistics

# Step 9
path_to_classification_folder<-file.path(path_to_processing_folder,"Classification") # Path to the fodler of output of the classification (including model)
path_to_output_classification_vector<-file.path(path_to_classification_folder,"classification_vector.gpkg")
# Step 10
path_to_output_classification_raster<-file.path(path_to_classification_folder,"classification_raster.tif")

## Create directories
directories<-list(path_to_dem_raw_folder,path_to_spot67_preprocessed_folder,path_to_sentinel2_raw_folder,path_to_dem_preprocessed_folder,path_to_sentinel2_preprocessed_folder)
lapply(directories, dir.create)

# Set GRASS environment and database location 
loc <- rgrass7::initGRASS(path_to_grassApplications_folder, home=getwd(), gisDbase="GRASS_TEMP", override=TRUE,mapset = "PERMANENT" )
execGRASS("g.proj",flags="c",parameters = list(proj4=proj_srs))

# Set saga environment
saga_work_env <- RSAGA::rsaga.env()

########################################################################################################################
############ Start Workflow ############
########################################################################################################################

########################################################################################################################
########################################################################################################################
############ Step 1 - Downloading the SRTM tiles  ############
########################################################################################################################
########################################################################################################################
### Uses otb applications: DownloadSRTMTiles

cat("Downloading the DEM SRTM tiles corresponding to the ROI...")
system(paste0(file.path(path_to_otbApplications_folder,"otbcli_DownloadSRTMTiles")," -vl ",file.path(path_to_processing_folder,path_to_roi_vector)," -mode download -tiledir ",file.path(path_to_processing_folder,path_to_dem_raw_folder)))
cat("OK")

########################################################################################################################
########################################################################################################################
############ Step 2 - Pre-processing the Spot6/7 products ############
########################################################################################################################
########################################################################################################################
### Uses otb applications: OpticalCalibration, BandMathX, ExtractROI, Mosaic, OrthoRectification, Superimpose, BundleToPerfectSensor

cat("Pre-processing the Spot6/7 products ")

## List the folders available in the folder (1 folder / couple of PAN + MS image)
products_to_preprocess<-dir(path_to_spot67_raw_folder,full.names = TRUE)

# Preprocess each product
for (i in 1:length(products_to_preprocess)){
  cat(paste0("Starting preprocessing of product ",products_to_preprocess[i]))
  
  spot67_preprocessing_output_folder_path=file.path(path_to_spot67_preprocessed_folder,dir(path_to_spot67_raw_folder,full.names = FALSE)[i])
  dir.create(spot67_preprocessing_output_folder_path)
  # Unzip both MS and PAN data
  products_to_unzip=list.files(products_to_preprocess[i],full.names = T)
  for (j in 1:length(products_to_unzip)){
  untar(products_to_unzip[j],exdir=products_to_preprocess[i])
  }
  # Identifiy MS and PAN folders
  MS_folder=dir(products_to_preprocess[i],pattern = "MS",full.names = T)
  MS_folder=MS_folder[!grepl('.tar', MS_folder)]
  PAN_folder=dir(products_to_preprocess[i],pattern = "PAN",full.names = T)
  PAN_folder=PAN_folder[!grepl('.tar', PAN_folder)]
  
  ##############################################################
  #### 2.1 - fusionning the tiles of the panchromatic image ####
  ##############################################################
  
  ## Identify folder containing the PAN tiles
  folders=list.dirs(path = PAN_folder, full.names = TRUE)
  PAN_tile_folder=folders[grepl('IMG', folders)][1]
  pan_tifs_paths=as.vector(list.files(path=PAN_tile_folder,pattern=".TIF",full.names=TRUE)) 
  ## Mosaic
  cat("Starting PAN tile fusionning ...")
  PAN_mosaic_path<-file.path(spot67_preprocessing_output_folder_path,"PAN.TIF")
  res<-mosaic(pan_tifs_paths,PAN_mosaic_path,path_to_otbApplications_folder)
  cat(res)
  
  ##############################################################
  #### 2.2 - converting multispectral and panchromatic images from digital numbers to TOA reflectance ####
  ##############################################################
  
  cat("Starting Convertion to TOA reflectance ...")
  ## PAN 
  PAN_toa_path<-file.path(spot67_preprocessing_output_folder_path,"PAN_TOA.TIF")
  res<-convert_dn_to_toa_reflectance(PAN_mosaic_path,PAN_toa_path,path_to_otbApplications_folder)
  cat(res)
  ## MS
  # Identify folder containing the MS tile
  folders=list.dirs(path = MS_folder, full.names = TRUE)
  MS_tile_folder=folders[grepl('IMG', folders)][1]
  ms_tif_path=as.vector(list.files(path=MS_tile_folder,pattern=".TIF",full.names=TRUE)) 
  MS_toa_path<-file.path(spot67_preprocessing_output_folder_path,"MS_TOA.TIF")
  res<-convert_dn_to_toa_reflectance(ms_tif_path,MS_toa_path,path_to_otbApplications_folder)
  cat(res)
  
  ##############################################################
  #### 2.3 - orthorectifying multispectral and panchromatic images ####
  ##############################################################
  
  cat("Starting orthorectification ...")
  ## PAN 
  PAN_orthorectified_path=file.path(spot67_preprocessing_output_folder_path,"PAN_ORTHO.TIF")
  res<-orthorectify_spot67(PAN_toa_path,PAN_orthorectified_path,path_to_otbApplications_folder,path_to_dem_raw_folder)
  cat(res)
  ## MS 
  MS_orthorectified_path=file.path(spot67_preprocessing_output_folder_path,"MS_ORTHO.TIF")
  res<-orthorectify_spot67(MS_toa_path,MS_orthorectified_path,path_to_otbApplications_folder,path_to_dem_raw_folder)
  cat(res)
  
  ##############################################################
  #### 2.4 - extracting the ROI ####
  ##############################################################
  
  cat("Extracting the ROI ...")
  ## PAN 
  PAN_roi_path=file.path(spot67_preprocessing_output_folder_path,"PAN.TIF")
  res<-extract_roi(PAN_orthorectified_path,PAN_roi_path,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  cat(res)
  ## MS
  MS_roi_path=file.path(spot67_preprocessing_output_folder_path,"MS.TIF")
  res<-extract_roi(MS_orthorectified_path,MS_roi_path,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  cat(res)
  
  ##############################################################
  #### 2.5 - pansharpening ####
  ##############################################################
  
  PANSHARPEN_roi_path<-file.path(spot67_preprocessing_output_folder_path,"PANSHARPEN.TIF")
  
    otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Superimpose")," -inr ",PAN_roi_path," -inm ",MS_roi_path," -out ",gsub(".TIF","_temp.TIF",PANSHARPEN_roi_path)," uint16")
    system(otb_appli)
    otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BundleToPerfectSensor")," -inp ",PAN_roi_path," -inxs ",gsub(".TIF","_temp.TIF",PANSHARPEN_roi_path)," -out ",PANSHARPEN_roi_path," uint16")
    system(otb_appli)
    file.remove(gsub(".TIF","_temp.TIF",PANSHARPEN_roi_path))
    file.remove(gsub(".TIF",".geom",gsub(".TIF","_temp.TIF",PANSHARPEN_roi_path)))
    
  ### Remove temporary files and folders
  file.remove(setdiff(list.files(spot67_preprocessing_output_folder_path,full.names = T),c(PAN_roi_path,MS_roi_path,PANSHARPEN_roi_path,gsub(".TIF",".geom",PAN_roi_path),gsub(".TIF",".geom",MS_roi_path),gsub(".TIF",".geom",PANSHARPEN_roi_path))))
  folder_to_remove=list.files(products_to_preprocess[i],full.names = T)[!grepl('.tar', list.files(products_to_preprocess[i],full.names = T))]
  for (i in 1:length(folder_to_remove)){
    system(paste0("rm -r ", folder_to_remove[i]))  
    }
}

##############################################################
#### 2.6 - mosaicing the various tiles if relevant ####
##############################################################

if(length(products_to_preprocess)>1){
  #dir.create(file.path(path_to_spot67_preprocessed_folder,"mosaic"))
  all_preprocessed_data<-list.files(path_to_spot67_preprocessed_folder,recursive = T,full.names = T)
  
  MS_to_mosaic_paths=rev(all_preprocessed_data[grepl('MS.TIF', all_preprocessed_data)])
  PAN_to_mosaic_paths=rev(all_preprocessed_data[grepl('PAN.TIF', all_preprocessed_data)])
  PANSHARPEN_to_mosaic_paths=rev(all_preprocessed_data[grepl('PANSHARPEN.TIF', all_preprocessed_data)])
  
  # We perform a very simple mosaic, by simply copying the last image over earlier ones in areas of overlap. We could perform much more advanced mosaicing operations with the otbcli_Mosaic application (for additional details: https://github.com/remicres/otb-mosaic)
  cat("Mosaicing the multiple images ...")
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(MS_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",path_to_spot67_preprocessed_ms," uint16")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(PAN_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",path_to_spot67_preprocessed_pan," uint16")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(PANSHARPEN_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",path_to_spot67_preprocessed_pansharpen," uint16")
  system(otb_appli)
  cat("Done")
} else {
  file.copy(list.files(spot67_preprocessing_output_folder_path,full.names = T),path_to_spot67_preprocessed_folder)
  system(paste0("rm -r ", spot67_preprocessing_output_folder_path)) 
}

cat("Pre-processing the Spot6/7 products OK")


########################################################################################################################
########################################################################################################################
############ Step 3 - Downloading ancillary data ############
########################################################################################################################
########################################################################################################################
### Uses Copernicus Scihub API

##############################################################
#### 3.1 - Download Sentinel 2 data ####
##############################################################

cat("Downloading the ancillary data: Sentinel 2 product(s) ...")
### NOT WORK
for (i in 1:length(Sentinel2_products_id)){
  # from: https://scihub.copernicus.eu/userguide/BatchScripting
  # wget --content-disposition --continue --user={USERNAME} --password={PASSWORD}"https://scihub.copernicus.eu/dhus/odata/v1/Products('22e7af63-07ad-4076-8541-f6655388dc5e')/\$value"
  system(paste0("wget --content-disposition --continue --user=",copernicus_scihub_username," --password=",copernicus_scihub_password," --directory-prefix ",path_to_sentinel2_raw_folder," \"https://scihub.copernicus.eu/dhus/odata/v1/Products('",Sentinel2_products_id[i],"')/\\$value\\"))
  }
cat("Downloading the ancillary data: Sentinel 2 product(s) OK")


########################################################################################################################
########################################################################################################################
############ Step 4 - Preprocessing the ancillary data ############
########################################################################################################################
########################################################################################################################
### Uses otb applications: Mosaic, ExtractROI

##############################################################
#### 4.1 - preprocessing the DEM : mosaicing the various tiles if relevant, and then extracting the ROI ####
##############################################################

cat("Pre-processing the DEM : if relevant, mosaicing the various tiles and extracting the ROI ...")

# List the products
products_to_preprocess<-list.files(path_to_dem_raw_folder,full.names = T)
# Unzip the files
for (i in 1:length(products_to_preprocess)){
unzip(products_to_preprocess[i],exdir = path_to_dem_raw_folder)
}
file.remove(products_to_preprocess)
products_to_preprocess<-list.files(path_to_dem_raw_folder,full.names = T)
# If there are multiple tiles, mosaic them and then extract the ROI, else only extract the ROI
if(length(products_to_preprocess)>1){
  res<-mosaic_and_extract_roi(products_to_preprocess,dem_output_path_file,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  } else {
  res<-extract_roi(products_to_preprocess,dem_output_path_file,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  }
cat(res)

# Convert from EPSG 4326 (default SRTM EPSG) to UTM EPSG
system(paste0("gdalwarp -t_srs '",proj_srs,"' -overwrite ",dem_output_path_file," ",gsub("DEM.TIF","DEM_temp.TIF",dem_output_path_file)))

file.remove(dem_output_path_file)
file.rename(gsub("DEM.TIF","DEM_temp.TIF",dem_output_path_file),dem_output_path_file)

##############################################################
#### 4.2 - preprocessing the Sentinel 2 product(s) : mosaicing the various tiles if relevant, and then extracting the ROI ####
##############################################################

cat("Pre-processing the S2 images : if relevant, mosaicing the various bands and extracting the ROI ...")

## List the products
products_to_preprocess<-dir(path_to_sentinel2_raw_folder,full.names = TRUE)

# Unzip all products
for (i in 1:length(products_to_preprocess)){
  unzip(products_to_preprocess[i],exdir=path_to_sentinel2_preprocessed_folder)
}

# If there are multiple products, mosaic them and then extract the ROI, else only extract the ROI
patterns=c("B01.jp2","B02.jp2","B03.jp2","B04.jp2","B05.jp2","B06.jp2","B07.jp2","B08.jp2","B8A.jp2","B09.jp2","B10.jp2","B11.jp2","B12.jp2")
  
for (i in 1:length(patterns)){
  paths_images_to_mosaic=list.files(path_to_sentinel2_preprocessed_folder,pattern = patterns[i],full.names = T,recursive = T)
  path_to_output_mosaiced_image=file.path(path_to_sentinel2_preprocessed_folder,gsub(".jp2",".TIF",patterns[i]))
  
  if(length(products_to_preprocess)>1){
  res<-mosaic_and_extract_roi(paths_images_to_mosaic,path_to_output_mosaiced_image,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  } else {
  res<-extract_roi(paths_images_to_mosaic,path_to_output_mosaiced_image,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  }
  cat(res)
   
}
  
########################################################################################################################
########################################################################################################################
############ Step 5 - Preparing the ancillary data for the classification ############
########################################################################################################################
########################################################################################################################
### Uses otb applications: HaralickTextureExtraction (or SelectiveHaralickTextures), RadiometricIndices, ConcatenateImages, SplitImages
### Uses grass applications: r.slope.aspect, r.terraflow, r.out.gdal
### Uses gdal application: gdal_translate

##############################################################
#### 5.1 - extract indices from the DEM : slope, aspect, flow accumulation, flow direction, topographic convergence index ####
##############################################################
## We use GRASS, calling it in R using the "rgrass7" package. We use two GRASS applications: r.slope.aspect and r.terraflow

# Set output paths
slope_output_path<-file.path(path_to_dem_preprocessed_folder,"slope.tif")
aspect_output_path<-file.path(path_to_dem_preprocessed_folder,"aspect.tif")
accumulation_output_path<-file.path(path_to_dem_preprocessed_folder,"accumulation.tif")
direction_output_path<-file.path(path_to_dem_preprocessed_folder,"direction.tif")
tci_output_path<-file.path(path_to_dem_preprocessed_folder,"tci.tif") 

# Import DEM to GRASS and set region
execGRASS("r.external", flags="o", parameters=list(input=dem_output_path_file, output="tmprast",band=1))
execGRASS("g.region", parameters=list(raster="tmprast")) 

# Compute slope and aspect and save to disk
execGRASS("r.slope.aspect", flags="overwrite", parameters=list(elevation="tmprast", slope="slope",aspect="aspect",format="percent", precision="FCELL",zscale=1,min_slope=0))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="slope", output=slope_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="aspect", output=aspect_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))

# Compute hydrograpy indices and save to disk
execGRASS("r.terraflow", flags="overwrite", parameters=list(elevation="tmprast", direction="direction",accumulation="accumulation",tci="tci"))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="direction", output=direction_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="accumulation", output=accumulation_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="tci", output=tci_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))

# Create accumulation vector file given the threshold. This dataset will be used afterwards for the computation of the zonal statistics (distance to hydographic network)
acc_raster<-raster(accumulation_output_path)
acc_raster[which(values(acc_raster)<threshold_accumulation_raster)]=0
acc_raster[which(values(acc_raster)>=threshold_accumulation_raster)]=1
accumulation_threshold_output_path<-gsub(".TIF","_treshold.tif",accumulation_output_path)
writeRaster(acc_raster,accumulation_threshold_output_path,overwrite=TRUE)
output_path=file.path(path_to_dem_preprocessed_folder,"accumulation_treshold.gpkg")
gdal_appli<-paste0("gdal_polygonize.py ",accumulation_threshold_output_path," ",output_path," -b 1 None DN")
system(gdal_appli)
acc_vect<-readOGR(output_path,stringsAsFactors = F)
acc_vect<-acc_vect[which(acc_vect$DN==1),]
writeOGR(acc_vect,path_to_output_accumulation_threshold,driver = "GPKG",layer="accumulation_vector")
file.remove(output_path)


# We disaggregate the resolution of the rasters because some small objets send back NA values with a 30 m resolution DEM
rasters_to_disaggregate<-c(dem_output_path_file,slope_output_path,aspect_output_path,accumulation_output_path,direction_output_path,tci_output_path)
cat("Disaggregating the resolutions of the rasters")
for (i in 1:length(rasters_to_disaggregate)){
rast<-raster(rasters_to_disaggregate[i])
rast<-disaggregate(rast,fact=3)
writeRaster(rast,rasters_to_disaggregate[i],overwrite=TRUE)
}

# To read and plot the output raster 
#r <- readRAST("slope")
#spplot(r)

##############################################################
#### 5.2 - extract textural indices from the Spot6/7 image ####
##############################################################

# Get image maximum and minimum
rast<-raster(path_to_spot67_preprocessed_pan)
min<-as.numeric(cellStats(rast,min))
max<-as.numeric(cellStats(rast,max))

## Compute textures at various window sizes (5×5, 9 × 9, 17 × 17, 21 × 21, and 35 × 35)
# Using the application available in the official release (does not enable to select the set of textures)
for (i in 1:length(xrad)){
  
path_to_simple_texture<-gsub(".TIF",paste0("_",xrad[i],"_",xrad[i],".TIF"),path_to_simple_textural_indices)
path_to_advanced_texture<-gsub(".TIF",paste0("_",xrad[i],"_",xrad[i],".TIF"),path_to_advanced_textural_indices)
  
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_HaralickTextureExtraction")," -in ",path_to_spot67_preprocessed_pan," -parameters.xrad ",xrad[i]," -parameters.yrad ",yrad[i]," -parameters.nbbin ",nbbin," -parameters.min ",min," -parameters.max ",max," -texture simple -out ",path_to_simple_texture)
system(otb_appli)

## Split the textures into n bands (1 / texture)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SplitImage")," -in ",path_to_simple_texture," -out ",path_to_simple_texture)
system(otb_appli)

## Remove useless textures (textures that we will not use for the classification)
# We keep : energy (HaralickTextures_simple_0), entropy (HaralickTextures_simple_1), correlation (HaralickTextures_simple_2), inertia (HaralickTextures_simple_4 - distinction sols nu / bati), Haralick correlation (HaralickTextures_simple_7)
file.remove(c(path_to_simple_texture,
          gsub(".TIF","_3.TIF",path_to_simple_texture),
          gsub(".TIF","_5.TIF",path_to_simple_texture),
          gsub(".TIF","_6.TIF",path_to_simple_texture)
          ))

## Compute advanced textures
# Using the application available in the official release (does not enable to select the set of textures)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_HaralickTextureExtraction")," -in ",path_to_spot67_preprocessed_pan," -parameters.xrad ",xrad[i]," -parameters.yrad ",yrad[i]," -parameters.nbbin ",nbbin," -parameters.min ",min," -parameters.max ",max," -texture advanced -out ",path_to_advanced_texture)
system(otb_appli)

## Split the textures into n bands (1 / texture)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SplitImage")," -in ",path_to_advanced_texture," -out ",path_to_advanced_texture)
system(otb_appli)

## Remove useless textures (textures that we will not use for the classification)
# We keep : mean (HaralickTextures_advanced_0)
file.remove(c(path_to_advanced_texture,
              gsub(".TIF","_1.TIF",path_to_advanced_texture),
              gsub(".TIF","_2.TIF",path_to_advanced_texture),
              gsub(".TIF","_3.TIF",path_to_advanced_texture),
              gsub(".TIF","_4.TIF",path_to_advanced_texture),
              gsub(".TIF","_5.TIF",path_to_advanced_texture),
              gsub(".TIF","_6.TIF",path_to_advanced_texture),
              gsub(".TIF","_7.TIF",path_to_advanced_texture),
              gsub(".TIF","_8.TIF",path_to_advanced_texture),
              gsub(".TIF","_9.TIF",path_to_advanced_texture)
))

}

# Using the application available in the personal release (enables to select the set of textures)
# indices_list<-c("enthropy","energy")
#otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SelectiveHaralickTextures")," -in ",file.path(path_to_processing_folder,PAN_roi_path)," -parameters.xrad 11 -parameters.yrad 11 -parameters.nbbin 32 -parameters.min ",minValue(rast)," -parameters.max ",maxValue(rast)," -indices ",paste(texture_indices_list, collapse = " ")," -out ",path_to_simple_textural_indices)
#system(otb_appli)

##############################################################
#### 5.3 - extract radiometric indices from the Spot6/7 image  ####
##############################################################

# The list of available radiometric indices automatically calculated with OTB can be found here: https://www.orfeo-toolbox.org/Applications/RadiometricIndices.html
for (i in 1:length(radiometric_indices_list_spot67)){
  
  indice<-sub(".*:","",radiometric_indices_list_spot67[i])
  path_to_output_indice<-file.path(path_to_spot67_preprocessed_folder,paste0(indice,".TIF"))
  
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_RadiometricIndices")," -in ",path_to_spot67_preprocessed_pansharpen," -channels.blue 1 -channels.green 2 -channels.red 3 -channels.nir 4 -list ",radiometric_indices_list_spot67[i]," -out ",path_to_output_indice)
  system(otb_appli)
}


##############################################################
#### 5.4 - extract radiometric indices from the Sentinel 2 image ####
##############################################################

b03<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B03.TIF"))
b04<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B04.TIF"))
b05<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B05.TIF"))
b06<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B06.TIF"))
b07<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B07.TIF"))
b08<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B08.TIF"))
b08A<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B8A.TIF"))
b11<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B11.TIF"))
b12<-raster(file.path(path_to_sentinel2_preprocessed_folder,"B12.TIF"))

# disagregate the resolutions of the 20m bands to 10m to be able to compute the indices. disagregate does not change the cells values.
b05<-disaggregate(b05,fact=2)
b06<-disaggregate(b06,fact=2)
b07<-disaggregate(b07,fact=2)
b08A<-disaggregate(b08A,fact=2)
b11<-disaggregate(b11,fact=2)
b12<-disaggregate(b12,fact=2)

# Compute NDVI
ndvi<-(b08-b04)/(b08+b04)
writeRaster(ndvi,file.path(path_to_sentinel2_preprocessed_folder,"NDVI.TIF"))

# Compute MNDVI
mndvi<-(b08-b11)/(b08+b11)
writeRaster(mndvi,file.path(path_to_sentinel2_preprocessed_folder,"MNDVI.TIF"))

# Compute RNDVI
rndvi<-(b08-b06)/(b08+b06)
writeRaster(rndvi,file.path(path_to_sentinel2_preprocessed_folder,"RNDVI.TIF"))

# Compute BRI
bri<-sqrt(b03^2+b04^2+b05^2+b06^2+b07^2+b08^2+b08A^2+b11^2+b12^2)
writeRaster(bri,file.path(path_to_sentinel2_preprocessed_folder,"BRI.TIF"))

# Compute NDWI
ndnwi<-(b03-b08)/(b03+b08)
writeRaster(ndnwi,file.path(path_to_sentinel2_preprocessed_folder,"NDWI.TIF"))

# Compute MNDWI
mndnwi<-(b03-b11)/(b03+b11)
writeRaster(mndnwi,file.path(path_to_sentinel2_preprocessed_folder,"MNDWI.TIF"))


##############################################################
#### 5.5 - Split the bands of the pansharpened Spot6/7 image ####
##############################################################

otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SplitImage")," -in ",path_to_spot67_preprocessed_pansharpen," -out ",path_to_spot67_preprocessed_pansharpen," uint16")
system(otb_appli)

########################################################################################################################
########################################################################################################################
############ Step 6 - Segmenting the Spot6/7 image ############
## IMPORTANT NOTE : USES AN OTB APPLICATION ONLY AVAILABLE IN THE PERSONAL RELEASE OF OTB 
########################################################################################################################
########################################################################################################################
## TODO test with https://github.com/RTOBIA/LSOBIA
### Uses otb applications: otbcli_LSGRM
### Uses gdal application: gdal_polygonize.py

cat("Segmenting the Spot6/7 image...")

if (length(products_to_preprocess)==1){
  path_to_image_to_segment<-PANSHARPEN_roi_path
} else {
  path_to_image_to_segment<-PANSHARPEN_mosaiced_path
}

cat("Starting segmentation...")

# Segment
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_LSGRM")," -in ",path_to_image_to_segment," -out ",gsub(".gpkg",".TIF",path_to_segmented_dataset)," int32 -threshold ",segmentation_threshold," -criterion.bs.cw ",segmentation_cw," -criterion.bs.sw ",segmentation_sw)
system(otb_appli)

# Output of the segmentation is a raster. Vectorize
#path_to_output_segmentation_vector=file.path(path_to_outputFiles_segmentation_folder,paste0(outputImageName,"_segmented_final.gpkg"))
gdal_appli<-paste0("gdal_polygonize.py ",gsub(".gpkg",".TIF",path_to_segmented_dataset)," ",path_to_segmented_dataset," -b 1 None DN")
system(gdal_appli)

cat("Segmenting the Spot6/7 OK")

########################################################################################################################
########################################################################################################################
############ Step 7 - Extract zonal statistics for each band / indice   ############
########################################################################################################################
########################################################################################################################
### Uses saga application: rsaga.import.gdal, rsaga.geoprocessor

# More info on the use of RSAGA here: https://www.r-bloggers.com/rsaga-getting-started/ and here: https://cran.r-project.org/web/packages/RSAGA/vignettes/RSAGA.html
# Saga must be installed in the computer to use RSAGA. To install on Unbuntu : sudo apt install saga

# raster data must be converted to .sgrd format to be used in saga
for (i in 1:length(indices_for_classif_paths)){
  cat(paste0("converting raster ",indices_for_classif_paths[i],"\n"))
  rsaga.import.gdal(indices_for_classif_paths[i])
}

### Function to compute the zonal statistics. It will be used to compute the zonal statistics of the ground truth DB + the segmented objects
function_compute_zonal_statistics<-function(path_to_input_gpkg_file,path_to_output_gpkg_file,indices_for_classif_paths,indices_for_classif_labels){

# First save an shp version of the geopackage datasets (ground truth and segmentation) since SAGA does not deal with geopackages
path_to_stats_shp<-gsub(".gpkg",".shp",path_to_output_gpkg_file)
system(paste0("ogr2ogr ",path_to_stats_shp," ",path_to_input_gpkg_file))


# Compute zonal statistics using the Saga Grid Statistics for Polygons function
# to see the function parameters : rsaga.get.usage('shapes_grid', 2) 
for (i in 1:length(indices_for_classif_paths)){
rsaga.geoprocessor('shapes_grid', module = 2, env = saga_work_env, param = list(
  PARALLELIZED = 'true',
  GRIDS = paste(gsub("TIF|tif","sgrd",indices_for_classif_paths[i]),collapse=';'), # Note that in theory we could integrate multiple raster at the same time, e.g. GRIDS = "/home/ptaconet/Documents/react/data_BF/VHR_SPOT6/processed_data/NDWI.sgrd;/home/ptaconet/Documents/react/data_BF/VHR_SPOT6/processed_data/NDVI.sgrd". However grids must have the same characteristics (extent)
  POLYGONS = path_to_stats_shp,
  RESULT = path_to_stats_shp,
  COUNT='false',
  MIN='false',
  MAX='false',
  RANGE='false',
  SUM='false',
  MEAN='true',
  VAR='false',
  STDDEV='true'))
}

## Compute shape statistics with Saga: rsaga.get.usage('shapes_polygons', 7) 
rsaga.geoprocessor('shapes_polygons', module = 7, env = saga_work_env, param = list(SHAPES = path_to_stats_shp))


# Open and rename columns (names have been cut since shp does not accept long names)
input_data<-sf::st_read(path_to_input_gpkg_file)
output_stats<-sf::st_read(path_to_stats_shp)
cols_to_rename<-setdiff(colnames(output_stats),c("geometry",colnames(input_data)))
new_colnames<-NULL
for (i in 1:length(indices_for_classif_labels)){
  new_colnames<-c(new_colnames,paste0(indices_for_classif_labels[i],"_avg"))
  new_colnames<-c(new_colnames,paste0(indices_for_classif_labels[i],"_stddev"))
}

new_colnames<-c(new_colnames,"shape_area","shape_perimeter","shape_interior_edge_ratio","shape_p_sqrt_a","shape_max_dist","shape_d_over_a","shape_d_sqrt_a","shape_index")

if (length(new_colnames)!=length(cols_to_rename)){
  cat("error")
} else {
  for (i in 1:length(new_colnames)){
  colnames(output_stats)[which(colnames(output_stats)==cols_to_rename[i])]<-new_colnames[i]
  }
}

# Write data as geopackage
sf::st_write(output_stats,path_to_output_gpkg_file,layer_options = "OVERWRITE=true")


## Compute shape statistics using the Schwartzberg and the reock indexes. It uses a function extracted from https://github.com/gerrymandr/compactr/blob/master/compactness.R
# Source the functions with : source("https://raw.githubusercontent.com/gerrymandr/compactr/master/compactness.R")
# We copied them here for offline processing

#' @rdname shape_index
polsby_popper = function(poly1) {
  require(sf, quietly = TRUE)
  require(units, quietly = TRUE)
  return(drop_units(4 * pi * st_area(poly1) / st_length(st_boundary(poly1))^2))
}
#' @rdname shape_index
schwartzberg = function(poly1) {
  return(polsby_popper(poly1)^-0.5)
}

#' @rdname area_compactness
reock = function(poly1, mbc = NULL) {
  require(sf, quietly = TRUE)
  require(units, quietly = TRUE)
  if (is.null(mbc)) {
    require(lwgeom, quietly = TRUE)
    mbc = st_minimum_bounding_circle(st_convex_hull(st_geometry(poly1)))
  }
  
  return(drop_units(st_area(poly1) / st_area(mbc)))
}

poly<-sf::st_read(path_to_output_gpkg_file)
#poly$shape_schwartzberg<-schwartzberg(poly) # equal to shape_index of Saga shapes_polygons index
poly$shape_reock<-reock(poly)
sf::st_write(poly,path_to_output_gpkg_file,layer_options = "OVERWRITE=true")

## Compute distance statistics (distance to hydrographic network)
# Distance to hydrographic network (using the accumulation dataset derived from the DEM)
execGRASS("v.in.ogr", flags=c("o","overwrite"), parameters=list(input=path_to_output_accumulation_threshold, output="accumulation",min_area=0.0001, snap=-1.0))
execGRASS("v.in.ogr", flags=c("o","overwrite"), parameters=list(input=path_to_output_gpkg_file, output="ground_truth",min_area=0.0001, snap=-1.0))
execGRASS("v.db.addcolumn", parameters=list(map="ground_truth", columns="dist_to_hydro double"))
execGRASS("v.distance", flags=c("overwrite"), parameters=list(from="ground_truth", from_type="point,line,area", to="accumulation",to_type="point,line,area",dmax=-1,dmin=-1,upload="dist",column="dist_to_hydro",output="gt_stats_updated"))

# Save file as geopackage and csv
writeOGR(readVECT("ground_truth"),path_to_output_gpkg_file,driver = "GPKG",layer="ground_truth_stats",overwrite_layer = TRUE)
write.csv(as.data.frame(readVECT("ground_truth")),gsub(".gpkg",".csv",path_to_output_gpkg_file))

if (file.exists(path_to_output_gpkg_file)){
  return("Done")
} else {
  return("Error")
}

}

##############################################################
#### 7.1 - Extract zonal statistics for the ground truth dataset ####
##############################################################

res<-function_compute_zonal_statistics(path_to_ground_truth_data,path_to_ground_truth_stats,indices_for_classif_paths,indices_for_classif_labels)

##############################################################
#### 7.2 - Extract zonal statistics for the segmented objects dataset ####
##############################################################

res<-function_compute_zonal_statistics(path_to_segmented_dataset,path_to_segmented_dataset_stats,indices_for_classif_paths,indices_for_classif_labels)


## Remove saga files
file.remove(list.files(path_to_processing_folder,pattern = ".mgrd|.sdat|.sgrd", recursive = TRUE))

########################################################################################################################
########################################################################################################################
############ Step 8 - Train classifier  ############
########################################################################################################################
########################################################################################################################

####### Prepare datasets for the classification

# To ensure that simulations or random objects can be reproduced (more info: http://rfunction.com/archives/62). It must be called each time a new simulation or random object is run
set.seed(1)

# Read the ground truth dataset with zonal statistics 
ground_truth_df<-read.csv(gsub(".gpkg",".csv",path_to_ground_truth_stats))
# Remove useless columns
ground_truth_df$X=NULL
ground_truth_df$cat=NULL

## Get primitives (features) paths for further calculation of the zonal statistics on the segmented objects 
column_names_primitives<-setdiff(colnames(ground_truth_df),column_names_lc_classes)
df_primitives_types_sources<-as.data.frame(column_names_primitives,stringsAsFactors=FALSE)
pattern<-strsplit(methods_to_compute, split=',')[[1]]
pattern<-paste(pattern,collapse = '|_')
pattern<-paste0("_",pattern)
df_primitives_types_sources$indice<-gsub(pattern,"",df_primitives_types_sources$column_names_primitives)
df_primitives_types_sources<-left_join(df_primitives_types_sources,data.frame(indices_for_classif_labels,indices_for_classif_paths,stringsAsFactors = F),by=c("indice"="indices_for_classif_labels"))

# Keep only useful columns in the GT dataset (i.e. column to classify + primitives) and rename column to classify (response)
ground_truth_df_model<-ground_truth_df[,c(column_name_lc_classification,column_names_primitives)]
colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==column_name_lc_classification)]<-"response"

# Recursive Feature Elimination does not accept NAs in the predictors. Fill the NAs in the ground truth dataset using the na.roughfix function. We do it this way since there are not many NAs in our datasets. Note: we could also use the randomForest::rfImpute function
ground_truth_df_model <- randomForest::na.roughfix(ground_truth_df_model)

##############################################################
#### 8.1 - Parameterize the model ####
##############################################################

# get_optimum_rf_model: Function to get the optimum parameters for the model.
  # - Selects the optimum features to keep in the model (by running the model with the default parameters). Feature selection uses the Recursive Feature Elimination (RFE) of the caret package. Useful information about the RFE can be found here: http://topepo.github.io/caret/recursive-feature-elimination.html . Examples of implentations of the RFE can be found here: https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/ and here: https://www.r-bloggers.com/feature-selection-using-the-caret-package/
  # - Selects the optimum number of trees (ntrees)
  # - Selects the optimum number of variables to test at each division (mtry)
# Input: 
  # - df_model: data.frame whose first column is the response vector (named "response") and the other columns are the classifiers (no NA allowed in the classifiers)
# Outputs: 
  # - [[1]] : optimum_primitives : vector of optimum primitives
  # - [[2]] : optimum_ntrees : integer optimum ntree
  # - [[3]] : optimum_mtry : integer optimum mtry
 
get_optimum_rf_model<-function(model_df){
  
  cat("\nGetting the optimum primitives (x) (might be long)...\n")
  
  ## Define the control using a random forest selection function (functions=rfFuncs). 
  # number: the number of folds to consider. we select the default values which is 10 folds
  control <- caret::rfeControl(functions=rfFuncs, method="cv", number=10, returnResamp = "final")
  ## Run the Recursive Feature Elimination algorithm
  # sizes = a vector of integers corresponding to the number of features that should be retained in the updated model. Here we choose the number of primitives divided by 30
  subsets <- seq(1,ncol(model_df)-1,round((ncol(model_df)-1)/30))
  results <- caret::rfe(model_df[,2:ncol(model_df)], model_df[,1], sizes=subsets, rfeControl=control)
  # to summarize the results : print(results) ; to list the chosen features : predictors(results) ; to plot the results : plot(results, type=c("g", "o"))
  
  optimum_primitives<-predictors(results)
  
  ground_truth_df_features_optimum_primitives<-model_df[,c("response",optimum_primitives)]
  # Get optimum ntrees
  cat("\nGetting the optimum number of trees, setting max to 500 (ntree)\n...")
  model<-randomForest::randomForest(response ~ ., data=ground_truth_df_features_optimum_primitives, ntree=500,  na.action = na.omit)
  # to print error OOB in function of the number of trees: plot(model$err.rate[, 1], type = "l", xlab = "nombre d'arbres", ylab = "erreur OOB")
  optimum_ntrees<-which(model$err.rate[, 1]==min(model$err.rate[, 1]))[1]
  # Get optimum mtry
  cat("\nGetting the optimum number of variables randomly sampled as candidates at each split (mtry)...")
  model_tuned<-randomForest::tuneRF(x=ground_truth_df_features_optimum_primitives[,predictors(results)],y=ground_truth_df_features_optimum_primitives$response,ntree=optimum_ntrees)
  optimum_mtry<-as.numeric(model_tuned[,1][which(model_tuned[,2]==min(model_tuned[,2]))])
  
  return(list(optimum_primitives,optimum_ntrees,optimum_mtry))
  
}

## Create the model with the optimum parameters
model_param<-get_optimum_rf_model(ground_truth_df_model)
ground_truth_df_model<-ground_truth_df_model[,c("response",model_param[[1]])]
model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model, ntree=model_param[[2]], mtry=model_param[[3]])

# Retrieve the paths of the features to keep
#columns_to_keep<-model_param[[1]]
#paths<-left_join(data.frame(columns_to_keep),df_primitives_types_sources,by = c("columns_to_keep"="column_names_primitives"))
#paths<-unique(data.frame(paths$indice,paths$indices_for_classif_paths))
#paths_indices_to_compute<-paths[which(!is.na(paths$paths.indices_for_classif_paths)),]

########################################################################################################################
########################################################################################################################
############ Step 9 - Classify  ############
########################################################################################################################
########################################################################################################################

# Open objects segmentation dataset
segmentation_gpkg<-sf::st_read(path_to_segmented_dataset_stats)
segmentation_df<-as.data.frame(segmentation_gpkg)
segmentation_df$geom<-NULL
# Fill-in NA values with median value of the column (using randomForest::na.roughfix)
segmentation_df<-randomForest::na.roughfix(segmentation_df)
# Predict land cover class using the predict fonction with the model extracted in step 8
segmentation_df$predicted<-predict(model,segmentation_df)
# Merge back to segmentation_gpkg and write geopackage with the predicted values of the objects
segmentation_df<-segmentation_df[,c("cat","predicted")]
segmentation_gpkg<-merge(segmentation_gpkg,segmentation_df,by="cat")
sf::st_write(segmentation_gpkg,path_to_output_classification_vector,layer_options = "OVERWRITE=true")


########################################################################################################################
########################################################################################################################
############ Step 10 - Post-processing  ############
########################################################################################################################
########################################################################################################################

##############################################################
#### 10.1 - get clean raster dataset  ####
##############################################################

## Rasterize the classification
# Set output raster characteristics
poly<-sf::st_read(path_to_output_classification_vector)
output_res<-res(raster(path_to_spot67_preprocessed_pan))[1]
r <- raster(poly, res = output_res)
# Rasterize using fasterize (fast version of rasterize)
predicted<-data.frame(predicted=unique(poly$predicted))
predicted$predicted_integer<-seq(1:nrow(predicted))
poly<-merge(poly,predicted,by="predicted")
r <- fasterize::fasterize(poly, r, field = "predicted_integer")
# Write the classification raster
writeRaster(r,path_to_output_classification_raster)

##############################################################
#### 10.2 - get clean vector dataset  ####
##############################################################

# Vectorize back (now that the adjacent polygons that have the same class have been gathered)
gdal_appli<-paste0("gdal_polygonize.py ",path_to_output_classification_raster," ",gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector)," -b 1 None DN")
system(gdal_appli)
# Get back the labels of the classes on the vector version of the classification
poly<-sf::st_read(gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector))
poly<-merge(poly,predicted,by.x="DN",by.y="predicted_integer")
sf::st_write(poly,path_to_output_classification_vector,layer_options = "OVERWRITE=true")
file.remove(gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector))


########################################################################################################################
############ Close workflow ############
########################################################################################################################

# Remove grass temorary folder
system(paste0("rm -r ", file.path(getwd(),"GRASS_TEMP")))  
file.remove(file.path(getwd(),".grassrc7"))









########################################################################################################################
############ Step 7 - Pre-processing and partition of the ground truth database  ############
########################################################################################################################

#ground_truth<-readOGR(file.path(path_to_processing_folder,path_to_ground_truth_data),stringsAsFactors = F)
#path_to_ground_truth_processed<-file.path(path_to_processing_folder,"Ground_truth","ground_truth.gpkg")

# Add an ID colum
#ground_truth$id<-seq(1:nrow(ground_truth))
# convert as data frame for future processings
#ground_truth_df<-as.data.frame(ground_truth)
# Set land cover classes as integers
#for (i in 1:length(column_names_lc_classes)){
#gt_lc_types<-data.frame(lc_type=unique(ground_truth_df[,column_names_lc_classes[i]]))
#gt_lc_types[,paste0(column_names_lc_classes[i],"_int")]<-seq(1:nrow(gt_lc_types))
#ground_truth<-merge(ground_truth,gt_lc_types,by.x=column_names_lc_classes[i],by.y="lc_type")
#ground_truth_df<-merge(ground_truth_df,gt_lc_types,by.x=column_names_lc_classes[i],by.y="lc_type")
#}

#writeOGR(ground_truth,path_to_ground_truth_processed,driver = "GPKG",layer="ground_truth",overwrite_layer = TRUE)

# Get training samples ids. We keep 70% of the data for training (size=0.7) and 30% for validation
#ground_truth_training<-ground_truth_df %>% group_by_(paste0(column_names_lc_classes,"_int")) %>% sample_frac(size = 0.7)
#ground_truth_training_ids<-ground_truth_training$id
# Get validation samples ids
#ground_truth_validation_ids<-setdiff(ground_truth_df$id,ground_truth_training_ids)

#ground_truth_df<-unique(ground_truth_df$type_1,ground_truth_df$type_1_int)
# Extract training dataset
#ground_truth_training_ogr<-ground_truth[which(ground_truth$id %in% ground_truth_training_ids),]

# Extract validation dataset
#ground_truth_validation_ogr<-ground_truth[which(ground_truth$id %in% ground_truth_validation_ids),]

# Save training and validation datasets as GPKG datasets
#file.remove(path_to_ground_truth_training)
#file.remove(path_to_ground_truth_validation)
#writeOGR(ground_truth_training_ogr,path_to_ground_truth_training,driver = "GPKG",layer="ground_truth")
#writeOGR(ground_truth_validation_ogr,path_to_ground_truth_validation,driver = "GPKG",layer="ground_truth")


###### Train classifier using OTB
#path_to_xml_feature_statistics<-file.path(path_to_processing_folder,path_to_groundtruth_folder,"feature_statistics.xml") # Path to the XML feature statistics used as input of the TrainVectorClassifier application
#path_to_output_model<-file.path(path_to_processing_folder,path_to_classification_folder,"model_file.txt") # Path to output model 
#path_to_confusion_matrix<-file.path(path_to_processing_folder,path_to_classification_folder,"confusion_matrix.csv") # Path to confusion matrix

# Set column names for the classification
#ground_truth_training<-as.data.frame(readOGR(path_to_ground_truth_training_stats,stringsAsFactors = F))
# All columns :
#columns_stats<-setdiff(colnames(ground_truth_training),c(column_names_lc_classes,"cat","cat_","id",paste0(column_names_lc_classes,"_int")))
# 

# Generate the input XML statistics file (used as input of the OTB TrainVectorClassifier application)
#otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_ComputeOGRLayersFeaturesStatistics")," -inshp ",path_to_ground_truth_training_stats," -outstats ",path_to_xml_feature_statistics," -feat ",paste(columns_stats, collapse = " "))
#system(otb_appli)

# Train vector classifier
#otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_TrainVectorClassifier")," -io.vd ",path_to_ground_truth_training_stats," -io.stats ",path_to_xml_feature_statistics," -io.out ",path_to_output_model," -io.confmatout ",path_to_confusion_matrix," -cfield ",paste0(column_names_lc_classification,"_int")," -feat ",paste(columns_stats, collapse = " ")," -valid.vd ",path_to_ground_truth_validation_stats)
#otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_TrainVectorClassifier")," -io.vd ",path_to_ground_truth_training_stats," -io.stats ",path_to_xml_feature_statistics," -io.out ",path_to_output_model," -io.confmatout ",path_to_confusion_matrix," -cfield ",paste0(column_names_lc_classification,"_int")," -feat ",paste(columns_stats, collapse = " ")," -valid.vd ",path_to_ground_truth_validation_stats," -classifier sharkrf -classifier.sharkrf.nbtrees 45 -classifier.sharkrf.nodesize 10 -classifier.sharkrf.mtry 0 -classifier.sharkrf.oobr 0.66")

#system(otb_appli)

# Open confusion matrix and pre-process
#conf_mat<-read.csv(path_to_confusion_matrix,skip = 1)
#class_labels<-unique(ground_truth_training[,c(column_names_lc_classification,paste0(column_names_lc_classification,"_int"))])
#class_labels<-class_labels[order(class_labels[,paste0(column_names_lc_classification,"_int")]),]
#colnames(conf_mat)<-paste0(class_labels[,column_names_lc_classification],"_reference")
#rownames(conf_mat)<-paste0(class_labels[,column_names_lc_classification],"_produced")

# Get values of the conf matrix in percentage instead of absolute numbers
#conf_mat_perc<-conf_mat
#for (i in 1:nrow(conf_mat_perc)){
#  conf_mat_perc[i,]<- conf_mat_perc[i,]/sum(conf_mat_perc[i,])*100
#  }
# for (i in 1:nrow(conf_mat_perc)){ print(paste(class_labels[,column_names_lc_classification][i],conf_mat_perc[i,i],"% (nb tot:",sum(conf_mat[i,]),")")) }


