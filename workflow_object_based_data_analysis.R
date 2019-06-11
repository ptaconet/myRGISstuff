######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = workflow_hierarchical_random_forest_and_obia_classification, title = Workflow for the generation of land cover map using a supervised object based image analysis combined with a hierarchical random forest classification approach , abstract = This script creates a land cover map on a region of interest using a supervised object based image analysis approach combined with a hierarchical random forest classification. The hierarchical classification enables to account for the class hierarchical structure. Additional information are provided below.
# wps.in: id = path_to_otbApplications_folder, type = string, title = Path to the folder containing the Orfeo Toolbox applications. , value = "/home/ptaconet/OTB-6.6.1-Linux64/bin";
# wps.in: id = path_to_grassApplications_folder, type = string, title = Path to the folder containing the GRASS applications. , value = "/usr/lib/grass74";
# wps.in: id = path_to_processing_folder, type = string, title = Path to working folder (i.e. where all input and ouput datasets are/will be stored) , value = "/home/ptaconet/Documents/react/data_CIV";
# wps.in: id = path_to_roi_vector, type = string, title = Relative path to the ROI in kml or shp format in the path_to_processing_folder folder, value = "ROI.kml";
# wps.in: id = path_to_spot67_raw_folder, type = string, title = Relative path to the folder containing the Spot 6/7 raw data (tar.gz files) in the path_to_processing_folder folder, value = "VHR_SPOT6/raw_data";
# wps.in: id = path_to_groundtruth_folder, type = string, title = Relative path to the folder containing the gound truth (i.e. training/validation) dataset in the path_to_processing_folder folder, value = "Ground_truth";
# wps.in: id = path_to_ground_truth_data, type = string, title = Name of the ground truth dataset in the vfolder, value = "groundtruth_bf.gpkg" ;
# wps.in: id = column_names_lc_classes_hierarchy, type = string, title = Vector of the columns names of land cover classes in the ground truth dataset, with the hierarchical structure. eg : c("type_1","type_2"). Typically type_1 is the most aggregated land cover, type_2 is a less aggregated classification, etc., value = c("level_1_fr","level_2_fr","level_3_fr","level_4_fr","level_5_fr") ;
# wps.in: id = column_name_lc_classification, type = string, title =  Desired hierarchy level for the final classification, value = "level_4_fr";
# wps.in: id = path_to_copernicus_scihub_credentials, type = string, title = Path to the text file containing the credentials to ESA Copernicus scihub, value = "credentials_copernicus.txt";
# wps.in: id = Sentinel2_products_uuid, type = string, title = Vector of Ids of the S2 products to download in the Copernicus Scihub, value = c("bc6bafd7-d44f-4d62-8754-3cc4ba4e8cc0","18895056-852f-4a4f-a3aa-ca7882fe79de") ;
# wps.in: id = proj_srs, type = string, title = Proj srs for the ROI, value = "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs" ;
# wps.in: id = threshold_accumulation_raster, type = integer, Threshold to use for the water accumulation raster file as hydrographic network. All cells above this threshold are considered as the hydrographic network., value = "800" ;
# wps.in: id = xrad, type = string, value = Vector of sizes of the moving windows in the x direction for the computation of the textures from the Spot 6/7 panchromatic images, value = c(5,9,17) ;
# wps.in: id = yrad, type = string, value = Vector of sizes of the moving windows in the y direction for the computation of the textures from the Spot 6/7 panchromatic images, value = c(5,9,17) ;
# wps.in: id = nnbin, type = integer, value = Number of bins to consider for the computation of the textures from the Spot 6/7 panchromatic images, value = 64 ;
# wps.in: id = radiometric_indices_list_spot67, type = string, value = Vector of radiometric indices to compute from the Spot 6/7 image for the classification, value = c("Vegetation:NDVI","Water:NDWI","Soil:BI") ;
# wps.in: id = segmentation_threshold, type = integer, value = Scale parameter for the VHRS segmentation, value = 100 ;
# wps.in: id = segmentation_cw, type = float, value = Shape parameter for the segmentation, value = 0.1 ;
# wps.in: id = segmentation_sw, type = float, value = Spectral parameter for the segmentation, value = 0.9 ;

## Author : Paul Taconet, French research institute for Sustainable Development (IRD)

## ideas for further improvements: run classification n times (e.g. 15 times), not only 1 time, to improve overall classif stats

#### Important information
# The script takes as input a Very high resolution image of the area of interesed (Spot 6/7) + the Copernicus scihub identifiers of a High resolution image (Sentinel 2) + a learning/validation ground truth dataset, enventually with a class hierarchical structure. 
# It performs all the pre-processing, data preparation, classification and post-processing steps. Details are provided below. 
# The script uses various R packages (see section "prepare workflow" for a list of packages used) and open-source libraries that are called through R : the Orfeo Toolbox v6.6.1 (https://www.orfeo-toolbox.org/), the GRASS library v7.4 (https://grass.osgeo.org), the GDAL library v2.2.1 (https://www.gdal.org/) and SAGA GIS v2.3.1 (http://www.saga-gis.org/en/index.html). In addition, the workflow uses a personal release of the Orfeo Toolbox for the segmentation process (since no relevant application was found in the official OTB release for the segmentation of very large images using the Baatz and Shape Generic Region Merging algorithm). The release was generated and kindly provided by Rafaelle Gaetano. It is available here: http://napoli.teledetection.fr/logiciels/otb_moringa_build_win_x64.zip
# The methodology used in this workflow was inspired from these two articles and uses an R package developed in the frame of the Gavish et al. article (unfortunataly neither available on a git repository nor on the CRAN) (with minor adaptations for this workflow): 
  # - Gavish et al., Comparing the performance of flat and hierarchical Habitat/Land-Cover classification models in a NATURA 2000 site, ISPRS Journal of Photogrammetry and Remote Sensing. Volume 136, February 2018, Pages 1-12 https://doi.org/10.1016/j.isprsjprs.2017.12.002
  # - Lebourgeois et al., A Combined Random Forest and OBIA Classification Scheme for Mapping Smallholder Agriculture at Different Nomenclature Levels Using Multisource Data (Simulated Sentinel-2 Time Series, VHRS and DEM). Remote Sens. 2017, 9, 259.   https://doi.org/10.3390/rs9030259

#### Workflow steps :
### Step 1 - Download the Digital Elevation Model (SRTM tiles) for the ROI
### Step 2 - Download the ancillary data :
  ## 2.1 - Download the HRS Sentinel 2 image(s) from the Copernicus scihub
### Step 3 - Pre-process the VHRS Spot6/7 image(s) :
  ## 3.1 - fusion the tiles of the panchromatic image
  ## 3.2 - convert the multispectral and panchromatic images from digital numbers to TOA reflectance
  ## 3.3 - orthorectify the multispectral and panchromatic images
  ## 3.4 - extract the ROI
  ## 3.5 - pansharpen the MS image using the PAN image
  ## 3.6 - mosaic the various tiles covering the ROI (if relevant)
### Step 4 - Preprocess the ancillary data :
  ## 4.1 - preprocess the DEM : mosaic the various tiles covering the ROI (if relevant), and then extract the ROI
  ## 4.2 - preprocess the Sentinel 2 image(s) : mosaic the various images covering the ROI (if relevant), and then extract the ROI
### Step 5 - Prepare the data for the classification : 
  ## 5.1 - extract ancillary indices from the DEM : slope, aspect, flow accumulation, flow direction, topographic convergence index
  ## 5.2 - extract textural indices from the Spot6/7 panchromatic image at various moving windows sizes
  ## 5.3 - extract radiometric indices from the Spot6/7 pansharpened image
  ## 5.4 - extract radiometric indices from the S2 image
  ## 5.5 - Split the bands of the Spot6/7 image
### Step 6 - Segment the Spot6/7 image using the Baatz and Shape Generic Region Merging algorithm
### Step 7 - Extract the zonal statistics (reflectance, spectral indices, textural indices, ancillary, shape, contextual)
  ## 7.1 - Extract zonal statistics for the ground truth dataset
  ## 7.2 - Extract zonal statistics for the segmented objects dataset
### Step 8 - Prepare classification : generate a set of random forest classifiers using the training dataset with two approaches: i) a flat classification approach, ii) a class hierarchical structure approach. Classify the objects output of the segmentation process using the approach that gives the best results
  ## 8.1 - Generate the RF classifiers at each class hierarchical level using i) a flat approach, ii) a hierarchical approach, and compare the results
  ## 8.2 - Get useful information on the classification (discriminant variables, etc.) considering the class hierarchical structure
### Step 9 - Classify
  ## 9.1 - Classify the objects output of the segmentation using the approach that gives the best results
  ## 9.2 - Save the classification as GIS data in various formats (vector and raster)
  ## 9.3 - Add user criterions to enhance the classification


#### Outputs: In the path_to_processing_folder, the following folders and files will be available : 
  ### folder "DEM_SRTM" : data regarding the SRTM Digital Elevation Model [GENERATED BY THE WF AS WELL AS ALL SUB-FOLDER AND FILES]
    ## DEM_SRTM/raw_data : SRTM tile(s) covering the ROI (step 1)
    ## DEM_SRTM/processed_data : a set of files derived from the DEM, output of the workflow :
      # - DEM.tif : SRTM DEM cut following the ROI and re-projected in the proj_srs (step 4) ;
      # - slope.tif : slope dataset derived from the DEM (step 5.1) ;
      # - aspect.tif : aspect dataset derived from the DEM (step 5.1);
      # - accumulation.tif : flow accumulation dataset derived from the DEM (step 5.1);
      # - direction.tif : flow direction dataset derived from the DEM (step 5.1);
      # - tci.tif : topographic convergence index derived from the DEM (note : not used in further analysis) (step 5.1);
      # - accumulation_treshold.tif : raster with two values : 0 if the flow accumulation is above the threshold_accumulation_raster , 1 if it is under. This information is used for the calculation of the distance from the objects to the flow accumulation network (used as primitive in the classification). (step 5.1) ;
      # - accumulation_treshold_vector.gpkg : vector version of accumulation_treshold.tif (step 5.1) ;
  ### folder "VHR_SPOT6" : data regarding the VHRS Spot 6/7 satellite images
    ## VHR_SPOT6/raw_data : the input Spot 6/7 data (as provided by the CNES). Must be provided by the user before execution of the workflow. There must be 1 sub-folder by image covering the ROI. Each sub-folder contains the two .tar.gz files as provided by the CNES : one for the panchromatic image and one for the multispectral image [PROVIDED BY THE USER]
    ## VHR_SPOT6/processed_data : a set of files derived from the Spot 6/7 datasets, output of the workflow:  [GENERATED BY THE WF AS WELL AS ALL SUB-FOLDER AND FILES]
      # - PAN.TIF : the panchromatic image, mosaiced, orthorectified and cut following the ROI (step 3.x) ;
      # - PANSHARPEN.TIF : the multispectral image pansharpened using the panshromatic image, orthorectified and cut following the ROI (step 3.x) ;
      # - PANSHARPEN_0.TIF, PANSHARPEN_1.TIF, PANSHARPEN_2.TIF, PANSHARPEN_3.TIF : respectively the blue, green, red and nir bands of the VHRS (step 3.x) ;
      # - HaralickTextures_xxxx.TIF : the set of textures generated using from the panchromatic image, at various moving windows sizes (step 5.2) ;
      # - NDVI.TIF, NDWI.TIF, BI.TIF : the radiometric indices derived from the VHRS (step 5.3).
  ### folder "HR_Sentinel2" : data regarding the HRS Sentinel 2 satellite images  [GENERATED BY THE WF AS WELL AS ALL SUB-FOLDER AND FILES]
    ## HR_Sentinel2/raw_data : the input Sentinel 2 data (as downloaded by the WF in the Copernicus Scihub) (step 2.1). 
    ## HR_Sentinel2/processed_data : a set of files derived from the Sentinel 2 dataset(s), output of the workflow 
      # B01.TIF, BO2.TIF, etc... B12.TIF : bands of the Sentinel 2 images cut following the ROI. If multiple images are provided as input (i.e. if the ROI is covered by multiple tiles), the images will be mosaiced first (step 4.2 and 5.5) ;
      # BRI.TIF, MNDVI.TIF, MNDWI.TIF, NDVI.TIF, NDWI.TIF, RNDVI.TIF : the radiometric indices derived from the HRS. Formulas are provided in the workflow, section 5.4 (step 5.4) ;
  ### folder "Segmentation" : [GENERATED BY THE WF AS WELL AS ALL SUB-FOLDER AND FILES]
      # segmentation_vector.gpkg : vector output of the segmentation process (i.e. objects that will further be classified) (step 6) ;
      # segmentation_dataset_stats.gpkg : segmented datasets with the zonal statistics for each object (step 7.2).
  ### folder "Ground_truth" : data regarding the ground truth (i.e. training/validation) dataset 
      # path_to_ground_truth_data : ground truth dataset provided by the user, including the class hierarchical structure ; 
      # ground_truth_stats.gpkg : ground truth dataset with the zonal statistics for each object (step 7.1)
  ### folder "Classification" : data regarding the classification [GENERATED BY THE WF AS WELL AS ALL SUB-FOLDER AND FILES]  ########## A FINIR 
      # classes_hierarchy.png : a figure of the class hierarchical structure (step 8.2) ; 
      # flat_classif_stats_xxx.png : a figure with the confusion matrix + plot of variable importance at the class hierarchical level xxx, using a flat classification approach  (step 8.1)
      # hierar_classif_stats_xx.png : a figure with the confusion matrix + plot of variable importance at the class level xxx (step 8.1)
      # var_importance_yyyy.png : a figure showing the variables importance at each classification level, sorted by yyy {variable source, variable stat type, variable type}
      # classification.gpkg : Vector with the objects classified + the zonal stats (step 9.2)
      # classification.tif : Raster version of the classification (step 9.2)
      # classification_group.gpkg : Vector with adjacent objects having the same class grouped (step 9.2)

#### Start Workflow
########################################################################################################################
############ Set Input parameters for the workflow ############
########################################################################################################################

rm(list = ls())

### Global variables used throughout the WF
path_to_otbApplications_folder<-"/home/ptaconet/OTB-6.6.1-Linux64/bin"
path_to_grassApplications_folder<-"/usr/lib/grass74" #<Can be retrieved with grass74 --config path . More info on the use of rgrass7 at https://grasswiki.osgeo.org/wiki/R_statistics/rgrass7
path_to_processing_folder<-"/home/ptaconet/Documents/react/data_CIV"  #<Path to the processing folder (i.e. where all the data produced by the workflow will be stored)>
path_to_roi_vector="ROI.kml" #<Path to the Region of interest in KML format>


### Parameters for step 3 : Downloading ancillary data
path_to_copernicus_scihub_credentials<-"credentials_copernicus.txt" # <path to the file containing the credential to the ESA Sentinel data server (Copenicus Scihub)>
Sentinel2_products_uuid<-c("6236fb46-41c1-4950-b6c8-602c48b90049","89bc8775-cc0c-4dc8-8ee3-6c9115d75fa3") #<Ids of the products to download in the Copernicus Scihub>
time_range <-  c("2018-12-01", "2018-12-10")
# BF: Sentinel2_products_uuid<-c("bc6bafd7-d44f-4d62-8754-3cc4ba4e8cc0","18895056-852f-4a4f-a3aa-ca7882fe79de")
# CIV: Sentinel2_products_uuid<-c("6236fb46-41c1-4950-b6c8-602c48b90049","89bc8775-cc0c-4dc8-8ee3-6c9115d75fa3)

### Parameters for step 4
proj_srs="+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs" #<proj srs for the ROI>
threshold_accumulation_raster<-800   #<Threshold for the water accumulation raster file. All cells above this threshold are considered as the hydrographic network.>
xrad=c(5,9,17)  #<size of the moving windows in the x direction for the computation of the textures>
yrad=c(5,9,17)  #<size of the moving windows in the y direction for the computation of the textures>
nbbin=64  #<number of bins in the output image for the computation of the textures>
# BF: threshold_accumulation_raster<-1000
# CIV: threshold_accumulation_raster<-800

### Parameters for step 5
radiometric_indices_list_spot67<-c("Vegetation:NDVI","Water:NDWI2","Soil:BI") #<Radiometric indices to compute for the Spot6/7 image>

### Parameters for step 6
segmentation_threshold=100       #<segmentation scale parameter> 
segmentation_cw=0.1              #<segmentation shape parameter> 
segmentation_sw=0.9             #<segmentation spectral parameter>
# BF: segmentation_threshold<-100
# CIV: segmentation_threshold<-160
# BF: segmentation_cw<-0.1
# CIV: segmentation_cw<-0.1
# BF: segmentation_sw<-0.9
# CIV: segmentation_sw<-0.8

### Parameters for step 7
path_to_groundtruth_folder<-file.path(path_to_processing_folder,"Ground_truth")
path_to_ground_truth_data<-file.path(path_to_groundtruth_folder,"civ_groundtruth_objects_segmentation_v_classes_update.gpkg") #file.path(path_to_groundtruth_folder,"groundtruth_bf.gpkg") #<Path to the ground truth dataset. The geometry column must be named "geom">
methods_to_compute<-"avg,stddev" #<methods_to_compute for the primitives. Available methods are: "minimum,maximum,range,average,stddev,variance,coeff_var,first_quartile,median,third_quartile,percentile">
indices_for_classif_labels<-c("DEM",
                              "slope",
                              "accumulation",
                              "NDVI_SPOT6",
                              "NDWI2_SPOT6",
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
                              "blue_SPOT6",
                              "green_SPOT6",
                              "red_SPOT6",
                              "nir_SPOT6",
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
                              ) # parameter indices_for_classif_paths sets the path to each variable
indices_for_classif_paths<-c(file.path(path_to_processing_folder,"DEM_SRTM/processed_data/DEM.tif"),
                             file.path(path_to_processing_folder,"DEM_SRTM/processed_data/slope.tif"),
                             file.path(path_to_processing_folder,"DEM_SRTM/processed_data/accumulation.tif"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/NDVI.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/NDWI2.TIF"),
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
column_names_lc_classes_hierarchy<-c("L1","L2","L3","L4","L5") #<Names of the columns of land cover classes in the ground truth dataset. eg : c("L1","L2"). Typically L1 is the most aggregated land cover, L2 is a less aggregated classification, etc.>
#column_name_lc_classification<-"L4" #<Name of the column of land cover class that we want to classify in the ground truth dataset. Column type must be character string>

########################################################################################################################
############ Prepare workflow ############
########################################################################################################################

#### Call useful packages
require(raster)
require(rgrass7)
require(RSAGA)
require(rgdal)
require(gdalUtils)
require(sf)
require(fasterize)
require(dplyr)
require(reshape)
require(reshape2)
require(ggplot2)
require(grid)
require(gridExtra)
require(cowplot)
require(randomForest)
require(caret)
require(HieRanFor) ## Downloaded in the supplementary material of the Gavish et al. article. To download here: https://ars.els-cdn.com/content/image/1-s2.0-S0924271617303696-mmc5.zip and manually install on the R packages folder of your system.
require(httr)


#### Provide useful functions

#all_functions<-list.files("/home/ptaconet/r_react/functions",full.names = T)
#lapply(all_functions, source)

### Functions for data pre-processing

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

## Function to mosaic multiple images. Details: https://github.com/remicres/otb-mosaic
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


### Functions for data classification

## Function to plot and save a set of information about a classification: confusion matrix and variables of inmportance
plot_classif_info<-function(conf_matrix,
                            var_importance,
                            classifier_name,
                            classifier_ID=NULL,
                            classif_overall_stats=NULL,
                            path_to_output_plot){
  
  conf_matrix_im<-as.data.frame(conf_matrix)
  conf_matrix_im$class.error<-NULL
  conf_matrix_im$sum=rowSums(conf_matrix_im)
  actual = data.frame(Actual=rownames(conf_matrix_im),ActualFreq=conf_matrix_im$sum)
  
  confusion<-reshape2::melt(conf_matrix)
  confusion <- confusion %>% filter (Var2!="class.error")
  names(confusion)<-c("Actual","Predicted","Freq")
  
  
  confusion = merge(confusion, actual, by=c("Actual"))
  confusion$Percent = round(confusion$Freq/confusion$ActualFreq*100)
  
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=Percent),data=confusion, color="black",size=0.1) +
    labs(x="Actual",y="Predicted")+ theme_grey(base_size = 15)
  tile = tile + 
    geom_text(aes(x=Actual,y=Predicted, label=Freq),data=confusion, size=4.5, colour="black") +
    #scale_fill_gradient(low="gray87",high="gray48") 
    scale_fill_gradientn(colours = rev(gray.colors(100))) + theme(axis.text.x=element_text(angle=40, hjust=1))
  # lastly we draw diagonal tiles. We use alpha = 0 so as not to hide previous layers but use size=0.3 to highlight border
  tile = tile + 
    geom_tile(aes(x=Actual,y=Predicted),data=subset(confusion, as.character(Actual)==as.character(Predicted)), color="black",size=0.3, fill="black", alpha=0) 
  
  # create an apporpriate viewport.  Modify the dimensions and coordinates as needed
  vp.BottomRight <- viewport(height=unit(0.45, "npc"), width=unit(1, "npc"), 
                             just=c("top"), 
                             y=0.45, x=0.5)
  
  png(path_to_output_plot,width = 1000, height = 1000, units = "px" )
  par(mfrow=c(2,1))
  
  # Plot of variable importance
  # Dotplot: Grouped Sorted and Colored and keep only the first 30 variables
  var_importance <- var_importance[order(-var_importance$MeanDecreaseGini)[1:30],] 
  var_importance$color[var_importance$type=="VHSR"] <- "red"
  var_importance$color[var_importance$type=="HSR"]  <- "blue"
  var_importance$color[var_importance$type=="ancillary"]  <- "forestgreen"
  
  if(!is.null(classif_overall_stats)){ 
    classif_stats<-paste0("\nKappa=",round(classif_overall_stats[2,1],2)," ; Accuracy=",round(classif_overall_stats[1,1],2))
  } else {
    classif_stats<-NULL
  }
  
  dotchart(rev(var_importance$MeanDecreaseGini),labels=rev(var_importance$column_names_primitives),cex=.7,#groups= var_importance$type,
           main=paste0(classifier_ID," - ",classifier_name," - \n",classif_stats),sub="red = VHSR ; blue = HSR ; green = ancillary",xlab="MeanDecAcc", col=rev(var_importance$color), xlim = c(0,max(var_importance$MeanDecreaseGini))) 
  
  # plot the ggplot using the print command
  print(tile, vp=vp.BottomRight)
  dev.off()
}


## Function to plot the importance of the primitives at each classification level on the hierarchical classification approach using random forest. This function is adapted from the HieRanFor::PlotImportanceHie function
PlotImportanceHie_v_taconet<-function (input.data, X.data = 2, Y.data = 3, imp.data = 4, explanatory.variable = NULL, plot.type = "Tile", 
                                       imp.title = colnames(input.data)[imp.data], X.Title = colnames(input.data)[X.data], 
                                       Y.Title = colnames(input.data)[Y.data], low.col = "blue", 
                                       high.col = "red", geom.tile.bor.col = "white", pos.col = "green", 
                                       zero.col = "white", neg.col = "red", supp.warn = TRUE, ...) 
{
  X.var <- Y.var <- NULL
  localenv <- environment()
  if (plot.type != "Tile" && plot.type != "Bubble") {
    stop("\nplot.type can be either 'Tile' or 'Bubble'\n")
  }
  if (!is.null(explanatory.variable)){
  work.data <- input.data[, c(X.data, Y.data, imp.data, explanatory.variable)]
  colnames(work.data) <- c("X.var", "Y.var", "imp.var", "type")
  } else {
    work.data <- input.data[, c(X.data, Y.data, imp.data)]
    colnames(work.data) <- c("X.var", "Y.var", "imp.var")
  }
  
  # BUG IN THE SOURCE FUNCTION WITH THE FOLLOWING LINES : 
  #work.data$X.var <- suppressWarnings(factor(work.data$X.var, 
  #                                          as.character(work.data$X.var)))
  #work.data$Y.var <- suppressWarnings(factor(work.data$Y.var, 
  #                                           as.character(work.data$Y.var)))
  
  # PTACONET FIX BUG :
  work.data$X.var <- suppressWarnings(as.factor(as.character(work.data$X.var)))
  work.data$Y.var <- suppressWarnings(as.factor(as.character(work.data$Y.var)))
  if (!is.null(explanatory.variable)){
  work.data$type <- suppressWarnings(as.factor(as.character(work.data$type)))
  work.data=work.data[order(work.data$type),]
  }
  if (plot.type == "Tile") {
    p <- ggplot(work.data, aes(X.var, Y.var), environment = localenv)
    if (!is.null(explanatory.variable)){
    p <- p + geom_tile(aes(fill = work.data$imp.var), color = geom.tile.bor.col) + facet_grid(type~.,scales="free_y",space="free_y")
    p <- p + labs(y="var. used",title = paste0("Mean decrease in accuracy of each explanatory variable in each local classifier by ",names(input.data[explanatory.variable]))) 
    } else {
      p <- p + geom_tile(aes(fill = work.data$imp.var), color = geom.tile.bor.col)
      p <- p + labs(y="var. used",title = paste0("Variables importance by classifier")) 
    }
    p <- p + xlab(X.Title)
    #p <- p + ylab(Y.Title)
    p <- p + scale_fill_continuous(low = low.col, high = high.col, 
                                   guide = "colourbar", guide_colourbar(title = imp.title))
  }
  if (plot.type == "Bubble") {
    #for (k3 in 1:dim(work.data)[1]) {
    #  if (work.data[k3, 3] > 0) {
    #    work.data[k3, 4] <- "Positive"
    #  }
    #  if (work.data[k3, 3] == 0) {
    #    work.data[k3, 4] <- "Zero"
    #  }
    #  if (work.data[k3, 3] <= 0) {
    #    work.data[k3, 4] <- "Negative"
    #  }
    #}
    #colnames(work.data)[4] <- "sign.imp"
    #fillScaleValues <- c(Positive = pos.col, Zero = zero.col, 
    #                     Negative = neg.col)
    p <- ggplot(work.data, aes(X.var, Y.var), environment = localenv)
    if (!is.null(explanatory.variable)){
    p <- p + geom_point(aes(size = abs(work.data$imp.var), 
                            fill = work.data[,4]), shape = 21) + facet_grid(type~.,scales="free_y",space="free_y")
    p <- p + labs(y="var. used",fill = names(input.data[explanatory.variable]),title = paste0("Mean decrease in accuracy of each explanatory variable in each local classifier by ",names(input.data[explanatory.variable]))) 
    } else {
      p <- p + geom_point(aes(size = abs(work.data$imp.var), 
                              fill = work.data[,4]), shape = 21)
      p <- p + labs(y="var. used",fill = names(input.data[explanatory.variable]),title = paste0("Variables importance by classifier"))
    }
    p <- p + xlab(X.Title)
    #p <- p + ylab(Y.Title)
    p <- p + scale_size_continuous(guide = "legend", guide_legend(title = imp.title))
    #p <- p + scale_fill_manual(values = fillScaleValues, 
    #                           name = explanatory.variable.name)
  }
  if (supp.warn) {
    return(suppressWarnings(print(p)))
  }
  else {
    return(print(p))
  }
}

# Function PerformanceHRF
path_to_PerformanceHRF_modif<-file.path(path_to_processing_folder,"PerformanceHRF_modif.R")
download.file("https://raw.githubusercontent.com/ptaconet/r_react/master/functions/PerformanceHRF_modif.R",path_to_PerformanceHRF_modif)
insertSource(path_to_PerformanceHRF_modif, package = "HieRanFor", functions="PerformanceHRF")
file.remove(path_to_PerformanceHRF_modif)


## Function to save the outputs of the classification to the disk
# - classification.gpkg : Vector with the objects classified + the zonal stats
# - classification.tif : Raster version of the classification
# - classification_group.gpkg : Vector with adjacent objects having the same class grouped
save_classif_to_disk_function<-function(dataset_classified, dataset_to_classify_sf, path_to_output_classification_vector, path_to_outputs_folder, save_objects_raster=FALSE, save_objects_vector_group_adj_polygons=FALSE ){
  
  cat("Saving the classification as vector gpkg...\n")
  #dataset_to_classify_sf<-merge(dataset_classified,dataset_to_classify_sf,by.x="cat",by.y="DN")  ## this is for hierarchical classif
  dataset_classified<-dataset_classified[,c("DN","predicted")]
  dataset_to_classify_sf<-merge(dataset_classified,dataset_to_classify_sf,by="DN")
  dataset_to_classify_sf<-st_as_sf(dataset_to_classify_sf)
  sf::st_write(dataset_to_classify_sf,path_to_output_classification_vector,layer_options = "OVERWRITE=true")
  
  if (save_objects_raster==TRUE){
    cat("Saving the classification as raster tif...\n")
    ## Rasterize the classification
    require(raster)
    require(fasterize)
    # Set output raster characteristics
    poly<-sf::st_read(file.path(path_to_outputs_folder,"classification.gpkg"))
    output_res<-1.633175
    r <- raster(poly, res = output_res)
    # Rasterize using fasterize (fast version of rasterize)
    predicted<-data.frame(predicted=unique(poly$predicted))
    predicted$predicted_integer<-seq(1:nrow(predicted))
    poly<-merge(poly,predicted,by="predicted")
    r <- fasterize::fasterize(poly, r, field = "predicted_integer")
    # Write the classification raster
    path_to_output_classification_raster<-file.path(path_to_outputs_folder,"classification.tif")
    writeRaster(r,path_to_output_classification_raster, overwrite=TRUE, datatype='INT2S')
  }
  
  if (save_objects_vector_group_adj_polygons==TRUE){
    cat("Saving the classification as vector gpkg with adjacent objects with the same class grouped...\n")
    # Vectorize back (now that the adjacent polygons that have the same class have been gathered)
    gdal_appli<-paste0("gdal_polygonize.py ",path_to_output_classification_raster," ",gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector)," -b 1 None DN")
    system(gdal_appli)
    # Get back the labels of the classes on the vector version of the classification
    poly<-sf::st_read(gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector))
    poly<-merge(poly,predicted,by.x="DN",by.y="predicted_integer")
    path_to_output_classification_vector_objects_grouped<-file.path(path_to_outputs_folder,"classification_group.gpkg")
    sf::st_write(poly,path_to_output_classification_vector_objects_grouped,layer_options = "OVERWRITE=true")
    file.remove(gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector))
  }
  
}


### Set working directory
setwd(path_to_processing_folder)

## Connection to Copernicus (Sentinel data)
copenicus_credentials<-readLines(path_to_copernicus_scihub_credentials)
copernicus_scihub_username<-strsplit(copenicus_credentials,"=")[[1]][2]
copernicus_scihub_password<-strsplit(copenicus_credentials,"=")[[2]][2]

### Set the paths of output folders / files
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
# Step 8
path_to_segmented_dataset_stats<-file.path(path_to_segmentation_folder,"segmented_dataset_stats.gpkg") # Path to the object segmented datasets with zonal statistics

# Step 9
path_to_classification_folder<-file.path(path_to_processing_folder,"Classification") # Path to the fodler of output of the classification (including model)
path_to_output_classification_vector<-file.path(path_to_classification_folder,"classification_vector.gpkg")
# Step 10
path_to_output_classification_raster<-file.path(path_to_classification_folder,"classification_raster.tif")

## Create the output folders
directories<-list(path_to_dem_raw_folder,path_to_spot67_preprocessed_folder,path_to_sentinel2_raw_folder,path_to_dem_preprocessed_folder,path_to_sentinel2_preprocessed_folder,path_to_classification_folder,path_to_segmentation_folder)
lapply(directories, dir.create)

## Set GRASS environment and database location 
loc <- rgrass7::initGRASS(path_to_grassApplications_folder, home=getwd(), gisDbase="GRASS_TEMP", override=TRUE,mapset = "PERMANENT" )
execGRASS("g.proj",flags="c",parameters = list(proj4=proj_srs))

## Set saga environment
saga_work_env <- RSAGA::rsaga.env()

########################################################################################################################
############ Start Workflow ############
########################################################################################################################
cat("Starting workflow")
########################################################################################################################
########################################################################################################################
############ Step 1 - Downloading the SRTM tiles  ############
########################################################################################################################
########################################################################################################################
### Uses otb applications: DownloadSRTMTiles

cat("Downloading the DEM SRTM tiles corresponding to the ROI...")
system(paste0(file.path(path_to_otbApplications_folder,"otbcli_DownloadSRTMTiles")," -vl ",file.path(path_to_processing_folder,path_to_roi_vector)," -mode download -tiledir ",file.path(path_to_dem_raw_folder)))
cat("OK")
products_to_preprocess<-list.files(path_to_dem_raw_folder,full.names = T)
# Unzip the files
for (i in 1:length(products_to_preprocess)){
  unzip(products_to_preprocess[i],exdir = path_to_dem_raw_folder)
}
file.remove(products_to_preprocess)

########################################################################################################################
########################################################################################################################
############ Step 2 - Downloading ancillary data ############
########################################################################################################################
########################################################################################################################
### Uses Copernicus Scihub API

##############################################################
#### 2.1 - Download Sentinel 2 data ####
##############################################################
####### BUG TO FIX

cat("Downloading the ancillary data: Sentinel 2 product(s) ...")

httr::set_config(authenticate(user=copernicus_scihub_username, password=copernicus_scihub_password, type = "basic"))

for (i in 1:length(Sentinel2_products_uuid)){
  url_s2_tile<-paste0("https://scihub.copernicus.eu/dhus/odata/v1/Products('",Sentinel2_products_uuid[i],"')/$value")
  httr::GET(url_s2_tile,write_disk(file.path(path_to_sentinel2_raw_folder,paste0(Sentinel2_products_uuid[i],".zip"))))
  # other option : from: https://scihub.copernicus.eu/userguide/BatchScripting
  # wget --content-disposition --continue --user={USERNAME} --password={PASSWORD}"https://scihub.copernicus.eu/dhus/odata/v1/Products('22e7af63-07ad-4076-8541-f6655388dc5e')/\$value"
  #system(paste0("wget --content-disposition --continue --user=",copernicus_scihub_username," --password=",copernicus_scihub_password," --directory-prefix ",path_to_sentinel2_raw_folder," \"https://scihub.copernicus.eu/dhus/odata/v1/Products('",Sentinel2_products_uuid[i],"')/\\$value\\"))
}

cat("Downloading the ancillary data: Sentinel 2 product(s) OK")


########################################################################################################################
########################################################################################################################
############ Step 3 - Pre-processing the Spot6/7 products ############
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
  #### 3.1 - fusionning the tiles of the panchromatic image ####
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
  #### 3.2 - converting multispectral and panchromatic images from digital numbers to TOA reflectance ####
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
  #### 3.3 - orthorectifying multispectral and panchromatic images ####
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
  #### 3.4 - extracting the ROI ####
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
  #### 3.5 - pansharpening ####
  ##############################################################
  
  PANSHARPEN_roi_path<-file.path(spot67_preprocessing_output_folder_path,"PANSHARPEN.TIF")
  
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BundleToPerfectSensor")," -inp ",PAN_roi_path," -inxs ",MS_roi_path," -out ",PANSHARPEN_roi_path," uint16")
  system(otb_appli)

  ### Remove temporary files and folders
  file.remove(setdiff(list.files(spot67_preprocessing_output_folder_path,full.names = T),c(PAN_roi_path,PANSHARPEN_roi_path,gsub(".TIF",".geom",PAN_roi_path),gsub(".TIF",".geom",PANSHARPEN_roi_path))))
  folder_to_remove=list.files(products_to_preprocess[i],full.names = T)[!grepl('.tar', list.files(products_to_preprocess[i],full.names = T))]
  for (i in 1:length(folder_to_remove)){
    system(paste0("rm -r ", folder_to_remove[i]))  
    }
}

##############################################################
#### 3.6 - mosaicing the various tiles if relevant ####
##############################################################

if(length(products_to_preprocess)>1){
  #dir.create(file.path(path_to_spot67_preprocessed_folder,"mosaic"))
  all_preprocessed_data<-list.files(path_to_spot67_preprocessed_folder,recursive = T,full.names = T)
  
  PAN_to_mosaic_paths=rev(all_preprocessed_data[grepl('PAN.TIF', all_preprocessed_data)])
  PANSHARPEN_to_mosaic_paths=rev(all_preprocessed_data[grepl('PANSHARPEN.TIF', all_preprocessed_data)])
  
  # We perform a very simple mosaic, by simply copying the last image over earlier ones in areas of overlap. We could perform much more advanced mosaicing operations with the otbcli_Mosaic application (for additional details: https://github.com/remicres/otb-mosaic)
  cat("Mosaicing the multiple images ...")
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(PAN_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",path_to_spot67_preprocessed_pan," uint16")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(PANSHARPEN_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",path_to_spot67_preprocessed_pansharpen," uint16")
  system(otb_appli)
  for (i in 1:length(products_to_preprocess)){
    system(paste0("rm -r ", products_to_preprocess[i]))  
  }
  cat("Done")
} else {
  file.copy(list.files(spot67_preprocessing_output_folder_path,full.names = T),path_to_spot67_preprocessed_folder)
  system(paste0("rm -r ", spot67_preprocessing_output_folder_path)) 
}

cat("Pre-processing the Spot6/7 products OK")


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
products_to_preprocess<-list.files(path_to_dem_raw_folder,pattern = ".hgt",full.names = T)
# If there are multiple tiles, mosaic them and then extract the ROI, else only extract the ROI
if(length(products_to_preprocess)>1){
  res<-mosaic_and_extract_roi(products_to_preprocess,dem_output_path_file,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  } else {
  res<-extract_roi(products_to_preprocess,dem_output_path_file,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_raw_folder)
  }
cat(res)

# Convert from EPSG 4326 (default SRTM EPSG) to UTM EPSG
gdalUtils::gdalwarp(srcfile=dem_output_path_file,dstfile=gsub("DEM.tif","DEM_temp.tif",dem_output_path_file),t_srs=proj_srs,overwrite=TRUE)

file.remove(dem_output_path_file)
file.rename(gsub("DEM.tif","DEM_temp.tif",dem_output_path_file),dem_output_path_file)

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

##############################################################
#### 5.1 - extract indices from the DEM : slope, aspect, flow accumulation, flow direction, topographic convergence index ####
##############################################################
## We use GRASS, calling it in R using the "rgrass7" package. We use two GRASS applications: r.slope.aspect and r.terraflow . Grass must be installed on the computer.

# Set output paths
dem_depressionless_output_path<-file.path(path_to_dem_preprocessed_folder,"DEM_depressionless.tif")
slope_output_path<-file.path(path_to_dem_preprocessed_folder,"slope.tif")
aspect_output_path<-file.path(path_to_dem_preprocessed_folder,"aspect.tif")
accumulation_output_path<-file.path(path_to_dem_preprocessed_folder,"accumulation.tif")
twi_output_path<-file.path(path_to_dem_preprocessed_folder,"twi.tif") 

# Import DEM to GRASS and set region
execGRASS("r.external", flags="o", parameters=list(input=dem_output_path_file, output="tmprast",band=1))
execGRASS("g.region", parameters=list(raster="tmprast")) 

# Filters and generates a depressionless elevation map
execGRASS("r.fill.dir", flags="overwrite", parameters=list(input="tmprast", output="DEM",direction="dir"))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="DEM", output=dem_depressionless_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))

# Compute slope and aspect and save to disk
execGRASS("r.slope.aspect", flags="overwrite", parameters=list(elevation="DEM", slope="slope",aspect="aspect",format="percent", precision="FCELL",zscale=1,min_slope=0))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="slope", output=slope_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="aspect", output=aspect_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))

# Compute hydrograpy indices and save to disk
execGRASS("r.terraflow", flags="overwrite", parameters=list(elevation="DEM", direction="direction",accumulation="accumulation",tci="tci"))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="accumulation", output=accumulation_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))

# Compute TWI indice
execGRASS("r.topidx", flags="overwrite", parameters=list(input="DEM", output="twi"))
execGRASS("r.out.gdal", flags=c("t","m","overwrite"), parameters=list(input="twi", output=twi_output_path, format="GTiff",  createopt="TFW=YES,COMPRESS=LZW" ))


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
rasters_to_disaggregate<-c(dem_depressionless_output_path,slope_output_path,aspect_output_path,accumulation_output_path,twi_output_path)
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
  
  cat(paste0("Computing simple texture indices for radius = ",xrad[i]))
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
  cat(paste0("Computing advanced texture indices for radius = ",xrad[i]))
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
ndwi<-(b03-b08)/(b03+b08)
writeRaster(ndwi,file.path(path_to_sentinel2_preprocessed_folder,"NDWI.TIF"))

# Compute MNDWI
mndwi<-(b03-b11)/(b03+b11)
writeRaster(mndwi,file.path(path_to_sentinel2_preprocessed_folder,"MNDWI.TIF"))


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
# An R function to polygonize (uses also gdal_polygonize.py) : https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
#path_to_output_segmentation_vector=file.path(path_to_outputFiles_segmentation_folder,paste0(outputImageName,"_segmented_final.gpkg"))
gdal_appli<-paste0("gdal_polygonize.py ",gsub(".gpkg",".TIF",path_to_segmented_dataset)," ",path_to_segmented_dataset," -b 1 None DN")
system(gdal_appli)

cat("Segmenting the Spot6/7 OK")

########################################################################################################################
########################################################################################################################
############ Step 7 - Extract zonal statistics for each band / indice   ############
########################################################################################################################
########################################################################################################################
### Uses saga application: Grid Statistics for Polygons (shapes_grid), shape statistics (shapes_polygons)
# Note: we tested with Grass, with the application 'v.rast.stats'. The results are equivalent to what Saga outputs, but the processing time is much longer... so we kept the Saga solution 
### Uses grass application: v.distance, v.in.ogr, v.db.addcolumn
# More info on the use of RSAGA here: https://www.r-bloggers.com/rsaga-getting-started/ and here: https://cran.r-project.org/web/packages/RSAGA/vignettes/RSAGA.html
# Saga must be installed in the computer to use RSAGA. To install on Unbuntu : sudo apt install saga

# raster data must be converted to .sgrd format to be used in saga
for (i in 1:length(indices_for_classif_paths)){
  cat(paste0("converting raster ",indices_for_classif_paths[i],"\n"))
  rsaga.import.gdal(indices_for_classif_paths[i])
}

### Function to compute the zonal statistics. It will be used to compute the zonal statistics of the ground truth dataset + the segmented objects
function_compute_zonal_statistics<-function(path_to_input_gpkg_file,path_to_output_gpkg_file,indices_for_classif_paths,indices_for_classif_labels){

# First create an shp version of the geopackage datasets (ground truth and segmentation) since SAGA does not deal with geopackages
path_to_stats_shp<-gsub(".gpkg",".shp",path_to_output_gpkg_file)
system(paste0("ogr2ogr ",path_to_stats_shp," ",path_to_input_gpkg_file))


# Compute zonal statistics using the Saga Grid Statistics for Polygons function
# to see the function parameters : rsaga.get.usage('shapes_grid', 2) 
for (i in 1:length(indices_for_classif_paths)){
rsaga.geoprocessor('shapes_grid', module = 2, env = saga_work_env, param = list(
  PARALLELIZED = 'true',
  GRIDS = paste(gsub("TIF|tif","sgrd",indices_for_classif_paths[i]),collapse=';'), # Note that in theory we could integrate multiple raster at the same time, e.g. GRIDS = "/home/ptaconet/Documents/react/data_BF/VHR_SPOT6/processed_data/NDWI.sgrd;/home/ptaconet/Documents/react/data_BF/VHR_SPOT6/processed_data/NDVI.sgrd". However all the grids must have the same characteristics (extent), which is not our case
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


# Open and rename columns (names have been automatically cut since shp does not accept long column names)
input_data<-sf::st_read(path_to_input_gpkg_file)
output_stats<-sf::st_read(path_to_stats_shp)
cols_to_rename<-setdiff(colnames(output_stats),c("geometry","landcover_",colnames(input_data)))
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

colnames(output_stats)[which(colnames(output_stats)=="landcover_")]<-"landcover_class"

# Write data as geopackage
sf::st_write(output_stats,path_to_output_gpkg_file,layer_options = "OVERWRITE=true")


## Compute shape statistics using reock indexes. It uses a function extracted from https://github.com/gerrymandr/compactr/blob/master/compactness.R
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

## Compute distance statistics (distance to hydrographic network using the accumulation dataset derived from the DEM)
execGRASS("v.in.ogr", flags=c("o","overwrite"), parameters=list(input=path_to_output_accumulation_threshold, output="accumulation",min_area=0.0001, snap=-1.0))
execGRASS("v.in.ogr", flags=c("o","overwrite"), parameters=list(input=path_to_output_gpkg_file, output="ground_truth",min_area=0.0001, snap=-1.0))
execGRASS("v.db.addcolumn", parameters=list(map="ground_truth", columns="dist_to_hydro double"))
execGRASS("v.distance", flags=c("overwrite"), parameters=list(from="ground_truth", from_type="point,line,area", to="accumulation",to_type="point,line,area",dmax=-1,dmin=-1,upload="dist",column="dist_to_hydro",output="gt_stats_updated"))

# Save file as geopackage
writeOGR(readVECT("ground_truth"),path_to_output_gpkg_file,driver = "GPKG",layer="ground_truth_stats",overwrite_layer = TRUE)
#write.csv(as.data.frame(readVECT("ground_truth")),gsub(".gpkg",".csv",path_to_output_gpkg_file))

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


## Remove saga files (very big files...)
file.remove(list.files(path_to_processing_folder,pattern = ".mgrd|.sdat|.sgrd", recursive = TRUE))

########################################################################################################################
########################################################################################################################
############ Step 8 - Prepare classification  ############
########################################################################################################################
########################################################################################################################
### For the hierarchical classification, uses the package HieRanFor developed by Dr Yoni Gavish. The package is neither available on a Git repository nor on the CRAN. It was downloaded from https://ars.els-cdn.com/content/image/1-s2.0-S0924271617303696-mmc5.zip and manually installed 
## This package uses the R package RandomForest to classify the dataset using Random Forest classifiers. 
## Additional info and a use cases here : Gavish et al., Comparing the performance of flat and hierarchical Habitat/Land-Cover classification models in a NATURA 2000 site, ISPRS Journal of Photogrammetry and Remote Sensing. Volume 136, February 2018, Pages 1-12 https://doi.org/10.1016/j.isprsjprs.2017.12.002

### Data preparation

# To ensure that simulations or random objects can be reproduced (more info: http://rfunction.com/archives/62). It must be called each time a new simulation or random object is run
set.seed(158)

# Read the ground truth dataset with zonal statistics 
ground_truth_df<-as.data.frame(sf::st_read(path_to_ground_truth_stats))
# Remove useless columns
ground_truth_df$X=NULL
ground_truth_df$landcover_class=NULL
ground_truth_df$geom=NULL

## Set dataframe of primitives types (HSR, VHSR, ancillary), sources (reflectance, spectral_indice, ancillary, texture), stat (average,stddev,other)
# Not necessary for the classification but interesting for additional information on the classification (variable importance by type/source/stat, etc.)
column_names_primitives<-setdiff(colnames(ground_truth_df),c(column_names_lc_classes_hierarchy,"cat"))
df_primitives_types_sources<-as.data.frame(column_names_primitives,stringsAsFactors=FALSE)

df_primitives_types_sources$type<-NA
df_primitives_types_sources$source<-NA
df_primitives_types_sources$stat<-NA

df_primitives_types_sources$type[which(grepl("S2",df_primitives_types_sources$column_names_primitives))]<-"HSR"
df_primitives_types_sources$type[which(grepl("SPOT6",df_primitives_types_sources$column_names_primitives))]<-"VHSR"
df_primitives_types_sources$type[which(grepl("text",df_primitives_types_sources$column_names_primitives))]<-"VHSR"
df_primitives_types_sources$type[which(is.na(df_primitives_types_sources$type))]<-"ancillary"

df_primitives_types_sources$source[which(grepl("text",df_primitives_types_sources$column_names_primitives))]<-"texture"
df_primitives_types_sources$source[which(grepl("NDVI|NDWI|BI|BRI|MNDVI|MNDWI|RNDVI",df_primitives_types_sources$column_names_primitives))]<-"spectral_indice"
df_primitives_types_sources$source[which(grepl("shape|dist_to_hydro|DEM|slope|accumulation",df_primitives_types_sources$column_names_primitives))]<-"ancillary"
df_primitives_types_sources$source[which(is.na(df_primitives_types_sources$source))]<-"reflectance"

df_primitives_types_sources$stat[which(grepl("avg",df_primitives_types_sources$column_names_primitives))]<-"average"
df_primitives_types_sources$stat[which(grepl("stddev",df_primitives_types_sources$column_names_primitives))]<-"stddev"
df_primitives_types_sources$stat[which(is.na(df_primitives_types_sources$stat))]<-"absolute"

pattern<-strsplit(methods_to_compute, split=',')[[1]]
pattern<-paste(pattern,collapse = '|_')
pattern<-paste0("_",pattern)
df_primitives_types_sources$indice<-gsub(pattern,"",df_primitives_types_sources$column_names_primitives)
df_primitives_types_sources<-left_join(df_primitives_types_sources,data.frame(indices_for_classif_labels,indices_for_classif_paths,stringsAsFactors = F),by=c("indice"="indices_for_classif_labels"))

# Keep only useful columns in the GT dataset (i.e. hierarchy of columns to classify + primitives + unique identifiers)
ground_truth_df_model<-ground_truth_df[,c("cat",column_names_lc_classes_hierarchy,column_names_primitives)]

# Fill the NAs in the ground truth dataset using the na.roughfix function. We do it this way since there are only few NAs in our datasets. Note: we could also use the randomForest::rfImpute function
ground_truth_df_model <- randomForest::na.roughfix(ground_truth_df_model)

## Get the dataset used for the hierarchical classif
ground_truth_df_hie_model<-ground_truth_df_model

## Get the dataset used for the flat classif
df_fill_missing_class_values <- ground_truth_df_hie_model[,column_names_lc_classes_hierarchy]
df_fill_missing_class_values <- df_fill_missing_class_values %>% mutate_all(as.character)

for (i in 1:ncol(df_fill_missing_class_values)-1){
  index_end_path<- which(df_fill_missing_class_values[,i+1]=="END.PATH")
  df_fill_missing_class_values[index_end_path,i+1]=df_fill_missing_class_values[index_end_path,i]
}

ground_truth_df_flat_model<-cbind(df_fill_missing_class_values,ground_truth_df_model[,column_names_primitives])


## To remove the variables linked to MS VHSR (but not the textures) :
#column_names_primitives<-df_primitives_types_sources$column_names_primitives[which((df_primitives_types_sources$type %in% c("ancillary","HSR")) | (df_primitives_types_sources$type=="VHSR" & df_primitives_types_sources$source=="texture"))]
## To remove the variables linked to MS VHSR ( :
#column_names_primitives<-df_primitives_types_sources$column_names_primitives[which(!(df_primitives_types_sources$type=="VHSR"))]
## To remove the variables linked to MS HSR :
#column_names_primitives<-df_primitives_types_sources$column_names_primitives[which(!(df_primitives_types_sources$type=="HSR"))]

#ground_truth_df_hie_model<-ground_truth_df_hie_model[,c("cat",column_names_lc_classes_hierarchy,column_names_primitives)]
#ground_truth_df_flat_model<-ground_truth_df_flat_model[,c(column_names_lc_classes_hierarchy,column_names_primitives)]

##############################################################
#### 8.1 - Generate the RF classifiers at each class hierarchical level using i) a flat approach, ii) a hierarchical approach, and compare the results ####
##############################################################

## Iterate on each class hierarchical level to produce for each classification level :
  # - a) a figure with the confusion matrix + the variable importance dotchart when using a flat classification
  # - b) a figure comparing the overall classification performance (Kappa, Accuracy, etc) when using a flat vs. a hierarchical classification

classification_perf_flat_vs_hie<-NULL

for (i in 1:length(column_names_lc_classes_hierarchy)){
  
  cat(paste0("Generating the classifiers and statistics al level ",column_names_lc_classes_hierarchy[i],"\n"))
  
  ### a) figure with the confusion matrix + the variable importance dotchart when using a flat classification
  
  ground_truth_df_flat_model_this_c<-ground_truth_df_flat_model[,c(column_names_lc_classes_hierarchy[i],column_names_primitives)]
  colnames(ground_truth_df_flat_model_this_c)[1]<-"response"
  ground_truth_df_flat_model_this_c$response<-as.factor(ground_truth_df_flat_model_this_c$response)
  
  ## Run RF flat classif
  model_tuned<-randomForest::tuneRF(x=ground_truth_df_flat_model_this_c[,2:ncol(ground_truth_df_flat_model_this_c)],y=ground_truth_df_flat_model_this_c$response,trace = FALSE)
  optimum_mtry<-as.numeric(model_tuned[,1][which(model_tuned[,2]==min(model_tuned[,2]))])
  model_flat_classif<-randomForest::randomForest(response ~ ., data=ground_truth_df_flat_model_this_c,mtry=optimum_mtry)
  conf <- caret::confusionMatrix(data = model_flat_classif$predicted, reference = ground_truth_df_flat_model_this_c$response)
  
  ## Get outputs of the flat classif 
  # Confusion matrix for flat classif
  confusion_matrix_flat<-conf$table
  # Overall stats for flat classif (Kappa, accuracy, etc)
  classif_overall_stats_flat<-as.data.frame(conf$overall)
  classif_overall_stats_flat$variable<-as.factor(rownames(classif_overall_stats_flat))
  colnames(classif_overall_stats_flat)=c("value","variable")
  # Stats by class for flat classif (f1 score, etc)
  classif_stats_by_class_flat<-conf$byClass
  classif_stats_by_class_flat<-as.data.frame(classif_stats_by_class_flat)
  if (ncol(classif_stats_by_class_flat)>1){
  classif_stats_by_class_flat$class<-sub('.*\\: ', '', row.names(classif_stats_by_class_flat))
  classif_stats_by_class_flat_f1<-classif_stats_by_class_flat[,c("class","F1")]
  classif_stats_by_class_flat_f1$classif_method<-"flat"
  }
  
  # Get variable importance
  var_importance<-as.data.frame(model_flat_classif$importance)
  var_importance$column_names_primitives<-rownames(var_importance)
  var_importance<-merge(var_importance,df_primitives_types_sources)
  var_importance$type<-as.factor(var_importance$type)
  var_importance$source<-as.factor(var_importance$source)
  var_importance$stat<-as.factor(var_importance$stat)
  var_importance$indice<-as.factor(var_importance$indice)
  
  
  plot_classif_info(conf_matrix=model_flat_classif$confusion,
                    var_importance=var_importance,
                    classifier_name=column_names_lc_classes_hierarchy[i],
                    classif_overall_stats=classif_overall_stats_flat,
                    path_to_output_plot=file.path(path_to_classification_folder,paste0("flat_classif_stats_",column_names_lc_classes_hierarchy[i],".png")))
  
  # b) a figure comparing the overall classification performance (Kappa, Accuracy, etc) when using a flat vs. a hierarchical classification
  
  ## Run RF hierarchical classif
  if (i>1){   # Works only if there is more than one hierarchy level
### Setup the classifiers at each hierarchical level.
column_names_lc_classes_hierarchy_classif<-column_names_lc_classes_hierarchy[1:which(column_names_lc_classes_hierarchy==column_names_lc_classes_hierarchy[i])]
ground_truth_df_hierarch_model_this_c<-ground_truth_df_hie_model[,c("cat",column_names_lc_classes_hierarchy_classif,column_names_primitives)]

if ("END.PATH" %in% unique(ground_truth_df_hierarch_model_this_c[,column_names_lc_classes_hierarchy[i]])) {
  internal.end.path = TRUE
} else {
  internal.end.path = FALSE
}


#Info on the function : help(RunHRF)
hie.RF <- HieRanFor::RunHRF(train.data = ground_truth_df_hierarch_model_this_c,
                            case.ID = "cat",
                            exp.var = column_names_primitives,
                            hie.levels = column_names_lc_classes_hierarchy_classif,
                            mtry = "tuneRF2",
                            internal.end.path = internal.end.path)


## Important note: The source PerformanceHRF and PerformanceFlatRF functions in the HieRanFor package havz a bug. To fix them, 
## if you run the script online do the following (done in the : 
# path_to_PerformanceHRF_modif<-file.path(path_to_processing_folder,"PerformanceHRF_modif.R")
# download.file("https://raw.githubusercontent.com/ptaconet/r_react/master/functions/PerformanceHRF_modif.R",path_to_PerformanceHRF_modif)
# insertSource(path_to_PerformanceHRF_modif, package = "HieRanFor")
# file.remove(path_to_PerformanceHRF_modif)
## If you run the script offline do the following : 
# run: fix("PerformanceHRF")
# replace the following line : nodes.acc.ind <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy")
# by :                         nodes.acc.ind <- c("Sensitivity", "Specificity","Precision", "Recall", "F1", "Pos Pred Value", "Neg Pred Value", "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy")
# Same for PerformanceFlatRF
# This trick must be done each time the script is sourced, i.e. each time the workspace is cleaned... (since no trick has been found -for now- to permanently save these changes )
# help(PerformanceHRF) for additional info on the function
perf.hRF <- PerformanceHRF (hie.RF = hie.RF,
                            per.index = c("flat.measures","hie.F.measure") , # Available options are:  c("flat.measures") ; c("hie.F.measure")
                            crisp.rule  =  c("stepwise.majority",  # Available options are:  c("stepwise.majority") ; c("multiplicative.majority") ; c("multiplicative.permutation")
                                             "multiplicative.majority"),
                            perm.num    = 10,
                            div.print   = 2  
)

crisp.case.class <- perf.hRF$crisp.case.class
hie.performance <- perf.hRF$hie.performance



## Compare flat and hierarchical classification performance measures

# Merge flat and hierarchical classification indicators
hie.performance<-as.data.frame(t(hie.performance),stringsAsFactors = F)
colnames(hie.performance)<-c("hier_stepwise_majority","hier_multiplicative_majority")
hie.performance<-hie.performance[2:nrow(hie.performance),]
hie.performance$hier_stepwise_majority<-as.numeric(hie.performance$hier_stepwise_majority)
hie.performance$hier_multiplicative_majority<-as.numeric(hie.performance$hier_multiplicative_majority)
hie.performance$variable=rownames(hie.performance)

hie.performance$class<-gsub(";.*$", '', hie.performance$variable)
hie.performance$perf_indic<-sub('.*\\;', '', hie.performance$variable)

hie.performance.f1<-hie.performance %>% filter(perf_indic=="F1")
hie.performance.f1$variable<- hie.performance.f1$perf_indic<- NULL
hie.performance.f1<-melt(hie.performance.f1,"class")
colnames(hie.performance.f1)<-c("class","classif_method","F1")

f1_scores_this_level<-rbind(classif_stats_by_class_flat_f1,hie.performance.f1)

## Compaare flat and hierarchical classifs
# f1 scores by class
png(file.path(path_to_classification_folder,paste0("f1_scores_",column_names_lc_classes_hierarchy[i],".png")),width = 1000, height = 1000, units = "px" )
print(ggplot(f1_scores_this_level, aes(fill=classif_method, y=F1, x=reorder(class,-F1))) + 
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Class", title = paste0("F1 scores by class and by classification method for ",column_names_lc_classes_hierarchy[i])) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)))
dev.off()

# Global classification performance indicatorrs (kappa, etc.)
classification_perf_flat_vs_hie_this_level<-merge(classif_overall_stats_flat,hie.performance,by="variable")
colnames(classification_perf_flat_vs_hie_this_level)[which(colnames(classification_perf_flat_vs_hie_this_level)=="value")]="flat"
classification_perf_flat_vs_hie_this_level<-melt(classification_perf_flat_vs_hie_this_level)
colnames(classification_perf_flat_vs_hie_this_level)[1]<-"perf_indicator"
colnames(classification_perf_flat_vs_hie_this_level)[2]<-"classif_method"
classification_perf_flat_vs_hie_this_level$classif_level<-column_names_lc_classes_hierarchy[i]
classification_perf_flat_vs_hie_this_level$color[classification_perf_flat_vs_hie_this_level[,2]=="flat"] <- "red"
classification_perf_flat_vs_hie_this_level$color[classification_perf_flat_vs_hie_this_level[,2]=="hier_stepwise_majority"]  <- "blue"
classification_perf_flat_vs_hie_this_level$color[classification_perf_flat_vs_hie_this_level[,2]=="hier_multiplicative_majority"]  <- "forestgreen"

#png(file.path(path_to_classification_folder,paste0("classification_perf_flat_vs_hie_",column_names_lc_classes_hierarchy[i],".png")),width = 1000, height = 1000, units = "px" )
#dotchart(classification_perf_flat_vs_hie_this_level$value,labels=classification_perf_flat_vs_hie_this_level[,2],cex=1.4,pch = 21,groups= classification_perf_flat_vs_hie_this_level[,1],main=paste0("Classification performance measures by type of classification\n at level : ",column_names_lc_classes_hierarchy[i]),xlab="Measure value", col = classification_perf_flat_vs_hie_this_level$color,xlim = c(0,1)) 
#dev.off()

classification_perf_flat_vs_hie<-rbind(classification_perf_flat_vs_hie,classification_perf_flat_vs_hie_this_level)

  }

  
}

# Plot evolution of classification performance indicators in function of the nomenclature level and method used
classif_perf_indicators<-unique(classification_perf_flat_vs_hie$perf_indicator)

list_plots<-list()
for (i in 1:length(classif_perf_indicators)){
  
  df_perf_indicators<-classification_perf_flat_vs_hie %>% filter(perf_indicator==classif_perf_indicators[i])
  p <- ggplot(df_perf_indicators, aes(x=classif_level, y=value, col=variable,group=variable)) + theme(legend.position = "none") + geom_line() + geom_point(size=2) + labs(x = "Nomenclature level", title = classif_perf_indicators[i])
  # Get plot
  list_plots[[i]]<-p
  # Get legend
  if (i == length(classif_perf_indicators)){
    p <- ggplot(df_perf_indicators, aes(x=classif_level, y=value, col=variable,group=variable)) + geom_line() + geom_point(size=2)
    legend <- get_legend(p)
    list_plots[[i+1]]<-legend
  }
  
}

png(file.path(path_to_classification_folder,"evolution_perf_indices.png"),width = 1000, height = 1000, units = "px" )
gridExtra::grid.arrange(list_plots[[1]], list_plots[[6]] , list_plots[[i+1]] , ncol = 2, nrow = 2 )
dev.off()

##############################################################
#### 8.2 - Get useful information on the classification (discriminant variables, etc.) considering the class hierarchical structure ####
##############################################################

### Generate a hierarchical classification on the whole class hierarchy 
cat("Generating useful information on the classification (discriminant variables, etc.)\n")

hie.RF <- HieRanFor::RunHRF(train.data = ground_truth_df_hie_model,
                            case.ID = "cat",
                            exp.var = column_names_primitives,
                            hie.levels = column_names_lc_classes_hierarchy,
                            internal.end.path = TRUE)


##Get and save info about the classif
### 1) A png of the hierarchy between the classes
png(file.path(path_to_classification_folder,"classes_hierarchy.png"),width = 1000, height = 1000, units = "px" )
plot(x = hie.RF, text.size = 9, split.text = 10)
dev.off()

### 2) Plots of the variable importance for each classifier, sorted by type, source and stat indicator
Importance.hie.RF <- ImportanceHie(hie.RF = hie.RF, format.out = c("col.4.out"))
Importance.hie.RF<-merge(Importance.hie.RF,df_primitives_types_sources,by.x="expl.var",by.y="column_names_primitives")
Importance.hie.RF$type<-as.factor(Importance.hie.RF$type)
Importance.hie.RF$source<-as.factor(Importance.hie.RF$source)
Importance.hie.RF$stat<-as.factor(Importance.hie.RF$stat)
Importance.hie.RF$indice<-as.factor(Importance.hie.RF$indice)

for (i in 5:7){
  for (j in c("Bubble","Tile")){
    png(file.path(path_to_classification_folder,paste0("var_importance_",names(Importance.hie.RF[i]),"_",j,".png")),width = 1000, height = 1000, units = "px" )
        PlotImportanceHie_v_taconet(input.data = Importance.hie.RF,
                                    X.data     = 2,
                                    Y.data     = 1,
                                    imp.data   = 4,
                                    explanatory.variable = i,
                                    plot.type  = j,
                                    supp.warn  = FALSE)
    dev.off()
  }
}
 

### 3) For each local classifier :  
  # a png with 2 infos : the local confusion matrix and the plot of variable importance 
  # the confusion matrix
for (i in 1:length(hie.RF$all.local.RF)){
  
  # Name of the class
  classifier.name_this_class<-hie.RF$hier.struc$lRF.info$par.name[i]
  classifier.ID_this_class<-hie.RF$hier.struc$lRF.info$classifier.ID[i]
  
  conf_matrix<-hie.RF$all.local.RF[[i]]$local.RF$confusion

  Importance.hie.RF_this_class<-Importance.hie.RF[which(Importance.hie.RF$classifier.ID==classifier.ID_this_class),]
  colnames(Importance.hie.RF_this_class)[which(colnames(Importance.hie.RF_this_class)=="mean.dec.accu")]<-"MeanDecreaseGini"
  colnames(Importance.hie.RF_this_class)[which(colnames(Importance.hie.RF_this_class)=="expl.var")]<-"column_names_primitives"
  
  plot_classif_info(conf_matrix,
                    Importance.hie.RF_this_class,
                    classifier_name = classifier.name_this_class,
                    classifier_ID = classifier.ID_this_class,
                    path_to_output_plot=file.path(path_to_classification_folder,paste0("hierar_classif_stats_",
                                                                                       classifier.ID_this_class,".png")))
  

}


########################################################################################################################
########################################################################################################################
############ Step 9 - Classify  ############
########################################################################################################################
########################################################################################################################

##############################################################
#### 9.1 - Classify the objects output of the segmentation using the approach that gives the best results, and save outputs to disk  ####
##############################################################

# Open objects segmentation dataset
segmentation_gpkg<-sf::st_read(path_to_segmented_dataset_stats)
segmentation_df<-as.data.frame(segmentation_gpkg)
segmentation_df$geom<-NULL
segmentation_df$cat=NULL
# Fill-in NA values with median value of the column (using randomForest::na.roughfix)
segmentation_df<-randomForest::na.roughfix(segmentation_df)

# Run the model that gives the best results at all the hierarchical levels (here for BF: flat at all the levels)

for (i in 1:length(column_names_lc_classes_hierarchy)){
  cat(paste0("Classifying the objects at level ",column_names_lc_classes_hierarchy[i],"\n"))
ground_truth_df_flat_model_this_class<-ground_truth_df_flat_model[,c(column_names_lc_classes_hierarchy[i],column_names_primitives)]
colnames(ground_truth_df_flat_model_this_class)[1]<-"response"
ground_truth_df_flat_model_this_class$response<-as.factor(ground_truth_df_flat_model_this_class$response)

model_tuned<-randomForest::tuneRF(x=ground_truth_df_flat_model_this_class[,2:ncol(ground_truth_df_flat_model_this_class)],y=ground_truth_df_flat_model_this_class$response,trace = FALSE)
optimum_mtry<-as.numeric(model_tuned[,1][which(model_tuned[,2]==min(model_tuned[,2]))])
model<-randomForest::randomForest(response ~ ., data=ground_truth_df_flat_model_this_class,mtry=optimum_mtry)

segmentation_df$predicted<-predict(model,segmentation_df)

cat("Saving outputs to the disk")
res<-save_classif_to_disk_function(dataset_classified = segmentation_df, 
                                   dataset_to_classify_sf = segmentation_gpkg, 
                                   path_to_output_classification_vector = file.path(path_to_classification_folder,paste0("classification_",column_names_lc_classes_hierarchy[i],".gpkg")), 
                                   save_objects_raster=FALSE, # turn to TRUE when script is OK
                                   save_objects_vector_group_adj_polygons=FALSE # turn to TRUE when script is OK
                                    )


}


### Here under is for a hierarchical classification
# Predict land cover class using the predict HieRanFor::PredictNewHRF function with the classifiers extracted in step 8
# Info on the function : help(PredictNewHRF)
# Outputs : 
# pred.segmentation.hRF$raw.votes : Data frame containing for each case, the proportion of votes for each node in each local classifier. 
# ppred.segmentation.hRF$crisp.case.class : Data frame containing the crisp class for each case based on all options defined by crisp.rule.

#pred.segmentation.hRF <- PredictNewHRF(hie.RF, 
#                              segmentation_df, 
#                              new.data.exp.var=c(2:ncol(segmentation_df)), 
#                              new.data.case.ID = 1,
#                              crisp.rule = c("stepwise.majority"), # available options: c("stepwise.majority","multiplicative.majority" ,"multiplicative.permutation"). See details in help(PerformanceHRF)
#                              perm.num   = 10,
#                              div.print  = 2)


#dataset_classified<-merge(pred.segmentation.hRF$crisp.case.class,pred.segmentation.hRF$raw.votes,by="cat")
#dataset_classified$train.or.test.x <- dataset_classified$train.or.test.y <-NULL

# change column names since sql and qgis do not like dots in the column names
#if ("stepwise.majority.rule" %in% colnames(dataset_classified)) {
#  colnames(dataset_classified)[which(colnames(dataset_classified)=="stepwise.majority.rule")]<-"predicted_class"
#}

########################################################################################################################
########################################################################################################################
############ Step 9 - Post-processing  ############
########################################################################################################################
########################################################################################################################


##############################################################
#### 9.2 - Enhance the classification with user-defined rules ####
##############################################################

## TODO

########################################################################################################################
############ Close the workflow ############
########################################################################################################################

# Remove grass temporary folder
system(paste0("rm -r ", file.path(getwd(),"GRASS_TEMP")))  
file.remove(file.path(getwd(),".grassrc7"))


########################################################################################################################
############ The end ############
########################################################################################################################
cat("End workflow")








########################### OLD TRIES not working but I keep them because they could be reused #########################

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
#for (i in 1:length(column_names_lc_classes_hierarchy)){
#gt_lc_types<-data.frame(lc_type=unique(ground_truth_df[,column_names_lc_classes_hierarchy[i]]))
#gt_lc_types[,paste0(column_names_lc_classes_hierarchy[i],"_int")]<-seq(1:nrow(gt_lc_types))
#ground_truth<-merge(ground_truth,gt_lc_types,by.x=column_names_lc_classes_hierarchy[i],by.y="lc_type")
#ground_truth_df<-merge(ground_truth_df,gt_lc_types,by.x=column_names_lc_classes_hierarchy[i],by.y="lc_type")
#}

#writeOGR(ground_truth,path_to_ground_truth_processed,driver = "GPKG",layer="ground_truth",overwrite_layer = TRUE)

# Get training samples ids. We keep 70% of the data for training (size=0.7) and 30% for validation
#ground_truth_training<-ground_truth_df %>% group_by_(paste0(column_names_lc_classes_hierarchy,"_int")) %>% sample_frac(size = 0.7)
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
#columns_stats<-setdiff(colnames(ground_truth_training),c(column_names_lc_classes_hierarchy,"cat","cat_","id",paste0(column_names_lc_classes_hierarchy,"_int")))
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


