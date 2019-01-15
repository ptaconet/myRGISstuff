######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = workflow_object_based_data_analysi, title = Workflow for the land cover classification using object based image analysis, abstract = This script is a workflow for the object based image analysis work used to create the land cover maps. It performs all the pre-processing, data preparation, classification and post-processing steps. The script uses applications of the Orfeo Toolbox v 6.6.1 (https://www.orfeo-toolbox.org/) and the GRASS library (https://grass.osgeo.org) v 7.4
# wps.in: id = path_to_otbApplications_folder, type = string, title = Path to the folder containing the OTB applications. , value = "/home/ptaconet/OTB-6.6.1-Linux64/bin";
# wps.in: id = path_to_grassApplications_folder, type = string, title = Path to the folder containing the GRASS applications. , value = "/usr/lib/grass74";
# wps.in: id = path_to_processing_folder, type = string, title = Path to working folder , value = "/home/ptaconet/Documents/react/data_CIV";
# wps.in: id = path_to_roi_vector, type = string, title = Path to the ROI in kml or shp format, value = "ROI.kml";
# wps.in: id = path_to_spot67_raw_folder, type = string, title = Path to the folder where input data (i.e. Spot6/7 products as tar.gz files) are stored , value = "VHR_SPOT6/raw_data";
# wps.in: id = path_to_dem_folder, type = string, title = Path to the folder containing the SRTM Digital Elevation Model files. value = "DEM_SRTM/raw_data";
# wps.out: id = output_zip, type = text/zip, title = Folder path_to_spot67_output_preprocessed_folder containing the data pre-processed: for each input product, the orthorectified multispectral image (extracted on the ROI), the orthorectified panchromatic image (extracted on the ROI), the orthorectified pansharpened image (extracted on the ROI). In addition, if there are multiple input products, the previously mentioned products will be mosaiced 

#### Workflow :
### Step 1 - Download the SRTM tiles for the ROI
### Step 2 - Pre-process the Spot6/7 image(s) :
  ## 2.1 - fusion the tiles of the panchromatic image
  ## 2.2 - convert the multispectral and panchromatic images from digital numbers to TOA reflectance
  ## 2.3 - orthorectifying the multispectral and panchromatic images
  ## 2.4 - extract the ROI
  ## 2.5 - pansharpen
  ## 2.6 - mosaic the various tiles (if relevant)
### Step 3 - Segment the Spot6/7 image 
### Step 4 - Download the ancillary data :
  ## 4.1 - Download Sentinel 2 image(s)
### Step 5 - Preprocess the ancillary data :
  ## 5.1 - preprocess the DEM : mosaic the various tiles (if relevant), and then extract the ROI
  ## 5.2 - preprocess the Sentinel 2 image(s) : mosaic the various images (if relevant), and then extract the ROI
#### Step 6 - Prepare the ancillary data for the classification : 
## 6.1 - extract indices from the DEM : slope, aspect, flow accumulation, flow direction, topographic convergence index
## 6.3 - extract textural indices from the Spot6/7 panchromatic image : 
## 6.4 - extract radiometric indices from the Spot6/7 pansharpened image : 
## 6.5 - extract each band from the Spot6/7 pansharpened image



### Set Input parameters for the workflow

### Global variables used throughout the WF
path_to_otbApplications_folder<-"/home/ptaconet/OTB-6.6.1-Linux64/bin"
path_to_grassApplications_folder<-"/usr/lib/grass74" ## Can be retrieved with grass74 --config path . More info on the use of rgrass7 at https://grasswiki.osgeo.org/wiki/R_statistics/rgrass7
path_to_processing_folder<-"/home/ptaconet/Documents/react/data_CIV"
path_to_roi_vector="ROI.kml"

### Parameters for step 1 : Downloading the SRTM tiles
path_to_dem_folder="DEM_SRTM/raw_data"

### Parameters for step 2 : Pre-processing the Spot6/7 products
path_to_spot67_raw_folder<-"VHR_SPOT6/raw_data"

### Parameters for step 3 : Segmenting the Spot6/7 image 
segmentation_threshold=300       #<paramètre d'échelle> 
segmentation_cw=0.3              # <poids radiométrie> 
segmentation_sw=0.75              # <poids compacité>

### Parameters for step 4 : Downloading ancillary data
copernicus_scihub_username="ptaconet"
copernicus_scihub_password="****"
Sentinel2_output_folder_path<-"HR_Sentinel2/raw_data"
Sentinel2_products_id<-c("bc6bafd7-d44f-4d62-8754-3cc4ba4e8cc0","18895056-852f-4a4f-a3aa-ca7882fe79de")

### Parameters for step 5 : Preprocessing the ancillary data
## None

### Parameters for step 6 : Preparing the ancillary data for the classification
radiometric_indices_list<-c("Vegetation:NDVI","Vegetation:RVI","Vegetation:SAVI","Vegetation:IPVI","Water:NDWI","Water:MNDVI","Water:NDPI","Water:NDTI","Soil:BI","Soil:BI2")


#### Prepare workflow

### Source useful functions
all_functions<-list.files("/home/ptaconet/r_react/functions",full.names = T)
lapply(all_functions, source)

### Call useful packages
require(raster)
require(rgrass7)


### Set working directory
setwd(path_to_processing_folder)


#### Start Workflow

#### Step 1 - Downloading the SRTM tiles (since they will be used throughout the WF)
cat("Downloading the DEM SRTM tiles corresponding to the ROI...")
system(paste0(file.path(path_to_otbApplications_folder,"otbcli_DownloadSRTMTiles")," -vl ",file.path(path_to_processing_folder,path_to_roi_vector)," -mode download -tiledir ",file.path(path_to_processing_folder,path_to_dem_folder)))
cat("OK")


#### Step 2 - Pre-processing the Spot6/7 products
## 2.1 - fusionning the tiles of the panchromatic image
## 2.2 - converting multispectral and panchromatic images from digital numbers to TOA reflectance
## 2.3 - orthorectifying multispectral and panchromatic images
## 2.4 - extracting the ROI
## 2.5 - pansharpening
## 2.6 - mosaicing the various tiles if relevant

cat("Pre-processing the Spot6/7 products ")

## List the folders available in the folder (1 folder / couple of PAN + MS image)
products_to_preprocess<-dir(path_to_spot67_raw_folder,full.names = TRUE)

# Create directory where output datasets will be stored
path_to_spot67_output_preprocessed_folder=gsub("raw","processed",path_to_spot67_raw_folder)
dir.create(path_to_spot67_output_preprocessed_folder)

# Preprocess each product
for (i in 1:length(products_to_preprocess)){
  cat(paste0("Starting preprocessing of product ",products_to_preprocess[i]))
  
  spot67_preprocessing_output_folder_path=file.path(path_to_spot67_output_preprocessed_folder,dir(path_to_spot67_raw_folder,full.names = FALSE)[i])
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
  
  ## 2.1 - fusionning the tiles of the panchromatic image
  
  ## Identify folder containing the PAN tiles
  folders=list.dirs(path = PAN_folder, full.names = TRUE)
  PAN_tile_folder=folders[grepl('IMG', folders)][1]
  pan_tifs_paths=as.vector(list.files(path=PAN_tile_folder,pattern=".TIF",full.names=TRUE)) 
  ## Mosaic
  cat("Starting PAN tile fusionning ...")
  PAN_mosaic_path<-file.path(spot67_preprocessing_output_folder_path,"PAN.TIF")
  res<-mosaic(pan_tifs_paths,PAN_mosaic_path,path_to_otbApplications_folder)
  cat(res)
  
  ## 2.2 - converting multispectral and panchromatic images from digital numbers to TOC reflectance
  
  cat("Starting Convertion to TOC reflectance ...")
  ## PAN 
  PAN_toc_path<-file.path(spot67_preprocessing_output_folder_path,"PAN_TOC.TIF")
  res<-convert_dn_to_toc_reflectance(PAN_mosaic_path,PAN_toc_path,path_to_otbApplications_folder)
  cat(res)
  ## MS
  # Identify folder containing the MS tile
  folders=list.dirs(path = MS_folder, full.names = TRUE)
  MS_tile_folder=folders[grepl('IMG', folders)][1]
  ms_tif_path=as.vector(list.files(path=MS_tile_folder,pattern=".TIF",full.names=TRUE)) 
  MS_toc_path<-file.path(spot67_preprocessing_output_folder_path,"MS_TOC.TIF")
  res<-convert_dn_to_toc_reflectance(ms_tif_path,MS_toc_path,path_to_otbApplications_folder)
  cat(res)
  
  ## 2.3 - orthorectifying multispectral and panchromatic images
  
  cat("Starting orthorectification ...")
  ## PAN 
  PAN_orthorectified_path=file.path(spot67_preprocessing_output_folder_path,"PAN_ORTHO.TIF")
  res<-orthorectify_spot67(PAN_toc_path,PAN_orthorectified_path,path_to_otbApplications_folder,path_to_dem_folder)
  cat(res)
  ## MS 
  MS_orthorectified_path=file.path(spot67_preprocessing_output_folder_path,"MS_ORTHO.TIF")
  res<-orthorectify_spot67(MS_toc_path,MS_orthorectified_path,path_to_otbApplications_folder,path_to_dem_folder)
  cat(res)
  
  ## 2.4 - extracting the ROI
  
  cat("Extracting the ROI ...")
  ## PAN 
  PAN_roi_path=file.path(spot67_preprocessing_output_folder_path,"PAN.TIF")
  res<-extract_roi(PAN_orthorectified_path,PAN_roi_path,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_folder)
  cat(res)
  ## MS
  MS_roi_path=file.path(spot67_preprocessing_output_folder_path,"MS.TIF")
  res<-extract_roi(MS_orthorectified_path,MS_roi_path,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_folder)
  cat(res)
  
  ## 2.5 - pansharpening
  
  PANSHARPEN_roi_path<-file.path(spot67_preprocessing_output_folder_path,"PANSHARPEN.TIF")
  res<-pansharpen(PAN_roi_path,MS_roi_path,PANSHARPEN_roi_path,path_to_otbApplications_folder,path_to_dem_folder)
  cat(res)
    
  ### Remove temporary files and folders
  file.remove(setdiff(list.files(spot67_preprocessing_output_folder_path,full.names = T),c(PAN_roi_path,MS_roi_path,PANSHARPEN_roi_path,gsub(".TIF",".geom",PAN_roi_path),gsub(".TIF",".geom",MS_roi_path),gsub(".TIF",".geom",PANSHARPEN_roi_path))))
  folder_to_remove=list.files(products_to_preprocess[i],full.names = T)[!grepl('.tar', list.files(products_to_preprocess[i],full.names = T))]
  for (i in 1:length(folder_to_remove)){
    system(paste0("rm -r ", folder_to_remove[i]))  
    }
}

## 2.6 - mosaicing the various tiles if relevant

if(length(products_to_preprocess)>1){
  #dir.create(file.path(path_to_spot67_output_preprocessed_folder,"mosaic"))
  all_preprocessed_data<-list.files(path_to_spot67_output_preprocessed_folder,recursive = T,full.names = T)
  
  MS_to_mosaic_paths=rev(all_preprocessed_data[grepl('MS.TIF', all_preprocessed_data)])
  PAN_to_mosaic_paths=rev(all_preprocessed_data[grepl('PAN.TIF', all_preprocessed_data)])
  PANSHARPEN_to_mosaic_paths=rev(all_preprocessed_data[grepl('PANSHARPEN.TIF', all_preprocessed_data)])
  
  MS_roi_path=file.path(path_to_spot67_output_preprocessed_folder,"MS.TIF")
  PAN_roi_path=file.path(path_to_spot67_output_preprocessed_folder,"PAN.TIF")
  PANSHARPEN_roi_path=file.path(path_to_spot67_output_preprocessed_folder,"PANSHARPEN.TIF")
  
  # We perform a very simple mosaic, by simply copying the last image over earlier ones in areas of overlap. We could perform much more advanced mosaicing operations with the otbcli_Mosaic application (for additional details: https://github.com/remicres/otb-mosaic)
  cat("Mosaicing the multiple images ...")
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(MS_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",MS_roi_path," uint16")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(PAN_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",PAN_roi_path," uint16")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(PANSHARPEN_to_mosaic_paths, collapse = " ")," -comp.feather none -harmo.method none -out ",PANSHARPEN_roi_path," uint16")
  system(otb_appli)
  cat("Done")
} else {
  file.copy(list.files(spot67_preprocessing_output_folder_path,full.names = T),path_to_spot67_output_preprocessed_folder)
  system(paste0("rm -r ", spot67_preprocessing_output_folder_path)) 
}

cat("Pre-processing the Spot6/7 products OK")

#### Step 3 - Segmenting the Spot6/7 image 
####### IMPORTANT NOTE : USES AN OTB APPLICATION ONLY AVAILABLE IN THE PERSONAL RELEASE OF OTB #########
cat("Segmenting the Spot6/7 image...")

if (length(products_to_preprocess)==1){
  path_to_image_to_segment<-PANSHARPEN_roi_path
  
} else {
  path_to_image_to_segment<-PANSHARPEN_mosaiced_path
}

dir.create("Segmentation")
path_to_output_segmented_data=file.path("Segmentation","segmentation_output.gpkg")
cat("Starting segmentation...")
res<-segment(path_to_image_to_segment,path_to_output_segmented_data,segmentation_threshold,segmentation_cw,segmentation_,otbapplications_folder_path)
cat(res)

cat("Segmenting the Spot6/7 OK")


#### Step 4 - Downloading ancillary data :
## 4.1 - Download Sentinel 2 data

cat("Downloading the ancillary data: Sentinel 2 product(s) ...")

dir.create(Sentinel2_output_folder_path)
for (i in 1:length(Sentinel2_products_id)){
  res<-download_copernicus_scihub_products(copernicus_scihub_username,copernicus_scihub_password,Sentinel2_products_id[i],Sentinel2_output_folder_path)
  cat(res)
}
cat("Downloading the ancillary data: Sentinel 2 product(s) OK")


#### Step 5 - Preprocessing the ancillary data :
## 5.1 - preprocessing the DEM : mosaicing the various tiles if relevant, and then extracting the ROI
## 5.2 - preprocessing the Sentinel 2 product(s) : mosaicing the various tiles if relevant, and then extracting the ROI

cat("Pre-processing the DEM : if relevant, mosaicing the various tiles and extracting the ROI ...")

# List the products
products_to_preprocess<-list.files(path_to_dem_folder,full.names = T)
# Create directory where output datasets will be stored
path_to_dem_output_preprocessed_folder=gsub("raw","processed",path_to_dem_folder)
# Set output DEM file path
dem_output_path_file<-file.path(path_to_dem_output_preprocessed_folder,"DEM.TIF")
# If there are multiple tiles, mosaic them and then extract the ROI, else only extract the ROI
if(length(products_to_preprocess)>1){
  res<-mosaic_and_extract_roi(products_to_preprocess,path_to_dem_output_preprocessed_folder,dem_output_path_file,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_folder)
  } else {
  res<-extract_roi(products_to_preprocess,dem_output_path_file,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_folder)
  }
cat(res)



## 5.2 - preprocessing the Sentinel 2 product(s) : mosaicing the various tiles if relevant, and then extracting the ROI

cat("Pre-processing the S2 images : if relevant, mosaicing the various bands and extracting the ROI ...")

## List the products
products_to_preprocess<-dir(Sentinel2_output_folder_path,full.names = TRUE)
# Create directory where output datasets will be stored
path_to_sentinel2_output_preprocessed_folder=gsub("raw","processed",Sentinel2_output_folder_path)
dir.create(path_to_sentinel2_output_preprocessed_folder)

# Unzip all products
for (i in 1:length(products_to_preprocess)){
  unzip(products_to_preprocess[i],exdir=path_to_sentinel2_output_preprocessed_folder)
}

# If there are multiple products, mosaic them and then extract the ROI, else only extract the ROI
patterns=c("B01.jp2","B02.jp2","B03.jp2","B04.jp2","B05.jp2","B06.jp2","B07.jp2","B08.jp2","B8A.jp2","B09.jp2","B10.jp2","B11.jp2","B12.jp2")
  
for (i in 1:length(patterns)){
  paths_images_to_mosaic=list.files(path_to_sentinel2_output_preprocessed_folder,pattern = patterns[i],full.names = T,recursive = T)
  path_to_output_mosaiced_image=file.path(path_to_sentinel2_output_preprocessed_folder,gsub(".jp2",".TIF",patterns[i]))
  
  if(length(products_to_preprocess)>1){
  res<-mosaic_and_extract_roi(paths_images_to_mosaic,path_to_sentinel2_output_preprocessed_folder,path_to_output_mosaiced_image,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_folder)
  } else {
  res<-extract_roi(paths_images_to_mosaic,path_to_output_mosaiced_image,path_to_roi_vector,path_to_otbApplications_folder,path_to_dem_folder)
  }
  cat(res)
}
  

#### Step 6 - Preparing the ancillary data for the classification : 
## 6.1 - extract indices from the DEM : slope, aspect, flow accumulation, flow direction, topographic convergence index
## 6.3 - extract textural indices from the Spot6/7 image : 
## 6.4 - extract radiometric indices from the Spot6/7 image : NDVI, BRI, MNDVI, NDWI, MNDWI, RNDVI (pour Sentinel 2 uniquement)
## 6.5 - extract each band from the Spot6/7 image



## 6.1 - extract indices from the DEM : slope, aspect, flow accumulation, flow direction, topographic convergence index

## We use GRASS, calling it in R using the "rgrass7" package. We use two GRASS applications: r.slope.aspect and r.terraflow

# Set output paths
slope_output_path<-file.path(path_to_dem_output_preprocessed_folder,"slope.TIF")
aspect_output_path<-file.path(path_to_dem_output_preprocessed_folder,"aspect.TIF")
accumulation_output_path<-file.path(path_to_dem_output_preprocessed_folder,"accumulation.TIF")
direction_output_path<-file.path(path_to_dem_output_preprocessed_folder,"direction.TIF")
tci_output_path<-file.path(path_to_dem_output_preprocessed_folder,"tci.TIF") 

# Set GRASS environment and database location 
loc <- rgrass7::initGRASS(path_to_grassApplications_folder, home=getwd(), gisDbase="GRASS_TEMP", override=TRUE,mapset = "PERMANENT" )

# Import DEM to GRASS and set region
execGRASS("g.proj",flags="c",parameters = list(proj4="+proj=longlat +datum=WGS84 +no_defs"))
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

# Remove temorary folder
system(paste0("rm -r ", file.path(getwd(),"GRASS_TEMP")))  
file.remove(file.path(getwd(),".grassrc7"))

# To read and plot the output raster :
#r <- readRAST("slope")
#spplot(r)

## 6.2 - extract textural indices from the Spot6/7 image : 
####### IMPORTANT NOTE : USES AN OTB APPLICATION ONLY AVAILABLE IN THE PERSONAL RELEASE OF OTB #########

# Get image maximum and minimum
rast<-raster(PAN_roi_path)

# Set output file path
path_to_textural_indices<-file.path(path_to_processing_folder,path_to_spot67_output_preprocessed_folder,"HaralickTextures.TIF")

## Compute textures
# Using the application available in the official release (does not enable to select the set of textures)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_HaralickTextureExtraction")," -in ",file.path(path_to_processing_folder,PAN_roi_path)," -parameters.xrad 11 -parameters.yrad 11 -parameters.nbbin 32 -parameters.min ",minValue(rast)," -parameters.max ",maxValue(rast)," -texture simple -out ",path_to_textural_indices)
system(otb_appli)

# Using the application available in the personal release (enables to select the set of textures)
# indices_list<-c("enthropy","energy")
#otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SelectiveHaralickTextures")," -in ",file.path(path_to_processing_folder,PAN_roi_path)," -parameters.xrad 11 -parameters.yrad 11 -parameters.nbbin 32 -parameters.min ",minValue(rast)," -parameters.max ",maxValue(rast)," -indices ",paste(indices_list, collapse = " ")," -out ",path_to_textural_indices)
#system(otb_appli)


## 6.3 - extract radiometric indices from the Spot6/7 image : NDVI, BRI, MNDVI, NDWI, MNDWI, RNDVI (pour Sentinel 2 uniquement)

# Set output file path

for (i in 1:length(radiometric_indices_list)){
  
  indice<-sub(".*:","",radiometric_indices_list[i])
  path_to_output_indice<-file.path(path_to_processing_folder,path_to_spot67_output_preprocessed_folder,paste0(indice,".TIF"))
  
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_RadiometricIndices")," -in ",file.path(path_to_processing_folder,PANSHARPEN_roi_path)," -channels.blue 1 -channels.green 2 -channels.red 3 -channels.nir 4 -list ",radiometric_indices_list[i]," -out ",path_to_output_indice)
  system(otb_appli)
}


## 6.5 - extract each band from the Spot6/7 image

otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SplitImage")," -in ",file.path(path_to_processing_folder,PANSHARPEN_roi_path)," -out ",file.path(path_to_processing_folder,PANSHARPEN_roi_path)," uint16")
system(otb_appli)




