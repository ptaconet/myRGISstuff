# R scripts developed in the frame of the REACT project   (« Gestion de la résistance aux insecticides au Burkina Faso et en Côte d’Ivoire : Recherche sur les stratégies de lutte anti-vectorielle »)

Description of the scripts :

## *object_based_image_analysis* : scripts for the processing of the satellite imagery (Spot 6/7) used for the generation of the land cover maps
 - *preprocessing* : This script is a workflow for the pre-processing SPOT6/7 satellite imagery data (tile fusionning + orthorectification + pansharpening + cloud mask extraction). The user can parameterize the operations he/she wants to be computed. The script uses applications coming from various libraries that must be installed : the Orfeo Toolbox (https://www.orfeo-toolbox.org/), R spatial packages ("sf" for vector and "raster" for raster) and the Geospatial Data Abstraction Library (https://www.gdal.org/). The user should be aware that some operations might work only on Linux OS.
 
 - *preprocessing_mosaicing* : This script performs a simple mosaicing of orthorectified TIF images using the Orfeo Toolbox otbcli_Mosaic application.
 
 - *segmentation* :  segmentation [...work in progress...]


## *outdated_scripts* : a set of scripts that are not used in the workflow
  - *preprocess_spot_6_7_data.R* : This script performs a set of pre-processing operations on SPOT6/7 satellite imagery data downloaded in the GEOSUD portal (mosaicing + pansharpening + georeferencing). The user can parameterize the operations he/she wants to be computed. Note: the script uses ssytem commands and spatial libraries that might work only on Linux OS (particularly for georeferencing).
  
  - *mosaic_tif_images.R* : Mosaic tif images that are split into many pieces.
  
  - *pansharp.R* : Pansharp. Adapted from https://www.r-bloggers.com/pan-sharpening-using-r/ . This function uses the same algorithm as the OTB Toolbox
  
  - *get_bingmaps_images/download_and_georeference_bingmaps_data.R* : Download and georeference Bing maps data covering an AOI (provided by the user). Outputs are georeferenced TIF files covering the AOI. This script uses the Bing API. More info here: https://msdn.microsoft.com/en-us/library/ff701724.aspx
    
  - *get_bingmaps_images/create_10km_square_tif_from_bingmaps_data* : Create a 10km² georeferenced TIF + a OGC Geopackage of a given AOI on Earth using Bing maps satellite imagery at best available zoom (19, which corresponds to a resolution of approx. 0.5 cm)
  