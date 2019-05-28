

require(raster)
require(RSAGA)

## Set saga environment
saga_work_env <- RSAGA::rsaga.env()

path_to_pedo_geotif<-"/home/ptaconet/Documents/react/data_CIV/pedology/pedo_civ2.tif"
num_clusters_unsuperv_classif<-13  # must be set by the user after tests  #for civ1 : 10
radius_majority_filter<-7  # must be set by the user after tests          #for civ1 : 7

path_to_output_unsuperv_classif<-gsub(".tif","_unsuperv_classif.tif",path_to_pedo_geotif)
path_to_output_unsuperv_classif_filter<-gsub(".tif","_unsuperv_classif_filter.tif",path_to_pedo_geotif)

## Unsupervised classification of the image with kmeans
#from https://rspatial.org/rs/4-unsupclassification.html
pedo <- stack(path_to_pedo_geotif)
nr <- getValues(pedo)
# It is important to set the seed generator because `kmeans` initiates the centers in random locations
set.seed(99)
# We want to create num_clusters_unsuperv_classif clusters, allow 500 iterations, start with 5 random sets using "Lloyd" method
kmncluster <- kmeans(na.omit(nr), centers = num_clusters_unsuperv_classif, iter.max = 500, nstart = 5, algorithm="Lloyd")

# First create a copy of the pedo layer
knr <- pedo
# Now replace raster cell values with kmncluster$cluster array
knr[] <- kmncluster$cluster

writeRaster(knr[[1]],path_to_output_unsuperv_classif)


## Apply majority filter on the output unsupervised_classif
rsaga.import.gdal(path_to_output_unsuperv_classif)

rsaga.geoprocessor('grid_filter', module = 6, env = saga_work_env, param = list(
  INPUT = gsub(".tif",".sgrd",path_to_output_unsuperv_classif), 
  RESULT = gsub(".tif",".sgrd",path_to_output_unsuperv_classif_filter),
  MODE = 1,
  RADIUS = radius_majority_filter
))

rsaga.geoprocessor('io_gdal', module = 2, env = saga_work_env, param = list(
  GRIDS = gsub(".tif",".sgrd",path_to_output_unsuperv_classif_filter),
  FILE = path_to_output_unsuperv_classif_filter
))


## Specific for our case, because we noticed that all the pixels that are not integers (ie pixel that are real) are the rivers.
  r<-raster(path_to_output_unsuperv_classif_filter)
  # r[r%%1!=0] => looks for all the pixels that are not integers 
  r[r%%1!=0]<-14
  # for civ2 (north part): pixel = 1 and pixel = 12 are the same classes
  r[r==12]<-1
  # for civ1 (south part): pixel = 14 and pixel = 2 are the same classes
  #r[r==14]<-2
  #r[r==6]<-2
  #r[r==1]<-10
  
  writeRaster(r,path_to_output_unsuperv_classif_filter,overwrite=T)

  
  ## Set classes of civ2 as classes of civ1
  r<-raster("/home/ptaconet/Documents/react/data_CIV/pedology/pedo_civ1_processed.tif")
  r[r==2]<-1
  r[r==3]<-NA   # cette classe est trop petite
  r[r==4]<-4
  r[r==5]<-NA   # cette classe est trop petite
  r[r==6]<-6

  
  r[r==8]<-5 
  r[r==9]<-NA   
  r[r==7]<-9 
  r[r==13]<-8 
  r[r==10]<-13 
  r[r==11]<-11 

  writeRaster(r,gsub(".tif","_corrected.tif",path_to_output_unsuperv_classif_filter),overwrite=T)
  
  
  r<-raster("/home/ptaconet/Documents/react/data_CIV/pedology/pedo_civ2_processed.tif")
  r[r==14]<-1 
  writeRaster(r,"/home/ptaconet/Documents/react/data_CIV/pedology/pedo_civ2_processed.tif",overwrite=T)
  
  # concatenate the 2 rasters
  
  source("/home/ptaconet/r_react/outdated_scripts/mosaic_tif_images.R")
  
  mosaic_tif_images(c("/home/ptaconet/Documents/react/data_CIV/pedology/pedo_civ1_unsuperv_classif_filter_corrected.tif","/home/ptaconet/Documents/react/data_CIV/pedology/pedo_civ2_processed.tif"),"/home/ptaconet/Documents/react/data_CIV/pedology/pedo_civ_final.tif")
