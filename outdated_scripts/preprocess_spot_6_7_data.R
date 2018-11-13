######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = preprocess_spot_6_7_data, title = Pre-processing of SPOT 6/7 images (mosaicing + pansharpening + georeferencing), abstract = This script performs a set of pre-processing operations on SPOT6/7 satellite imagery data downloaded in the GEOSUD portal (mosaicing + pansharpening + georeferencing). The user can parameterize the operations he/she wants to be computed. Note: the script uses ssytem commands and spatial libraries that might work only on Linux OS (particularly for georeferencing);
# wps.in: id = mosaic_panchromatic_images, type = boolean, title = Panchromatic images (PAN) are split into many pieces (i.e. TIF files). Mosaic panchromatic TIF files so that the panchromatic image is available into 1 single TIF file? , value="TRUE|FALSE";
# wps.in: id = resize_multispectral_image, type = boolean, title = The multispectral image (MS) is not georeferenced. Hence the MS resolution being coarser than the PAN one it makes PAN and MS images incomparables. Resize the MS image to reach the spatial extent of the PAN image? , value="TRUE|FALSE"; 
# wps.in: id = preprocess_cloud_mask, type = boolean, title = Some masks are included in the products delivered however there coordinates are not aligned with the images coordinates: the latitude coordinates are inverted (i.e. positive while the images one are negative). Process cloud mask so that it can be used with the images? , value="TRUE|FALSE"; 
# wps.in: id = pansharpen, type = boolean, title = Pansharpen? Pansharpening uses a function available here: https://www.r-bloggers.com/pan-sharpening-using-r/ . This function uses the same algorithm as the OTB Toolbox , value="TRUE|FALSE"; 
# wps.in: id = georeference, type = boolean, title = Images as delivered by GEOSUD are not georeferenced. georeference images using Ground Control Points (GCP)? , value="TRUE|FALSE"; 
# wps.in: id = pan_tifs_folder, type = string, title = Path to the folder where the PAN .TIF images are stored, value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A";
# wps.in: id = pan_cloudmask_path, type = string, title = Path to the cloud mask file, value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A/MASKS/CLD_SPOT7_P_201710111019265_SEN_SPOT7_20171012_1350221m2hw3hzfkke0_1_MSK.GML";
# wps.in: id = ms_tif_file, type = string, title = Path to the MS .TIF file, value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_MS_001_A/IMG_SPOT7_MS_201710111019265_SEN_SPOT7_20171012_1350131l9pjds0m4wzt_1_R1C1.TIF";
# wps.in: id = ground_control_points_path, type = string, title = Path to the Ground Control Points (GCP) file. The file must be structured as presented here: [link], value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/ground_control_points.csv";
# wps.in: id = georeference_resampling_method, type = string, title = Resampling method to use for georeferencing with GCP. , value = "near|bilinear|cubic|cubicspline|lanczos|average|mode|max|min|med|Q1|Q3"
# wps.in: id = output_zip_name, type = string, title = Path to the folder where the new data will be stored. The folder will be created if not already existing in the system. , value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_108";
# wps.out: id = output_zip, type = text/zip, title = ZIP file containing the following datasets : PAN_l1.tif : PAN image as one single tif file (output of operation mosaic_panchromatic_images) / MS_l1.tif : MS image resized to fit PAN extent (output of operation resize_multispectral_image) / PANSHARP_l1.tif : Image pansharpened (output of operation pansharpen) / PAN_l2.tif : PAN image georeferenced (output of operation georeference) / MS_l2.tif : MS image georeferenced (output of operation georeference) / PANSHARP_l2.tif : Pansharpened image georeferenced (output of operation georeference) / CLD_MSK.gpkg : (output of operation preprocess_cloud_mask) / GCP.csv : Ground control points used for georeferencing .

mosaic_panchromatic_images=TRUE
resize_multispectral_image=TRUE
preprocess_cloud_mask=TRUE
pansharpen=TRUE
georeference=TRUE
pan_tifs_folder="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A"
pan_cloudmask_path="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A/MASKS/CLD_SPOT7_P_201710111019265_SEN_SPOT7_20171012_1350221m2hw3hzfkke0_1_MSK.GML"
ms_tif_file="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_MS_001_A/IMG_SPOT7_MS_201710111019265_SEN_SPOT7_20171012_1350131l9pjds0m4wzt_1_R1C1.TIF"
ground_control_points_path="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/ground_control_points.csv"
georeference_resampling_method="bilinear"
output_zip_name="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_108"

# Call useful libraries (+ install them if not already installed)
if(!require(raster)){
  install.packages("raster")
}
if(!require(XML)){
  install.packages("XML")
}
if(!require(gdalUtils)){
  install.packages("gdalUtils")
}
if(!require(sf)){
  install.packages("sf")
}
if(!require(dplyr)){
  install.packages("dplyr")
}

require(raster)
require(XML)
require(gdalUtils)
require(sf)
require(dplyr)


## global variables. To be set as local variables in future developments 
pansharpen_cell_size=5000

output_data_epsg=32630 # Corresponds to the UTM EPSG in our AOI
output_data_resolution=1.63313  # Can be found in the iso-metadata.xml file provided by Geosud

# Source useful functions
source("/home/ptaconet/Documents/r_react/mosaic_tif_images.R")
source("/home/ptaconet/Documents/r_react/pansharp.R")

# Create folder to store data that will be generated
dir.create(output_zip_name)

path_to_output_pan_tif=file.path(output_zip_name,"PAN_l1.tif") # Path to the output PAN tif
path_to_output_ms_tif=file.path(output_zip_name,"MS_l1.tif") # Path to the output MS tif
path_to_output_pansharpen_tif=file.path(output_zip_name,"PANSHARPEN_l1.tif") # Path to the output pansharpened tif

path_to_output_pan_tif_georeferenced=file.path(output_zip_name,"PAN_l2.tif") # Path to the output PAN tif
path_to_output_ms_tif_georeferenced=file.path(output_zip_name,"MS_l2.tif") # Path to the output MS tif
path_to_output_pansharpen_tif_georeferenced=file.path(output_zip_name,"PANSHARPEN_l2.tif") # Path to the output pansharpened tif

dir.create(file.path(output_zip_name,"PANSHARPENED_PIECES"))
folder_pansharpened_pieces=file.path(file.path(output_zip_name,"PANSHARPENED_PIECES"))


####### Start workflow
cat(paste0("Starting workflow"))

if (mosaic_panchromatic_images==TRUE){
## 1) Panchromatic images are split into many pieces. Step 1: Mosaic Panchromatic images

cat("Mosaic PAN images...\n")
pan_tifs_paths=as.vector(list.files(path=pan_tifs_folder,pattern=".TIF",full.names=TRUE)) # panchromatic tifs names
mosaic_tif_images(pan_tifs_paths,path_to_output_pan_tif)
cat("Mosaic PAN images OK\n")
} 

if (resize_multispectral_image==TRUE){
## 2) Multispectral image is not georeferenced. Step 2: Resize Multispectral images so that it reaches the size of PAN images  
cat("Resizing MS images...\n")

# Get MS image
ms=brick(ms_tif_file)
# set extent of ms image as the extent of the pan image
extent_pan=extent(raster(path_to_output_pan_tif))
extent(ms)=extent_pan
writeRaster(ms, file=path_to_output_ms_tif, format="GTiff",overwrite=TRUE)

cat("Resizing MS images OK\n")
}

if (preprocess_cloud_mask==TRUE){
## 3) Masks delivered (including cloud mask) are not properly georeferenced. Step 3: georeference cloud masks

cat("Pre-processing cloud masks ...\n")

mask_pan=st_read(pan_cloudmask_path)
coordinates_mask=as.data.frame(st_coordinates(mask_pan))
mask_corrected=st_sfc()
for (i in 1:length(unique(coordinates_mask[,"L2"]))){
  this_poly<-coordinates_mask %>% filter (L2==i)
  seq_coord_this_poly=NULL
  for (j in 1:nrow(this_poly)){
  seq_coord_this_poly=c(seq_coord_this_poly,this_poly[j,"X"],-this_poly[j,"Y"])
  }
  mask_corrected=c(mask_corrected,st_sfc(st_polygon(list(matrix(seq_coord_this_poly,ncol=2, byrow=TRUE)))))
}

st_write(mask_corrected,file.path(output_zip_name,"CLD_MSK.gpkg"))

cat("Pre-processing cloud masks OK\n")
}

if (pansharpen==TRUE){
## 4) Pan sharpening. We use a function developed by some guy, located here: https://www.r-bloggers.com/pan-sharpening-using-r/ .  This function uses the same algorithm as the OTB Toolbox

cat("Pan sharpening ...\n")

# Read PAN
pan=raster(path_to_output_pan_tif)
# Read MS band 1
ms=brick(path_to_output_ms_tif)

# Pansharpen by pieces of 5000 px * 5000 px, to avoid memory problems. With our images, 1 MS parsharpened image is approx. 10 GB
seq_coord_bbox=c(bbox(pan)["s1","min"],bbox(pan)["s2","max"],  
                 bbox(pan)["s1","max"],bbox(pan)["s2","max"],  
                 bbox(pan)["s1","max"],bbox(pan)["s2","min"],  
                 bbox(pan)["s1","min"],bbox(pan)["s2","min"],  
                 bbox(pan)["s1","min"],bbox(pan)["s2","max"])

bbox_sf=sf::st_sfc(st_polygon(list(matrix(seq_coord_bbox,ncol=2, byrow=TRUE))))

grid_crop_pansharpen=sf::st_make_grid(bbox_sf,what="polygons",cellsize = pansharpen_cell_size)

for (i in 1:length(grid_crop_pansharpen)){
  cat(paste0("pansharpening piece n° ",i," over ",length(grid_crop_pansharpen)," pieces to pansharpen (",round((i-1)*100/length(grid_crop_pansharpen))," %)...\n"))
grid_crop_pansharpen_this_box=grid_crop_pansharpen[i]
extent_crop=extent(trunc(as.numeric(st_bbox(grid_crop_pansharpen_this_box)$xmin))-10,
                   trunc(as.numeric(st_bbox(grid_crop_pansharpen_this_box)$xmax)),
                   trunc(as.numeric(st_bbox(grid_crop_pansharpen_this_box)$ymin))-10,
                   trunc(as.numeric(st_bbox(grid_crop_pansharpen_this_box)$ymax))
                   )
pan_crop=crop(pan,extent_crop)
ms_crop=crop(ms,extent_crop)

# Pansharpen
pansharp <- processingPansharp(pan_crop, ms_crop, filter = 'auto', fun = mean)

writeRaster(pansharp,file.path(folder_pansharpened_pieces,paste0(i,".tif")),format = "GTiff",datatype = 'INT2S',overwrite = T)
}

# Assemble pansharpened images:
mosaic_tif_images(list.files(folder_pansharpened_pieces,full.names = TRUE),path_to_output_pansharpen_tif)
# Remove the pieces of pansharpened images
file.remove()
cat("Pansharpening OK\n")
}

if (georeference==TRUE){
## 5) Images as delivered by GEOSUD are not georefereced. Step 5: georeference images using Ground Control Points (GCP)

# 4 GPC, corresponding to the 4 corners of the image, are provided in the INDEX.HTM file
# The pixel resolution is provided in the metadata-iso.xml file
GCP=read.csv(ground_control_points_path)

# convert GCP to output epsg
gcp_project=sf_project("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs", cbind(GCP$longitude,GCP$latitude))
GCP$longitude=gcp_project[,1]
GCP$latitude=gcp_project[,2]

path_images_to_georeference=c(path_to_output_pan_tif)
path_georeferenced_images=c(path_to_output_pan_tif_georeferenced)

gcp_system_command<-NULL
for (i in 1:nrow(GCP)){
  gcp_system_command=paste(gcp_system_command,"-gcp",GCP$col[i],GCP$row[i],GCP$longitude[i],GCP$latitude[i],sep=" ")
}

for (i in 1:length(path_images_to_georeference)){
system(paste0("gdal_translate -of GTiff ",gcp_system_command," -a_srs epsg:",output_data_epsg," ",path_images_to_georeference[i]," ",file.path(output_zip_name,"tmp_l2.tif")))
system(paste0("gdalwarp -of GTiff -r ",georeference_resampling_method," -tr ",output_data_resolution," ",output_data_resolution," ",file.path(output_zip_name,"tmp_l2.tif")," ",path_georeferenced_images[i]))
file.remove(file.path(output_zip_name,"tmp_l2.tif"))
}

}


## Zip all the files generated


#metadata=readHTMLTable(paste0(folder,"/",image_folder,"/",gsub("MD_","",image_folder),"/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A/INDEX.HTM"),useInternalNodes = TRUE)
#coordinates<-metadata[[14]]

## Retrieve ground control points and then execute the bash script to georeference the image