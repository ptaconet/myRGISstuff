######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = preprocessing, title = Pre-processing of SPOT 6/7 images (tile fusionning + orthorectification + pansharpening + cloud mask extraction), abstract = This script is a workflow for the pre-processing SPOT6/7 satellite imagery data (tile fusionning + orthorectification + pansharpening + cloud mask extraction). The user can parameterize the operations he/she wants to be computed. The script uses applications coming from various libraries that must be installed : the Orfeo Toolbox (https://www.orfeo-toolbox.org/), R spatial packages ("sf" for vector and "raster" for raster) and the Geospatial Data Abstraction Library (https://www.gdal.org/).
# wps.in: id = pan_tile_fusionning, type = boolean, title = Panchromatic images (PAN) are split into many tiles. Fusion tiles into 1 single TIF file? , value="TRUE|FALSE" ;
# wps.in: id = toa_reflectance_conversion, type = boolean, title = Convert MS and PAN to Top-of-atmosphere reflectance? , value="TRUE|FALSE" ;
# wps.in: id = pan_orthorectification, type = boolean, title = Orthorectify the panchromatic tiled image? , value="TRUE|FALSE" ;
# wps.in: id = ms_orthorectification, type = boolean, title = Orthorectify the multispectral image? , value="TRUE|FALSE" ;
# wps.in: id = ms_pansharpening, type = boolean, title = Pansharpen MS image ? , value="TRUE|FALSE"; 
# wps.in: id = extract_roi, type = boolean, title = Extract following a region of interest ? , value="TRUE|FALSE"; 
# wps.in: id = cloud_mask_generation, type = boolean, title = Some masks are included in the products delivered however they are not geo-refererenced. Generate a geo-referenced cloud mask ? , value="TRUE|FALSE"; 
# wps.in: id = path_to_panImages_folder, type = string, title = Path to the folder where the PAN .TIF images are stored, value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A";
# wps.in: id = path_to_msImage_folder, type = string, title = Path to the MS .TIF file, value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_MS_001_A/IMG_SPOT7_MS_201710111019265_SEN_SPOT7_20171012_1350131l9pjds0m4wzt_1_R1C1.TIF";
# wps.in: id = path_to_msCloudMask_file, type = string, title = Path to the cloud mask file, value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A/MASKS/CLD_SPOT7_P_201710111019265_SEN_SPOT7_20171012_1350221m2hw3hzfkke0_1_MSK.GML";
# wps.in: id = path_to_outputFiles_folder, type = string, title = Path to the folder where the new data will be stored. The folder will be created if not already existing in the system. , value = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_108";
# wps.in: id = outputImageName, type = string, title = Prefix use for the output files (i.e. files that will be generated in the workflow) , value = "SEN_SPOT7_20171012";
# wps.in: id = path_to_roi_polygon, type = string, title = Path to the ROI in kml or shp format, value = "/home/ptaconet/Documents/react/REACT_BF.kml";
# wps.in: id = path_to_path_to_dem_folder, type = string, title = Path to the folder containing the SRTM Digital Elevation Model files (the DEM is used for orthorectification). The DEM should be available in .hgt format. SRTM DEM can be easily downloaded here: http://dwtkns.com/srtm30m/, value = "/home/ptaconet/react/SRTM_dem";

# wps.in: id = path_to_otbApplications_folder, type = string, title = Path to the folder containing the OTB applications. , value = "/home/ptaconet/OTB-6.6.0-Linux64/bin";
# wps.out: id = output_zip, type = text/zip, title = ZIP file containing the following datasets : PAN_mosaic.tif : PAN mosaiced tif file / PAN_ortho.tif : PAN orthorectified tif file / MS_ortho.tif : MS orthorectified tif file / PANSHARPEN_ortho.TIF : MS pansherpened orthorectified tif file / CLDMASK_ortho.TIF : cloud mask georeferenced tif file / CLDMASK_ortho.shp : cloud mask georeferenced shapefile

#### Additional informations on how Spot 6/7 images are delivered can be found on the SPOT 6/7 user guide available here: https://www.spaceoffice.nl/blobs/Dataportaal/User_Guide_SPOT6_V1.0.pdf ;
#### Some operations might work only on Linux OS

### Variables to be set by the user
# Workflow parameterization
pan_tile_fusionning<-FALSE
toa_reflectance_conversion<-TRUE
pan_orthorectification<-TRUE
ms_orthorectification<-TRUE
ms_pansharpening<-TRUE
extract_roi<-TRUE
cloud_mask_generation<-TRUE

# Paths to input data folders / files
path_to_panImages_folder="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_113/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_113/PROD_SPOT6_001/VOL_SPOT6_001_A/IMG_SPOT6_P_001_A"
path_to_msImage_folder="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_113/SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_113/PROD_SPOT6_001/VOL_SPOT6_001_A/IMG_SPOT6_MS_001_A"
path_to_msCloudMask_file="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_113/SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_113/PROD_SPOT6_001/VOL_SPOT6_001_A/IMG_SPOT6_MS_001_A/MASKS/CLD_SPOT6_MS_201710171025039_SEN_SPOT6_20171018_1304181fuxjrlyazjr4_1_MSK.GML"
path_to_outputFiles_folder="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_113"
outputImageName="SPOT6_113_20171018_CIV"
path_to_roi_polygon="/home/ptaconet/Documents/react/REACT_CIV.kml"
path_to_dem_folder="/home/ptaconet/react/SRTM_dem"
path_to_otbApplications_folder<-"/home/ptaconet/OTB-6.6.1-Linux64/bin"


## Start of the workflow

## Call useful packages
library(sf)
library(raster)
library(dplyr)

# Create folder to store data that will be generated
dir.create(path_to_outputFiles_folder)

# Get or set useful paths
path_to_output_pan_mosaic_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_PAN_mosaic.TIF")) # Path to the output mosaic PAN tif
path_to_output_pan_toa_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_PAN_TOAx10000.TIF")) # Path to the output PAN TOAtif
path_to_output_ms_toa_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_MS_TOAx10000.TIF")) # Path to the output MS TOA tif
path_to_output_pan_ortho_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_PAN_ortho.TIF")) # Path to the output mosaic PAN tif orthorectified 
path_to_output_ms_ortho_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_MS_ortho.TIF")) # Path to the output MS tif orthorectified 
path_to_output_pansharpen_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_PANSHARPEN_ortho.TIF")) # Path to the output pansharpened tif orthorectified
path_to_output_pansharpen_tif_roi=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_PANSHARPEN_ortho_roi.TIF")) # Path to the region of interest - output pansharpened tif orthorectified 
path_to_output_cloudmask_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_CLDMASK_ortho.TIF")) # Path to the output cloud mask tif orthorectified
path_to_output_cloudmask_shp=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_CLDMASK_ortho.shp")) # Path to the output cloud mask shapefile orthorectified

pan_tifs_paths=as.vector(list.files(path=path_to_panImages_folder,pattern=".TIF",full.names=TRUE)) 
ms_tif_path=as.vector(list.files(path=path_to_msImage_folder,pattern=".TIF",full.names=TRUE))
ms_tif_outputImageName=as.vector(list.files(path=path_to_msImage_folder,pattern=".TIF"))


### Start workflow

######## 1) PAN tile fusionning ###########
if (pan_tile_fusionning==TRUE){
cat("Starting PAN tile fusionning ...")
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(pan_tifs_paths, collapse = ' ')," -out ",path_to_output_pan_mosaic_tif," uint16")
system(otb_appli)
cat("PAN tile fusionning OK")
}

######## 2) Convertion to TOA reflectance ########
if (toa_reflectance_conversion==TRUE){
cat("Starting Convertion to TOA reflectance ...")
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OpticalCalibration")," -in ",path_to_output_pan_mosaic_tif," -out ",gsub(".TIF","_temp.TIF",path_to_output_pan_toa_tif)," float")
system(otb_appli)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BandMathX")," -il ",gsub(".TIF","_temp.TIF",path_to_output_pan_toa_tif)," -exp \"10000*im1\" -out ",path_to_output_pan_toa_tif," uint16")
system(otb_appli)
file.remove(gsub(".TIF","_temp.TIF",path_to_output_pan_toa_tif))

otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OpticalCalibration")," -in ",ms_tif_path," -out ",gsub(".TIF","_temp.TIF",path_to_output_ms_toa_tif)," float")
system(otb_appli)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BandMathX")," -il ",gsub(".TIF","_temp.TIF",path_to_output_ms_toa_tif)," -exp \"10000*im1\" -out ",path_to_output_ms_toa_tif," uint16")
system(otb_appli)
file.remove(gsub(".TIF","_temp.TIF",path_to_output_ms_toa_tif))

ms_tif_path = path_to_output_ms_toa_tif
path_to_output_pan_mosaic_tif = path_to_output_pan_toa_tif

cat("Convertion to TOA reflectance OK")
}


######## 2) PAN orthorectification ###########
if (pan_orthorectification==TRUE){
cat("Starting PAN orthorectification ...")
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OrthoRectification")," -io.in \"",path_to_output_pan_mosaic_tif,"?&skipcarto=true\""," -io.out ",path_to_output_pan_ortho_tif," uint16 -elev.dem ",path_to_dem_folder)
system(otb_appli)
cat("PAN orthorectification OK")
}

########### 3) MS orthorectification ###########
if (ms_orthorectification==TRUE){
cat("Starting MS orthorectification ...")
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OrthoRectification")," -io.in \"",ms_tif_path,"?&skipcarto=true\""," -io.out ",path_to_output_ms_ortho_tif," uint16 -elev.dem ",path_to_dem_folder)
system(otb_appli)
cat("MS orthorectification OK")
}

########### 4) ROI Extraction ###########
if (extract_roi==TRUE){ 
  cat("Extracting the ROI for the MS image...")
  
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_ExtractROI")," -in ",path_to_output_ms_ortho_tif," -out ",path_to_output_ms_ortho_tif_roi," uint16 -mode fit -mode.fit.vect ",path_to_roi_polygon," -elev.dem ",path_to_dem_folder)
  system(otb_appli)
  
  # rename files
  file.remove(path_to_output_ms_ortho_tif)
  file.remove(gsub(".TIF",".geom",path_to_output_ms_ortho_tif))
  file.rename(path_to_output_ms_ortho_tif_roi,path_to_output_ms_ortho_tif)
  file.rename(gsub(".TIF",".geom",path_to_output_ms_ortho_tif_roi),gsub(".TIF",".geom",path_to_output_ms_ortho_tif))
  
  cat("Extracting the ROI for the PAN image...")
  
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_ExtractROI")," -in ",path_to_output_pan_ortho_tif," -out ",path_to_output_pan_ortho_tif_roi," uint16 -mode fit -mode.fit.vect ",path_to_roi_polygon," -elev.dem ",path_to_dem_folder)
  system(otb_appli)
  
  # rename files
  file.remove(path_to_output_pan_ortho_tif)
  file.remove(gsub(".TIF",".geom",path_to_output_pan_ortho_tif))
  file.rename(path_to_output_pan_ortho_tif_roi,path_to_output_pan_ortho_tif)
  file.rename(gsub(".TIF",".geom",path_to_output_pan_ortho_tif_roi),gsub(".TIF",".geom",path_to_output_pan_ortho_tif))
  
  cat("Extracting the ROI OK")
  
}

########### 5) MS Pansharpening ###########
if (ms_pansharpening==TRUE){
cat("Starting MS Pansharpening ...")
path_to_output_ms_ortho_superimpose_tif=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_MS_L1_superimpose.TIF"))
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Superimpose")," -inr ",path_to_output_pan_ortho_tif," -inm ",path_to_output_ms_ortho_tif," -out ",path_to_output_ms_ortho_superimpose_tif," uint16")
system(otb_appli)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BundleToPerfectSensor")," -inp ",path_to_output_pan_ortho_tif," -inxs ",path_to_output_ms_ortho_superimpose_tif," -out ",path_to_output_pansharpen_tif," uint16")
system(otb_appli)
file.remove(path_to_output_ms_ortho_superimpose_tif)
file.remove(gsub(".TIF",".geom",path_to_output_ms_ortho_superimpose_tif))
cat("MS Pansharpening OK")
}

########### 6) Cloud mask creation (both in raster and vector formats) ###########
if (cloud_mask_generation==TRUE){ 
cat("Starting cloud mask creation ...")
dimap_path=list.files(path=path_to_msImage_folder,pattern="DIM",full.names=TRUE)
tfw_path=list.files(path=path_to_msImage_folder,pattern="TFW",full.names=TRUE)
rpc_path=list.files(path=path_to_msImage_folder,pattern="RPC",full.names=TRUE)
file.copy(from = dimap_path, to = path_to_outputFiles_folder)
file.copy(from = tfw_path, to = path_to_outputFiles_folder)
file.copy(from = rpc_path, to = path_to_outputFiles_folder)
dimap_path_in_new_folder=list.files(path=path_to_outputFiles_folder,pattern="DIM",full.names=TRUE)
tfw_path_in_new_folder=list.files(path=path_to_outputFiles_folder,pattern="TFW",full.names=TRUE)
rpc_path_in_new_folder=list.files(path=path_to_outputFiles_folder,pattern="RPC",full.names=TRUE)

# Align with raster grid (i.e. negative latitudes)
mask=sf::st_read(path_to_msCloudMask_file)
coordinates_mask=as.data.frame(st_coordinates(mask))
mask_corrected=sf::st_sfc()
for (i in 1:length(unique(coordinates_mask[,"L2"]))){
  this_poly<-coordinates_mask %>% filter (L2==i) %>% filter (L1==1)
  seq_coord_this_poly=NULL
  for (j in 1:nrow(this_poly)){
    seq_coord_this_poly=c(seq_coord_this_poly,this_poly[j,"X"],-this_poly[j,"Y"])
  }
  mask_corrected=c(mask_corrected,st_sfc(st_polygon(list(matrix(seq_coord_this_poly,ncol=2, byrow=TRUE)))))
}  

# Convert to raster, reading extent, ncols, nrows and resolution from tif
tif=raster(ms_tif_path)
r <- raster(as(mask_corrected, "Spatial"), ncols = ncol(tif) , nrows = nrow(tif), ext=extent(tif),res=c(1,1))
rr <- rasterize(as(mask_corrected, "Spatial"), r, progress = "text")
rr[!is.na(rr)]=1
writeRaster(rr,gsub(".TFW",".tif",tfw_path_in_new_folder), format="GTiff",overwrite=TRUE,datatype='INT4U')

# Orthorectify
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OrthoRectification")," -io.in \"",gsub(".TFW",".tif",tfw_path_in_new_folder),"?&skipcarto=true\""," -io.out ",path_to_output_cloudmask_tif," int16 -elev.dem ",path_to_dem_folder)
system(otb_appli)

# Pre-process before rasterization
cld_mask_rast<-raster(path_to_output_cloudmask_tif)
cld_mask_rast[cld_mask_rast != 1] <- NA
writeRaster(cld_mask_rast,path_to_output_cloudmask_tif, format="GTiff",overwrite=TRUE,datatype='INT4U')
file.rename(gsub(".TIF",".tif",path_to_output_cloudmask_tif),path_to_output_cloudmask_tif)

# Convert to vector
gdal_appli<-paste0("gdal_polygonize.py ",path_to_output_cloudmask_tif," ",path_to_output_cloudmask_shp," -b 1 -f \"ESRI Shapefile\" None DN")
system(gdal_appli)

# Remove temporary files
file.remove(c(tfw_path_in_new_folder,dimap_path_in_new_folder,rpc_path_in_new_folder))

cat("Cloud mask creation OK")

}
