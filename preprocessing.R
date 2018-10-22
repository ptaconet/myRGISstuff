dem_folder="/home/ptaconet/react/SRTM_dem"
otb_applications_path<-"/home/ptaconet/OTB-6.6.0-Linux64/bin"

pan_tifs_folder="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A"
pan_cloudmask_path="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A/MASKS/CLD_SPOT7_P_201710111019265_SEN_SPOT7_20171012_1350221m2hw3hzfkke0_1_MSK.GML"
ms_tif_file="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_MS_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_MS_001_A/IMG_SPOT7_MS_201710111019265_SEN_SPOT7_20171012_1350131l9pjds0m4wzt_1_R1C1.TIF"
output_zip_name="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_108"
image_name="SEN_SPOT7_20171012"

# Create folder to store data that will be generated
dir.create(output_zip_name)

path_to_output_pan_mosaic_tif=file.path(output_zip_name,paste0(image_name,"_PAN_mosaic.tif")) # Path to the output mosaic PAN tif
path_to_output_pan_ortho_tif=file.path(output_zip_name,paste0(image_name,"_PAN_l1.tif")) # Path to the output mosaic PAN tif
path_to_output_ms_ortho_tif=file.path(output_zip_name,paste0(image_name,"_MS_l1.tif")) # Path to the output MS tif
path_to_output_ms_ortho_superimpose_tif=file.path(output_zip_name,paste0(image_name,"_MS_L1_superimpose.tif"))
path_to_output_pansharpen_tif=file.path(output_zip_name,paste0(image_name,"_PANSHARPEN_l1.tif")) # Path to the output pansharpened tif
path_to_output_pansharpen_gpkg=file.path(output_zip_name,paste0(image_name,"_PANSHARPEN.gpkg"))

path_to_output_pansharpen_ROI_tif=file.path(output_zip_name,paste0(image_name,"_PANSHARPEN_ROI.tif"))

# 1.2) Mosaiquage des PAN
pan_tifs_paths=as.vector(list.files(path=pan_tifs_folder,pattern=".TIF",full.names=TRUE)) # panchromatic tifs names
output_file_name=list.files(path=pan_tifs_folder,pattern="DIM")
output_file_name=gsub("DIM","IMG",output_file_name)
output_file_name=gsub(".XML",".TIF",output_file_name)
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_Mosaic")," -il ",paste(pan_tifs_paths, collapse = ' ')," -out ",path_to_output_pan_mosaic_tif," uint16")
system(otb_appli)

# 1.2) Ortorectification de la PAN
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_OrthoRectification")," -io.in \"",path_to_output_pan_mosaic_tif,"?&skipcarto=true\""," -io.out ",path_to_output_pan_ortho_tif," uint16 -elev.dem ",dem_folder)
system(otb_appli)

# 2) Orthorectificaiton de la MS
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_OrthoRectification")," -io.in \"",ms_tif_file,"?&skipcarto=true\""," -io.out ",path_to_output_ms_ortho_tif," uint16 -elev.dem ",dem_folder)
system(otb_appli)

# 3) Pansharpening
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_Superimpose")," -inr ",path_to_output_pan_ortho_tif," -inm ",path_to_output_ms_ortho_tif," -out ",path_to_output_ms_ortho_superimpose_tif," uint16")
system(otb_appli)
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_BundleToPerfectSensor")," -inp ",path_to_output_pan_ortho_tif," -inxs ",path_to_output_ms_ortho_superimpose_tif," -out ",path_to_output_pansharpen_tif," uint16")
system(otb_appli)

## Extract ROI to test segmentation
mode.extent.unit="pxl" #[pxl/phy/lonlat]
ulx=1000
uly=1000
lrx=3000
lry=3000
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_ExtractROI")," -in ",path_to_output_pansharpen_tif," -mode extent -mode.extent.ulx ",ulx," -mode.extent.uly ",uly," -mode.extent.lrx ",lrx," -mode.extent.lry ",lry," -mode.extent.unit ",mode.extent.unit," -elev.dem ",dem_folder," -out ",path_to_output_pansharpen_ROI_tif," uint16")
system(otb_appli)


# 4) NOT WORKING Convert to geopackage raster and build pyramids (more compressed than tif = better for visualization in qfield)
#system(paste0("gdal_translate -ot Byte -of GPKG ",path_to_output_pansharpen_ROI_tif," ",path_to_output_pansharpen_gpkg))
#system(paste0("gdaladdo --config OGR_SQLITE_SYNCHRONOUS OFF -r AVERAGE ",path_to_output_pansharpen_gpkg," 2 4 8 16 32 64 128 256"))


### SEGMENTATION using Baatz and Schape algorithm. Other available algorithms are Euclidian Distance and Full Lambda Schedule
# First test on the ROI. Then when parameters seem to be ok, apply on all the image (very long process...)
threshold=25       #<paramètre d’échelle> 
cw=0.7              # <poids radiométrie> 
sw=0.3              # <poids compacité>
path_to_output_segmentation=file.path(output_zip_name,paste0(image_name,"_segmented_",threshold,"_",gsub("..","",cw),"_",gsub("..","",sw),".tif"))
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_GenericRegionMerging")," -in ",path_to_output_pansharpen_ROI_tif," -out ",path_to_output_segmentation," int32 -threshold ",threshold," -cw ",cw," -sw ",sw)
system(otb_appli)

# Then convert to vector to visualize the results on QGIS and visualy compare them with the segmented image
path_to_output_segmentation_vector=file.path(output_zip_name,paste0(image_name,"_segmented_",threshold,"_",gsub("..","",cw),"_",gsub("..","",sw),".gpkg"))
gdal_appli<-paste0("gdal_polygonize.py ",path_to_output_segmentation," ",path_to_output_segmentation_vector," -b 1 None DN")
system(gdal_appli)

# Segment the whole image when the segmentation paramters are valid and then convert to vector (note that the application used - otbcli_LSGRM - is not available in the official release of OTB. We use a personalized release available here: http://napoli.teledetection.fr/logiciels/otb_moringa_build_win_x64.zip)
# Segment
path_to_output_segmentation=file.path(output_zip_name,paste0(image_name,"_segmented_final.tif"))
otb_appli<-paste0(file.path(otb_applications_path,"otbcli_LSGRM")," -in ",path_to_output_pansharpen_tif," -out ",path_to_output_segmentation," int32 -threshold ",threshold," -criterion.bs.cw ",cw," -criterion.bs.sw ",sw)
system(gdal_appli)
# Vectorize
path_to_output_segmentation_vector=file.path(output_zip_name,paste0(image_name,"_segmented_final.gpkg"))
gdal_appli<-paste0("gdal_polygonize.py ",path_to_output_segmentation," ",path_to_output_segmentation_vector," -b 1 None DN")
system(gdal_appli)


