######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = segmentation, title = Workflow to segment a satellite image, abstract = This scripts enables to segment a satellite image using the Large Scale Generic Region Merging algorithm of the Orfeo Toolbox. It first segment the image using the Baatz and Schape algorithm and then vectorize the output. Note: the application used - otbcli_LSGRM - is not available in the official release of OTB. We use a personalized release available here: http://napoli.teledetection.fr/logiciels/otb_moringa_build_win_x64.zip). At the time being the personal release is available only in Windows OS, hence, this script works only in Windows
# wps.in: id = path_to_image_to_segment, type = string, title = Path to the image to segment , value = "/home/ptaconet/Documents/react/data_BF/VHR_SPOT6/processed_data/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_108/PANSHARPEN_ROI.TIF";
# wps.in: id = path_to_output_segmented_data, type = string, title = Path to the output vector segmented dataset. , value = "/home/ptaconet/Documents/react/data_BF/Segmentation/segmentation_result.gpkg";
# wps.in: id = path_to_otbApplications_folder, type = string, title = Path to the folder containing the OTB applications. , value = "/home/ptaconet/OTB-6.6.1-Linux64/bin";
# wps.in: id = threshold, type = numeric, title = Threshold for the Generic Region Merging algorithm. , value = 300 ;
# wps.in: id = cw, type = numeric, title = poids radiométrie , value = 0.3 ;
# wps.in: id = sw, type = numeric, title = poids compacité, value = 0.75 ;
# wps.out: id = output_zip, type = text/zip, title = Result of the segmentation process : a TIF raster file and the related GPKG vector file 


path_to_image_to_segment<-"/home/ptaconet/Documents/react/data_BF/VHR_SPOT6/processed_data/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_108/PANSHARPEN_ROI.TIF"
path_to_output_segmented_data="/home/ptaconet/Documents/react/data_BF/Segmentation/segmentation_result.gpkg"
path_to_otbApplications_folder<-"/home/ptaconet/OTB-6.6.1-Linux64/bin"
threshold=300       #<paramètre d'échelle> 
cw=0.3              # <poids radiométrie> 
sw=0.75              # <poids compacité>

## Source segmentation function
source("functions_image_segmentation/segment.R")

### Start workflow
cat("Starting segmentation...")
res<-segment(path_to_image_to_segment,path_to_output_segmented_data,threshold,cw,sw,otbapplications_folder_path)
cat(res)


######## Test the segmentation on a ROI ########

## Extract a ROI
#mode.extent.unit="pxl" #[pxl/phy/lonlat]
#ulx=1000
#uly=1000
#lrx=3000
#lry=3000

#path_to_output_pansharpen_ROI_tif=file.path(path_to_outputFiles_segmentation_folder,paste0("ROI_",ulx,"_",uly,"_",lrx,"_",lry,".TIF"))
#otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_ExtractROI")," -in ",path_to_output_pansharpen_tif," -mode extent -mode.extent.ulx ",ulx," -mode.extent.uly ",uly," -mode.extent.lrx ",lrx," -mode.extent.lry ",lry," -mode.extent.unit ",mode.extent.unit," -elev.dem ",path_to_dem_folder," -out ",path_to_output_pansharpen_ROI_tif," uint16")
#system(otb_appli)

## Segment using Baatz and Schape algorithm. Other available algorithms are Euclidian Distance and Full Lambda Schedule
#path_to_output_segmentation=file.path(path_to_outputFiles_segmentation_folder,paste0("ROI_segmented_",ulx,"_",uly,"_",lrx,"_",lry,"_",threshold,"_",gsub("..","",cw),"_",gsub("..","",sw),".tif"))
#otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_GenericRegionMerging")," -in ",path_to_output_pansharpen_ROI_tif," -out ",path_to_output_segmentation," int32 -threshold ",threshold," -cw ",cw," -sw ",sw)
#system(otb_appli)

## Convert output raster to vector to visualize the results on QGIS and visualy compare them with the segmented image
#path_to_output_segmentation_vector=file.path(path_to_outputFiles_segmentation_folder,paste0("ROI_segmented_",ulx,"_",uly,"_",lrx,"_",lry,"_",threshold,"_",gsub("..","",cw),"_",gsub("..","",sw),".gpkg"))
#gdal_appli<-paste0(file.path(path_to_gdal_folder,"gdal_polygonize.py "),path_to_output_segmentation," ",path_to_output_segmentation_vector," -b 1 None DN")
#system(gdal_appli)


