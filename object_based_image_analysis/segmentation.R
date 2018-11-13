

path_to_outputFiles_folder="C:\\Users\\taconet\\Documents\\data\\SPOT6_7_processed\\MD_SPOT6_2017_HC_BRUT_GEOSUD_108"
path_to_output_pansharpen_tif="C:\\Users\\taconet\\Documents\\data\\SPOT6_7_processed\\MD_SPOT6_2017_HC_BRUT_GEOSUD_108\\SPOT7_108_20171012_BF_PANSHARPEN_ortho.TIF"

path_to_dem_folder="C:\\Users\\taconet\\Documents\\data\\DEM_SRTM"
path_to_otbApplications_folder<-"C:\\otb\\bin"
gdal_polygonize_path<-"C:\\OSGeo4W64\\bin"



# Create directory to store output of segmentation process
path_to_outputFiles_segmentation_folder=file.path(path_to_outputFiles_folder,"segmentation")
dir.create(path_to_outputFiles_segmentation_folder)





## Extract ROI to test segmentation
mode.extent.unit="pxl" #[pxl/phy/lonlat]
ulx=1000
uly=1000
lrx=3000
lry=3000

path_to_output_pansharpen_ROI_tif=file.path(path_to_outputFiles_segmentation_folder,paste0("ROI_",ulx,"_",uly,"_",lrx,"_",lry,".TIF"))
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_ExtractROI")," -in ",path_to_output_pansharpen_tif," -mode extent -mode.extent.ulx ",ulx," -mode.extent.uly ",uly," -mode.extent.lrx ",lrx," -mode.extent.lry ",lry," -mode.extent.unit ",mode.extent.unit," -elev.dem ",path_to_dem_folder," -out ",path_to_output_pansharpen_ROI_tif," uint16")
system(otb_appli)



### SEGMENTATION using Baatz and Schape algorithm. Other available algorithms are Euclidian Distance and Full Lambda Schedule
# First test on the ROI. Then when parameters seem to be ok, apply on all the image (very long process...)


# Loop to test several combinations 
threshold=300       #<param??tre d?????chelle> 
cw=0.3              # <poids radiom??trie> 
sw=0.75              # <poids compacit??>

path_to_output_segmentation=file.path(path_to_outputFiles_segmentation_folder,paste0("ROI_segmented_",ulx,"_",uly,"_",lrx,"_",lry,"_",threshold,"_",gsub("..","",cw),"_",gsub("..","",sw),".tif"))
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_GenericRegionMerging")," -in ",path_to_output_pansharpen_ROI_tif," -out ",path_to_output_segmentation," int32 -threshold ",threshold," -cw ",cw," -sw ",sw)
system(otb_appli)


# Then convert to vector to visualize the results on QGIS and visualy compare them with the segmented image
path_to_output_segmentation_vector=file.path(path_to_outputFiles_segmentation_folder,paste0("ROI_segmented_",ulx,"_",uly,"_",lrx,"_",lry,"_",threshold,"_",gsub("..","",cw),"_",gsub("..","",sw),".gpkg"))
gdal_appli<-paste0(file.path(gdal_polygonize_path,"gdal_polygonize.py "),path_to_output_segmentation," ",path_to_output_segmentation_vector," -b 1 None DN")
system(gdal_appli)





# Segment the whole image when the segmentation paramters are valid and then convert to vector (note that the application used - otbcli_LSGRM - is not available in the official release of OTB. We use a personalized release available here: http://napoli.teledetection.fr/logiciels/otb_moringa_build_win_x64.zip)
# Segment
path_to_output_segmentation=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_segmented_final.tif"))
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_LSGRM")," -in ",path_to_output_pansharpen_tif," -out ",path_to_output_segmentation," int32 -threshold ",threshold," -criterion.bs.cw ",cw," -criterion.bs.sw ",sw)
system(gdal_appli)
# Vectorize
path_to_output_segmentation_vector=file.path(path_to_outputFiles_folder,paste0(outputImageName,"_segmented_final.gpkg"))
gdal_appli<-paste0("gdal_polygonize.py ",path_to_output_segmentation," ",path_to_output_segmentation_vector," -b 1 None DN")
system(gdal_appli)
