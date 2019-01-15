## Function to segment an image using the Baatz and shape algorithm (and then vectorize the result). Note: the application used - otbcli_LSGRM - is not available in the official release of OTB. We use a personalized release available here: http://napoli.teledetection.fr/logiciels/otb_moringa_build_win_x64.zip).
segment<-function(input_path,output_path,threshold,cw,sw,otbapplications_folder_path){
  
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_LSGRM")," -in ",input_path," -out ",gsub(".gpkg",".TIF",output_path)," int32 -threshold ",threshold," -criterion.bs.cw ",cw," -criterion.bs.sw ",sw)
system(gdal_appli)

# Vectorize
path_to_output_segmentation_vector=file.path(path_to_outputFiles_segmentation_folder,paste0(outputImageName,"_segmented_final.gpkg"))
gdal_appli<-paste0("gdal_polygonize.py ",gsub(".gpkg",".TIF",output_path)," ",output_path," -b 1 None DN")
system(gdal_appli)

if (file.exists(output_path)){
  return("Done")
} else {
  return("Error")
}
}
