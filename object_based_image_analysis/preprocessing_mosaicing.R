## Mosaic images to get one single image that covers the whole area

path_to_otbApplications_folder<-"/home/ptaconet/OTB-6.6.0-Linux64/bin"
paths_images_to_mosaic<-c("/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_97/SPOT6_97_20170914_CIV_PANSHARPEN_ortho.TIF","/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_113/SPOT6_113_20171018_CIV_PANSHARPEN_ortho.TIF")
path_output_mosaic_image<-"/home/ptaconet/react/resultat_mosaique_97_113.TIF"

# We perform a very simple mosaic, by simply copying the last image over earlier ones in areas of overlap. We could perform much more advanced mosaicing operations with the otbcli_Mosaic application (for additional details: https://github.com/remicres/otb-mosaic)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(paths_images_to_mosaic, collapse = " ")," -comp.feather none -harmo.method none -out ",path_output_mosaic_image," uint16")
system(otb_appli)

