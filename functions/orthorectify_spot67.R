## Function to orthorectify images. Details: https://www.orfeo-toolbox.org/Applications/OrthoRectification.html
orthorectify_spot67<-function(input_path,output_path,otbapplications_folder_path,dem_folder_path){
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OrthoRectification")," -io.in \"",input_path,"?&skipcarto=true\""," -io.out ",output_path," uint16 -elev.dem ",dem_folder_path)
  system(otb_appli)
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}