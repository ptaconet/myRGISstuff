## Function to mosaic multiple remote sensing images. Details: https://github.com/remicres/otb-mosaic
mosaic<-function(input_path,output_path,otbapplications_folder_path){
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Mosaic")," -il ",paste(input_path, collapse = ' ')," -out ",output_path," uint16")
  system(otb_appli)
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}
