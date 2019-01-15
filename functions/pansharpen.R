## Function to Pansharpen. Details: https://www.orfeo-toolbox.org/Applications/Superimpose.html and https://www.orfeo-toolbox.org/Applications/BundleToPerfectSensor.html
pansharpen<-function(pan_input_path,ms_input_path,output_path,otbapplications_folder_path,dem_path){
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_Superimpose")," -inr ",pan_input_path," -inm ",ms_input_path," -out ",gsub(".TIF","_temp.TIF",output_path)," uint16")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BundleToPerfectSensor")," -inp ",pan_input_path," -inxs ",gsub(".TIF","_temp.TIF",output_path)," -out ",output_path," uint16")
  system(otb_appli)
  file.remove(gsub(".TIF","_temp.TIF",output_path))
  file.remove(gsub(".TIF",".geom",gsub(".TIF","_temp.TIF",output_path)))
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}