## Function to convert digital numbers to Top of Atmosphere reflectance. Details: https://www.orfeo-toolbox.org/Applications/OpticalCalibration.html
convert_dn_to_toa_reflectance<-function(input_path,output_path,otbapplications_folder_path){
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_OpticalCalibration")," -in ",input_path," -out ",gsub(".TIF","_temp.TIF",output_path)," float -level toa")
  system(otb_appli)
  otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_BandMathX")," -il ",gsub(".TIF","_temp.TIF",output_path)," -exp \"10000*im1\" -out ",output_path," uint16")
  system(otb_appli)
  file.remove(gsub(".TIF","_temp.TIF",output_path))
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}