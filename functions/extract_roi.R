## Function to extract a region of interest defined by a user. Details: https://www.orfeo-toolbox.org/Applications/ExtractROI.html 
extract_roi<-function(input_path,output_path,roi_vector_path,otbapplications_folder_path,dem_folder_path){
  otb_appli<-paste0(file.path(otbapplications_folder_path,"otbcli_ExtractROI")," -in ",input_path," -out ",output_path," uint16 -mode fit -mode.fit.vect ",roi_vector_path," -elev.dem ",dem_folder_path)
  system(otb_appli)
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
}
