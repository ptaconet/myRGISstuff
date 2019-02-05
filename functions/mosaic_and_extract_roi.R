mosaic_and_extract_roi<-function(input_path_files,output_path_file,roi_vector_path,otbapplications_folder_path,dem_folder_path){
  
  cat("Mosaicing the tiles...")
  output_path_temp_file=gsub("\\.","_temporary.",output_path_file)
  res<-mosaic(input_path_files,output_path_temp_file,otbapplications_folder_path)
  if (res=="Done"){
    cat ("Mosaic OK")
  }
  
  cat("Extracting the ROI...")
  res<-extract_roi(output_path_temp_file,output_path_file,roi_vector_path,otbapplications_folder_path,dem_folder_path)
  if (res=="Done"){
    cat ("ROI Extraction OK")
  }
  
  file.remove(output_path_temp_file)
  
  if (file.exists(output_path_file)){
    return("Done")
  } else {
    return("Error")
  }
}