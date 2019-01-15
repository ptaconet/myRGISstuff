download_copernicus_scihub_products<-function(username,password,product_id,output_path){
  
  system(paste0("wget --content-disposition --continue --user=",username," --password=",password," --directory-prefix ",output_path," \"https://scihub.copernicus.eu/dhus/odata/v1/Products('",product_id,"')/\\$value\\"))
  
  if (file.exists(output_path)){
    return("Done")
  } else {
    return("Error")
  }
  
}


