######################################################################
##### 52North WPS annotations ##########
######################################################################
# wps.des: id = mosaic_tif_images, title = Mosaic (i.e. assemble into 1 single file) many tif images that are split into many images, abstract = Mosaic tif images that are split into many pieces. ;
# wps.in: id = path_to_input_tifs, type = string, title = Path to the tif images to mosaic (as vector - each element of the output vector is the path to one input image to mosaic). All tifs must be with the same projection., value = c("/home/documents/tif_1.tif","/home/documents/tif_2.tif","/home/documents/tif_3.tif","/home/documents/tif_4.tif");
# wps.in: id = path_to_output_tif, type = string, title = Path to the output mosaic tif generated., value = "/home/documents/tif_mosaic.tif";


mosaic_tif_images<-function(path_to_input_tifs,path_to_output_tif){
cat("in case there is no projection associated with the input tifs, R sends an error ERROR: failed to load SRS definition from ... . However, the mosaic tif is created (without any projection)")
require(raster)
#require(XML)
require(gdalUtils)

# Open all the rasters
rasts<-lapply(path_to_input_tifs,raster)

# Get extents of the rasters
extent_list<-lapply(rasts,extent)

# make a matrix out of it, each column represents a raster, rows the values
extent_list<-lapply(extent_list,as.matrix)
matrix_extent<-matrix(unlist(extent_list),ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin","ymin","xmax","ymax")

# create an extent with the extrem values of your extent
best_extent<-extent(min(matrix_extent[1,]),max(matrix_extent[3,]),min(matrix_extent[2,]),max(matrix_extent[4,]))

template <- raster(best_extent)
projection(template) <- projection(rasts[[1]])
writeRaster(template, file=path_to_output_tif, format="GTiff",overwrite=TRUE)
mosaic_rasters(gdalfile=path_to_input_tifs,dst_dataset=path_to_output_tif,of="GTiff")
return(gdalinfo(path_to_output_tif))
}