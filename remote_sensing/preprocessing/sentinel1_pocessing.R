require(stars)
require(gdalUtils)
require(sf)
require(getRemoteData)

path_to_otbApplications_folder<-"/home/ptaconet/OTB-6.6.1-Linux64/bin"
path_to_s1_folder<-"/home/ptaconet/react/datasets/data_BF/Sentinel_1/data_test_SAR_Sentinel1/raw_data"
path_to_s1_zip<-"S1B_IW_GRDH_1SDV_20181110T182640_20181110T182716_013546_01911D_32A2.zip" #S1B_IW_GRDH_1SDV_20181110T182640_20181110T182716_013546_01911D_32A2.zip     S1B_IW_GRDH_1SDV_20190427T182639_20190427T182714_015996_01E10D_6E59.zip

roi<-sf::st_read("/home/ptaconet/react/datasets/data_BF/ROI.kml")
utm_zone<-getRemoteData::getUTMepsg(roi)
roi_utm <- roi %>% sf::st_transform(utm_zone)


#unzip(file.path(path_to_s1_folder,path_to_s1_zip),exdir=path_to_s1_folder)

path_to_vh<-list.files(file.path(path_to_s1_folder,gsub(".zip",".SAFE",path_to_s1_zip),"measurement"),pattern = "s1b-iw-grd-vh",recursive = T,full.names = T)
path_to_vv<-list.files(file.path(path_to_s1_folder,gsub(".zip",".SAFE",path_to_s1_zip),"measurement"),pattern = "s1b-iw-grd-vv",recursive = T,full.names = T)

# Radiometric correction / calibration (https://www.orfeo-toolbox.org/CookBook/Applications/app_SARCalibration.html and https://www.asf.alaska.edu/asf-tutorials/data-recipes/correct-sentinel-data/esa-toolbox/)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SARCalibration")," -in ",path_to_vh," -out ",gsub(".tiff","_radcal.tif",path_to_vh))
system(otb_appli)
otb_appli<-paste0(file.path(path_to_otbApplications_folder,"otbcli_SARCalibration")," -in ",path_to_vv," -out ",gsub(".tiff","_radcal.tif",path_to_vv))
system(otb_appli)

# create vv/vh dataset
vv<-raster(gsub(".tiff","_radcal.tif",path_to_vv))
vh<-raster(gsub(".tiff","_radcal.tif",path_to_vh))
ratio_vv_vh<-vv/vh


# Geocoding  (tutorial available at : https://media.asf.alaska.edu/uploads/pdf/current_data_recipe_pdfs/GeocodingSentinelData_GDAL_v7.1.pdf)
gdalUtils::gdalwarp(srcfile=gsub(".tiff","_radcal.tif",path_to_vh),dstfile=gsub(".tiff","_radcal_geocoded.tif",path_to_vh),tps=T,r="bilinear",tr=c(10,10),srcnodata=0,dstnodata=0,t_srs=paste0("EPSG:",utm_zone))
gdalUtils::gdalwarp(srcfile=gsub(".tiff","_radcal.tif",path_to_vv),dstfile=gsub(".tiff","_radcal_geocoded.tif",path_to_vv),tps=T,r="bilinear",tr=c(10,10),srcnodata=0,dstnodata=0,t_srs=paste0("EPSG:",utm_zone))

# Update paths
path_to_vv<-gsub(".tiff","_radcal_geocoded.tif",path_to_vv)
path_to_vh<-gsub(".tiff","_radcal_geocoded.tif",path_to_vh)

# Open as a star object and crop to ROI
vh <- stars::read_stars(path_to_vh)[roi_utm]
vv <- stars::read_stars(path_to_vv)[roi_utm]

# Compute raster ratio : vv/vh
ratio_vv_vh<-vv/vh

