# article: https://www.mdpi.com/2072-4292/8/4/354

library(raster)

path_to_sentinel2_folder="/home/ptaconet/Documents/react/data_BF/HR_Sentinel2/processed_data/S2B_MSIL1C_20181116T103309_N0207_R108_T30PVS_20181116T160025.SAFE/GRANULE/L1C_T30PVS_A008857_20181116T104132/IMG_DATA"

b11<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B11.jp2"))
b2<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B02.jp2"))
b3<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B03.jp2"))
b4<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B04.jp2"))
b8<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B08.jp2"))

# upscale resolution of b2,b3,b4,b8
system(paste0('gdal_translate -tr 20 20 "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B02.jp2"),'" "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B02_20m.jp2"),'"'))
b2_20m<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B02_20m.jp2"))

system(paste0('gdal_translate -tr 20 20 "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B03.jp2"),'" "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B03_20m.jp2"),'"'))
b3_20m<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B03_20m.jp2"))

system(paste0('gdal_translate -tr 20 20 "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B04.jp2"),'" "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B04_20m.jp2"),'"'))
b4_20m<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B04_20m.jp2"))

system(paste0('gdal_translate -tr 20 20 "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B08.jp2"),'" "',file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B08_20m.jp2"),'"'))
b8_20m<-raster(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B08_20m.jp2"))

b<-values(b11)


plot(b2_20m,b11)
a<-values(b2_20m)
abline(fit<-lm((b ~ a),col='red')
legend("topright",bty="n",legend=paste("R2 is",format(summary(fit)$adj.r.squared,digits=4)))

plot(b3_20m,b11)
a<-values(b3_20m)
abline(fit<-lm((b ~ a),col='red'))
legend("topright",bty="n",legend=paste("R2 is",format(summary(fit)$adj.r.squared,digits=4)))

plot(b4_20m,b11)
a<-values(b4_20m)
abline(fit<-lm((b ~ a),col='red'))
legend("topright",bty="n",legend=paste("R2 is",format(summary(fit)$adj.r.squared,digits=4)))

plot(b8_20m,b11)
a<-values(b8_20m)
abline(fit<-lm((b ~ a),col='red'))
legend("topright",bty="n",legend=paste("R2 is",format(summary(fit)$adj.r.squared,digits=4)))



# Panshaperning of B11 using B4

pansharpen(file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B04.jp2"),
           file.path(path_to_sentinel2_folder,"T30PVS_20181116T103309_B11.jp2"),
           file.path(path_to_sentinel2_folder,"B11_pansharpened.jp2"),
           "/home/ptaconet/OTB-6.6.1-Linux64/bin",
           "/home/ptaconet/Documents/react/data_BF/DEM_SRTM/raw_data")

# Calculate MDNWI

