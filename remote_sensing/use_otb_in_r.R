# https://gis.stackexchange.com/questions/261607/use-orfeo-toolbox-in-r
meanshift.segm <- function(otb.path = "", raster.in = "", out.path = "", name ="", filter.meanshift.spatialr = "5",
                           filter.meanshift.ranger = "15", filter.meanshift.thres = "0.1",
                           filter.meanshift.maxiter = "100", filter.meanshift.minsize = "100",
                           mode.vector.outmode = "ovw", mode.vector.inmask = "", mode.vector.neighbor = "false",
                           mode.vector.stitch = "true", mode.vector.minsize = 1, mode.vector.simplify = 0.1,
                           mode.vector.layername = "layer", mode.vector.fieldname = "DN", mode.vector.tilesize = 1024,
                           mode.vector.startlabel = 1){
  # Set configuration      
  conf <- paste("-in",raster.in,"-filter meanshift","-filter.meanshift.spatialr",filter.meanshift.spatialr,
                "-filter.meanshift.ranger",filter.meanshift.ranger,"-filter.meanshift.thres",filter.meanshift.thres,
                "-filter.meanshift.maxiter",filter.meanshift.maxiter,"-filter.meanshift.minsize",filter.meanshift.minsize,
                "-mode vector","-mode.vector.out",paste(out.path,"/",name,".shp",sep=""),"-mode.vector.outmode",mode.vector.outmode,
                ifelse(missingArg(mode.vector.inmask),"",paste("-mode.vector.inmask",mode.vector.inmask)),
                "-mode.vector.neighbor", mode.vector.neighbor,
                "-mode.vector.stitch",mode.vector.stitch,
                "-mode.vector.minsize",mode.vector.minsize,
                "-mode.vector.simplify",mode.vector.simplify,
                "-mode.vector.layername",mode.vector.layername,
                "-mode.vector.fieldname",mode.vector.fieldname,
                "-mode.vector.tilesize",mode.vector.tilesize,
                "-mode.vector.startlabel",mode.vector.startlabel)
  # apply function in command line
  system(paste(otb.path,"/otbcli_Segmentation"," ",conf,sep=""))
  # save configuration for futher use
  write.table(x = conf,file = paste(out.path,"/",name,"_conf.txt",sep=""),row.names = F, col.names = F)
}

# usage, you can set any option listed above
meanshift.segm(otb.path = "/home/ptaconet/OTB-6.6.0-Linux64/bin", raster.in = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A/IMG_SPOT7_P_201710111019265_SEN_SPOT7_20171012_1350221m2hw3hzfkke0_1_R1C3.TIF", out.path="/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_108", name = "test")

# or, to read the output into R
out_path = '/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_GEOSUD_108'
lyr_name = 'test'

# usage
meanshift.segm(otb.path = "/home/ptaconet/OTB-6.6.0-Linux64/bin", raster.in = "/home/ptaconet/react/MD_SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/SPOT6_2017_HC_BRUT_NC_GEOSUD_PAN_108/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_P_001_A/IMG_SPOT7_P_201710111019265_SEN_SPOT7_20171012_1350221m2hw3hzfkke0_1_R1C3.TIF", out.path=out_path, name = lyr_name)

shp <- readOGR(dsn=out_path, layer = lyr_name, driver = 'ESRI Shapefile')