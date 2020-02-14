# myRGISstuff

Many scripts and functions for GIS stuff with R. Mostly tests


  - *miscellaneous/get_bingmaps_images/download_and_georeference_bingmaps_data.R* : Download and georeference Bing maps data covering an AOI (provided by the user). Outputs are georeferenced TIF files covering the AOI. This script uses the Bing API. More info here: https://msdn.microsoft.com/en-us/library/ff701724.aspx
    
    - *miscellaneous/get_bingmaps_images/create_10km_square_tif_from_bingmaps_data* : Create a 10km² georeferenced TIF + a OGC Geopackage of a given AOI on Earth using Bing maps satellite imagery at best available zoom (19, which corresponds to a resolution of approx. 0.5 cm)
  
