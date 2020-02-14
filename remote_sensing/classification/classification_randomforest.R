require(randomForest)
require(rgdal)
library(caret)

path_to_processing_folder<-"/home/ptaconet/Documents/react/data_BF"  #<Path to the processing folder (i.e. where all the data produced by the workflow will be stored)>
indices_for_classif_labels<-c("DEM",
                              "slope",
                              "accumulation",
                              "NDVI_SPOT6",
                              "NDWI_SPOT6",
                              "BI_SPOT6",
                              "B02_S2",
                              "B03_S2",
                              "B04_S2",
                              "B05_S2",
                              "B06_S2",
                              "B07_S2",
                              "B08_S2",
                              "B8A_S2",
                              "B11_S2",
                              "B12_S2",
                              "B0_SPOT6",
                              "B1_SPOT6",
                              "B2_SPOT6",
                              "B3_SPOT6",
                              "PAN_SPOT6",
                              "text_energy_5",
                              "text_entropy_5",
                              "text_correlation_5",
                              "text_inertia_5",
                              "text_haralickcorellation_5",
                              "text_mean_5",
                              "text_energy_9",
                              "text_entropy_9",
                              "text_correlation_9",
                              "text_inertia_9",
                              "text_haralickcorellation_9",
                              "text_mean_9",
                              "text_energy_17",
                              "text_entropy_17",
                              "text_correlation_17",
                              "text_inertia_17",
                              "text_haralickcorellation_17",
                              "text_mean_17",
                              "NDVI_S2",
                              "NDWI_S2",
                              "BRI_S2",
                              "MNDWI_S2",
                              "MNDVI_S2",
                              "RNDVI_S2"
) # Go to parameter indices_for_classif_paths to set the path to each band
indices_for_classif_paths<-c(file.path(path_to_processing_folder,"DEM_SRTM/processed_data/DEM.TIF"),
                             file.path(path_to_processing_folder,"DEM_SRTM/processed_data/slope.TIF"),
                             file.path(path_to_processing_folder,"DEM_SRTM/processed_data/accumulation.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/NDVI.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/NDWI.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/BI.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B02.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B03.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B04.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B05.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B06.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B07.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B08.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B8A.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B11.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/B12.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PANSHARPEN_3.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/PAN.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_4.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_5_5_7.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_advanced_5_5_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_4.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_9_9_7.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_advanced_9_9_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_0.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_1.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_2.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_4.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_simple_17_17_7.TIF"),
                             file.path(path_to_processing_folder,"VHR_SPOT6/processed_data/HaralickTextures_advanced_17_17_0.TIF"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/NDVI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/NDWI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/BRI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/MNDWI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/MNDVI.tif"),
                             file.path(path_to_processing_folder,"HR_Sentinel2/processed_data/RNDVI.tif")
                             
)

segmented_dataset_path<-"/home/ptaconet/Documents/react/data_BF/Segmentation/segmentation_vector.gpkg"
ground_truth_path<-"/home/ptaconet/Documents/react/data_BF/Ground_truth/ground_truth_stats.csv"
column_names_lc_classes<-c("type_1","type_2","type_3") #<Names of the columns of land cover classes in the ground truth dataset. eg : c("type_1","type_2"). Typically type_1 is the most detailed land cover, type_2 is a more aggregated classification, etc.>
column_name_lc_classification<-"type_1"

# To ensure that simulations or random objects can be reproduced (more info: http://rfunction.com/archives/62). It must be called each time a new simulation or random object is run
set.seed(1)
# Read ground truth dataset with zonal statistics 
ground_truth_df<-read.csv(ground_truth_path)

# Remove useless columns
ground_truth_df$X=NULL
ground_truth_df$cat=NULL


## Set dataframe of primitives types (HSR, VHSR, ancillary), sources (reflectance, spectral_indice, ancillary, texture), stat (average,stddev,other)
column_names_primitives<-setdiff(colnames(ground_truth_df),column_names_lc_classes)

df_primitives_types_sources<-as.data.frame(column_names_primitives,stringsAsFactors=FALSE)
df_primitives_types_sources$type<-NA
df_primitives_types_sources$source<-NA
df_primitives_types_sources$stat<-NA

df_primitives_types_sources$type[which(grepl("S2",df_primitives_types_sources$column_names_primitives))]<-"HSR"
df_primitives_types_sources$type[which(grepl("SPOT6",df_primitives_types_sources$column_names_primitives))]<-"VHSR"
df_primitives_types_sources$type[which(grepl("text",df_primitives_types_sources$column_names_primitives))]<-"VHSR"
df_primitives_types_sources$type[which(is.na(df_primitives_types_sources$type))]<-"ancillary"

df_primitives_types_sources$source[which(grepl("text",df_primitives_types_sources$column_names_primitives))]<-"texture"
df_primitives_types_sources$source[which(grepl("NDVI|NDWI|BI|BRI|MNDVI|MNDWI|RNDVI",df_primitives_types_sources$column_names_primitives))]<-"spectral_indice"
df_primitives_types_sources$source[which(grepl("shape|dist_to_hydro|DEM|slope|accumulation",df_primitives_types_sources$column_names_primitives))]<-"ancillary"
df_primitives_types_sources$source[which(is.na(df_primitives_types_sources$source))]<-"reflectance"

df_primitives_types_sources$stat[which(grepl("average",df_primitives_types_sources$column_names_primitives))]<-"average"
df_primitives_types_sources$stat[which(grepl("stddev",df_primitives_types_sources$column_names_primitives))]<-"stddev"
df_primitives_types_sources$stat[which(is.na(df_primitives_types_sources$stat))]<-"other"

df_primitives_types_sources$indice<-gsub("_average|_stddev","",df_primitives_types_sources$column_names_primitives)

df_primitives_types_sources<-left_join(df_primitives_types_sources,data.frame(indices_for_classif_labels,indices_for_classif_paths,stringsAsFactors = F),by=c("indice"="indices_for_classif_labels"))

# Plots of correlations between the variables:
# variables_indexes<-c(5:10)
# plot (ground_truth_df[variables_indexes])


######### Model construction workflow (make a function out of it) : 
##### Step 1.1: Select the features to keep in the model (by running the model with the default parameters)
##### Step 1.2: Fine-tune the parameters of the model with the features selected (number of trees (ntrees), number of variables to test at each division (mtry))
## Output of the function: object of class randomForest


######### Hierarchical approach classification : We split the ground truth database and then we run the model construction function at these various levels

########### Classification workflow 
##### Step 2.1: Get the zonal statistics of the objects to classify (outputs of the segmentation)
##### Step 2.2: Run the model on the segmented objects 

########## Outputs : 
# Objects classified
# Plots of the importance of the variables 
# Confusion matrix
# Global information about the model (print(model))

######## Other business to test :
# Differences in the model error rates in function of the use of :
  # the source of the primitive (reflectance, spectral_indice, ancillary, texture)  ?
  # the nature of the primitives (HRS, VHRS, ancillary) (e.g. when the VHRS are used vs. when they are not) ?






######### Hierachical classif
## 1) check which variables are discriminative at each level of classif
## 2) Extract f-scores with a hierarchical method
classif_hierarchy<-c("type_2","type_3","type_1")

for (i in 1:(length(classif_hierarchy)-1)){
  # Create classifier to discriminate the whole hierarchy
  #ground_truth_df_model<-ground_truth_df[,c(classif_hierarchy[i],column_names_primitives)]
  #colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==classif_hierarchy[i])]<-"response"
  #ground_truth_df_model <- randomForest::rfImpute(response ~ ., ground_truth_df_model)
  #model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model)
  #varImpPlot(model)
  
  classes<-as.character(unique(ground_truth_df[,classif_hierarchy[i]]))
  
  # Take each class of the hierarchy 
  for (j in 1:length(classes)){
    ground_truth_df_model<-ground_truth_df[which(ground_truth_df[,classif_hierarchy[i]]==as.character(classes[j])),]
    
    if (length(as.character(unique(ground_truth_df_model[,classif_hierarchy[i+1]]))) > 1){ # i.e. If there is more than one class to classify
      ground_truth_df_model<-ground_truth_df_model[,c(classif_hierarchy[i+1],column_names_primitives)]
      colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==classif_hierarchy[i+1])]<-"response"
      ground_truth_df_model$response<-as.factor(as.character(ground_truth_df_model$response))
      if (TRUE %in% unique(apply(ground_truth_df_model, 2, function(x) any(is.na(x))))){
        ground_truth_df_model <- randomForest::rfImpute(response ~ ., ground_truth_df_model)
      }
      model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model)
      varImpPlot(model,main=paste0(classes[j]," : ",paste(as.character(unique(ground_truth_df_model$response)),collapse = ' / ')))
    }
    
  }
  
}





# Model df: data frame whose first column is the response vector (named "response") and the other columns are the classifiers
get_optimum_rf_model<-function(model_df){

## Step 1.1
# We first run the model with all the primitives and with the randomForest default values

cat("\nGetting the model with all the primitives and default values...")
model<-randomForest::randomForest(response ~ ., data=model_df, na.action = na.omit)

## Get optimum primitives
####################################################
##### Feature Selection (Selection of the primitives to keep) #####
####################################################


## Solution 1: Look for the variables that are correlated and remove them
# Get correlation matrix of the variables (i.e. how the variables are correlated). Argument use="complete.obs" tells the correlation to ignore the NAs
#correlationMatrix <- stats::cor(ground_truth_df_model[,2:ncol(ground_truth_df_model)],use = "complete.obs")
# Find attributes that are highly corrected (correlation >0.8). Use names = TRUE to get column names instead of column indices
#highlyCorrelated <- caret::findCorrelation(correlationMatrix, cutoff=0.8, names = TRUE)
# Remove the higly correlated columns from the dataset, re-run the model and compare error rate
#ground_truth_df_model<-ground_truth_df_model[-which(colnames(ground_truth_df_model) %in% highlyCorrelated )]
#model2<-randomForest(type_1 ~ ., data=ground_truth_df_model, ntree = 500, na.action = na.omit)


### Solution 2: Recursive Feature Elimination (RFE) using the caret::rfe function
## Useful information about the RFE can be found here: http://topepo.github.io/caret/recursive-feature-elimination.html
## Examples of implentations of the RFE can be found here: https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/ and here: https://www.r-bloggers.com/feature-selection-using-the-caret-package/

cat("\nGetting the optimum primitives (x) (might be long)...\n")

## RFE does not accept NAs in the predictors. Fill the NAs in the ground truth dataset using the rfImpute function. Note: we could also use the randomForest::na.roughfix function
ground_truth_df.imputed <- randomForest::rfImpute(response ~ ., ground_truth_df_model)
## Define the control using a random forest selection function (functions=rfFuncs)
# Method: boot or cv
# Number: the number of folds to consider. we select the default values which is 10 folds
control <- caret::rfeControl(functions=rfFuncs, method="cv", number=10, returnResamp = "final")
## Parallelisation
# control$workers <- numCores-1
# control$computeFunction <- mclapply
# control$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
## Run the RFE algorithm
# sizes = a vector of integers corresponding to the number of features that should be retained in the updated model.
subsets <- seq(1,ncol(ground_truth_df.imputed)-1,round((ncol(ground_truth_df.imputed)-1)/30))
results <- caret::rfe(ground_truth_df.imputed[,2:ncol(ground_truth_df.imputed)], ground_truth_df.imputed[,1], sizes=subsets, rfeControl=control)
# summarize the results
#print(results)
# list the chosen features
#predictors(results)
# plot the results
#plot(results, type=c("g", "o"))

optimum_primitives<-predictors(results)

#### Step 1.2: Fine-tune model parameters with features retained
ground_truth_df_features_optimum_primitives<-model_df[,c("response",optimum_primitives)]
# Get optimum ntrees
cat("\nGetting the optimum number of trees, setting max to 500 (ntree)\n...")
model<-randomForest::randomForest(response ~ ., data=ground_truth_df_features_optimum_primitives, ntree=500,  na.action = na.omit)
# to print error OOB in function of the number of trees: plot(model$err.rate[, 1], type = "l", xlab = "nombre d'arbres", ylab = "erreur OOB")
optimum_ntrees<-which(model$err.rate[, 1]==min(model$err.rate[, 1]))[1]
# Get optimum mtry
cat("\nGetting the optimum number of variables randomly sampled as candidates at each split (mtry)...")

ground_truth_df.imputed <- randomForest::rfImpute(response ~ ., ground_truth_df_features_optimum_primitives)
model_tuned<-randomForest::tuneRF(x=ground_truth_df.imputed[,predictors(results)],y=ground_truth_df.imputed$response,ntree=optimum_ntrees,na.action = na.omit)
optimum_mtry<-as.numeric(model_tuned[,1][which(model_tuned[,2]==min(model_tuned[,2]))])

return(list(optimum_primitives,optimum_ntrees,optimum_mtry))

}

## Implement standard approach :
## We create a table df_error_rates where we put the error rates 
df_error_rates<-data.frame(approach_type=character(),class=character(),error_rate=numeric())
for (i in 1:length(column_names_lc_classes)){
  ground_truth_df_model<-ground_truth_df[,c(column_names_lc_classes[i],column_names_primitives)]
  colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==column_names_lc_classes[i])]<-"response"
  model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model, na.action = na.omit)
  
  ## Get min error_rate
  error_rate <- min(model$err.rate[, 1])[1]
  
  df_error_rates<-rbind(df_error_rates,data.frame(approach_type="standard",class=column_names_lc_classes[i],error_rate=error_rate))
}



## On fait la classif au niveau 1
column_name_lc_classification="type_1"
ground_truth_df_model<-ground_truth_df[,c(column_name_lc_classification,column_names_primitives)]
colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==column_name_lc_classification)]<-"response"
model_param<-get_optimum_rf_model(ground_truth_df_model)

# We retrieve the paths of the features to keep
columns_to_keep<-model_param[[1]]
paths<-left_join(data.frame(columns_to_keep),df_primitives_types_sources,by = c("columns_to_keep"="column_names_primitives"))
paths<-unique(data.frame(paths$indice,paths$indices_for_classif_paths))
paths_indices_to_compute<-paths[which(!is.na(paths$paths.indices_for_classif_paths)),]

# We compute the zonal indices on the objects segmented

poly<-sf::st_read(segmented_dataset_path)

if ("shape_schwartzberg" %in% columns_to_keep){
  poly$shape_schwartzberg<-schwartzberg(poly)
}
if ("shape_reock" %in% columns_to_keep){
  poly$shape_reock<-reock(poly)
}
sf::st_write(poly,segmented_dataset_path,layer_options = "OVERWRITE=true")

execGRASS("v.in.ogr", flags=c("o","overwrite"), parameters=list(input=segmented_dataset_path, output="segmentation",min_area=0.0001, snap=-1.0))
if ("dist_to_hydro" %in% columns_to_keep){
  execGRASS("v.in.ogr", flags=c("o","overwrite"), parameters=list(input=path_to_output_accumulation_threshold, output="accumulation",min_area=0.0001, snap=-1.0))
  execGRASS("v.db.addcolumn", parameters=list(map="segmentation", columns="dist_to_hydro double"))
  execGRASS("v.distance", flags=c("overwrite"), parameters=list(from="segmentation", from_type="point,line,area", to="accumulation",to_type="point,line,area",dmax=-1,dmin=-1,upload="dist",column="dist_to_hydro",output="gt_stats_updated"))
}

for (i in 1:length(paths_indices_to_compute)){
  cat(paste0("Computing zonal statistics for indice ",paths_indices_to_compute$paths.indice[i],"\n"))
  execGRASS("r.external", flags="overwrite", parameters=list(input=paths_indices_to_compute$paths.indices_for_classif_paths[i], output="tmprast",band=1))
  execGRASS("g.region", parameters=list(raster="tmprast")) 
  execGRASS("v.rast.stats", flags=c("c","verbose"), parameters=list(map="segmentation", raster="tmprast",column_prefix=paths_indices_to_compute$paths.indice[i],method=methods_to_compute,percentile=90))
  ### ligne ci dessous a décocher quand on calculera les stats sur les objets issus de la segmentation : 
  #execGRASS("v.rast.stats", flags=c("c","verbose"), parameters=list(map="segmented_dataset", raster="tmprast",column_prefix=indices_for_classif_labels[i],method=methods_to_compute,percentile=90))
  execGRASS("g.remove", flags="f", parameters=list(type="raster",name="tmprast"))
}


# Save file as geopackage and csv
writeOGR(readVECT("segmentation"),path_to_ground_truth_stats,driver = "GPKG",layer="ground_truth_stats",overwrite_layer = TRUE)
write.csv(as.data.frame(readVECT("segmentation")),gsub(".gpkg",".csv",path_to_ground_truth_stats))

# Predict the class
predict(model,)


## Implement hierarchical approach : TODO
#classify first level
#column_name_lc_classification="type_1"
#ground_truth_df_model<-ground_truth_df[,c(column_name_lc_classification,column_names_primitives)]

#colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==column_name_lc_classification)]<-"response"

#train <- as.data.frame(ground_truth_df_model %>% group_by(response) %>% sample_frac(size = 0.75))
#test <- anti_join(ground_truth_df_model, train)

#model_param<-get_optimum_rf_model(train)

#train<-train[,c("response",model_param[[1]])]
#train<-randomForest::rfImpute(response ~ ., train)

#model<-randomForest::randomForest(response ~ ., data=train, ntree=model_param[[2]] , mtry=model_param[[3]] )

#test$predicted<-predict(model,test)
#table(test$predicted, test$response) 
#conf <- caret::confusionMatrix(data = test$predicted, reference = test$Titulaire)






# global info about the model
print(model)

# Taux d'erreur: Le taux d'erreur correspond à la proportion de cas où la prédiction est incorrecte
error_rate=1-sum(diag(model$confusion))/sum(model$confusion)
print(error_rate)


# Importance of the variables: L'importance d'une variable dans la classification correspond à la diminution moyenne de l'impureté qu'elle apporte. Pour chaque arbre, la diminution totale de l'impureté liée à une variable correspond à la diminution de l'impureté cumulée sur l'ensemble des noeuds qu'elle régit. Cette diminution est ensuite moyennée sur l'ensemble des arbres. Autrement dit: elle est calculée par l’index de Gini : la diminution pour chaque noeud est cumulée, puis une moyenne sur l’ensemble des arbres est effectuée. 
variables_importance_df<-model$importance
randomForest::varImpPlot(model)

# Confusion matrix
conf_matrix<-model$confusion

# Votes: Répartition des votes pour chaque individu.
model$votes

# Marge: Différence entre la proportion de votes pour la classe correcte (i) et la proportion de votes pour la classe sortie majoritaire parmi les autres classes (j ≠ i).
# plus la marge est proche de 1 et plus la confiance accordée à la prédiction est grande... Au contraire, quand la marge est faible ou même négative, la confiance à accorder à la classification pour l'individu considéré est faible.
m=randomForest::margin(model)
print(m)    


#library(caret)
#conf <- confusionMatrix(data = model$predicted, reference = model$type_1)
#conf$byClass["Sensitivity"]
#conf$byClass["Specificity"]
