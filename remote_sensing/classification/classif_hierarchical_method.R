# To ensure that simulations or random objects can be reproduced (more info: http://rfunction.com/archives/62). It must be called each time a new simulation or random object is run
set.seed(1)

# Read the ground truth dataset with zonal statistics 
#ground_truth_df<-read.csv(gsub(".gpkg",".csv",path_to_ground_truth_stats))  
ground_truth_df<-as.data.frame(sf::st_read(path_to_ground_truth_stats))

# Remove useless columns
ground_truth_df$X=NULL
ground_truth_df$cat=NULL
ground_truth_df$geom=NULL


## Get primitives (features) paths for further calculation of the zonal statistics on the segmented objects 
classification_hierarchy_colnames<-c("level_1_fr","level_2_fr","level_3_fr","level_4_fr","level_5_fr")
primitives_column_names<-setdiff(colnames(ground_truth_df),column_names_lc_classes_hierarchy)
#df_primitives_types_sources<-as.data.frame(primitives_column_names,stringsAsFactors=FALSE)
#pattern<-strsplit(methods_to_compute, split=',')[[1]]
#pattern<-paste(pattern,collapse = '|_')
#pattern<-paste0("_",pattern)
#df_primitives_types_sources$indice<-gsub(pattern,"",df_primitives_types_sources$primitives_column_names)
#df_primitives_types_sources<-left_join(df_primitives_types_sources,data.frame(indices_for_classif_labels,indices_for_classif_paths,stringsAsFactors = F),by=c("indice"="indices_for_classif_labels"))

#  Fill the NAs in the ground truth dataset using the na.roughfix function. We do it this way since there are not many NAs in our datasets. Note: we could also use the randomForest::rfImpute function
ground_truth_df <- randomForest::na.roughfix(ground_truth_df)

# Open objects segmentation dataset
dataset_to_classify_sf<-sf::st_read(path_to_segmented_dataset_stats)
segmentation_df<-as.data.frame(dataset_to_classify_sf)
segmentation_df$geom<-NULL
# Fill-in NA values with median value of the column (using randomForest::na.roughfix)
segmentation_df<-randomForest::na.roughfix(segmentation_df)


response_column_name="level_4_fr"

###### Standard approach classif
path_to_outputs_folder=file.path(path_to_classification_folder,paste0("standard__",response_column_name))

## Classify 
classif_standard_rf<-rf_standard_classif_function(ground_truth_df=ground_truth_df,
                                            response_column_name=response_column_name,
                                            primitives_column_names=primitives_column_names,
                                            model_with_optimum_parameters=FALSE,
                                            classify=TRUE,
                                            dataset_to_classify=segmentation_df,
                                            save_output_stats_to_disk=TRUE,
                                            path_to_outputs_folder=path_to_outputs_folder)


## Save output of the classification on the disk
res<-save_classif_to_disk_function(dataset_classified=classif_standard_rf[[7]],
                                   dataset_to_classify_sf=dataset_to_classify_sf,
                                   path_to_outputs_folder=path_to_outputs_folder,
                                   save_objects_raster=FALSE,
                                   save_objects_vector_group_adj_polygons=FALSE)


###### Hierarchical approach classif

path_to_outputs_folder=file.path(path_to_classification_folder,paste0("hierarchical__",response_column_name))

classif_hierarchical_rf<-rf_hierarchical_classif_function(ground_truth_df=ground_truth_df,
                                           classification_hierarchy_colnames=classification_hierarchy_colnames,
                                           classification_level_colname=response_column_name,
                                           primitives_column_names=primitives_column_names,
                                           model_with_optimum_parameters=FALSE,
                                           classify=TRUE,
                                           dataset_to_classify=segmentation_df,
                                           save_output_stats_to_disk=TRUE,
                                           path_to_outputs_folder=path_to_outputs_folder)

## Save output of the classification on the disk
res<-save_classif_to_disk_function(dataset_classified=classif_hierarchical_rf[[5]],
                                   dataset_to_classify_sf=dataset_to_classify_sf,
                                   path_to_outputs_folder=path_to_outputs_folder,
                                   save_objects_raster=FALSE,
                                   save_objects_vector_group_adj_polygons=FALSE)


## Compare results of standard and hierarchical approaches
# TODO

# Functions

###### get_optimum_rf_model: Function to get the optimum parameters for the model.
# - Selects the optimum features to keep in the model (by running the model with the defafult parameters). Feature selection uses the Recursive Feature Elimination (RFE) of the caret package. Useful information about the RFE can be found here: http://topepo.github.io/caret/recursive-feature-elimination.html . Examples of implentations of the RFE can be found here: https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/ and here: https://www.r-bloggers.com/feature-selection-using-the-caret-package/
# - Selects the optimum number of trees (ntrees)
# - Selects the optimum number of variables to test at each division (mtry)
# Input: 
# - ground_truth_df : learning/validation dataset (data.frame). No NA allowed in the classifiers columns
# - response_column_name : column name in the learning/validation dataset to use as response vector
# - get_optimum_primitives : boolean to choose whether to get (TRUE) or not (FALSE) the optimum primitives 
# - verbose : get infos (cat)
# Outputs: 
# - [[1]] : optimum_ntrees : integer optimum ntree
# - [[2]] : optimum_mtry : integer optimum mtry
# - [[3]] : optimum_primitives : vector of optimum primitives  or NULL is get_optimum_primitives is FALSE


get_optimum_rf_model<-function(ground_truth_df, get_optimum_primitives=FALSE,verbose=TRUE){
  
  if (get_optimum_primitives==TRUE){
  if (verbose) { cat("\nGetting the optimum primitives (x) (might be long)...\n") }
  
  ## Define the control using a random forest selection function (functions=rfFuncs). 
  # number: the number of folds to consider. we select the default values which is 10 folds
  control <- caret::rfeControl(functions=rfFuncs, method="cv", number=10, returnResamp = "final")
  ## Run the Recursive Feature Elimination algorithm
  # sizes = a vector of integers corresponding to the number of features that should be retained in the updated model. Here we choose the number of primitives divided by 30
  subsets <- seq(1,ncol(ground_truth_df)-1,round((ncol(ground_truth_df)-1)/30))
  results <- caret::rfe(ground_truth_df[,2:ncol(ground_truth_df)], ground_truth_df[,1], sizes=subsets, rfeControl=control)
  # to summarize the results : print(results) ; to list the chosen features : predictors(results) ; to plot the results : plot(results, type=c("g", "o"))
  
  optimum_primitives<-predictors(results)
  ground_truth_df<-ground_truth_df[,c("response",optimum_primitives)]
  
  } else {
    optimum_primitives<-NULL
  }
  
  # Get optimum ntrees
  if (verbose) { cat("\nGetting the optimum number of trees, setting max to 500 (ntree)\n...") }
  model<-randomForest::randomForest(response ~ ., data=ground_truth_df, ntree=500)
  # to print error OOB in function of the number of trees: plot(model$err.rate[, 1], type = "l", xlab = "nombre d'arbres", ylab = "erreur OOB")
  optimum_ntrees<-which(model$err.rate[, 1]==min(model$err.rate[, 1]))[1]
  # Get optimum mtry
  if(optimum_ntrees>=3){
  if (verbose) { cat("\nGetting the optimum number of variables randomly sampled as candidates at each split (mtry)...") }
  model_tuned<-randomForest::tuneRF(x=ground_truth_df[,2:ncol(ground_truth_df)],y=ground_truth_df$response,ntree=optimum_ntrees,trace = FALSE)
  optimum_mtry<-as.numeric(model_tuned[,1][which(model_tuned[,2]==min(model_tuned[,2]))])
  } else {
  optimum_mtry=floor(sqrt(ncol(ground_truth_df[,2:ncol(ground_truth_df)])))  # default value
  }
  return(list(optimum_ntrees,optimum_mtry,optimum_primitives))
  
}


######## rf_standard_classif_function: Function to perform a classification using random forest with standard (flat) approach

## ground_truth_df : learning/validation dataset (data.frame). No NA allowed in the classifiers columns
## response_column_name : column name in the learning/validation dataset to use as response vector
## primitives_column_names : vector of primitives column names 
## model_with_optimum_parameters : boolean. Create the model with the optimum parameters (mtry and ntrees) ? default FALSE
## classify : boolean. Classify a dataset with the classifier generated? default FALSE
## dataset_to_classify : dataset to classify with the classifier generated (data.frame). Mandatory if dataset_to_classify is set to TRUE
## save_output_stats_to_disk : boolean. Save statistics (confusion matrix, etc.) on the classification to disk (TRUE) or not (FALSE). Default FALSE
## path_to_outputs_folder : path to the folder where outputs will be stored (if save_outputs_to_disk is set to TRUE). Outputs are:
  ## - the data classified (classif.gpkg)
  ## - the confusion matrix (conf_matric.csv)
  ## - a table with kappa and accuracy (classif_stats.csv)
  ## - the dotchart of variable importance as measured by the Random Forest (varImpPlot.png)
  ## - a typical decision tree made with the partykit::ctree function (might be unreadable if too many classes)

rf_standard_classif_function<-function(ground_truth_df,
                                       response_column_name,
                                       primitives_column_names,
                                       model_with_optimum_parameters=FALSE,
                                       classify=FALSE,
                                       dataset_to_classify=NULL,
                                       save_output_stats_to_disk=FALSE,
                                       path_to_outputs_folder=NULL){
  
  # Keep only useful columns in the GT dataset (i.e. column to classify + primitives) and rename column to classify (response)
  ground_truth_df<-ground_truth_df[,c(response_column_name,primitives_column_names)]
  colnames(ground_truth_df)[which(colnames(ground_truth_df)==response_column_name)]<-"response"
  
  if(model_with_optimum_parameters==TRUE){
  # Get optimized parameters for the model
  rf_optimum_parameters<-get_optimum_rf_model(ground_truth_df,verbose = FALSE)
  # Create the model
  model<-randomForest::randomForest(response ~ ., data=ground_truth_df, ntree=rf_optimum_parameters[[1]][1], mtry=rf_optimum_parameters[[2]][1])
  } else {
    model<-randomForest::randomForest(response ~ ., data=ground_truth_df)
    rf_optimum_parameters<-NULL
  }
  # confusion matrix, VarImpPlot, classif statistics, typical decision tree
  conf <- caret::confusionMatrix(data = model$predicted, reference = ground_truth_df$response)
  
  # Confusion matrix
  confusion_matrix<-conf$table
  # Overall stats for the classif (Kappa, accuracy, etc)
  classif_overall_stats<-as.data.frame(conf$overall)
  # Stats by class (f1 score, etc)
  classif_classes_stats<-conf$byClass
  classif_classes_stats<-as.data.frame(classif_classes_stats)
  classif_classes_stats$class<-sub('.*\\: ', '', row.names(classif_classes_stats))
  
  # Variables importance plot
  VarImportancePlot<-randomForest::varImpPlot(model)
  
  # Decision tree 
  model.ctree = partykit::ctree(ground_truth_df[,1] ~ ., data=ground_truth_df[,2:ncol(ground_truth_df)])
  
  if (save_output_stats_to_disk==TRUE){
    cat("Saving outputs to disk...\n")
  # Create directory
  dir.create(path_to_outputs_folder)
  # Save outputs
  write.csv(confusion_matrix,file.path(path_to_outputs_folder,"confusion_matrix.csv"))
  write.csv(classif_overall_stats,file.path(path_to_outputs_folder,"classif_overall_stats.csv"))
  write.csv(classif_classes_stats,file.path(path_to_outputs_folder,"classif_classes_stats.csv"))
  
  png(file.path(path_to_outputs_folder,"varImpPlot.png"))
  randomForest::varImpPlot(model)
  dev.off()
  
  png(file.path(path_to_outputs_folder,"ctree.png"))
  plot(model.ctree)
  dev.off()
  }
  
  if (classify==TRUE){
    cat("Classifing the dataset...\n")
    dataset_to_classify$predicted<-predict(model,dataset_to_classify)
  } else {
    dataset_to_classify<-NULL
  }
  
  return(list(model,rf_optimum_parameters,confusion_matrix,classif_overall_stats,classif_classes_stats,VarImportancePlot,dataset_to_classify))
}


######## rf_hierarchical_classif_function: Function to perform a classification using random forest with hierachical approach

## ground_truth_df : learning/validation dataset (data.frame). No NA allowed in the classifiers columns
## classification_hierarchy_column_names : vector of the column names of the hierarchy, in the order of the classification hierarchy (ie from the most aggregated classification to the less aggregated)
## classification_level_index : index in classification_hierarchy_column_names of the level desired for the classification
## primitives_column_names : vector of primitives column names 
## model_with_optimum_parameters : boolean. Create the model with the optimum parameters (mtry and ntrees) ? default FALSE
## classify : boolean. Classify a dataset with the classifier generated? default FALSE
## dataset_to_classify : dataset to classify with the classifier generated (data.frame). Mandatory if dataset_to_classify is set to TRUE
## save_output_stats_to_disk : boolean. Save statistics (confusion matrix, etc.) on the classification to disk (TRUE) or not (FALSE). Default FALSE
## path_to_outputs_folder : path to the folder where outputs will be stored (if save_outputs_to_disk is set to TRUE).



rf_hierarchical_classif_function<-function(ground_truth_df,
                                           classification_hierarchy_colnames,
                                           classification_level_colname,
                                           primitives_column_names,
                                           model_with_optimum_parameters=FALSE,
                                           classify=FALSE,
                                           dataset_to_classify=NULL,
                                           save_output_stats_to_disk=FALSE,
                                           path_to_outputs_folder=NULL){
 
  list_model<-list()
  list_confusion_matrix<-list()
  classif_overall_stats<-list()
  classif_classes_stats<-NULL
  
  classification_hierarchy_colnames<-classification_hierarchy_colnames[1:which(classification_hierarchy_colnames==classification_level_colname)]
  
  for (i in 1:(length(classification_hierarchy_colnames)-1)){
    cat(paste0("Classifing at level ",i," (",classification_hierarchy_colnames[i],") ...\n"))
    # Model and classif at level 1
    if (i==1){
      ground_truth_df_model<-ground_truth_df[,c(classification_hierarchy_colnames[i],primitives_column_names)]
      colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==classification_hierarchy_colnames[i])]<-"response"
      if (model_with_optimum_parameters==TRUE){
        # Get optimized parameters for the model
        rf_optimum_parameters<-get_optimum_rf_model(ground_truth_df_model,verbose = FALSE)
        # Create the model
        model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model, ntree=rf_optimum_parameters[[1]][1], mtry=rf_optimum_parameters[[2]][1])
      } else {
        model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model)
      }
      if (classify){
      dataset_to_classify$predicted<-predict(model,dataset_to_classify)
      }
    }
    
    # Get child classes 
    classes<-as.character(unique(ground_truth_df[,classification_hierarchy_colnames[i]]))

    if (classify){
    df_classified<-NULL
    }
    
    # Create child classifiers for each mother class (j)
    for (j in 1:length(classes)){
      # filter the GT dataset
      ground_truth_df_model<-ground_truth_df[which(ground_truth_df[,classification_hierarchy_colnames[i]]==as.character(classes[j])),]
      # Filter the dataset to classify
      if (classify){
      df_to_classify<-dataset_to_classify %>% filter (dataset_to_classify$predicted==classes[j])
      }
      # Create the classifier for the child classes of the mother class 
      child_classes<-as.character(unique(ground_truth_df_model[,classification_hierarchy_colnames[i+1]]))
      number_of_child_classes<-length(child_classes)
      if (number_of_child_classes > 1){ # i.e. classify only if there is more than one class to classify
        cat(paste0("       Classifying children classes of mother class : ",classes[j]," (",paste(child_classes,collapse = " / "),") ...\n"))
        ground_truth_df_model<-ground_truth_df_model[,c(classification_hierarchy_colnames[i+1],primitives_column_names)]
        colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==classification_hierarchy_colnames[i+1])]<-"response"
        ground_truth_df_model$response<-as.factor(as.character(ground_truth_df_model$response))
        #if (TRUE %in% unique(apply(ground_truth_df_model, 2, function(x) any(is.na(x))))){
        #  ground_truth_df_model <- randomForest::rfImpute(response ~ ., ground_truth_df_model)
        #}
        
        if (model_with_optimum_parameters==TRUE){
          # Get optimized parameters for the model
          rf_optimum_parameters<-get_optimum_rf_model(ground_truth_df_model,verbose = FALSE)
          # Create the model
          model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model, ntree=rf_optimum_parameters[[1]][1], mtry=rf_optimum_parameters[[2]][1])
        } else {
          model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model)
        }
        
        ## Generate the classif statistics

        # Output of randomForest::model
        list_model<-c(list_model,list(model))
        names(list_model)[length(list_model)]=paste0(classification_hierarchy_colnames[i],"_",classes[j])
        
        conf <- caret::confusionMatrix(data = model$predicted, reference = ground_truth_df_model$response)
        
        # Confusion Matrix
        confusion_matrix<-conf$table
        list_confusion_matrix<-c(list_confusion_matrix,list(confusion_matrix))
        names(list_confusion_matrix)[length(list_confusion_matrix)]=paste0(classification_hierarchy_colnames[i],"_",classes[j])
        # Overall stats for the classif (Kappa, accuracy, etc)
        classif_overall_stats_this_class<-as.data.frame(t(conf$overall))
        classif_overall_stats_this_class$mother_class_name<-classes[j]
        classif_overall_stats_this_class$mother_class_level<-i
        classif_overall_stats_this_class$children_class_names<-paste(as.character(unique(ground_truth_df_model$response)),collapse = ' / ')
        classif_overall_stats_this_class$children_class_level<-i+1
        classif_overall_stats<-rbind(classif_overall_stats,classif_overall_stats_this_class)
        # Stats by class (f1 score, etc)
        conf<-as.data.frame(conf$byClass)
        if (!("conf$byClass" %in% colnames(conf))){
          conf$class<-sub('.*\\: ', '', row.names(conf))
          conf$level<-i+1
          classif_classes_stats<-rbind(classif_classes_stats,conf)
        } else {
          conf<-rbind(t(conf),t(conf))
          conf<-as.data.frame(conf)
          conf<-cbind(conf,class=as.character(unique(ground_truth_df_model$response)))
          conf$level<-i+1
          classif_classes_stats<-rbind(classif_classes_stats,conf)
        }
        
        if (classify){
          df_to_classify$predicted<-NULL
          df_to_classify$predicted<-predict(model,df_to_classify)
        }
        
      } else {
        cat(paste0("       Only one class in the children classes of mother class : ",classes[j],". Hence only transferring the whole mother data to child data without classifying\n"))

        if (classify){
        df_to_classify <- dataset_to_classify %>% filter (dataset_to_classify$predicted==classes[j])
        }
      }
      
      if (classify){
      df_classified<-rbind(df_classified,df_to_classify)
      }
    }
    
    if (classify){
    dataset_to_classify<-df_classified
    }
    
  }

  if (save_output_stats_to_disk==TRUE){
    dir.create(path_to_outputs_folder)
    for (i in 1:length(list_confusion_matrix)){
    write.csv(list_confusion_matrix[[i]],file.path(path_to_outputs_folder,paste0(names(list_confusion_matrix[i]),".csv")))
    png(file.path(path_to_outputs_folder,paste0(names(list_model[i]),".png")))
    varImpPlot(list_model[[i]],main = paste0(names(list_model[i])," \n ",paste(list_model[i][[1]]$classes, collapse = " / ")))
    dev.off()
    }
    write.csv(classif_overall_stats,file.path(path_to_outputs_folder,"classif_overall_stats.csv"))
    write.csv(classif_classes_stats,file.path(path_to_outputs_folder,"classif_classes_stats.csv"))
  }
  
  return(list(list_model,list_confusion_matrix,classif_overall_stats,classif_classes_stats,dataset_to_classify))
}



##### save_classif_to_disk_function: function to save outputs of classification in the disk
## dataset_classified : dataset that has been classified through the predict() function. The predicted column must be named 'predicted'
## dataset_to_classify_sf : dataset to classify (sf class)
## path_to_outputs_folder : path to the folder where outputs will be stored (if save_outputs_to_disk is set to TRUE). Output are the objects save as gpkg file (classification.gpkg). Optional outputs are : 
## save_objects_raster : save the objects as raster (tif) file. default FALSE
## save_objects_vector_group_adj_polygons : save the objects as vector (gpkg) file, where adjacent objects with the same class are grouped. default FALSE. Cannot be TRUE if save_objects_raster is FALSE 

save_classif_to_disk_function<-function(dataset_classified, dataset_to_classify_sf, path_to_outputs_folder, save_objects_raster=FALSE, save_objects_vector_group_adj_polygons=FALSE ){
  
  cat("Saving the classification as vector gpkg...\n")
  dataset_classified<-dataset_classified[,c("cat","predicted")]
  dataset_to_classify_sf<-merge(dataset_to_classify_sf,dataset_classified,by="cat")
  path_to_output_classification_vector<-file.path(path_to_outputs_folder,"classification.gpkg")
  sf::st_write(dataset_to_classify_sf,path_to_output_classification_vector,layer_options = "OVERWRITE=true")
  
  if (save_objects_raster==TRUE){
  cat("Saving the classification as raster tif...\n")
  ## Rasterize the classification
  require(raster)
  require(fasterize)
  # Set output raster characteristics
  poly<-sf::st_read(file.path(path_to_outputs_folder,"classification.gpkg"))
  output_res<-1.633175
  r <- raster(poly, res = output_res)
  # Rasterize using fasterize (fast version of rasterize)
  predicted<-data.frame(predicted=unique(poly$predicted))
  predicted$predicted_integer<-seq(1:nrow(predicted))
  poly<-merge(poly,predicted,by="predicted")
  r <- fasterize::fasterize(poly, r, field = "predicted_integer")
  # Write the classification raster
  path_to_output_classification_raster<-file.path(path_to_outputs_folder,"classification.tif")
  writeRaster(r,path_to_output_classification_raster, overwrite=TRUE, datatype='INT2S')
  }
  
if (save_objects_vector_group_adj_polygons==TRUE){
  cat("Saving the classification as vector gpkg with adjacent objects with the same class grouped...\n")
  # Vectorize back (now that the adjacent polygons that have the same class have been gathered)
  gdal_appli<-paste0("gdal_polygonize.py ",path_to_output_classification_raster," ",gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector)," -b 1 None DN")
  system(gdal_appli)
  # Get back the labels of the classes on the vector version of the classification
  poly<-sf::st_read(gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector))
  poly<-merge(poly,predicted,by.x="DN",by.y="predicted_integer")
  sf::st_write(poly,path_to_output_classification_vector,layer_options = "OVERWRITE=true")
  file.remove(gsub(".gpkg","_temp.gpkg",path_to_output_classification_vector))
  }
  
}






















####### Statistics: standard classif vs hierarchical classif

stats_classif_by_class<-NULL
# Get statistics by class for standard (flat) classif
for (i in 2:(length(classif_hierarchy))){
  ground_truth_df_model<-ground_truth_df[,c(classif_hierarchy[i],column_names_primitives)]
  colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==classif_hierarchy[i])]<-"response"
  model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model)
  
  conf <- caret::confusionMatrix(data = model$predicted, reference = ground_truth_df_model$response)
  conf<-as.data.frame(conf$byClass)
  conf$class<-sub('.*\\: ', '', row.names(conf))
  conf$level<-i
  conf$classif_type<-"standard"
  stats_classif_by_class<-rbind(stats_classif_by_class,conf)
}



# Get statitistics by class for hierarchical classif
list_model<-list()
list_confusion_matrix<-list()
list_conf_stats<-list()
list_conf<-list()

for (i in 1:(length(classif_hierarchy)-1)){
  # Create classifier to discriminate the whole hierarchy using the whole DB ('flat' or classical approach)
  ground_truth_df_model<-ground_truth_df[,c(classif_hierarchy[i],column_names_primitives)]
  colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==classif_hierarchy[i])]<-"response"
  #ground_truth_df_model <- randomForest::rfImpute(response ~ ., ground_truth_df_model)
  model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model)
  varImpPlot(model,main = classif_hierarchy[i])
  
  classes<-as.character(unique(ground_truth_df[,classif_hierarchy[i]]))
  
  # Create classifier for each sub-category (hierarchical approach)
  for (j in 1:length(classes)){
    ground_truth_df_model<-ground_truth_df[which(ground_truth_df[,classif_hierarchy[i]]==as.character(classes[j])),]
    
    if (length(as.character(unique(ground_truth_df_model[,classif_hierarchy[i+1]]))) > 1){ # i.e. classify only if there is more than one class to classify
      ground_truth_df_model<-ground_truth_df_model[,c(classif_hierarchy[i+1],column_names_primitives)]
      colnames(ground_truth_df_model)[which(colnames(ground_truth_df_model)==classif_hierarchy[i+1])]<-"response"
      ground_truth_df_model$response<-as.factor(as.character(ground_truth_df_model$response))
      if (TRUE %in% unique(apply(ground_truth_df_model, 2, function(x) any(is.na(x))))){
        ground_truth_df_model <- randomForest::rfImpute(response ~ ., ground_truth_df_model)
      }
      model<-randomForest::randomForest(response ~ ., data=ground_truth_df_model)
      # Plot the most discrimating primitives
      varImpPlot(model,main=paste0(classif_hierarchy[i+1]," - ",classes[j]," : ",paste(as.character(unique(ground_truth_df_model$response)),collapse = ' / ')))
      # save model in list
      list_model<-c(list_model,list(model))
      names(list_model)[length(names(list_model))]=paste0(classif_hierarchy[i],"_",classes[j])
      # save confusion matrix in list
      list_confusion_matrix<-c(list_confusion_matrix,list(model$confusion))
      names(list_confusion_matrix)[length(names(list_confusion_matrix))]=paste0(classif_hierarchy[i],"_",classes[j])
      # save confusion matrix stats in list
      conf <- caret::confusionMatrix(data = model$predicted, reference = ground_truth_df_model$response)
      list_conf<-c(list_conf,conf)
      list_conf_stats<-c(list_conf_stats,list(as.data.frame(conf$byClass)))
      names(list_conf_stats)[length(names(list_conf_stats))]=paste0(classif_hierarchy[i],"_",classes[j])
      conf<-as.data.frame(conf$byClass)
      if (!("conf$byClass" %in% colnames(conf))){
        conf$class<-sub('.*\\: ', '', row.names(conf))
        conf$level<-i+1
        conf$classif_type<-"hierarchical"
        stats_classif_by_class<-rbind(stats_classif_by_class,conf)
      } else {
        conf<-rbind(t(conf),t(conf))
        conf<-as.data.frame(conf)
        conf<-cbind(conf,class=as.character(unique(ground_truth_df_model$response)))
        conf$level<-i+1
        conf$classif_type<-"hierarchical"
        stats_classif_by_class<-rbind(stats_classif_by_class,conf)
      }
      
    }
    
  }
  
}

row.names(stats_classif_by_class)<-NULL

classes<-unique(stats_classif_by_class$class)
levels<-unique(stats_classif_by_class$level)

stats_diff_f1<-NULL

for (i in 1:length(classes)){
  for (j in 1:length(levels)){
    df<- stats_classif_by_class %>% filter(class==classes[i] & level==levels[j])
    if (nrow(df)>1){
      diff_rf1<-df$F1[which(df$classif_type=="hierarchical")] - df$F1[which(df$classif_type=="standard")]
      stats_diff_f1<-rbind(stats_diff_f1,c(classes[i],levels[j],diff_rf1))
    }
  }
}

stats_diff_f1<-as.data.frame(stats_diff_f1)
colnames(stats_diff_f1)<-c("class","hierarchy_level","diff_f1")

