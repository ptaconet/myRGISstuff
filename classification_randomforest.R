require(randomForest)
require(rgdal)
library(caret)
# Read ground truth dataset with zonal statistics 
ground_truth_df<-as.data.frame(readOGR("/home/ptaconet/Documents/react/data_BF/Ground_truth/ground_truth_stats.gpkg"))

# primitives column indices
#primitives_column_indices<-which(colnames(ground_truth_df)) # TODO

# Plots of correlations between variables:
# variables_indexes<-c(5:10)
# plot (ground_truth_df[variables_indexes])


# Create random forest model

# subset the dataset with columns to use for the model
ground_truth_df_model<-ground_truth_df[,c("type_1",colnames(ground_truth_df[9:ncol(ground_truth_df)]))]

# Run the model with all the columns
model<-randomForest(type_1 ~ ., data=ground_truth_df_model, ntree = 500, na.action = na.omit)

####################################################
##### Feature Selection (variables to keep) #####
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

## RFE does not accept NAs in the predictors. Fill the NAs in the ground truth dataset using the rfImpute function. Note: we could also use the randomForest::na.roughfix function
ground_truth_df.imputed <- randomForest::rfImpute(type_1 ~ ., ground_truth_df_model)
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
subsets <- seq(1,ncol(ground_truth_df.imputed)-1,round((ncol(ground_truth_df.imputed)-1)/20))
results <- caret::rfe(ground_truth_df.imputed[,2:ncol(ground_truth_df.imputed)], ground_truth_df.imputed[,1], sizes=subsets, rfeControl=control)
# summarize the results
#print(results)
# list the chosen features
#predictors(results)
# plot the results
#plot(results, type=c("g", "o"))




# global info about the model
print(model)

# Taux d'erreur: Le taux d'erreur correspond à la proportion de cas où la prédiction est incorrecte
error_rate=1-sum(diag(model$confusion))/sum(model$confusion)
print(error_rate)


# Importance of the variables: L'importance d'une variable dans la classification correspond à la diminution moyenne de l'impureté qu'elle apporte. Pour chaque arbre, la diminution totale de l'impureté liée à une variable correspond à la diminution de l'impureté cumulée sur l'ensemble des noeuds qu'elle régit. Cette diminution est ensuite moyennée sur l'ensemble des arbres. Autrement dit: elle est calculée par l’index de Gini : la diminution pour chaque noeud est cumulée, puis une moyenne sur l’ensemble des arbres est effectuée. 
variables_importance_df<-model$importance
varImpPlot(model)

# Confusion matrix
conf_matrix<-model$confusion

# Votes: Répartition des votes pour chaque individu.
model$votes

# Marge: Différence entre la proportion de votes pour la classe correcte (i) et la proportion de votes pour la classe sortie majoritaire parmi les autres classes (j ≠ i).
# plus la marge est proche de 1 et plus la confiance accordée à la prédiction est grande... Au contraire, quand la marge est faible ou même négative, la confiance à accorder à la classification pour l'individu considéré est faible.
m=margin(model)
print(m)    


#library(caret)
#conf <- confusionMatrix(data = model$predicted, reference = model$type_1)
#conf$byClass["Sensitivity"]
#conf$byClass["Specificity"]
