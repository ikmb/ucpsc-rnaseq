library(data.table)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(tidyverse)
library(ranger)
library(Information)
library(InformationValue)

#Perform full analysis of relabComp.R first, you need the R objects "shannon_table" and "samples" from this script.


abtable <- as.data.table(t(assay(vsd)))
abmeta <- as.data.table(metadata[!(metadata$severity %in% "Unknown")])

abmeta$Diagnose <- factor(abmeta$Diagnose)
abmeta <- droplevels(abmeta)



abtable <- as.data.table(sapply(abtable, as.numeric))
abtable <- apply(abtable,2,function(x){x-mean(x)})

##########FOR VALIDATION COHORT
vsd_validation <- getCohortsVSD(cohortname = "Mo_UC")
# dds_validation <- getCohortsDDS(cohortname = "Mo_UC")
# pcaExplorer(dds_validation,vsd_validation)
validation_counts <- data.table(t(assay(vsd_validation)))
abtable <- data.table(abtable[,colnames(abtable)%in% colnames(validation_counts)])
validation_counts <- validation_counts[,colnames(validation_counts) %in% colnames(abtable),with=FALSE]
validation_counts <- validation_counts[,colnames(abtable),with=FALSE]

validation_counts <- data.table(apply(validation_counts,2,function(x){x-mean(x)}))#[!as.logical(as.integer(vsd_validation$Diagnose)-1)]
##########


abtable[,UCpositive:=((as.integer(abmeta$Diagnose)-1))]

seednr <- 222
set.seed(seednr)
splitted <- rsample::initial_split(abtable[,c(unique(rfeResults$variables$var),"UCpositive"),with=FALSE], 0.5)
rftrain <- training(splitted)
rftest <- testing(splitted)

#infotables_training <- create_infotables(data = rftrain[,], valid = NULL, y="UCpositive")
#infotables_test <- create_infotables(data = rftest[,], valid = NULL, y="UCpositive")

rftrain$UCpositive <- as.factor(rftrain$UCpositive)
rftest$UCpositive <- as.factor(rftest$UCpositive)


rfmodel <- ranger(
  formula         = UCpositive ~ ., 
  data            = rftrain, 
  num.trees       = 1000,
  seed            = seednr,
  probability     = TRUE,
  importance      = 'impurity'
)
#######
validation_counts[,UCpositive:=((as.integer(vsd_validation$Diagnose)-1))]
validation_counts$UCpositive <- factor(validation_counts$UCpositive)
predictedmodel_validation <- predict(rfmodel, validation_counts[,c(unique(rfeResults$variables$var),"UCpositive"),with=FALSE])
predicted_validation <- as.data.table(predictedmodel_validation$predictions)[,2]%>%unlist()%>%as.vector()

optCutOff <- optimalCutoff(validation_counts[["UCpositive"]], predicted_validation, optimiseFor = "Both")[1]

misClassErrorResult <- misClassError(validation_counts[["UCpositive"]], predicted_validation, threshold = optCutOff)

ROCplotResult <- plotROC(validation_counts[["UCpositive"]], predicted_validation, returnSensitivityMat=T)

concordanceResult <- Concordance(validation_counts[["UCpositive"]], predicted_validation)
sensitivityResult <- sensitivity(validation_counts[["UCpositive"]], predicted_validation, threshold = optCutOff)
specificityResult <- specificity(validation_counts[["UCpositive"]], predicted_validation, threshold = optCutOff)
confusionMatrixResult <- confusionMatrix(validation_counts[["UCpositive"]], predicted_validation, threshold = optCutOff)

variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
top_importance <- head(variable_importance[order(variable_importance$importance, decreasing=T),],50)
top_importance$feature[1:25]


mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]
#######


predictedmodel <- predict(rfmodel, rftest)
predicted <- as.data.table(predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()


#evaluate
optCutOff <- optimalCutoff(rftest[["UCpositive"]], predicted, optimiseFor = "Both")[1]

misClassErrorResult <- misClassError(rftest[["UCpositive"]], predicted, threshold = optCutOff)

ROCplotResult <- plotROC(rftest[["UCpositive"]], predicted, returnSensitivityMat=T)


concordanceResult <- Concordance(rftest[["UCpositive"]], predicted)
sensitivityResult <- sensitivity(rftest[["UCpositive"]], predicted, threshold = optCutOff)
specificityResult <- specificity(rftest[["UCpositive"]], predicted, threshold = optCutOff)
confusionMatrixResult <- confusionMatrix(rftest[["UCpositive"]], predicted, threshold = optCutOff)

variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
head(variable_importance[order(variable_importance$importance, decreasing=T),],25)

mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]

#t.test(abtable[,variable_importance[variable_importance$importance==max(variable_importance$importance),taxa],with=F],abtable$UCpositive)


ggplot(abtable, aes(x=UCpositive,y=ZFP36L2,fill=as.factor(UCpositive)))+geom_boxplot()+geom_dotplot(binaxis="y",stackdir = "center",dotsize=0.1,binwidth = 1/100) + scale_y_continuous(trans="log1p")

confusionMatrixResult
#list(model=rfmodel, misClassErrorResult=misClassErrorResult,ROCplotResult=ROCplotResult, concordanceResult=concordanceResult, sensitivityResult=sensitivityResult, specificityResult=specificityResult, confusionMatrixResult=confusionMatrixResult)






# vsd_validation$severity <- vsd_validation$endoscopic_mayo
# vsd_validation$supervised <- vsd_validation$endoscopic_mayo
# vsd_validation$severity <- ordered(replace(vsd_validation$supervised, vsd_validation$supervised%in%NA,"Control"), levels=c("Control","0","1","2","3"))