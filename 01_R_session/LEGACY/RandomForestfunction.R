RandomForestFeatureSelection <- function(table=NULL,set_seed=123, splittestratio=0.5,featuretablesize=50){
library(data.table)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(tidyverse)
library(ranger)
library(Information)
library(InformationValue)

#Perform full analysis of relabComp.R first, you need the R objects "shannon_table" and "samples" from this script.

normalisedtable <- table
#normalisedtable <- as.data.table(t(assay(vsd)))
RFmetadata <- as.data.table(metadata[!(metadata$severity %in% "Unknown")])

RFmetadata$Diagnose <- factor(RFmetadata$Diagnose)
RFmetadata <- droplevels(RFmetadata)



normalisedtable <- as.data.table(sapply(normalisedtable, as.numeric))
# normalisedtable <- apply(normalisedtable,2,function(x){x-mean(x)})


normalisedtable[,UCpositive:=((as.integer(RFmetadata$Diagnose)-1))]
seednr <- set_seed
#seednr <- 222
set.seed(seednr)
splitted <- rsample::initial_split(normalisedtable[,], splittestratio)
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

predictedmodel <- predict(rfmodel, rftest)
predicted <- as.data.table(predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()


#evaluate
optCutOff <- optimalCutoff(rftest[["UCpositive"]], predicted, optimiseFor = "Both")[1]

misClassErrorResult <- misClassError(rftest[["UCpositive"]], predicted, threshold = optCutOff)

ROCplotResult <- plotROC(rftest[["UCpositive"]], predicted, returnSensitivityMat=T)

rftest[["UCpositive"]]
concordanceResult <- InformationValue::Concordance(rftest[["UCpositive"]], predicted)
sensitivityResult <- InformationValue::sensitivity(rftest[["UCpositive"]],predicted , threshold = optCutOff)#factor(as.integer(predicted>optCutOff), levels=c(0,1))
specificityResult <- InformationValue::specificity(rftest[["UCpositive"]], predicted, threshold = optCutOff)
confusionMatrixResult <- InformationValue::confusionMatrix(rftest[["UCpositive"]], predicted)

variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
return_features <- head(variable_importance[order(variable_importance$importance, decreasing=T),],featuretablesize)
return(return_features)
}
# 
# mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]
# 
# #t.test(normalisedtable[,variable_importance[variable_importance$importance==max(variable_importance$importance),taxa],with=F],normalisedtable$UCpositive)
# 
# 
# ggplot(normalisedtable, aes(x=UCpositive,y=ZFP36L2,fill=as.factor(UCpositive)))+geom_boxplot()+geom_dotplot(binaxis="y",stackdir = "center",dotsize=0.1,binwidth = 1/100) + scale_y_continuous(trans="log1p")
# 
# confusionMatrixResult
# #list(model=rfmodel, misClassErrorResult=misClassErrorResult,ROCplotResult=ROCplotResult, concordanceResult=concordanceResult, sensitivityResult=sensitivityResult, specificityResult=specificityResult, confusionMatrixResult=confusionMatrixResult)
# 
# 
# 
# 
# 
# 
# # vsd_validation$severity <- vsd_validation$endoscopic_mayo
# # vsd_validation$supervised <- vsd_validation$endoscopic_mayo
# # vsd_validation$severity <- ordered(replace(vsd_validation$supervised, vsd_validation$supervised%in%NA,"Control"), levels=c("Control","0","1","2","3"))