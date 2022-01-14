setwd(paste0("/ssdpool/ewacker/RNAseq","/01_R_session"))

library(data.table)
library(plyr)
library(tidyverse)
library(DESeq2)
#library(rvcheck)
# library(CEMiTool)
# library(rmarkdown)
library(RColorBrewer)
# library(gridExtra)
# library(ggpubr)
# library(ggridges)
# library(ggplotify)
# library(extrafont)
# library(reshape2)
# library(viridis)
# library(ggrepel)
# library(ggvenn)
library(rsample)      # data splitting 
# library(randomForest) # basic implementation
library(ranger)
library(Information)
library(InformationValue)
source("functions.R")

seednr <- 555555

# projectdir <- dirname(getwd())
projectdir <- "/ssdpool/ewacker/RNAseq"

#####Inputs#####PSC
corhorttablelocation <- paste0(projectdir,"/00_RawData/cohorts.csv")

cohorttable <- fread(corhorttablelocation, header=T)
cohortname <- "PSC"
Disease <- unlist(cohorttable[Cohort %in% cohortname, 5])

gmt_file <- paste0(projectdir,"/data_rescources/GO_Biological_Process_2018.gmt")
output_location <- paste0("output")

PSC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="PSC")
# vsd <- getnormalizedDESeq2object(corhorttablelocation=corhorttablelocation,cohortname=cohortname)
PSC_metadata <- getmetadata(corhorttablelocation=corhorttablelocation,cohortname="PSC")
PSC_metadata <- PSC_metadata[PSC_metadata$Diagnose %in% c("PSC","Control"),]


PSC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="PSC")
PSC_rawcounts <- PSC_rawcounts[,colnames(PSC_rawcounts) %in% PSC_metadata$SampleID]
UC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="Our")
UC_metadata <- getmetadata(corhorttablelocation=corhorttablelocation,cohortname="Our")


SampleIDs <- c(colnames(PSC_rawcounts),colnames(UC_rawcounts))
merged_rawcounts <- rbindlist(list(data.frame(t(PSC_rawcounts)),data.frame(t(UC_rawcounts))),fill=TRUE)
intersectedgenes <- intersect(colnames(data.frame(t(PSC_rawcounts))), colnames(data.frame(t(UC_rawcounts))))
merged_rawcounts <- merged_rawcounts[,..intersectedgenes]

#Hard filtering to exclude zero expression in any sample genes:
merged_rawcounts <- merged_rawcounts[,!(apply(merged_rawcounts <= 1, 2, any)),with=FALSE]


merged_rawcounts <- t(merged_rawcounts)
colnames(merged_rawcounts) <- SampleIDs

merged_metadata <- rbind(PSC_metadata[,c("SampleID","Diagnose","PlateNr")],UC_metadata[,c("SampleID","Diagnose","PlateNr")])
merged_metadata$Cohort <- factor(c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))

dds <- DESeqDataSetFromMatrix(countData=merged_rawcounts, colData=merged_metadata, design= ~Diagnose + PlateNr)
dds <- DESeq(dds)
#res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
#resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
vsd <- vst(dds, blind=FALSE)

merged_vst_counts <- as.data.table(sapply(data.table(assay(vsd)), as.numeric))
#merged_controls <- ifelse(merged_metadata$Diagnose == "Control",TRUE,FALSE)


##### Validation VST-Data
test_vst_counts <- as.data.table(t(merged_vst_counts))
colnames(test_vst_counts) <- rownames(assay(vsd))
rownames(test_vst_counts) <- SampleIDs
ggplot(test_vst_counts)+geom_jitter(aes(x=merged_metadata$Diagnose,y=HIF1A))+facet_wrap(~merged_metadata$Cohort)
#####

merged_controls_PSC <- ifelse(merged_metadata[grepl("^I",merged_metadata$SampleID)]$Diagnose == "Control",TRUE,FALSE)
merged_controls_UC <- ifelse(merged_metadata[grepl("^DE",merged_metadata$SampleID)]$Diagnose == "Control",TRUE,FALSE)

##Perform outlierfiltered transformation based on controls:
# Legacy version
# UC_z_con_stabilised <- data.table(apply(merged_vst_counts[,grepl("^DE",colnames(merged_vst_counts)),with=FALSE],
#                                          1,function(x){y <- x[merged_controls_UC];z <- y[remove_outlier_filter(y)];
#                                          (x-mean(y))/sd(z)}),
#                                    keep.rownames=TRUE)
PSC_z_con_stabilised <- outlierfiltered_control_transformation(dataset = merged_vst_counts[,grepl("^I",colnames(merged_vst_counts)),with=FALSE],
                                                               controlsvector = merged_controls_PSC)
UC_z_con_stabilised <- outlierfiltered_control_transformation(dataset = merged_vst_counts[,grepl("^DE",colnames(merged_vst_counts)),with=FALSE],
                                                               controlsvector = merged_controls_UC)

merged_z_con_stabilised <- data.table(rn=rownames(assay(vsd)),PSC_z_con_stabilised,UC_z_con_stabilised)
merged_z_con_stabilised <- transpose_datatable(merged_z_con_stabilised)

merged_RF <- data.table(merged_z_con_stabilised, Diagnose=merged_metadata$Diagnose, Cohort=c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))
# ggplot(merged_RF[Diagnose == "Control",])+geom_jitter(aes(x=Cohort, y=MYL12A))

#Convert Diagnose into a factor: 0=Control, 1=PSC, 2=UC
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
levels(merged_RF$Diagnose) <- c("0","1","2")

merged_RF$Cohort  <- factor(merged_RF$Cohort )

###### Check for Controls in different Cohorts####
set.seed(seednr)
splitted <- rsample::initial_split(merged_RF[merged_RF$Diagnose == 0,-c("rn","Diagnose")], 0.5)#,"DiagnoseCohort"
rftrain <- training(splitted)
rftest <- testing(splitted)

#check whether RF can distinguish between Cohorts in our merged dataset (it should not!)
Resultitem <- "Cohort"

rfmodel <- ranger(
  formula         = Cohort ~ ., 
  data            = rftrain, 
  num.trees       = 10000,
  seed            = seednr,
  probability     = TRUE,
  importance      = 'impurity'
)

predictedmodel <- predict(rfmodel, rftest, type="response")
predicted <- as.data.table(predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()

# #evaluate
optCutOff <- optimalCutoff(rftest[[Resultitem]], predicted, optimiseFor = "Both")[1]
misClassErrorResult <- misClassError(rftest[[Resultitem]], predicted, threshold = optCutOff)
ROCplotResult <- plotROC(rftest[[Resultitem]], predicted, returnSensitivityMat=T)

concordanceResult <- Concordance(rftest[[Resultitem]], predicted)
sensitivityResult <- sensitivity(rftest[[Resultitem]], predicted, threshold = optCutOff)
specificityResult <- specificity(rftest[[Resultitem]], predicted, threshold = optCutOff)
confusionMatrixResult <- confusionMatrix(rftest[[Resultitem]], predicted, threshold = optCutOff)
confusionMatrixResult

variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
head(variable_importance[order(variable_importance$importance, decreasing=T),],25)
mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]

#ggplot(merged_RF[merged_RF$Diagnose ==0,])+geom_jitter(aes(x=Cohort,y=RPS28),height=0)

#####DIAGNOSE####
#Convert Diagnose into a factor: 0=Control, 1=PSC, 2=UC
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
levels(merged_RF$Diagnose) <- c("0","1","2")

merged_RF$Cohort  <- factor(merged_RF$Cohort )

filtered_merged_RF <- merged_RF[merged_RF$Cohort == 1,]
levels(filtered_merged_RF$Diagnose) <- droplevels(factor(filtered_merged_RF$Diagnose))
levels(filtered_merged_RF$Diagnose) <- c("0", "1")#0 = Control, 1= UC

# filtered_merged_RF[filtered_merged_RF$Diagnose == 0 & filtered_merged_RF$Cohort ==0,DiagnoseCohort:=0]
# filtered_merged_RF[filtered_merged_RF$Diagnose == 0 & filtered_merged_RF$Cohort ==1,DiagnoseCohort:=1]
# filtered_merged_RF[filtered_merged_RF$Diagnose == 1 & filtered_merged_RF$Cohort ==0,DiagnoseCohort:=2]
# filtered_merged_RF[filtered_merged_RF$Diagnose == 2 & filtered_merged_RF$Cohort ==1,DiagnoseCohort:=3]
# filtered_merged_RF$DiagnoseCohort <- factor(filtered_merged_RF$DiagnoseCohort)
#levels(rftest$DiagnoseCohort) <- c("0","1","1")
# table(filtered_merged_RF$Diagnose, filtered_merged_RF$DiagnoseCohort)

set.seed(seednr)
splitted <- rsample::initial_split(filtered_merged_RF[,-c("Cohort","rn")], 0.5)
rftrain <- training(splitted)
rftest <- testing(splitted)
# levels(rftrain$Diagnose) <- c("0","1","1")
# levels(rftest$Diagnose) <- c("0","1","1")
#rftest$DiagnoseCohort <- factor(rftest$DiagnoseCohort)

#case.weights calculation
training_sample_weights <- case_when(rftrain$Diagnose == "0" ~ 1/sum(rftrain$Diagnose=="0"),
                                     rftrain$Diagnose == "1" ~ 1/sum(rftrain$Diagnose=="1"),
                                     rftrain$Diagnose == "2" ~ 1/sum(rftrain$Diagnose=="2"),
                                     rftrain$Diagnose == "3" ~ 1/sum(rftrain$Diagnose=="3"))

rfmodel <- ranger(
  formula         = Diagnose ~ ., 
  data            = rftrain, 
  num.trees       = 10000,
  seed            = seednr,
  probability     = TRUE,
  importance      = 'impurity',
  classification = TRUE,
  #class.weights = c(1/45,1/114,1/108,1/256),
  case.weights = training_sample_weights
)

predictedmodel <- predict(rfmodel, rftest, type="response")
predicted <- as.data.table(predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()
# predicted <- predictions(predictedmodel,type="response")

# colnames(predicted) <- c("0","1","2","3")
# colnames(predicted) <- c("1","2","3")
# predicted_table<-colnames(predicted)[apply(predicted,1,which.max)]
#  table(predicted_table,rftest$Diagnose)
#compare uc to others
# predicted[predicted=="2"] <- "0"
#evaluate
optCutOff <- optimalCutoff(rftest[["Diagnose"]], predicted, optimiseFor = "Both")[1]

misClassErrorResult <- misClassError(rftest[["Diagnose"]], predicted)#, threshold = optCutOff)

ROCplotResult <- plotROC(rftest[["Diagnose"]], predicted, returnSensitivityMat=T)


concordanceResult <- Concordance(rftest[["Diagnose"]], predicted)
sensitivityResult <- sensitivity(rftest[["Diagnose"]], predicted, threshold = optCutOff)
specificityResult <- specificity(rftest[["Diagnose"]], predicted, threshold = optCutOff)
confusionMatrixResult <- confusionMatrix(rftest[["Diagnose"]], predicted, threshold = optCutOff)
confusionMatrixResult
variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
top25 <- head(variable_importance[order(variable_importance$importance, decreasing=T),],25)

mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]

# ggplot(filtered_merged_RF)+geom_jitter(aes(x=Diagnose,y=ZFP36L2),height=0)
# ggplot(test_vst_counts)+geom_jitter(aes(x=merged_metadata$Diagnose,y=ZFP36L2))+facet_wrap(~merged_metadata$Cohort)


####Validation####
validation_cohortname<- "Mo_UC"
validationcohorts <- list(Ostrowski ="UCAI Based\nDisease Severity",
                          Mo_UC="Diagnose",
                          Planell="Endoscopic\nMayo Score")
vsd_validation <- getCohortsVSD(cohortname=validation_cohortname)
severityscale <- validationcohorts[[validation_cohortname]]
dds_validation <- getCohortsDDS(cohortname=validation_cohortname)
vsd_validation$severity <- vsd_validation$Diagnose
vsd_validation$severity[is.na(vsd_validation$severity)] <- "Control"  


res_validation <- results(dds_validation, c("Diagnose","Control","UC"),tidy=TRUE)
res_validation[res_validation$row %in% top25$feature,]

validation_controls <- ifelse(dds_validation$Diagnose == "Control",TRUE,FALSE)
validation_vst_counts <- data.table(assay(vsd_validation)[,])
# colnames(validation_vst_counts)
# rownames(validation_vst_counts)

# validation_z_con_stabilised <- data.table(apply(validation_vst_counts,
#                  1,function(x){y <- x[validation_controls];z <- y[remove_outlier_filter(y)];
#                  (x-mean(y))/sd(z)}),
#                 keep.rownames=TRUE)

validation_z_con_stabilised <- outlierfiltered_control_transformation(validation_vst_counts, validation_controls)

validation_z_con_stabilised <- data.table(rn=rownames(assay(vsd_validation)), validation_z_con_stabilised)
validation_z_con_stabilised <- transpose_datatable(validation_z_con_stabilised)


#colnames(validation_z_con_stabilised) <- c("rn",rownames(assay(vsd_validation)))

# ggplot(as.data.table(t(assay(vsd_validation)[,])),aes(x=colData(vsd_validation)$severity,y=S100A6))+geom_boxplot()+geom_jitter(width=0.05)


####RF for validation####
Mo_genefeatures <- c("rn",rownames(assay(vsd_validation)))
sum(colnames(merged_RF) %in% Mo_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo
Mo_dataset <- data.table(validation_z_con_stabilised, Diagnose=vsd_validation$Diagnose)


reduced_merged_RF <- merged_RF[,intersect(colnames(Mo_dataset), colnames(merged_RF)),with=F]
reduced_Mo_RF <- Mo_dataset[,intersect(colnames(Mo_dataset), colnames(merged_RF)),with=F]


set.seed(seednr)
reduced_splitted <- rsample::initial_split(reduced_merged_RF[,-"rn"], 0.5)
reduced_rftrain <- training(reduced_splitted)
reduced_rftest <- testing(reduced_splitted)

levels(reduced_Mo_RF$Diagnose) <- c("0","1")

reduced_Resultitem <- "Diagnose"

#case.weights calculation
reduced_training_sample_weights <- case_when(reduced_rftrain$Diagnose == "0" ~ 1/sum(reduced_rftrain$Diagnose=="0"),
                                             reduced_rftrain$Diagnose == "1" ~ 1/sum(reduced_rftrain$Diagnose=="1"),
                                             reduced_rftrain$Diagnose == "2" ~ 1/sum(reduced_rftrain$Diagnose=="2"),
                                             reduced_rftrain$Diagnose == "3" ~ 1/sum(reduced_rftrain$Diagnose=="3"))

reduced_rfmodel <- ranger(
  formula         = Diagnose ~ ., 
  data            = reduced_rftrain, 
  num.trees       = 10000,
  seed            = seednr,
  probability     = TRUE,
  importance      = 'impurity',
  classification = TRUE,
  #class.weights = c(1/45,1/114,1/108,1/256),
  case.weights = reduced_training_sample_weights
)
#predict in test
reduced_predictedmodel <- predict(reduced_rfmodel, reduced_rftest, type="response")
reduced_predicted <- as.data.table(reduced_predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()

reduced_optCutOff <- optimalCutoff(reduced_rftest[[reduced_Resultitem]], reduced_predicted, optimiseFor = "Both")[1]
reduced_misClassErrorResult <- misClassError(reduced_rftest[[reduced_Resultitem]], reduced_predicted)#, threshold = optCutOff)
ROCplotResult <- plotROC(reduced_rftest[[reduced_Resultitem]], reduced_predicted, returnSensitivityMat=T)
#repeat for Mo
reduced_predictedmodel <- predict(reduced_rfmodel, reduced_Mo_RF, type="response")
reduced_predicted <- as.data.table(reduced_predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()
reduced_optCutOff <- optimalCutoff(reduced_Mo_RF[[reduced_Resultitem]], reduced_predicted, optimiseFor = "Both")[1]
ROCplotResult <- plotROC(reduced_Mo_RF[[reduced_Resultitem]], reduced_predicted, returnSensitivityMat=T)
Mo_variable_importance <- data.table(feature=names(reduced_rfmodel$variable.importance), importance=reduced_rfmodel$variable.importance)
head(Mo_variable_importance[order(Mo_variable_importance$importance, decreasing=T),],25)

# ggplot(validation_z_con_stabilised,aes(x=colData(vsd_validation)$severity,y=ZFP36L2))+geom_boxplot()+geom_jitter(width=0.05)


####validation2####
validation2_cohortname<- "Planell"
validation2cohorts <- list(Ostrowski ="UCAI Based\nDisease Severity",
                           Mo_UC="Diagnose",
                           Planell="Endoscopic\nMayo Score")
vsd_validation2 <- getCohortsVSD(cohortname=validation2_cohortname)
severityscale <- validation2cohorts[[validation2_cohortname]]
dds_validation2 <- getCohortsDDS(cohortname=validation2_cohortname)
vsd_validation2$severity <- vsd_validation2$Diagnose
vsd_validation2$severity[is.na(vsd_validation2$severity)] <- "Control"  

# ggplot(as.data.table(t(assay(vsd_validation2)[,])),aes(x=colData(vsd_validation2)$severity,y=S100A6))+geom_boxplot()+geom_jitter(width=0.05)

# res_validation2 <- results(dds_validation2, c("Diagnose","Control","UC"),tidy=TRUE)
# res_validation2[res_validation2$row %in% top25$feature,]

validation2_controls <- ifelse(dds_validation2$Diagnose == "Control",TRUE,FALSE)
validation2_vst_counts <- data.table((assay(vsd_validation2)[,]))


# validation2_z_con_stabilised <- data.table(apply(validation2_vst_counts,
#                                                  1,function(x){y <- x[validation2_controls];z <- y[remove_outlier_filter(y)];
#                                                  (x-mean(y))/sd(z)}),
#                                            keep.rownames=TRUE)
validation2_z_con_stabilised <- outlierfiltered_control_transformation(validation2_vst_counts, controlsvector = validation2_controls)

validation2_z_con_stabilised <- data.table(rn=rownames(assay(vsd_validation2)), validation2_z_con_stabilised)
validation2_z_con_stabilised <- transpose_datatable(validation2_z_con_stabilised)

# colnames(validation2_z_con_stabilised) <- c("rn",rownames(assay(vsd_validation2)))

# ggplot(validation2_z_con_stabilised,aes(x=colData(vsd_validation2)$severity,y=S100A12))+geom_boxplot()+geom_jitter(width=0.05)



####RF for validation2####
Planell_genefeatures <- c("rn",rownames(assay(vsd_validation2)))
sum(colnames(merged_RF) %in% Planell_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo

Planell_dataset <- data.table(validation2_z_con_stabilised, Diagnose=vsd_validation2$Diagnose)

reduced2_merged_RF <- merged_RF[,intersect(colnames(Planell_dataset), colnames(merged_RF)),with=F]
reduced2_Planell_RF <- Planell_dataset[,intersect(colnames(Planell_dataset), colnames(merged_RF)),with=F]

set.seed(seednr)
reduced2_splitted <- rsample::initial_split(reduced2_merged_RF[,-"rn"], 0.5)
reduced2_rftrain <- training(reduced2_splitted)
reduced2_rftest <- testing(reduced2_splitted)
# levels(rftrain$Diagnose) <- c("0","1","1")
# levels(rftest$Diagnose) <- c("0","1","1")
levels(reduced2_Planell_RF$Diagnose) <- c("0","1")
#rftest$DiagnoseCohort <- factor(rftest$DiagnoseCohort)
reduced2_Resultitem <- "Diagnose"
#case.weights calculation
reduced2_training_sample_weights <- case_when(reduced2_rftrain$Diagnose == "0" ~ 1/sum(reduced2_rftrain$Diagnose=="0"),
                                              reduced2_rftrain$Diagnose == "1" ~ 1/sum(reduced2_rftrain$Diagnose=="1"),
                                              reduced2_rftrain$Diagnose == "2" ~ 1/sum(reduced2_rftrain$Diagnose=="2"),
                                              reduced2_rftrain$Diagnose == "3" ~ 1/sum(reduced2_rftrain$Diagnose=="3"))

reduced2_rfmodel <- ranger(
  formula         = Diagnose ~ ., 
  data            = reduced2_rftrain, 
  num.trees       = 10000,
  seed            = seednr,
  probability     = TRUE,
  importance      = 'impurity',
  classification = TRUE,
  #class.weights = c(1/45,1/114,1/108,1/256),
  case.weights = reduced2_training_sample_weights
)
#predict in test
reduced2_predictedmodel <- predict(reduced2_rfmodel, reduced2_rftest, type="response")
reduced2_predicted <- as.data.table(reduced2_predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()

reduced2_optCutOff <- optimalCutoff(reduced2_rftest[[reduced2_Resultitem]], reduced2_predicted, optimiseFor = "Both")[1]
reduced2_misClassErrorResult <- misClassError(reduced2_rftest[[reduced2_Resultitem]], reduced2_predicted)#, threshold = optCutOff)
ROCplotResult <- plotROC(reduced2_rftest[[reduced2_Resultitem]], reduced2_predicted, returnSensitivityMat=T)

#repeat for Planell
reduced2_predictedmodel <- predict(reduced2_rfmodel, reduced2_Planell_RF, type="response")
reduced2_predicted <- as.data.table(reduced2_predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()
reduced2_optCutOff <- optimalCutoff(reduced2_Planell_RF[[reduced2_Resultitem]], reduced2_predicted, optimiseFor = "Both")[1]
ROCplotResult <- plotROC(reduced2_Planell_RF[[reduced2_Resultitem]], reduced2_predicted, returnSensitivityMat=T)
Planell_variable_importance <- data.table(feature=names(reduced2_rfmodel$variable.importance), importance=reduced2_rfmodel$variable.importance)
head(Planell_variable_importance[order(Planell_variable_importance$importance, decreasing=T),],25)

# ggplot(validation2_z_con_stabilised,aes(x=colData(vsd_validation2)$severity,y=ZFP36L2))+geom_boxplot()+geom_jitter(width=0.05)
# ggplot(validation2_z_con_stabilised,aes(x=colData(vsd_validation2)$severity,y=HNRNPA0))+geom_boxplot()+geom_jitter(width=0.05)

confusionMatrixResult <- confusionMatrix(reduced2_Planell_RF[[reduced2_Resultitem]], reduced2_predicted, threshold = reduced2_optCutOff)
confusionMatrixResult


####LightGMB####
# set.seed(seednr)
# reduced_splitted <- rsample::initial_split(reduced_merged_RF[,-"rn"], 0.5)
# reduced_rftrain <- training(reduced_splitted)
# reduced_rftest <- testing(reduced_splitted)
# # levels(rftrain$Diagnose) <- c("0","1","1")
# # levels(rftest$Diagnose) <- c("0","1","1")
# levels(reduced_Mo_RF$Diagnose) <- c("0","1")
# #rftest$DiagnoseCohort <- factor(rftest$DiagnoseCohort)
# reduced_Resultitem <- "Diagnose"
# #case.weights calculation
# reduced_training_sample_weights <- case_when(reduced_rftrain$Diagnose == "0" ~ 1/sum(reduced_rftrain$Diagnose=="0"),
#                                              reduced_rftrain$Diagnose == "1" ~ 1/sum(reduced_rftrain$Diagnose=="1"),
#                                              reduced_rftrain$Diagnose == "2" ~ 1/sum(reduced_rftrain$Diagnose=="2"),
#                                              reduced_rftrain$Diagnose == "3" ~ 1/sum(reduced_rftrain$Diagnose=="3"))
# 
# reduced_dtrain <- lgb.Dataset(as.matrix(reduced_rftrain[,-"Diagnose"]), label = as.integer(as.factor(reduced_rftrain$Diagnose))-1)
# reduced_model <- lgb.train(
#   params = list(
#     objective = "binary", 
#     metric = "auc",
#     boosting = "goss",
#     #feature_fraction = 0.5,
#     feature_fraction_bynode = 0.5
#   )
#   , data = reduced_dtrain
# )
# 
# reduced_dtest <- lgb.Dataset(as.matrix(reduced_rftest[,-"Diagnose"]), label = as.integer(as.factor(rftest$Diagnose))-1)
# 
# reduced_predictions <- predict(reduced_model,as.matrix(reduced_rftest[,-"Diagnose"]))
# 
# table(reduced_rftest$Diagnose,round(reduced_predictions)-1)
# plotROC(reduced_rftest[[reduced_Resultitem]], reduced_predictions, returnSensitivityMat=T)
# lgb.importance(reduced_model, percentage = TRUE)
# # cumsum(lgb.importance(reduced_model, percentage = TRUE)$Gain)
# # cumsum(lgb.importance(reduced_model, percentage = TRUE)$Frequency)
# #repeat for Mo
# reduced_Mo_predictions <- predict(reduced_model,as.matrix(reduced_Mo_RF[,-c("rn","Diagnose")]))
# reduced_Mo_predictions <- round(reduced_Mo_predictions)
# #reduced_predictedmodel <- predict(reduced_rfmodel, reduced_Mo_RF, type="response")
# #reduced_Mo_predicted <- as.data.table(reduced_predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()
# reduced_optCutOff <- optimalCutoff(reduced_Mo_RF[[reduced_Resultitem]], reduced_Mo_predictions, optimiseFor = "Both")[1]
# ROCplotResult <- plotROC(reduced_Mo_RF[[reduced_Resultitem]], reduced_Mo_predictions, returnSensitivityMat=T)
# table(reduced_Mo_RF$Diagnose,round(reduced_Mo_predictions))
# 
# 
# #Compare top5 features in boxplots: cohort and validation cohort
# plot_list <- list()
# for(gene in top25$feature){
#   require(ggpubr)
#   try(
#     p1 <- ggplot(
#       data=rbindlist(
#               list(data.table("expression"=unlist(merged_RF[Cohort == 1,..gene]),
#                               Diagnose=merged_metadata[merged_metadata$Diagnose %in% c("Control","UC") & merged_metadata$Cohort == 1]$Diagnose,
#                               Cohort="Our"),
#                     data.table("expression"=unlist(validation_z_con_stabilised[,..gene]),
#                                Diagnose=vsd_validation$Diagnose,
#                                Cohort="Validation")
#                    )
#               ),
#       aes(x=Diagnose,y=expression)
#       )+geom_boxplot()+facet_wrap(~Cohort)+ ylab(gene)+ stat_compare_means( aes(label = ..p.signif..), 
#                                                                 label.x = 1.5)
#   )
#   
#   plot_list[[gene]]<- p1
#   
#   #ggsave(file=p1,filename = paste0("comparison_",gene,".jpg"))
# }
# plot_list