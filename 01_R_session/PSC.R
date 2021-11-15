setwd(paste0("/ssdpool/ewacker/RNAseq","/01_R_session"))

library(data.table)
library(plyr)
library(tidyverse)
library(DESeq2)
library(rvcheck)
library(CEMiTool)
library(rmarkdown)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggridges)
library(ggplotify)
library(extrafont)
library(reshape2)
library(viridis)
library(ggrepel)
library(ggvenn)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)
library(Information)
library(InformationValue)
source("functions.R")
source("RandomForestfunction.R")

remove_outlier_filter <- function(x=NULL){
  calcedSD <- sd(x)
  calcedmean <- mean(x)
  booleanvector <- (x > (calcedmean - 3*calcedSD)) & (x < (calcedmean + 3*calcedSD))
  return(booleanvector)
}

# projectdir <- dirname(getwd())
projectdir <- "/ssdpool/ewacker/RNAseq"
setwd(paste0("/ssdpool/ewacker/RNAseq","/01_R_session"))

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
vsd <- vst(dds, blind=FALSE)

#merged_vst_counts <- t(assay(vsd))
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


#merged_z_con_stabilised <- data.table(apply(merged_vst_counts,1,function(x){(x-mean(x[merged_controls_PSC]))/sd(x[merged_controls_PSC])}),keep.rownames=TRUE)

# PSC_z_con_stabilised <- data.table(apply(merged_vst_counts[,grepl("^I",colnames(merged_vst_counts)),with=FALSE],
#                                             1,function(x){
#                                               (x-mean(x[merged_controls_PSC]))/sd(x[merged_controls_PSC])}),
#                                       keep.rownames=TRUE)

PSC_z_con_stabilised <- data.table(apply(merged_vst_counts[,grepl("^I",colnames(merged_vst_counts)),with=FALSE],
                                         1,function(x){y <- x[merged_controls_PSC];z <- y[remove_outlier_filter(y)];
                                         (x-mean(y))/sd(z)}),
                                   keep.rownames=TRUE)

# UC_z_con_stabilised <- data.table(apply(merged_vst_counts[,grepl("^DE",colnames(merged_vst_counts)),with=FALSE],
#                                          1,function(x){
#                                            (x-mean(x[merged_controls_UC]))/sd(x[merged_controls_UC])}),
#                                    keep.rownames=TRUE)

UC_z_con_stabilised <- data.table(apply(merged_vst_counts[,grepl("^DE",colnames(merged_vst_counts)),with=FALSE],
                                        1,function(x){y <- x[merged_controls_UC];z <- y[remove_outlier_filter(y)];
                                        (x-mean(y))/sd(z)}),
                                  keep.rownames=TRUE)

merged_z_con_stabilised <- rbindlist(list(PSC_z_con_stabilised,UC_z_con_stabilised))
colnames(merged_z_con_stabilised) <- c("rn",rownames(assay(vsd)))



merged_RF <- data.table(merged_z_con_stabilised, Diagnose=merged_metadata$Diagnose, Cohort=c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))
ggplot(merged_RF[Diagnose == "Control",])+geom_jitter(aes(x=Cohort, y=MYL12A))





#merged_z_con_stabilised <- as.data.table(sapply(merged_z_con_stabilised, as.numeric))


# merged_z_con_stabilised[450,merged_controls, with=FALSE]%>%unlist()%>%mean()
# ggplot(merged_z_con_stabilised)+geom_jitter(aes(x=vsd$Cohort,y=unlist(merged_z_con_stabilised[55,]),height=0))




# mat <- as.matrix(merged_z_con_stabilised[,-1])
# #rownames(mat) <- rownames(assay(vsd))
# # mat <- assay(vsd)
# mat <- t(mat)
# mm <- model.matrix(~Diagnose, colData(vsd))
# mat <- sva::ComBat(mat, batch=vsd$Cohort, mod=mm)#, group=as.factor(merged_metadata$Diagnose
# #mat <- limma::removeBatchEffect(mat, batch=vsd$PlateNr, design=mm)#, batch2=c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID)))
# merged_vst_counts_batched <- t(mat)
# merged_vst_counts_batched <- as.data.table(merged_vst_counts_batched)


# ggplot(merged_vst_counts)+geom_boxplot(aes(x=merged_metadata$Cohort,y=S100A12))
# ggplot(merged_z_con_stabilised)+geom_boxplot(aes(x=merged_metadata$Diagnose,y=S100A12))


#scaled_obj <- scale(merged_RF[Cohort==0,-c("Cohort","Diagnose")],center=FALSE)




# merged_RF <- data.table(merged_vst_counts_batched, Diagnose=merged_metadata$Diagnose, Cohort=c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))
# 
# merged_RF_scaled <- data.table(rbindlist(list(data.table(scale(as.matrix(merged_RF[Cohort==0,-c("Cohort","Diagnose")]),center=FALSE)),
#                                     data.table(scale(as.matrix(merged_RF[Cohort==1,-c("Cohort","Diagnose")]),center=FALSE)))),
#                                     Cohort=merged_RF$Cohort,Diagnose=merged_RF$Diagnose)
# 
# 
# ggplot(merged_RF_scaled[Diagnose=="Control",])+geom_jitter(aes(x=Cohort, y=ZFP36L2))
# ggplot(merged_RF_scaled[Diagnose!="Control",])+geom_jitter(aes(x=Cohort, y=ZFP36L2))
# ggplot(merged_RF_scaled)+geom_jitter(aes(x=Cohort, y=ZFP36L2))
# sd(merged_RF_scaled[Diagnose=="Control"&Cohort==0,ZFP36L2])
# sd(merged_RF_scaled[Diagnose=="Control"&Cohort==1,ZFP36L2])
# sd(merged_RF_scaled[Diagnose=="Control"&Cohort==0,ZFP36L2])
# sd(merged_RF_scaled[Diagnose=="Control"&Cohort==1,ZFP36L2])
# sd(merged_RF_scaled[Cohort==0,ZFP36L2])
# sd(merged_RF_scaled[Cohort==1,ZFP36L2])
# merged_RF[Cohort==1,-c("Cohort","Diagnose")]


merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
merged_RF$Cohort  <- factor(merged_RF$Cohort )
levels(merged_RF$Diagnose) <- c("0","1","2")


seednr <- 5495
set.seed(seednr)
splitted <- rsample::initial_split(merged_RF[merged_RF$Diagnose == 0,-c("rn","Diagnose","DiagnoseCohort")], 0.5)
rftrain <- training(splitted)
rftest <- testing(splitted)

# rftrain$Cohort <- droplevels(as.factor(rftrain$Cohort))
# rftest$Cohort <- droplevels(as.factor(rftest$Cohort))


# rftrain$Diagnose <- droplevels(as.factor(rftrain$Diagnose))
# rftest$Diagnose <- droplevels(as.factor(rftest$Diagnose))

####COHORT####

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
# predicted <- predictions(predictedmodel,type="response")
# predicted<-colnames(predicted)[apply(predicted,1,which.max)]
# table(predicted,rftest$Diagnose)
# #compare uc to others
# predicted[predicted=="2"] <- "0"
# #evaluate
optCutOff <- optimalCutoff(rftest[["Cohort"]], predicted, optimiseFor = "Both")[1]

misClassErrorResult <- misClassError(rftest[["Cohort"]], predicted, threshold = optCutOff)

ROCplotResult <- plotROC(rftest[["Cohort"]], predicted, returnSensitivityMat=T)


concordanceResult <- Concordance(rftest[["Cohort"]], predicted)
sensitivityResult <- sensitivity(rftest[["Cohort"]], predicted, threshold = optCutOff)
specificityResult <- specificity(rftest[["Cohort"]], predicted, threshold = optCutOff)
confusionMatrixResult <- confusionMatrix(rftest[["Cohort"]], predicted, threshold = optCutOff)
confusionMatrixResult
variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
head(variable_importance[order(variable_importance$importance, decreasing=T),],25)

mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]

ggplot(merged_RF[merged_RF$Diagnose ==0,])+geom_jitter(aes(x=Cohort,y=PABPC1),height=0)



#####DIAGNOSE####

merged_RF[merged_RF$Diagnose == 0 & merged_RF$Cohort ==0,DiagnoseCohort:=0]
merged_RF[merged_RF$Diagnose == 0 & merged_RF$Cohort ==1,DiagnoseCohort:=1]
merged_RF[merged_RF$Diagnose == 1 & merged_RF$Cohort ==0,DiagnoseCohort:=2]
merged_RF[merged_RF$Diagnose == 2 & merged_RF$Cohort ==1,DiagnoseCohort:=3]
merged_RF$DiagnoseCohort <- factor(merged_RF$DiagnoseCohort)
#levels(rftest$DiagnoseCohort) <- c("0","1","1")
# table(merged_RF$Diagnose, merged_RF$DiagnoseCohort)



seednr <- 22505
set.seed(seednr)
splitted <- rsample::initial_split(merged_RF[merged_RF$Diagnose %in% c("1","2"),-c("Cohort","Diagnose","rn")], 0.5)
rftrain <- training(splitted)
rftest <- testing(splitted)
# levels(rftrain$Diagnose) <- c("0","1","1")
# levels(rftest$Diagnose) <- c("0","1","1")
#rftest$DiagnoseCohort <- factor(rftest$DiagnoseCohort)

#case.weights calculation
training_sample_weights <- case_when(rftrain$DiagnoseCohort == "0" ~ 1/sum(rftrain$DiagnoseCohort=="0"),
                                     rftrain$DiagnoseCohort == "1" ~ 1/sum(rftrain$DiagnoseCohort=="1"),
                                     rftrain$DiagnoseCohort == "2" ~ 1/sum(rftrain$DiagnoseCohort=="2"),
                                     rftrain$DiagnoseCohort == "3" ~ 1/sum(rftrain$DiagnoseCohort=="3"))

rfmodel <- ranger(
  formula         = DiagnoseCohort ~ ., 
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
#predicted <- as.data.table(predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()
predicted <- predictions(predictedmodel,type="response")

# colnames(predicted) <- c("0","1","2","3")
colnames(predicted) <- c("1","2","3")
predicted<-colnames(predicted)[apply(predicted,1,which.max)]
table(predicted,rftest$DiagnoseCohort)
#compare uc to others
# predicted[predicted=="2"] <- "0"
#evaluate
optCutOff <- optimalCutoff(rftest[["DiagnoseCohort"]], predicted, optimiseFor = "Both")[1]

misClassErrorResult <- misClassError(rftest[["DiagnoseCohort"]], predicted)#, threshold = optCutOff)

ROCplotResult <- plotROC(rftest[["DiagnoseCohort"]], predicted, returnSensitivityMat=T)


concordanceResult <- Concordance(rftest[["DiagnoseCohort"]], predicted)
sensitivityResult <- sensitivity(rftest[["DiagnoseCohort"]], predicted, threshold = optCutOff)
specificityResult <- specificity(rftest[["DiagnoseCohort"]], predicted, threshold = optCutOff)
confusionMatrixResult <- confusionMatrix(rftest[["DiagnoseCohort"]], predicted, threshold = optCutOff)
confusionMatrixResult
variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
head(variable_importance[order(variable_importance$importance, decreasing=T),],25)

mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]

ggplot(merged_RF)+geom_jitter(aes(x=DiagnoseCohort,y=HIF1A),height=0)

grep("HIF1A",rownames(assay(vsd)))

ggplot(test_vst_counts)+geom_jitter(aes(x=merged_metadata$Diagnose,y=HIF1A))+facet_wrap(~merged_metadata$Cohort)

