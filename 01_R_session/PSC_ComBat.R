library(data.table)
library(plyr)
library(tidyverse)
library(DESeq2)
# library(pheatmap)
library(rvcheck)
# library(CEMiTool)
# library(M3C)
library(rmarkdown)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggridges)
library(ggplotify)
# library(VennDiagram)
library(extrafont)
# library(SCORPIUS)
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


projectdir <- dirname(getwd())


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
# UC_metadata$Diagnose


SampleIDs <- c(colnames(PSC_rawcounts),colnames(UC_rawcounts))
merged_rawcounts <- rbindlist(list(data.frame(t(PSC_rawcounts)),data.frame(t(UC_rawcounts))),fill=TRUE)
intersectedgenes <- intersect(colnames(data.frame(t(PSC_rawcounts))), colnames(data.frame(t(UC_rawcounts))))
merged_rawcounts <- merged_rawcounts[,..intersectedgenes]

#Hard filtering to exclude zero expression in any sample genes:
merged_rawcounts <- merged_rawcounts[,!(apply(merged_rawcounts == 0, 2, any)),with=FALSE]


merged_rawcounts <- t(merged_rawcounts)
# rownames(merged_rawcounts)
colnames(merged_rawcounts) <- SampleIDs


merged_metadata <- rbind(PSC_metadata[,c("SampleID","Diagnose","PlateNr")],UC_metadata[,c("SampleID","Diagnose","PlateNr")])

mat <- merged_rawcounts
#mm <- model.matrix(~Diagnose, colData(vsd))
#mat <- limma::removeBatchEffect(mat, batch=vsd$PlateNr, design=mm)
mat <- sva::ComBat_seq(mat, batch=merged_metadata$PlateNr, group=as.factor(merged_metadata$Diagnose))
# merged_rawcounts_debatched <- t(mat)

dds <- DESeqDataSetFromMatrix(countData=mat, colData=merged_metadata, design= ~Diagnose)
dds <- DESeq(dds)
#res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
vsd <- vst(dds, blind=FALSE)
# vsd <- rlog(dds, blind=FALSE)

# mat <- assay(vsd)
# mm <- model.matrix(~Diagnose, colData(vsd))
# #mat <- limma::removeBatchEffect(mat, batch=vsd$PlateNr, design=mm)
# mat <- sva::ComBat_seq(mat, batch=vsd$PlateNr, group=as.factor(merged_metadata$Diagnose))
# merged_vst_counts <- t(mat)

merged_vst_counts <- t(assay(vsd))
merged_vst_counts <- as.data.table(merged_vst_counts)
merged_vst_counts <- as.data.table(sapply(merged_vst_counts, as.numeric))

merged_controls <- ifelse(merged_metadata$Diagnose == "Control",TRUE,FALSE)
merged_z_con_stabilised <- data.table(apply(merged_vst_counts[,],2,function(x){x-median(x[merged_controls])}),keep.rownames=TRUE)


ggplot(merged_vst_counts)+geom_boxplot(aes(x=merged_metadata$Diagnose,y=S100A12))
ggplot(merged_z_con_stabilised)+geom_boxplot(aes(x=merged_metadata$Diagnose,y=S100A12))



merged_RF <- data.table(merged_z_con_stabilised, Diagnose=merged_metadata$Diagnose, Cohort=c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
merged_RF$Cohort  <- factor(merged_RF$Cohort )
levels(merged_RF$Diagnose) <- c("0","1","2")




seednr <- 22505
set.seed(seednr)
splitted <- rsample::initial_split(merged_RF[merged_RF$Diagnose == 0,-"Diagnose"], 0.5)
rftrain <- training(splitted)
rftest <- testing(splitted)

# rftest$Cohort

# rftrain$Diagnose <- droplevels(as.factor(rftrain$Diagnose))
# rftest$Diagnose <- droplevels(as.factor(rftest$Diagnose))

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
#compare uc to others
# predicted[predicted=="2"] <- "0"
#evaluate
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


ggplot(merged_RF[merged_RF$Diagnose == 0,])+geom_jitter(aes(x=Cohort,y=CHIC2),height=0)


































# PSC_vsd <- getnormalizedDESeq2object(corhorttablelocation=corhorttablelocation,cohortname="PSC")
# PSC_vst_counts <- t(assay(PSC_vsd))
# PSC_vst_counts <- as.data.table(PSC_vst_counts, keep.rownames = TRUE)
# PSC_vst_counts <- PSC_vst_counts[PSC_vst_counts$rn %in% PSC_metadata$SampleID,]
# ggplot(PSC_vst_counts)+geom_boxplot(aes(x=PSC_metadata$Diagnose,y=S100A12))
# 
# PSCcontrols <- ifelse(PSC_metadata$Diagnose == "Control",TRUE,FALSE)
# PSC_z_con_stabilised <- data.table(apply(PSC_vst_counts[,-"rn"],2,function(x){x-median(x[PSCcontrols])}),keep.rownames=TRUE)
# ggplot(PSC_z_con_stabilised)+geom_boxplot(aes(x=PSC_metadata$Diagnose,y=S100A12))
# 
# # merged_RF <- data.table(PSC_z_con_stabilised, Diagnose=PSC_metadata$Diagnose)
# 
# 
# 
# #####Inputs#####UC
# corhorttablelocation <- paste0(projectdir,"/00_RawData/cohorts.csv")
# cohorttable <- fread(corhorttablelocation, header=T)
# cohortname <- "Our"
# Disease <- unlist(cohorttable[Cohort %in% cohortname, 5])
# 
# gmt_file <- paste0(projectdir,"/data_rescources/GO_Biological_Process_2018.gmt")
# output_location <- paste0("output")
# 
# UC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="Our")
# UC_metadata <- getmetadata(corhorttablelocation=corhorttablelocation,cohortname="Our")
# UC_vsd <- getnormalizedDESeq2object(corhorttablelocation=corhorttablelocation,cohortname="Our")
# UC_vst_counts <- t(assay(UC_vsd))
# UC_vst_counts <- as.data.table(UC_vst_counts, keep.rownames = TRUE)
# UC_vst_counts <- UC_vst_counts[UC_vst_counts$rn %in% UC_metadata$SampleID,]
# UCcontrols <- ifelse(UC_metadata$Diagnose == "Control",TRUE,FALSE)
# #UC_z_con_stabilised <- data.table(apply(assay(UC_vsd),2,function(x){x-median(x)}),keep.rownames=TRUE)
# UC_z_con_stabilised <- data.table(apply(UC_vst_counts[,-"rn"],2,function(x){x-median(x[UCcontrols])}),keep.rownames=TRUE)
# 
# ggplot(UC_vst_counts)+geom_boxplot(aes(x=UC_metadata$Diagnose,y=S100A12))
# ggplot(UC_z_con_stabilised)+geom_boxplot(aes(x=UC_metadata$Diagnose,y=S100A12))
# 
# # merged_RF <- data.table(UC_z_con_stabilised, Diagnose=UC_metadata$Diagnose)
# 
# 
# 
# merged_RF <- rbindlist(list(PSC_z_con_stabilised, UC_z_con_stabilised),fill=TRUE, idcol = TRUE)
# 
# intersectedgenes <- intersect(colnames(UC_z_con_stabilised), colnames(PSC_z_con_stabilised))
# 
# merged_RF <- merged_RF[,c(intersectedgenes,".id"),with=FALSE]

merged_RF$Diagnose <- c(PSC_metadata$Diagnose,UC_metadata$Diagnose)
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
levels(merged_RF$Diagnose) <- c("0","1","2")

merged_RF$.id <- factor(merged_RF$.id)

# UCcontrols <- !as.logical(as.integer(UC_vsd$Diagnose)-1)
# UC_z_con_stabilised <- data.table(apply(assay(UC_vsd),2,function(x){x-median(x)}),keep.rownames=TRUE)
# PSCcontrols <- ifelse(PSC_vsd[,PSC_vsd$Diagnose %in% c("PSC","Control")]$Diagnose == "Control",TRUE,FALSE)
# PSC_z_con_stabilised <- data.table(apply(assay(PSC_vsd)[,PSC_vsd$Diagnose %in% c("PSC","Control")],2,function(x){x-median(x)}),keep.rownames=TRUE)




#View(data.table(colnames(assay(PSC_vsd)[,PSC_vsd$Diagnose %in% c("PSC","Control")]),ifelse(PSC_vsd[,PSC_vsd$Diagnose %in% c("PSC","Control")]$Diagnose == "Control",TRUE,FALSE),PSC_vsd[,PSC_vsd$Diagnose %in% c("PSC","Control")]$Diagnose))

# merged_z_con_stabilised <- merge(UC_z_con_stabilised,PSC_z_con_stabilised,by="rn")
# merged_metadata <- data.table(rbind(PSC_metadata[PSC_metadata$Diagnose %in% c("PSC","Control"),c("SampleID","Diagnose")],UC_metadata[,c("SampleID","Diagnose")]))
# 
# transposed_merged_z_con_stabilised <- t(merged_z_con_stabilised)
# transposed_merged_z_con_stabilised <- as.data.table(transposed_merged_z_con_stabilised)
# colnames(transposed_merged_z_con_stabilised) <- unlist(transposed_merged_z_con_stabilised[1,])
# transposed_merged_z_con_stabilised <- transposed_merged_z_con_stabilised[-1,]
# 
# merged_RF <- data.table(as.data.table(sapply(transposed_merged_z_con_stabilised, as.numeric)), Diagnose=merged_metadata$Diagnose)
# merged_RF$Diagnose <- factor(merged_RF$Diagnose)
# levels(merged_RF$Diagnose) <- c("0", "1")
levels(merged_RF$.id) <- c("0","1")


seednr <- 223
set.seed(seednr)
splitted <- rsample::initial_split(merged_RF[merged_RF$Diagnose == 0,-"Diagnose"], 0.5)
rftrain <- training(splitted)
rftest <- testing(splitted)

rftrain$Diagnose <- droplevels(as.factor(rftrain$Diagnose))
rftest$Diagnose <- droplevels(as.factor(rftest$Diagnose))

rfmodel <- ranger(
  formula         = .id ~ ., 
  data            = rftrain[,c("SGIP1",".id")], 
  num.trees       = 100,
  seed            = seednr,
  probability     = TRUE,
  importance      = 'impurity'
)

ggplot(data=rftrain[,c("SGIP1",".id")])+geom_jitter(aes(x=.id,y=SGIP1))
ggplot(data=rftest[,c("SGIP1",".id")])+geom_jitter(aes(x=.id,y=SGIP1))
predictedmodel <- predict(rfmodel, rftest, type="response")
#predicted <- as.data.table(predictedmodel$predictions)[,2]%>%unlist()%>%as.vector()
predicted <- predictions(predictedmodel,type="response")
# predicted<-colnames(predicted)[apply(predicted,1,which.max)]
# table(predicted,rftest$Diagnose)
#compare uc to others
# predicted[predicted=="2"] <- "0"
#evaluate
optCutOff <- optimalCutoff(rftest[[".id"]], predicted, optimiseFor = "Both")[1]

misClassErrorResult <- misClassError(rftest[[".id"]], predicted, threshold = optCutOff)

ROCplotResult <- plotROC(rftest[[".id"]], predicted, returnSensitivityMat=T)


concordanceResult <- Concordance(rftest[["Diagnose"]], predicted)
sensitivityResult <- sensitivity(rftest[["Diagnose"]], predicted, threshold = optCutOff)
specificityResult <- specificity(rftest[["Diagnose"]], predicted, threshold = optCutOff)
confusionMatrixResult <- confusionMatrix(rftest[["Diagnose"]], predicted, threshold = optCutOff)

variable_importance <- data.table(feature=names(rfmodel$variable.importance), importance=rfmodel$variable.importance)
head(variable_importance[order(variable_importance$importance, decreasing=T),],25)

mostimpactfeature <- variable_importance[variable_importance$importance==max(variable_importance$importance),feature]
