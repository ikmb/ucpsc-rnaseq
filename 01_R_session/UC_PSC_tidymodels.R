setwd(paste0("/ssdpool/ewacker/RNAseq","/01_R_session"))
# setwd(paste0("/home/sukmb465/Documents/Eike/Projects/uc-rnaseq/01_R_session"))

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
library(ComplexHeatmap)
source("functions.R")

seednr <- 555555

projectdir <- dirname(getwd())
# projectdir <- "/ssdpool/ewacker/RNAseq"

# Inputs PSC+UC#####
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


PSC_metadata$Maindiagnose <- case_when(PSC_metadata$`IBD diagnose` == "Ulcerøs colitt" ~ "PSCUC",
                                       PSC_metadata$`IBD diagnose` == "Ingen" & PSC_metadata$Diagnose=="PSC" ~ "PSC",
                                       PSC_metadata$`IBD diagnose` == "Ingen" & PSC_metadata$Diagnose=="Control" ~ "Control",
                                       PSC_metadata$`IBD diagnose` %in% c('Crohn',"Indeterminate","Uavklart","Ukjent") ~ "REMOVE",
)
PSC_metadata <- PSC_metadata[!PSC_metadata$Maindiagnose == "REMOVE",]
PSC_metadata$Diagnose <- PSC_metadata$Maindiagnose

PSC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="PSC")
PSC_rawcounts <- PSC_rawcounts[,colnames(PSC_rawcounts) %in% PSC_metadata$SampleID]
UC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="Our")
UC_metadata <- getmetadata(corhorttablelocation=corhorttablelocation,cohortname="Our")


SampleIDs <- c(colnames(PSC_rawcounts),colnames(UC_rawcounts))
merged_rawcounts <- rbindlist(list(data.frame(t(PSC_rawcounts)),data.frame(t(UC_rawcounts))),fill=TRUE)
intersectedgenes <- intersect(colnames(data.frame(t(PSC_rawcounts))), colnames(data.frame(t(UC_rawcounts))))
merged_rawcounts <- merged_rawcounts[,..intersectedgenes]

#Hard filtering to exclude zero expression in any sample genes:
# merged_rawcounts <- merged_rawcounts[,!(apply(merged_rawcounts <= 1, 2, any)),with=FALSE]

#Filtering for allowing only 10% of samples with zero values per gene
merged_rawcounts <- merged_rawcounts[,as.vector(colSums(merged_rawcounts == 0)<(dim(merged_rawcounts)[1]/10)),with=FALSE]

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
# Transformation PSC+UC####

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

#Convert Diagnose into a factor: 0=Control, 1=PSC, 2=PSCUC, 3=UC
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
levels(merged_RF$Diagnose) <- c("0","1","2", "3")

merged_RF$Cohort  <- factor(merged_RF$Cohort )


PSC_UC_COHORT_results <- performML(dataset = merged_RF[merged_RF$Diagnose == 0,-c("rn","Diagnose")], splitfactor = 0.5, outcome_column = "Cohort",seed = seednr)

pROC::auc(PSC_UC_COHORT_results$pROC_object)
#ggplot(merged_RF[merged_RF$Diagnose ==0,])+geom_jitter(aes(x=Cohort,y=RPS28),height=0)

# UC vs CON RF Model####
#Convert Diagnose into a factor: 0=Control, 1=PSC, 2=PSCUC, 3=UC
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
levels(merged_RF$Diagnose) <- c("0","1","2", "3")

merged_RF$Cohort  <- factor(merged_RF$Cohort )

filtered_merged_RF <- merged_RF[merged_RF$Cohort == 1,]
levels(filtered_merged_RF$Diagnose) <- droplevels(factor(filtered_merged_RF$Diagnose))
#0 = Control, 1= UC
levels(filtered_merged_RF$Diagnose) <- c("0", "1")#0 = Control, 1= UC

# filtered_merged_RF[filtered_merged_RF$Diagnose == 0 & filtered_merged_RF$Cohort ==0,DiagnoseCohort:=0]
# filtered_merged_RF[filtered_merged_RF$Diagnose == 0 & filtered_merged_RF$Cohort ==1,DiagnoseCohort:=1]
# filtered_merged_RF[filtered_merged_RF$Diagnose == 1 & filtered_merged_RF$Cohort ==0,DiagnoseCohort:=2]
# filtered_merged_RF[filtered_merged_RF$Diagnose == 2 & filtered_merged_RF$Cohort ==1,DiagnoseCohort:=3]
# filtered_merged_RF$DiagnoseCohort <- factor(filtered_merged_RF$DiagnoseCohort)
#levels(rftest$DiagnoseCohort) <- c("0","1","1")
# table(filtered_merged_RF$Diagnose, filtered_merged_RF$DiagnoseCohort)


CON_UC_Diagnose_results <- performML(dataset = filtered_merged_RF[,-c("Cohort","rn")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(CON_UC_Diagnose_results$pROC_object)


# ggplot(test_vst_counts)+geom_jitter(aes(x=merged_metadata$Diagnose,y=ZFP36L2))+facet_wrap(~merged_metadata$Cohort)


# Mo Validation Preparation####
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
# res_validation[res_validation$row %in% top25$feature,]

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


# Mo RF for validation Model####
Mo_genefeatures <- c("rn",rownames(assay(vsd_validation)))
sum(colnames(filtered_merged_RF) %in% Mo_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo
Mo_dataset <- data.table(validation_z_con_stabilised, Diagnose=vsd_validation$Diagnose)


reduced_merged_RF <- filtered_merged_RF[,intersect(colnames(Mo_dataset), colnames(filtered_merged_RF)),with=F]
reduced_Mo_RF <- Mo_dataset[,intersect(colnames(Mo_dataset), colnames(filtered_merged_RF)),with=F]


Mo_CON_UC_Diagnose_results <- performML(dataset = reduced_merged_RF[,-c("rn")],testing_dataset = reduced_Mo_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(CON_UC_Diagnose_results$pROC_object)


# ggplot(validation_z_con_stabilised,aes(x=colData(vsd_validation)$severity,y=ZFP36L2))+geom_boxplot()+geom_jitter(width=0.05)


# Planell validation Preparation####
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


validation2_z_con_stabilised <- data.table(apply(validation2_vst_counts,
                                                 1,function(x){y <- x[validation2_controls];z <- y[remove_outlier_filter(y)];
                                                 (x-mean(y))/sd(z)}),
                                           keep.rownames=TRUE)
# validation2_z_con_stabilised <- outlierfiltered_control_transformation(validation2_vst_counts, controlsvector = validation2_controls)
# 
# validation2_z_con_stabilised <- data.table(rn=rownames(assay(vsd_validation2)), validation2_z_con_stabilised)
# validation2_z_con_stabilised <- transpose_datatable(validation2_z_con_stabilised)

colnames(validation2_z_con_stabilised) <- c("rn",rownames(assay(vsd_validation2)))

# ggplot(validation2_z_con_stabilised,aes(x=colData(vsd_validation2)$severity,y=S100A12))+geom_boxplot()+geom_jitter(width=0.05)



# Planell RF for validation2 Model####
Planell_genefeatures <- c("rn",rownames(assay(vsd_validation2)))
sum(colnames(filtered_merged_RF) %in% Planell_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo

Planell_dataset <- data.table(validation2_z_con_stabilised, Diagnose=vsd_validation2$Diagnose)

reduced2_merged_RF <- filtered_merged_RF[,intersect(colnames(Planell_dataset), colnames(filtered_merged_RF)),with=F]
reduced2_Planell_RF <- Planell_dataset[,intersect(colnames(Planell_dataset), colnames(filtered_merged_RF)),with=F]

levels(reduced2_Planell_RF[["Diagnose"]]) <- c(0,1)

Planell_CON_UC_Diagnose_results <- performML(dataset = reduced2_merged_RF[,-c("rn")],testing_dataset = reduced2_Planell_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)
reducedforplanell_CON_UC_Diagnose_results <- predictwithatunedmodel(tuned_model = Planell_CON_UC_Diagnose_results$tuned_model, testing_dataset = Planell_CON_UC_Diagnose_results$test_df, outcome_column = "Diagnose",seed = seednr)
#Auroc for external validation dataset
pROC::auc(Planell_CON_UC_Diagnose_results$pROC_object)
# pROC::ci.auc(Planell_CON_UC_Diagnose_results$pROC_object)
#validation with own test dataset
pROC::auc(reducedforplanell_CON_UC_Diagnose_results$pROC_object)



# Ostrowski validation Preparation####
validation3_cohortname<- "Ostrowski"
validation3cohorts <- list(Ostrowski ="UCAI Based\nDisease Severity",
                           Mo_UC="Diagnose",
                           Planell="Endoscopic\nMayo Score")
vsd_validation3 <- getCohortsVSD(cohortname=validation3_cohortname)
severityscale <- validation3cohorts[[validation3_cohortname]]
dds_validation3 <- getCohortsDDS(cohortname=validation3_cohortname)
vsd_validation3$severity <- vsd_validation3$Diagnose
vsd_validation3$severity[is.na(vsd_validation3$severity)] <- "Control"  

# ggplot(as.data.table(t(assay(vsd_validation3)[,])),aes(x=colData(vsd_validation3)$severity,y=S100A6))+geom_boxplot()+geom_jitter(width=0.05)

# res_validation3 <- results(dds_validation3, c("Diagnose","Control","UC"),tidy=TRUE)
# res_validation3[res_validation3$row %in% top25$feature,]

validation3_controls <- ifelse(dds_validation3$Diagnose == "Control",TRUE,FALSE)
validation3_vst_counts <- data.table((assay(vsd_validation3)[,]))


validation3_z_con_stabilised <- data.table(apply(validation3_vst_counts,
                                                 1,function(x){y <- x[validation3_controls];z <- y[remove_outlier_filter(y)];
                                                 (x-mean(y))/sd(z)}),
                                           keep.rownames=TRUE)
# validation3_z_con_stabilised <- outlierfiltered_control_transformation(validation3_vst_counts, controlsvector = validation3_controls)
# 
# validation3_z_con_stabilised <- data.table(rn=rownames(assay(vsd_validation3)), validation3_z_con_stabilised)
# validation3_z_con_stabilised <- transpose_datatable(validation3_z_con_stabilised)

colnames(validation3_z_con_stabilised) <- c("rn",rownames(assay(vsd_validation3)))

# ggplot(validation3_z_con_stabilised,aes(x=colData(vsd_validation3)$severity,y=S100A12))+geom_boxplot()+geom_jitter(width=0.05)



# Ostrowski RF for validation3 Model####
Ostrowski_genefeatures <- c("rn",rownames(assay(vsd_validation3)))
sum(colnames(filtered_merged_RF) %in% Ostrowski_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo

Ostrowski_dataset <- data.table(validation3_z_con_stabilised, Diagnose=vsd_validation3$Diagnose)

reduced3_merged_RF <- filtered_merged_RF[,intersect(colnames(Ostrowski_dataset), colnames(filtered_merged_RF)),with=F]
reduced3_Ostrowski_RF <- Ostrowski_dataset[,intersect(colnames(Ostrowski_dataset), colnames(filtered_merged_RF)),with=F]

levels(reduced3_Ostrowski_RF[["Diagnose"]]) <- c(0L,1L)

Ostrowski_CON_UC_Diagnose_results <- performML(dataset = reduced3_merged_RF[,-c("rn")],testing_dataset = reduced3_Ostrowski_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)
reducedforOstrowski_CON_UC_Diagnose_results <- predictwithatunedmodel(tuned_model = Ostrowski_CON_UC_Diagnose_results$tuned_model, testing_dataset = Ostrowski_CON_UC_Diagnose_results$test_df, outcome_column = "Diagnose",seed = seednr)
#Auroc for external validation dataset
pROC::auc(Ostrowski_CON_UC_Diagnose_results$pROC_object)
# pROC::ci.auc(Ostrowski_CON_UC_Diagnose_results$pROC_object)
#validation with own test dataset
pROC::auc(reducedforOstrowski_CON_UC_Diagnose_results$pROC_object)
# LightGMB (NOT USED)####
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


# PSC vs UC model, no controls, no PSCUC ####

PSC_UC_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 2,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSC_UC_Diagnose_results$pROC_object)

ggplot(merged_RF[,-c("rn","Cohort")])+geom_jitter(aes(x=Diagnose,y=ACLY))+geom_boxplot(aes(x=Diagnose,y=ACLY))#+facet_wrap(~merged_metadata$Cohort)


# PSC vs Controls model, no UC, no PSCUC ####

PSC_CONTROL_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 2 & merged_RF$Diagnose != 3,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSC_CONTROL_Diagnose_results$pROC_object)

#0=Control, 1=PSC, 2=PSCUC, 3=UC
ggplot(merged_RF[,-c("rn","Cohort")],aes(x=Diagnose,y=BACH2))+geom_jitter()+geom_boxplot()+scale_x_discrete(labels = c('Control',"PSC","PSC/UC","UC"))#+facet_wrap(~merged_metadata$Cohort)


# PSC vs PSCUC, no CON, no UC ###########

PSC_PSCUC_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 3,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSC_PSCUC_Diagnose_results$pROC_object)


# ggplot(merged_RF[,-c("rn","Cohort")],aes(x=Diagnose,y=UBASH3A))+
#   geom_jitter()+#geom_boxplot()+
#   scale_x_discrete(labels = c('Control',"PSC","PSC/UC","UC"))#+facet_wrap(~merged_metadata$Cohort)

# PSCUC vs UC, no CON, no PSC ###########
PSCUC_UC_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 1,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSCUC_UC_Diagnose_results$pROC_object)

# ggplot(merged_RF[,-c("rn","Cohort")],aes(x=Diagnose,y=CD14))+
#   geom_jitter()+
#   geom_boxplot()+
#   scale_x_discrete(labels = c('Control',"PSC","PSC/UC","UC"))#+facet_wrap(~merged_metadata$Cohort)

# PSCUC+PSC vs UC; no CON###########
#Diagnose: 0=Control, 1=PSC, 2=PSCUC, 3=UC
PSCUCPSC_merged_RF <- merged_RF[merged_RF$Diagnose != 0,-c("rn","Cohort")]
PSCUCPSC_merged_RF$Diagnose <- droplevels(PSCUCPSC_merged_RF$Diagnose)
#merge PSC and PSCUC to one level; PSC+PSCUC = 0; UC = 1:
levels(PSCUCPSC_merged_RF$Diagnose) <- c("0","0","1")
#PSCUCPSC_merged_RF$Diagnose <- droplevels(PSCUCPSC_merged_RF$Diagnose)

PSCUCPSC_UC_Diagnose_results <- performML(dataset = PSCUCPSC_merged_RF,
                                          splitfactor = 0.5, 
                                          outcome_column = "Diagnose",
                                          seed = seednr)

pROC::auc(PSCUCPSC_UC_Diagnose_results$pROC_object)

# ggplot(merged_RF[,-c("rn","Cohort")],aes(x=Diagnose,y=CD14))+
#   geom_jitter()+
#   geom_boxplot()+
#   scale_x_discrete(labels = c('Control',"PSC","PSC/UC","UC"))#+facet_wrap(~merged_metadata$Cohort)

# PSCUC+PSC vs CON; no UC###########
#Diagnose: 0=Control, 1=PSC, 2=PSCUC, 3=UC
PSCUCPSC_CON_merged_RF <- merged_RF[merged_RF$Diagnose != 3,-c("rn","Cohort")]
PSCUCPSC_CON_merged_RF$Diagnose <- droplevels(PSCUCPSC_CON_merged_RF$Diagnose)
#merge PSC and PSCUC to one level; PSC+PSCUC = 1; Con = 0:
levels(PSCUCPSC_CON_merged_RF$Diagnose) <- c("0","1","1")
#PSCUCPSC_CON_merged_RF$Diagnose <- droplevels(PSCUCPSC_merged_RF$Diagnose)

PSCUCPSC_CON_Diagnose_results <- performML(dataset = PSCUCPSC_CON_merged_RF,
                                          splitfactor = 0.5, 
                                          outcome_column = "Diagnose",
                                          seed = seednr)

pROC::auc(PSCUCPSC_CON_Diagnose_results$pROC_object)

# ggplot(merged_RF[,-c("rn","Cohort")],aes(x=Diagnose,y=CD14))+
#   geom_jitter()+
#   geom_boxplot()+
#   scale_x_discrete(labels = c('Control',"PSC","PSC/UC","UC"))#+facet_wrap(~merged_metadata$Cohort)
# CEMiTool #####
#create an expressionmatrix dataframe for CEMItool:
merged_RF_expressionmatrix <- transpose_datatable(merged_RF[,-c("Diagnose", "Cohort")])
merged_RF_expressionmatrix <- data.frame(merged_RF_expressionmatrix,row.names = merged_RF_expressionmatrix$rn)[,-1]

#perform CEMItool analysis
cem_merged_RF <- CEMiwrapper(expressionmatrix=merged_RF_expressionmatrix, ID=merged_metadata$SampleID, Groups=merged_metadata$Diagnose, reportname=paste0("PSC_PSCUC_UC","_Diagnose"), applyfiltering = FALSE)
genesinput <- data.table(cem_merged_RF@module)[modules == "M2",]$genes
# save.image(file="performedanalysis.RData")
# Cell Blueprints ####
#Cell Blueprints for every module created from cemitools:
enrichmentslist <- list()
for (x in unique(cem_merged_RF@module$modules)){
  print(x)
  genesinput <- data.table(cem_merged_RF@module)[modules == x,]$genes
  result <- blueprintenrichments(genes=genesinput)
  # result <- result[result$enrichment >= 0.1,]
  enrichmentslist[[x]] <- result
}

# lapply(enrichmentslist, function(x) print(x[order(x$enrichment, decreasing = T),]))

enrichment_df <- data.table(bind_rows(enrichmentslist, .id = "column_label")) %>% pivot_wider(names_from = celltype, values_from = enrichment) %>% as.data.table() %>% as.matrix(., rownames="column_label")

#M3C::M3C(dcast(melt(data.table(rn=rownames(enrichment_df),enrichment_df),id.vars="rn"), variable ~ rn)[,-"variable"], clusteralg = "km",cores = 4)



htmp <- Heatmap(enrichment_df, 
                name="Z-score",
                heatmap_legend_param = list(
                  legend_direction = "horizontal", 
                  legend_width = unit(6, "cm")),
                height = unit(10, "cm"),
                row_km = 3,
                row_km_repeats = 10,
                row_names_side = "left",
                width = unit(14, "cm"),
                column_names_rot = 45, 
                # column_names_max_height = unit(2, "cm"),
                column_names_side = "bottom",
                column_names_gp = gpar(fontsize = 10)
)
draw(htmp, heatmap_legend_side="top")

# topGO function results #####
#use topgo function to evaluate the biological function of features, chosen by RF:
res <- results(dds,  contrast=c("Diagnose","UC","Control"),tidy=T)

topGo_object <- topGO_enrichment(genelist = PSCUC_UC_Diagnose_results$variable_importance[1:1000,]$feature, statistical_test = "ks.ties", weights = PSCUC_UC_Diagnose_results$variable_importance[1:1000,]$importance, #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                              DESeq_result = res, match_by_expression = TRUE)
topGo_object <- topGO_enrichment(genelist = PSCUC_UC_Diagnose_results$variable_importance[1:100,]$feature, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                              DESeq_result = res, match_by_expression = TRUE, ontology_type = "MF")




#
# Analysis Results ####

# DEseq2 results:
#PSC vs CON
res_PSC <- results(dds,  contrast=c("Diagnose","PSC","Control"))
resLFC_PSC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

#UC vs CON
res_UC <- results(dds,  contrast=c("Diagnose","UC","Control"))
resLFC_UC <- lfcShrink(dds, coef="Diagnose_UC_vs_Control", type="apeglm")

#UC vs PSC
res_UC_PSC <- results(dds,  contrast=c("Diagnose","UC","PSC"))
resLFC_UC_PSC <- lfcShrink(dds, coef="Diagnose_UC_vs_PSC", type="apeglm")

#PSCUC vs PSC
res_PSCUC_PSC <- results(dds,  contrast=c("Diagnose","PSCUC","PSC"))
resLFC_PSCUC_PSC <- lfcShrink(dds, coef="Diagnose_PSCUC_vs_PSC", type="apeglm")

# Random Forest results:
#Validations:
Planell_CON_UC_Diagnose_results
Mo_CON_UC_Diagnose_results
Ostrowski_CON_UC_Diagnose_results

#Our Dataset:
CON_UC_Diagnose_results

PSC_UC_COHORT_results
PSC_UC_Diagnose_results
PSCUC_UC_Diagnose_results
PSC_PSCUC_Diagnose_results
PSC_CONTROL_Diagnose_results
PSCUCPSC_UC_Diagnose_results
PSCUCPSC_CON_Diagnose_results

#Cemitool
cem_merged_RF

#Blueprints
draw(htmp, heatmap_legend_side="top")

#TopGO
topGo_object