## Author: Florian Uellendahl-Werth & Eike Matthias Wacker
##
## IKMB Kiel, 2023
## Email: e.wacker at ikmb.uni-kiel.de

library(data.table)
library(plyr)
library(tidyverse)
library(tidymodels)
library(DESeq2)
library(RColorBrewer)
library(rsample)      # data splitting 
library(ranger)
library(Information)
library(InformationValue)
library(ComplexHeatmap)
library(ggpubr)
library(CEMiTool)
library(rmarkdown)
library(EnhancedVolcano)

source("functions.r")

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

PSC_metadata <- getmetadata(corhorttablelocation=corhorttablelocation,cohortname="PSC")
PSC_metadata <- PSC_metadata[PSC_metadata$Diagnose %in% c("PSC","Control"),]

PSC_metadata$yearssincediagnose <- PSC_metadata$`Age at samplingdate` - PSC_metadata$`Years at liverdiagnose`
as.Date(PSC_metadata$`Date liverdiagnose`,tryFormats = c("%m/%d/%Y")) - as.Date(PSC_metadata$`Samplingdate`,tryFormats = c("%m/%d/%Y"))# - PSC_metadata$`Years at liverdiagnose`

PSC_metadata$Maindiagnose <- case_when(PSC_metadata$`IBD diagnose` == "Ulcerøs colitt" ~ "PSCUC",
                                       PSC_metadata$`IBD diagnose` == "Ingen" & PSC_metadata$Diagnose=="PSC" ~ "PSC",
                                       PSC_metadata$`IBD diagnose` == "Ingen" & PSC_metadata$Diagnose=="Control" ~ "Control",
                                       PSC_metadata$`IBD diagnose` %in% c('Crohn',"Indeterminate","Uavklart","Ukjent") ~ "REMOVE",
)
PSC_metadata <- PSC_metadata[!PSC_metadata$Maindiagnose == "REMOVE",]
PSC_metadata$Diagnose <- PSC_metadata$Maindiagnose


PSC_metadata$yearssincediagnose <- PSC_metadata$`Age at samplingdate` - PSC_metadata$`Years at liverdiagnose`

PSC_metadata$yearssincediagnose[PSC_metadata$yearssincediagnose < 0] <- NA_real_



PSC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="PSC")
PSC_rawcounts <- PSC_rawcounts[,colnames(PSC_rawcounts) %in% PSC_metadata$SampleID]
UC_rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname="Our")
UC_metadata <- getmetadata(corhorttablelocation=corhorttablelocation,cohortname="Our")
UC_metadata$yearssincediagnose <- NA_real_

SampleIDs <- c(colnames(PSC_rawcounts),colnames(UC_rawcounts))
merged_rawcounts <- rbindlist(list(data.frame(t(PSC_rawcounts)),data.frame(t(UC_rawcounts))),fill=TRUE)
intersectedgenes <- intersect(colnames(data.frame(t(PSC_rawcounts))), colnames(data.frame(t(UC_rawcounts))))
merged_rawcounts <- merged_rawcounts[,..intersectedgenes]

#Filtering for allowing only 10% of samples with zero values per gene
merged_rawcounts <- merged_rawcounts[,as.vector(colSums(merged_rawcounts == 0)<(dim(merged_rawcounts)[1]/10)),with=FALSE]

merged_rawcounts <- t(merged_rawcounts)
colnames(merged_rawcounts) <- SampleIDs

merged_metadata <- rbind(PSC_metadata[,c("SampleID","Diagnose","PlateNr", "yearssincediagnose")],UC_metadata[,c("SampleID","Diagnose","PlateNr", "yearssincediagnose")])
merged_metadata$Cohort <- factor(c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))

merged_metadata$yearssincebinned <- case_when(merged_metadata$yearssincediagnose == 0 ~ "newly",
          merged_metadata$yearssincediagnose > 0 & merged_metadata$yearssincediagnose < 4 ~ "newly",
          merged_metadata$yearssincediagnose >= 4 & merged_metadata$yearssincediagnose < 10 ~ "medium",
          merged_metadata$yearssincediagnose >= 10 & merged_metadata$yearssincediagnose < 30 ~ "long",
          merged_metadata$yearssincediagnose >= 30 & merged_metadata$yearssincediagnose < 40 ~ "long",
          merged_metadata$yearssincediagnose >= 50 & merged_metadata$yearssincediagnose < 60 ~ "SixthDecade",
          is.na(merged_metadata$yearssincediagnose) ~ "Missing")

dds <- DESeqDataSetFromMatrix(countData=merged_rawcounts, colData=merged_metadata, design= ~Diagnose + PlateNr)# + yearssincebinned
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

merged_vst_counts <- as.data.table(sapply(data.table(assay(vsd)), as.numeric))

#Second DESEQ2 object for combined PSC+PSC/UC  vs UC (or CON) contrast:
merged_metadata2 <- merged_metadata
merged_metadata2$Diagnose <- gsub("PSC$|PSCUC$","PSC+PSCUC", merged_metadata2$Diagnose)
dds2 <- DESeqDataSetFromMatrix(countData=merged_rawcounts, colData=merged_metadata2, design= ~Diagnose + PlateNr)
dds2 <- DESeq(dds2)
#results(dds2,  contrast=c("Diagnose","PSC+PSCUC","Control"),tidy=T)

#####Optional: Validation VST-Data
# test_vst_counts <- as.data.table(t(merged_vst_counts))
# colnames(test_vst_counts) <- rownames(assay(vsd))
# rownames(test_vst_counts) <- SampleIDs
##Change y to any gene you wish to look up:
# ggplot(test_vst_counts)+geom_jitter(aes(x=merged_metadata$Diagnose,y=HIF1A))+facet_wrap(~merged_metadata$Cohort)

# Transformation PSC+UC####
#Two logical vectors to filter for Controls of each Cohort in dataset:
merged_controls_PSC <- ifelse(merged_metadata[grepl("^I",merged_metadata$SampleID)]$Diagnose == "Control",TRUE,FALSE)
merged_controls_UC <- ifelse(merged_metadata[grepl("^DE",merged_metadata$SampleID)]$Diagnose == "Control",TRUE,FALSE)

##Perform outlierfiltered transformation based on controls for each Cohort:
PSC_z_con_stabilised <- outlierfiltered_control_transformation(dataset = merged_vst_counts[,grepl("^I",colnames(merged_vst_counts)),with=FALSE],
                                                               controlsvector = merged_controls_PSC)
UC_z_con_stabilised <- outlierfiltered_control_transformation(dataset = merged_vst_counts[,grepl("^DE",colnames(merged_vst_counts)),with=FALSE],
                                                              controlsvector = merged_controls_UC)
#Combine transformated result tables to one merged dataset:
merged_z_con_stabilised <- data.table(rn=rownames(assay(vsd)),PSC_z_con_stabilised,UC_z_con_stabilised)
merged_z_con_stabilised <- transpose_datatable(merged_z_con_stabilised)
#Add Metadata to merged count-table:
merged_RF <- data.table(merged_z_con_stabilised, Diagnose=merged_metadata$Diagnose, Cohort=c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))

#Convert Diagnose into a factor: 0=Control, 1=PSC, 2=PSCUC, 3=UC
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
levels(merged_RF$Diagnose) <- c("0","1","2", "3")
#Cohort: 0=PSC, 1=UC
merged_RF$Cohort  <- factor(merged_RF$Cohort )

PSC_UC_COHORT_results <- performML(dataset = merged_RF[merged_RF$Diagnose == 0,-c("rn","Diagnose")], splitfactor = 0.5, outcome_column = "Cohort",seed = seednr)

pROC::auc(PSC_UC_COHORT_results$pROC_object)

# UC vs CON RF Model####
#Convert Diagnose into a factor: 0=Control, 1=PSC, 2=PSCUC, 3=UC
merged_RF$Diagnose  <- factor(merged_RF$Diagnose )
levels(merged_RF$Diagnose) <- c("0","1","2", "3")

merged_RF$Cohort  <- factor(merged_RF$Cohort )

filtered_merged_RF <- merged_RF[merged_RF$Cohort == 1,]
save(filtered_merged_RF, file="filtered_merged_RF.Rdata")
levels(filtered_merged_RF$Diagnose) <- droplevels(factor(filtered_merged_RF$Diagnose))
#0 = Control, 1= UC
levels(filtered_merged_RF$Diagnose) <- c("0", "1")#0 = Control, 1= UC

CON_UC_Diagnose_results <- performML(dataset = filtered_merged_RF[,-c("Cohort","rn")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(CON_UC_Diagnose_results$pROC_object)



# #glmnet performance:
# pROC::auc(performML(dataset = filtered_merged_RF[,-c("Cohort","rn")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr, mlmethod="glmnet")$pROC_object)
# pROC::ci.auc(performML(dataset = filtered_merged_RF[,-c("Cohort","rn")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr, mlmethod="glmnet")$pROC_object)


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

validation_controls <- ifelse(dds_validation$Diagnose == "Control",TRUE,FALSE)
validation_vst_counts <- data.table(assay(vsd_validation)[,])

validation_z_con_stabilised <- outlierfiltered_control_transformation(validation_vst_counts, validation_controls)

validation_z_con_stabilised <- data.table(rn=rownames(assay(vsd_validation)), validation_z_con_stabilised)
validation_z_con_stabilised <- transpose_datatable(validation_z_con_stabilised)


# Mo RF for validation Model####
Mo_genefeatures <- c("rn",rownames(assay(vsd_validation)))
#sum(colnames(filtered_merged_RF) %in% Mo_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo
Mo_dataset <- data.table(validation_z_con_stabilised, Diagnose=vsd_validation$Diagnose)


reduced_merged_RF_formo <- filtered_merged_RF[,intersect(colnames(Mo_dataset), colnames(filtered_merged_RF)),with=F]
reduced_Mo_RF <- Mo_dataset[,intersect(colnames(Mo_dataset), colnames(filtered_merged_RF)),with=F]


Mo_CON_UC_Diagnose_results <- performML(dataset = reduced_merged_RF_formo[,-c("rn")], testing_dataset = reduced_Mo_RF,  splitfactor = 0.5, outcome_column = "Diagnose", seed = seednr)

pROC::auc(Mo_CON_UC_Diagnose_results$pROC_object)
pROC::ci.auc(Mo_CON_UC_Diagnose_results$pROC_object)


# #glmnet performance:
# pROC::auc(performML(dataset = reduced_merged_RF_formo[,-c("rn")], testing_dataset = reduced_Mo_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr, mlmethod="glmnet")$pROC_object)
# pROC::ci.auc(performML(dataset = reduced_merged_RF_formo[,-c("rn")], testing_dataset = reduced_Mo_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr, mlmethod="glmnet")$pROC_object)
# # perforamnce on our dataset:
# performML(dataset = reduced_merged_RF_formo[,-c("rn")], testing_dataset = reduced_Mo_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr, mlmethod="glmnet") %>% 
#   pluck("tuned_model") %>%
#   predictwithatunedmodel(tuned_model = ., testing_dataset = Mo_CON_UC_Diagnose_results$test_df, outcome_column = "Diagnose",seed = seednr) %>% 
#   pluck("pROC_object") %>%
#   pROC::auc(.)


# Planell validation Preparation####
validation_planell_cohortname<- "Planell"
validation_planellcohorts <- list(Ostrowski ="UCAI Based\nDisease Severity",
                                  Mo_UC="Diagnose",
                                  Planell="Endoscopic\nMayo Score")
vsd_validation_planell <- getCohortsVSD(cohortname=validation_planell_cohortname)
severityscale <- validation_planellcohorts[[validation_planell_cohortname]]
dds_validation_planell <- getCohortsDDS(cohortname=validation_planell_cohortname)
vsd_validation_planell$severity <- vsd_validation_planell$Diagnose
vsd_validation_planell$severity[is.na(vsd_validation_planell$severity)] <- "Control"  

validation_planell_controls <- ifelse(dds_validation_planell$Diagnose == "Control",TRUE,FALSE)
validation_planell_vst_counts <- data.table((assay(vsd_validation_planell)[,]))

validation_planell_z_con_stabilised <- data.table(rn=c(rownames(assay(vsd_validation_planell))),
                                                  outlierfiltered_control_transformation(validation_planell_vst_counts, controlsvector = validation_planell_controls)
)

validation_planell_z_con_stabilised <- transpose_datatable(validation_planell_z_con_stabilised)

# Planell RF for validation_planell Model####
Planell_genefeatures <- c("rn",rownames(assay(vsd_validation_planell)))
sum(colnames(filtered_merged_RF) %in% Planell_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our reduced to intersect-genes merged dataset and then try to predict outcome in Mo

Planell_dataset <- data.table(validation_planell_z_con_stabilised, Diagnose=vsd_validation_planell$Diagnose)

reduced2_merged_RF <- filtered_merged_RF[,intersect(colnames(Planell_dataset), colnames(filtered_merged_RF)),with=F]
reduced2_Planell_RF <- Planell_dataset[,intersect(colnames(Planell_dataset), colnames(filtered_merged_RF)),with=F]

levels(reduced2_Planell_RF[["Diagnose"]]) <- c(0L,1L)

Planell_CON_UC_Diagnose_results <- performML(dataset = reduced2_merged_RF[,-c("rn")],testing_dataset = reduced2_Planell_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)
reducedforplanell_CON_UC_Diagnose_results <- predictwithatunedmodel(tuned_model = Planell_CON_UC_Diagnose_results$tuned_model, testing_dataset = Planell_CON_UC_Diagnose_results$test_df, outcome_column = "Diagnose",seed = seednr)
#Auroc for external validation dataset
pROC::auc(Planell_CON_UC_Diagnose_results$pROC_object)
pROC::ci.auc(Planell_CON_UC_Diagnose_results$pROC_object)
#validation with own test dataset
pROC::auc(reducedforplanell_CON_UC_Diagnose_results$pROC_object)



# #glmnet performance:
# pROC::auc(performML(dataset = reduced2_merged_RF[,-c("rn")],testing_dataset = reduced2_Planell_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr,mlmethod="glmnet")$pROC_object)
# pROC::ci.auc(performML(dataset = reduced2_merged_RF[,-c("rn")],testing_dataset = reduced2_Planell_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr,mlmethod="glmnet")$pROC_object)
# 
# 


# Ostrowski UC validation Preparation####
validation_ostrowski_cohortname<- "Ostrowski"
validation_ostrowskicohorts <- list(Ostrowski ="UCAI Based\nDisease Severity",
                                    Mo_UC="Diagnose",
                                    Planell="Endoscopic\nMayo Score")
vsd_validation_ostrowski <- getCohortsVSD(cohortname=validation_ostrowski_cohortname)

# ostrowski_res <- getCohortsRES(cohortname = "Ostrowski")

severityscale <- validation_ostrowskicohorts[[validation_ostrowski_cohortname]]
dds_validation_ostrowski <- getCohortsDDS(cohortname=validation_ostrowski_cohortname)
vsd_validation_ostrowski$severity <- vsd_validation_ostrowski$Diagnose
vsd_validation_ostrowski$severity[is.na(vsd_validation_ostrowski$severity)] <- "Control"  

validation_ostrowski_controls <- ifelse(dds_validation_ostrowski$Diagnose == "Control",TRUE,FALSE)
validation_ostrowski_vst_counts <- data.table((assay(vsd_validation_ostrowski)[,]))

validation_ostrowski_z_con_stabilised <- data.table(rn=c(rownames(assay(vsd_validation_ostrowski))),
                                                    outlierfiltered_control_transformation(validation_ostrowski_vst_counts, controlsvector = validation_ostrowski_controls)
)
validation_ostrowski_z_con_stabilised <- transpose_datatable(validation_ostrowski_z_con_stabilised)

# Ostrowski UC RF for validation_ostrowski Model####
Ostrowski_genefeatures <- c("rn",rownames(assay(vsd_validation_ostrowski)))
sum(colnames(filtered_merged_RF) %in% Ostrowski_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo

Ostrowski_dataset <- data.table(validation_ostrowski_z_con_stabilised, Diagnose=vsd_validation_ostrowski$Diagnose)

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





# #glmnet performance:
# pROC::auc(performML(dataset = reduced3_merged_RF[,-c("rn")],testing_dataset = reduced3_Ostrowski_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr,mlmethod="glmnet")$pROC_object)
# pROC::ci.auc(performML(dataset = reduced3_merged_RF[,-c("rn")],testing_dataset = reduced3_Ostrowski_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr,mlmethod="glmnet")$pROC_object)
# 
# performML(dataset = reduced3_merged_RF[,-c("rn")],testing_dataset = reduced3_Ostrowski_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr,mlmethod="glmnet") %>% 
#   pluck("tuned_model") %>%
#   predictwithatunedmodel(tuned_model = ., testing_dataset = Ostrowski_CON_UC_Diagnose_results$test_df, outcome_column = "Diagnose",seed = seednr) %>% 
#   pluck("pROC_object") %>%
#   pROC::auc(.)





# PSC vs UC model, no controls, no PSCUC ####

PSC_UC_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 2,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSC_UC_Diagnose_results$pROC_object)

ggplot(merged_RF[,-c("rn","Cohort")])+geom_jitter(aes(x=Diagnose,y=ACLY))+geom_boxplot(aes(x=Diagnose,y=ACLY))#+facet_wrap(~merged_metadata$Cohort)

PSC_UC_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 2,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)
#glmnet performance:
pROC::auc(performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 2,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr,mlmethod="glmnet")$pROC_object)
pROC::ci.auc(performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 2,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr,mlmethod="glmnet")$pROC_object)


# PSC vs Controls model, no UC, no PSCUC ####

PSC_CONTROL_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 2 & merged_RF$Diagnose != 3,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSC_CONTROL_Diagnose_results$pROC_object)

#0=Control, 1=PSC, 2=PSCUC, 3=UC
ggplot(merged_RF[,-c("rn","Cohort")],aes(x=Diagnose,y=BACH2))+geom_jitter()+geom_boxplot()+scale_x_discrete(labels = c('Control',"PSC","PSC/UC","UC"))#+facet_wrap(~merged_metadata$Cohort)


# PSC vs PSCUC, no CON, no UC ###########

PSC_PSCUC_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 3,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSC_PSCUC_Diagnose_results$pROC_object)


# PSCUC vs UC, no CON, no PSC ###########
PSCUC_UC_Diagnose_results <- performML(dataset = merged_RF[merged_RF$Diagnose != 0 & merged_RF$Diagnose != 1,-c("rn","Cohort")], splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)

pROC::auc(PSCUC_UC_Diagnose_results$pROC_object)


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

# PSCUC vs CON; no UC; no PSC###########
#Diagnose: 0=Control, 1=PSC, 2=PSCUC, 3=UC
PSCUC_CON_merged_RF <- merged_RF[!(merged_RF$Diagnose %in% c(3,1)),-c("rn","Cohort")]
PSCUC_CON_merged_RF$Diagnose <- droplevels(PSCUC_CON_merged_RF$Diagnose)
#merge PSC and PSCUC to one level; PSC+PSCUC = 1; Con = 0:
levels(PSCUC_CON_merged_RF$Diagnose) <- c("0","1")
#PSCUC_CON_merged_RF$Diagnose <- droplevels(PSCUCPSC_merged_RF$Diagnose)

PSCUC_CON_Diagnose_results <- performML(dataset = PSCUC_CON_merged_RF,
                                        splitfactor = 0.5, 
                                        outcome_column = "Diagnose",
                                        seed = seednr)

pROC::auc(PSCUC_CON_Diagnose_results$pROC_object)


# CEMiTool #####
#create an expressionmatrix dataframe for CEMItool:
merged_RF_expressionmatrix <- transpose_datatable(merged_RF[,-c("Diagnose", "Cohort")])
merged_RF_expressionmatrix <- data.frame(merged_RF_expressionmatrix,row.names = merged_RF_expressionmatrix$rn)[,-1]

#perform CEMItool analysis
cem_merged_RF <- CEMiwrapper(expressionmatrix=merged_RF_expressionmatrix, ID=merged_metadata$SampleID, Groups=merged_metadata$Diagnose, reportname=paste0("PSC_PSCUC_UC","_Diagnose"), applyfiltering = FALSE)
# save.image(file="performedanalysis.RData")

# # with KEGG
# cem_merged_RF_KEGG <- CEMiwrapper(expressionmatrix=merged_RF_expressionmatrix, ID=merged_metadata$SampleID, Groups=merged_metadata$Diagnose, reportname=paste0("PSC_PSCUC_UC","_Diagnose","_KEGG"), applyfiltering = FALSE, gmt_location=paste0(projectdir,"/data_rescources/KEGG_2021_Human.txt"))
# # With Reactome
# cem_merged_RF_REACTOME <- CEMiwrapper(expressionmatrix=merged_RF_expressionmatrix, ID=merged_metadata$SampleID, Groups=merged_metadata$Diagnose, reportname=paste0("PSC_PSCUC_UC","_Diagnose","_REACTOME"), applyfiltering = FALSE, gmt_location=paste0(projectdir,"/data_rescources/Reactome_2022.txt"))


# Cell Blueprints ####
#Cell Blueprints LM22 for every module created from cemitools:
enrichmentslist <- list()
for (x in unique(cem_merged_RF@module$modules)){
  print(x)
  genesinput <- data.table(cem_merged_RF@module)[modules == x,]$genes
  result <- blueprintenrichments(genes=genesinput)[[1]]
  # result <- result[result$enrichment >= 0.1,]
  enrichmentslist[[x]] <- result
}

enrichment_df <- data.table(bind_rows(enrichmentslist, .id = "column_label")) %>% pivot_wider(names_from = celltype, values_from = enrichment) %>% as.data.table() %>% as.matrix(., rownames="column_label")

htmp <- Heatmap(enrichment_df, 
                name="Z-score",
                heatmap_legend_param = list(
                  legend_direction = "horizontal", 
                  legend_width = unit(6, "cm")),
                height = unit(10, "cm"),
                row_km = 1,
                row_km_repeats = 10,
                row_names_side = "left",
                width = unit(14, "cm"),
                column_names_rot = 45, 
                # column_names_max_height = unit(2, "cm"),
                column_names_side = "bottom",
                column_names_gp = gpar(fontsize = 10)
)
draw(htmp, heatmap_legend_side="top")
fwrite(data.table(module=rownames(enrichment_df),enrichment_df), file = "output/cell_type_enrichment_modules_zscores.csv")
# topGO function results #####
#use topgo function to evaluate the biological function of features, chosen by RF:
res <- results(dds,  contrast=c("Diagnose","UC","Control"),tidy=T)

topGo_object <- topGO_enrichment(genelist = PSCUC_UC_Diagnose_results$variable_importance[1:1000,]$feature, statistical_test = "ks.ties", weights = PSCUC_UC_Diagnose_results$variable_importance[1:1000,]$importance, #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                 DESeq_result = res, match_by_expression = TRUE)
# topGo_object <- topGO_enrichment(genelist = PSCUC_UC_Diagnose_results$variable_importance[1:100,]$feature, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
#                                  DESeq_result = res, match_by_expression = TRUE, ontology_type = "MF")



#
# Analysis Results ####

# DEseq2 results:
#PSC vs CON
res_PSC <- results(dds,  contrast=c("Diagnose","PSC","Control"),tidy=T)
resLFC_PSC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

#PSC/PSCUC vs CON
res_PSCUCPSC_CON <- results(dds2, contrast=c("Diagnose", "PSC+PSCUC", "Control"),tidy=T)
fwrite(res_PSCUCPSC_CON, file="output/deseq_PSCUCPSC_CON.csv")
#UC vs PSC/PSCUC
res_UC_PSCUCPSC <- results(dds2, contrast=c("Diagnose", "UC","PSC+PSCUC"),tidy=T)
fwrite(res_UC_PSCUCPSC, file="output/deseq_PSCUCPSC_UC.csv")

#UC vs CON
res_UC <- results(dds,  contrast=c("Diagnose","UC","Control"),tidy=T)
# resLFC_UC <- lfcShrink(dds, coef="Diagnose_UC_vs_Control", type="apeglm")
fwrite(res_UC, file="output/deseq_UC_CON.csv")
#UC vs PSC
res_UC_PSC <- results(dds,  contrast=c("Diagnose","UC","PSC"),tidy=T)
# resLFC_UC_PSC <- lfcShrink(dds, coef="Diagnose_UC_vs_PSC", type="apeglm")
fwrite(res_UC_PSC, file="output/deseq_UC_PSC.csv")
#PSCUC vs PSC
res_PSCUC_PSC <- results(dds,  contrast=c("Diagnose","PSCUC","PSC"),tidy=T)
# resLFC_PSCUC_PSC <- lfcShrink(dds, coef="Diagnose_PSCUC_vs_PSC", type="apeglm")
fwrite(res_PSCUC_PSC, file="output/deseq_PSCUC_PSC.csv")

res_PSCUC_CON <- results(dds,  contrast=c("Diagnose","PSCUC","Control"),tidy=T)
fwrite(res_PSCUC_CON  %>% arrange(padj), file="output/deseq_PSCUC_CON.csv")
res_PSC_CON <- results(dds,  contrast=c("Diagnose","PSC","Control"),tidy=T)
fwrite(res_PSC_CON %>% arrange(padj), file="output/deseq_PSC_CON.csv")



#Overlap index between PSCUC and PSC, both against CON:
intersect(res_PSC_CON[res_PSC_CON$padj < 0.01,"row"],
          res_PSCUC_CON[res_PSCUC_CON$padj < 0.01,"row"]) %>%
  length() / min(res_PSC_CON[res_PSC_CON$padj < 0.01,"row"] %>% length(),
               res_PSCUC_CON[res_PSCUC_CON$padj < 0.01,"row"] %>% length())


# Random Forest results:
#Validations:
Planell_CON_UC_Diagnose_results
Mo_CON_UC_Diagnose_results
Ostrowski_CON_UC_Diagnose_results
# Ostrowski_CON_PSC_Diagnose_results


#Our Dataset:
CON_UC_Diagnose_results

PSC_UC_COHORT_results
PSC_UC_Diagnose_results
PSCUC_UC_Diagnose_results
PSC_PSCUC_Diagnose_results
PSC_CONTROL_Diagnose_results
PSCUCPSC_UC_Diagnose_results
PSCUCPSC_CON_Diagnose_results
PSCUC_CON_Diagnose_results

pROC::auc(PSCUC_CON_Diagnose_results$pROC_object)
pROC::ggroc(PSCUC_CON_Diagnose_results$pROC_object)
pROC::ci.auc(PSCUC_CON_Diagnose_results$pROC_object)


#Cemitool
cem_merged_RF

#Blueprints

png("output/celltype_modules_heatmap.png",width=7.5,height=7.5,units="in",res=1200)
draw(htmp, heatmap_legend_side="top")
dev.off()
#TopGO
topGo_object


# NES Heatmap ####
cem_merged_RF@enrichment$nes
nes_matrix <- as.matrix(cem_merged_RF@enrichment$nes[,-1])
rownames(nes_matrix) <- cem_merged_RF@enrichment$nes$pathway 


nes_heatmap <- Heatmap(nes_matrix, 
                       name="NES",
                       heatmap_legend_param = list(
                         legend_direction = "horizontal", 
                         legend_width = unit(6, "cm")),
                       height = unit(10, "cm"),
                       row_km = 1,
                       row_km_repeats = 10,
                       row_names_side = "left",
                       width = unit(14, "cm"),
                       column_names_rot = 0, 
                       # column_names_max_height = unit(2, "cm"),
                       column_names_side = "bottom",
                       column_names_gp = gpar(fontsize = 10)
)


png("output/figure3_nesheatmap.png",width=7.5,height=7.5,units="in",res=1200)
draw(nes_heatmap, heatmap_legend_side="top")
dev.off()



# TOPGO RF enrichment sets ####
#Supplemmentary Tables 8 to 11:
CON_UC_features <- mostimportantfeatures(variable_importance_table = CON_UC_Diagnose_results$variable_importance, limit = 50)

CON_UC_Diagnose_results_topgo_50_results <- topGO_enrichment(genelist = CON_UC_features$feature,
                                                             statistical_test = "Fisher",#"ks.ties", 
                                                             ontology_type = "BP",
                                                             weights = mostimportantfeatures(variable_importance_table = CON_UC_Diagnose_results$variable_importance, limit = 50)$importance, 
                                                             #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                                             DESeq_result = res_UC,
                                                             match_by_expression = TRUE)
CON_UC_Diagnose_results_topgo_50_results2 <- topGO_enrichment(genelist = CON_UC_features$feature,
                                                              statistical_test = "ks.ties", 
                                                              ontology_type = "BP",
                                                              weights = mostimportantfeatures(variable_importance_table = CON_UC_Diagnose_results$variable_importance, limit = 50)$importance, 
                                                              #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                                              DESeq_result = res_UC,
                                                              match_by_expression = TRUE)
fwrite(CON_UC_Diagnose_results_topgo_50_results, file = "output/ST8_CON_UC_Diagnose_results_topgo_50_results.csv")


PSCUCPSC_CON_Diagnose_results_topgo_50_results <- topGO_enrichment(genelist = mostimportantfeatures(variable_importance_table = PSCUCPSC_CON_Diagnose_results$variable_importance, limit = 50)$feature,
                                                                   statistical_test = "ks.ties",
                                                                   weights = mostimportantfeatures(variable_importance_table = PSCUCPSC_CON_Diagnose_results$variable_importance, limit = 50)$importance,
                                                                   gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature,
                                                                   #DESeq_result = res_UC,
                                                                   ontology_type = "BP",
                                                                   match_by_expression = TRUE)
fwrite(PSCUCPSC_CON_Diagnose_results_topgo_50_results, file = "output/PSCUCPSC_CONTROL_Diagnose_results_topgo_50_results.csv")
fwrite(PSCUCPSC_CON_Diagnose_results_topgo_50_results, file = "output/ST9_PSCUCPSC_CONTROL_Diagnose_results_topgo_50_results.csv")

PSC_CONTROL_Diagnose_results_topgo_50_results <- topGO_enrichment(genelist = mostimportantfeatures(variable_importance_table = PSC_CONTROL_Diagnose_results$variable_importance, limit = 50)$feature,
                                                                  statistical_test = "ks.ties", 
                                                                  weights = mostimportantfeatures(variable_importance_table = PSC_CONTROL_Diagnose_results$variable_importance, limit = 50)$importance, 
                                                                  #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                                                  DESeq_result = res_UC, 
                                                                  ontology_type = "BP",
                                                                  match_by_expression = TRUE)
fwrite(CON_UC_Diagnose_results_topgo_50_results, file = "output/PSC_CONTROL_Diagnose_results_topgo_50_results.csv")



PSCUCPSC_UC_Diagnose_results_topgo_50_results <- topGO_enrichment(genelist = mostimportantfeatures(variable_importance_table = PSCUCPSC_UC_Diagnose_results$variable_importance, limit = 50)$feature,
                                                                  statistical_test = "ks.ties", 
                                                                  weights = mostimportantfeatures(variable_importance_table = PSCUCPSC_UC_Diagnose_results$variable_importance, limit = 50)$importance, 
                                                                  #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                                                  DESeq_result = res_UC, 
                                                                  ontology_type = "BP",
                                                                  match_by_expression = TRUE)
fwrite(PSCUCPSC_UC_Diagnose_results_topgo_50_results, file = "output/PSCUCPSC_UC_Diagnose_results_topgo_50_results.csv")
fwrite(PSCUCPSC_UC_Diagnose_results_topgo_50_results, file = "output//ST10_PSCUCPSC_UC_Diagnose_results_topgo_50_results.csv")


#FEATURE IMPORTANCES SUPPLEMENTARY TABLE:
ST13 <- left_join(CON_UC_Diagnose_results$variable_importance, PSCUCPSC_CON_Diagnose_results$variable_importance, by = "feature") %>%
  left_join(., PSCUCPSC_UC_Diagnose_results$variable_importance, by = "feature") %>%
  left_join(., Planell_CON_UC_Diagnose_results$variable_importance, by = "feature") %>%
  left_join(., Mo_CON_UC_Diagnose_results$variable_importance, by = "feature") %>%
  left_join(., Ostrowski_CON_UC_Diagnose_results$variable_importance, by = "feature") 
colnames(ST13) <- c("features", "CON_UC_Diagnose_results", "PSCUCPSC_CON_Diagnose_results", "PSCUCPSC_UC_Diagnose_results", "Planell_CON_UC_Diagnose_results", "Mo_CON_UC_Diagnose_results", "Ostrowski_CON_UC_Diagnose_results" )
fwrite(ST13, file = "output/ST13_RF_feature_importances.csv")


# SUPPLEMENTARY TABLE, CEMITOOLS Modules
fwrite(cem_merged_RF@module[order(cem_merged_RF@module$modules),], file = "output/ST14_RF_feature_importances.csv")


# SUPPLEMENTARY TABLE, CEMITOOLS Modules



fwrite(cem_merged_RF@ora, file = "output/ST15_cemitool_ora.csv")


# FUW ######

# sex imputation
data.table(females=merged_RF[XIST>0,rn]) %>% fwrite(file = "output/females.txt")
data.table(males=merged_RF[XIST<=0,rn]) %>% fwrite(file = "output/males.txt")



#merged_counts export

merged_rawcounts_DT <- as.data.table(merged_rawcounts,keep.rownames = T)
setnames(merged_rawcounts_DT,old="rn",new = "Gene")
#subset to smaller tables because cibersort is crashing with to many data due to memory quota 500mb
#PSC 2-298
fwrite(merged_rawcounts_DT[,c(1:298)],file = "merged_rawcounts_PSC.txt",sep = "\t")
#PSC first section
#UC 299-1036 split in 299-675 and 676-1036
fwrite(merged_rawcounts_DT[,c(1,299:675)],file = "merged_rawcounts_UC1.txt",sep = "\t")
fwrite(merged_rawcounts_DT[,c(1,676:1036)],file = "merged_rawcounts_UC2.txt",sep = "\t")
#PSC controls
psccontrols <- merged_metadata[Diagnose=="Control"&Cohort==0,SampleID] %>% c("Gene",.)
fwrite(merged_rawcounts_DT[,..psccontrols],file = "merged_rawcounts_PSC_Controls.txt",sep = "\t")
#PSC cases
psccases <- merged_metadata[Diagnose=="PSC"&Cohort==0,SampleID] %>% c("Gene",.)
fwrite(merged_rawcounts_DT[,..psccases],file = "merged_rawcounts_PSC_Cases.txt",sep = "\t")
#PSCUC cases
pscuccases <- merged_metadata[Diagnose=="PSCUC"&Cohort==0,SampleID] %>% c("Gene",.)
fwrite(merged_rawcounts_DT[,..pscuccases],file = "merged_rawcounts_PSCUC_Cases.txt",sep = "\t")
#UC controls
uccontrols <- merged_metadata[Diagnose=="Control"&Cohort==1,SampleID] %>% c("Gene",.)
fwrite(merged_rawcounts_DT[,..uccontrols],file = "merged_rawcounts_UC_Controls.txt",sep = "\t")
#UC cases
uccases <- merged_metadata[Diagnose=="UC"&Cohort==1,SampleID] %>% c("Gene",.)
fwrite(merged_rawcounts_DT[,..uccases],file = "merged_rawcounts_UC_Cases.txt",sep = "\t")


fwrite(merged_rawcounts_DT,file = "merged_rawcounts.txt",sep = "\t")
# DESEq2 results
merged_metadata <- as.data.table(merged_metadata)
merged_metadata[,Diagnose2:=case_when(Diagnose=="Control" ~ "Control",
                                      Diagnose=="UC" ~ "UC",
                                      Diagnose=="PSC" ~ "PSC",
                                      Diagnose=="PSCUC" ~ "PSC",
                                      TRUE ~ NA_character_)]
dds_PSC_as_one <- DESeqDataSetFromMatrix(countData=merged_rawcounts, colData=merged_metadata, design= ~Diagnose2 + PlateNr)
dds_PSC_as_one <- DESeq(dds_PSC_as_one)
#res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
#resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
vsd_PSC_as_one <- vst(dds_PSC_as_one, blind=FALSE)
merged_vst_counts_PSC_as_one <- as.data.table(sapply(data.table(assay(vsd_PSC_as_one)), as.numeric))
results(dds, contrast = c("Diagnose","PSC","UC"),tidy=T)
results(dds_PSC_as_one, contrast = c("Diagnose2","PSC","UC"),tidy=T) #should be different


##### FIGURE 2 (if you exculde graphical abstract, it's 1)####
##Lookup single gene expressions:
#ggplot(test_vst_counts)+geom_jitter(aes(x=merged_metadata$Diagnose,y=MAN2A))+facet_wrap(~merged_metadata$Cohort)
size.rel = 1L

#UC vs CON
res_UC_untidy <- results(dds,  contrast=c("Diagnose","UC","Control"),tidy=F)
#resLFC_UC <- lfcShrink(dds, coef="Diagnose_UC_vs_Control", type="apeglm")
res_UC_tidy <- results(dds,  contrast=c("Diagnose","UC","Control"),tidy=T)
#PSC vs CON
res_PSC_CON_untidy <- results(dds,  contrast=c("Diagnose","PSC","Control"),tidy=F)
#resLFC_UC_PSC <- lfcShrink(dds, coef="Diagnose_UC_vs_PSC", type="apeglm")
res_PSC_CON_tidy <- results(dds,  contrast=c("Diagnose","PSC","Control"),tidy=T)
#UC vs PSC
res_UC_PSC_untidy <- results(dds,  contrast=c("Diagnose","UC","PSC"),tidy=F)
#resLFC_UC_PSC <- lfcShrink(dds, coef="Diagnose_UC_vs_PSC", type="apeglm")
res_UC_PSC_tidy <- results(dds,  contrast=c("Diagnose","UC","PSC"),tidy=T)
#PSCUC vs PSC
res_PSCUC_PSC_untidy <- results(dds,  contrast=c("Diagnose","PSCUC","PSC"),tidy=F)
#resLFC_PSCUC_PSC <- lfcShrink(dds, coef="Diagnose_PSCUC_vs_PSC", type="apeglm")
res_PSCUC_PSC_tidy <- results(dds,  contrast=c("Diagnose","PSCUC","PSC"),tidy=T)

res_PSCandPSCUC_CON_tidy <- results(dds2,  contrast=c("Diagnose","PSC+PSCUC","Control"),tidy=T)
res_PSCandPSCUC_CON_untidy <- results(dds2,  contrast=c("Diagnose","PSC+PSCUC","Control"),tidy=F)

res_UC_PSCandPSCUC_tidy <- results(dds2,  contrast=c("Diagnose","UC","PSC+PSCUC"),tidy=T)
res_UC_PSCandPSCUC_untidy <- results(dds2,  contrast=c("Diagnose","UC","PSC+PSCUC"),tidy=F)

volcano_UCvsCON <- EnhancedVolcano_function(res_UC_untidy, plottitle="UC vs CON",plotnumber=element_blank())+ 
  theme_linedraw()+
  theme(legend.position = "none", 
        plot.margin = margin(5.5,20,5.5,27, "pt"),
        axis.title.x = element_blank(),        
        plot.tag.position = c(-.01, .97),
        #plot.tag = element_text(size = rel(1.5 * size.rel)),
        axis.text = element_text(size=unit(17,"points")),
        axis.title = element_text(size=unit(17,"points")),
        plot.tag =element_text(size=unit(22,"points")),
        plot.subtitle = element_text(size=unit(22,"points"))
        #axis.title.y = element_blank()
  )+
  labs(tag="a)")# axis.title=element_text(size=12,face="bold"),

volcano_PSCvsCON <- EnhancedVolcano_function(res_PSCandPSCUC_CON_untidy, plottitle="PSC alone + PSC/UC vs CON",plotnumber=element_blank())+  
  theme_linedraw()+
  theme(legend.position = "none", 
        plot.margin = margin(5.5,20,5.5,27, "pt"),
        axis.title.x = element_blank(),        
        plot.tag.position = c(-.01, .97),
        #plot.tag = element_text(size = rel(1.5 * size.rel)),
        axis.text = element_text(size=unit(17,"points")),
        axis.title = element_text(size=unit(17,"points")),
        plot.tag =element_text(size=unit(22,"points")),
        plot.subtitle = element_text(size=unit(22,"points"))
        #axis.title.y = element_blank())
  )+
  
  labs(tag="b)")
volcano_UCvsPSC <- EnhancedVolcano_function(res_UC_PSCandPSCUC_untidy, plottitle="UC vs PSC alone + PSC/UC",plotnumber=element_blank())+ 
  theme_linedraw()+
  theme(legend.position = "bottom", 
        #legend.position = "none", 
        plot.margin = margin(5.5,20,5.5,27, "pt"),
        #axis.title.x = element_blank(),        
        legend.title = element_blank(),
        plot.tag.position = c(-.01, .97),
        #plot.tag = element_text(size = rel(1.5 * size.rel)),
        legend.text=element_text(size=rel(1.5 * size.rel)),
        axis.text = element_text(size=unit(17,"points")),
        axis.title = element_text(size=unit(17,"points")),
        plot.tag =element_text(size=unit(22,"points")),
        plot.subtitle = element_text(size=unit(22,"points"))
        #axis.title.y = element_blank()
  )+
  labs(tag="c)")

#MAN2A1 is heavily influenced by one sample outlier, so we remove it from the plot
volcano_PSCUCvsPSC <- EnhancedVolcano_function(res_PSCUC_PSC_untidy[!(rownames(res_PSCUC_PSC_untidy) %in% "MAN2A1"),],
                                               plottitle="PSC/UC vs PSC alone",
                                               plotnumber=element_blank(),
                                               #labsize=6
)+ 
  theme_linedraw()+
  theme(legend.position = "none", 
        plot.margin = margin(5.5,20,5.5,27, "pt"),
        #axis.title.x = element_blank(),        
        legend.title = element_blank(),
        plot.tag.position = c(-.01, .97),
        #plot.tag = element_text(size = rel(1.5 * size.rel)),
        legend.text=element_text(size=rel(1.5 * size.rel)),
        axis.text = element_text(size=unit(17,"points")),
        axis.title = element_text(size=unit(17,"points")),
        plot.tag =element_text(size=unit(22,"points")),
        plot.subtitle = element_text(size=unit(22,"points"))
        #axis.title.y = element_blank()
  )+
  labs(tag="d)")


#TopGO enrichments for sign genesets from volcano pvalue and log2 FC:
sign_UC <- res_UC_tidy[res_UC_tidy$padj < 0.01 & abs(res_UC_tidy$log2FoldChange)>1,]$row
sign_UC_PSC <-res_UC_PSC_tidy[res_UC_PSC_tidy$padj < 0.01 & abs(res_UC_PSC_tidy$log2FoldChange)>1,]$row %>% .[!is.na(.)]
sign_PSC_CON <-res_PSC_CON_tidy[res_PSC_CON_tidy$padj < 0.01 & abs(res_PSC_CON_tidy$log2FoldChange)>1,]$row %>% .[!is.na(.)]
sign_PSCUC_PSC <-res_PSCUC_PSC_tidy[res_PSCUC_PSC_tidy$padj < 0.01 & abs(res_PSCUC_PSC_tidy$log2FoldChange)>1,]$row #is empty



res_UC_CON_topgo <- topGO_enrichment(genelist = sign_UC, statistical_test = "Fisher", gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature,
                                     DESeq_result = res_UC_tidy, match_by_expression = TRUE, ontology_type = "BP")
res_UC_PSC_topgo <- topGO_enrichment(genelist = sign_UC_PSC, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                     DESeq_result = res_UC_tidy, match_by_expression = TRUE, ontology_type = "BP")
res_PSC_CON_topgo <- topGO_enrichment(genelist = sign_PSC_CON, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                      DESeq_result = res_UC_tidy, match_by_expression = TRUE, ontology_type = "BP")

fwrite(res_UC_CON_topgo, file = "output/topgo_UC_CON.csv")
fwrite(res_UC_PSC_topgo, file = "output/topgo_UC_PSC.csv")
fwrite(res_PSC_CON_topgo, file = "output/topgo_PSC_CON.csv")

sign_UC_PSCUCPSC <- res_UC_PSCandPSCUC_tidy[res_UC_PSCandPSCUC_tidy$padj < 0.01 & abs(res_UC_PSCandPSCUC_tidy$log2FoldChange)>1,]$row

res_UC_PSCUCPSC_topgo <- topGO_enrichment(genelist = sign_UC_PSCUCPSC, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                          DESeq_result = res_UC_tidy, match_by_expression = TRUE, ontology_type = "BP")
fwrite(res_UC_PSCUCPSC_topgo, file = "output/topgo_UC_PSCUCPSC.csv")

sign_PSCUCPSC_CON <- res_PSCandPSCUC_CON_tidy[res_PSCandPSCUC_CON_tidy$padj < 0.01 & abs(res_PSCandPSCUC_CON_tidy$log2FoldChange)>1,]$row

sign_PSCUCPSC_CON_topgo <- topGO_enrichment(genelist = sign_PSCUCPSC_CON, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                                            DESeq_result = res_PSCandPSCUC_CON_tidy, match_by_expression = TRUE, ontology_type = "BP")
fwrite(sign_PSCUCPSC_CON_topgo, file = "output/topgo_PSCUCPSC_CON.csv")

#GO term plots:

GO_plot_UC_CON <- ggplot(data=gotermpreview(res_UC_CON_topgo), aes(y=Term, x=padj,#as.numeric(test.pvalue),#
                                                                   size=Significant/Annotated, color=as.numeric(padj)))+
  geom_point()+
  scale_colour_distiller(name = "P-Adjusted",
                         palette = "Reds",
                         direction="1",
                         trans = reverselog_trans(base=10),
                         limits=c(1, min(c(res_UC_CON_topgo$padj, sign_PSCUCPSC_CON_topgo$padj, res_UC_PSCUCPSC_topgo$padj))),
                         guide="colourbar"
  )+
  #theme_bw()+
  theme_linedraw()+
  scale_x_continuous(trans=reverselog_trans(base=10),
                     limits = c(1, min(c(res_UC_CON_topgo$padj, sign_PSCUCPSC_CON_topgo$padj, res_UC_PSCUCPSC_topgo$padj)))
  )+
  theme(
    plot.margin = margin(5.5,20,5.5,27, "pt"),
    legend.box="vertical",
    text = element_text(size=unit(17,"points")),
    #axis.text.y = element_text(size=unit(12,"points")),
    legend.position = "none", 
    legend.justification = "right",
    axis.title.x = element_blank()
  )+
  scale_size_continuous(name="Significant/Annotated",
                        limits =c(0,1),
                        breaks=c(0,0.2,0.4,0.6,0.8,1)
  )+
  guides(size=guide_legend(nrow = 2),#guide_legend("Significant/Annotated"),
         #colour=guide_legend("P-Value")
         colour=guide_colorbar(direction ="horizontal",title.position = "left",barwidth = unit(2, "inch"))
  )+
  labs(x="-log10(P-Adjusted)",y=element_blank())

GO_plot_PSCUCPSC_CON <- ggplot(data= gotermpreview(sign_PSCUCPSC_CON_topgo), aes(y=Term, x=padj,#as.numeric(test.pvalue),#Significant/Annotated,
                                                                                 size=Significant/Annotated, color=as.numeric(padj)))+
  geom_point()+
  scale_colour_distiller(name = "P-Adjusted",
                         palette = "Reds",
                         direction="1",
                         trans = reverselog_trans(base=10),
                         limits=c(1, min(c(res_UC_CON_topgo$padj, sign_PSCUCPSC_CON_topgo$padj, res_UC_PSCUCPSC_topgo$padj))),
                         guide="colourbar"
  )+
  #theme_bw()+
  theme_linedraw()+
  scale_x_continuous(trans=reverselog_trans(base=10),
                     limits = c(1, min(c(res_UC_CON_topgo$padj, sign_PSCUCPSC_CON_topgo$padj, res_UC_PSCUCPSC_topgo$padj)))
  )+
  theme(
    plot.margin = margin(5.5,20,5.5,27, "pt"),
    legend.box="vertical",
    text = element_text(size=unit(17,"points")),
    #axis.text.y = element_text(size=unit(12,"points")),
    legend.position = "none", 
    legend.justification = "right",
    axis.title.x = element_blank()
  )+
  scale_size_continuous(name="Significant/Annotated",
                        limits =c(0,1),
                        breaks=c(0,0.2,0.4,0.6,0.8,1)
  )+
  guides(size=guide_legend(nrow = 2),#guide_legend("Significant/Annotated"),
         #colour=guide_legend("P-Value")
         colour=guide_colorbar(direction ="horizontal",title.position = "left",barwidth = unit(2, "inch"))
  )+
  labs(x="- log10(P-Adjusted)",y=element_blank())




GO_plot_UC_PSCUCPSC <- ggplot(data=gotermpreview(res_UC_PSCUCPSC_topgo), aes(y=Term, x=as.double(padj),#as.numeric(test.pvalue),#Significant/Annotated,
                                                                             size=Significant/Annotated, colour=as.double(padj)))+
  geom_point()+
  scale_colour_distiller(name = "P-Adjusted",
                         palette = "Reds",
                         direction="1",
                         trans = reverselog_trans(base=10),
                         limits=c(1, min(c(res_UC_CON_topgo$padj, sign_PSCUCPSC_CON_topgo$padj, res_UC_PSCUCPSC_topgo$padj))),
                         guide="colourbar"
  )+
  #theme_bw()+
  theme_linedraw()+
  scale_x_continuous(trans=reverselog_trans(base=10),
                     limits = c(1, min(c(res_UC_CON_topgo$padj, sign_PSCUCPSC_CON_topgo$padj, res_UC_PSCUCPSC_topgo$padj)))
  )+
  theme(
    plot.margin = margin(5.5,20,5.5,27, "pt"),
    legend.box="vertical",
    text = element_text(size=unit(17,"points")),
    #axis.text.y = element_text(size=unit(12,"points")),
    legend.position = "bottom", 
    legend.justification = "right"
  )+
  scale_size_continuous(name=bquote(frac("Significant","Annotated")),
                        limits =c(0,1),
                        breaks=c(0,0.2,0.4,0.6,0.8,1)
  )+
  guides(size=guide_legend(nrow = 2),#guide_legend("Significant/Annotated"),
         #colour=guide_legend("P-Value")
         colour=guide_colorbar(direction ="horizontal",title.position = "left",barwidth = unit(2, "inch"))
  )+
  labs(x=bquote(~ - ~Log[10] ~ italic(P- Adjusted)),y=element_blank())


# 
arranged_volcano <-ggarrange(ncol = 2, nrow=1,
                             ggarrange( ncol = 1,nrow=5, heights=c(1,1,1,1,0.2), align ="v",
                                        volcano_UCvsCON,volcano_PSCvsCON,volcano_UCvsPSC+theme(legend.position = "none"),volcano_PSCUCvsPSC+theme(legend.position= "none"), as_ggplot(get_legend(volcano_UCvsPSC))+theme(panel.background = element_rect(fill="white", color = NA))
                             ),
                             ggarrange(ncol=1, nrow=5,heights=c(1,1,1,1,0.2), GO_plot_UC_CON, GO_plot_PSCUCPSC_CON, GO_plot_UC_PSCUCPSC+theme(legend.position= "none"), as_ggplot(get_legend(GO_plot_UC_PSCUCPSC))+theme(panel.background = element_rect(fill="white", color = NA)), align ="v")
)+theme(panel.background = element_rect(fill="white", color = NA))


arranged_volcano
plot_width=16
#plot_height=18
plot_height=24
#plotname <- "output/VolcanoplotRNAseq"
plotname <- "output/figure2"



ggsave(arranged_volcano, filename = paste0(plotname, ".pdf"),device = "pdf",width = plot_width, height = plot_height)
ggsave(arranged_volcano, filename = paste0(plotname, ".png"),device = "png",width = plot_width, height = plot_height)
ggsave(arranged_volcano, filename = paste0(plotname, ".svg"),device = "svg",width = plot_width, height = plot_height)



#### FIGURE 4 (if you exclude graphical abstract, it's 3)#####

#umap_plot_marginal <- as.ggplot(ggMarginal(umap_plot, type="boxplot", groupColour = TRUE#, groupFill = TRUE))#+labs(tag="c)")
roclist <- list(
  Planell_CON_UC_Diagnose_results$pROC_object,
  Mo_CON_UC_Diagnose_results$pROC_object,
  Ostrowski_CON_UC_Diagnose_results$pROC_object,
  #Ostrowski_CON_PSC_Diagnose_results$pROC_object,
  
  
  #Our Dataset:
  CON_UC_Diagnose_results$pROC_object,
  
  PSC_UC_COHORT_results$pROC_object,
  PSC_UC_Diagnose_results$pROC_object,
  PSCUC_UC_Diagnose_results$pROC_object,
  PSC_PSCUC_Diagnose_results$pROC_object,
  PSC_CONTROL_Diagnose_results$pROC_object,
  PSCUCPSC_UC_Diagnose_results$pROC_object,
  PSCUCPSC_CON_Diagnose_results$pROC_object,
  PSCUC_CON_Diagnose_results$pROC_object)

roclist_names <- c(
  #"NULL",
  "Planell_CON_UC",
  "Mo_CON_UC",
  "Ostrowski_CON_UC",
#  "Ostrowski_CON_PSC",
  #Our Dataset:
  "CON_UC",
  "PSC_UC_COHORT",
  "PSC_UC",
  "PSCUC_UC",
  "PSC_PSCUC",
  "PSC_CONTROL",
  "PSCUCPSC_UC",
  "PSCUCPSC_CON",
  "PSCUC_CON")

roclist_UC_CON <- list(
  CON_UC_Diagnose_results$pROC_object,
  Planell_CON_UC_Diagnose_results$pROC_object,
  Mo_CON_UC_Diagnose_results$pROC_object,
  Ostrowski_CON_UC_Diagnose_results$pROC_object
  # Ostrowski_CON_PSC_Diagnose_results$pROC_object,
)

roclist_names_UC_CON <- c(
  #"NULL",
  "(i) This study: UC vs Control",
  "GSE94648: UC vs Control",
  "GSE112057: UC vs Control",
  "PRJEB28822: UC vs Control"
)




auroc_plot_PSC_UC_COHORT_controls <- plot_pROC_rocs(proclist = list(PSC_UC_COHORT_results$pROC_object), procnames = c("UC Control vs PSC Control"), plot_tag = "", plot_title=c("UC Control vs PSC Control"),setcolors = "#6FBBA1")#brewer.pal(12, "Set3")[1])

UC_CON_plot <- plot_pROC_rocs(proclist = roclist_UC_CON, procnames = roclist_names_UC_CON,plot_tag = element_blank(), plot_title="UC vs Control",setcolors = c("#BEBADA","#9B8EA9","#51497A","#D6D0D4"))#c("#bebada","#9894ae","#726f82","#4c4a57"))#c("darkblue","blue","lightblue","steelblue"))
UC_CON_plot

#brewer.pal(12, "Set3")[3]
PSC_CON_plot <- plot_pROC_rocs(proclist = list(PSC_CONTROL_Diagnose_results$pROC_object,
                                               PSCUC_CON_Diagnose_results$pROC_object,
                                               PSCUCPSC_CON_Diagnose_results$pROC_object
), 
procnames = c("PSC alone vs Control",
              "PSC/UC vs Control",
              "(ii)  PSC alone + PSC/UC vs Control"),
plot_tag = element_blank(), 
plot_title="PSC alone + PSC/UC vs Control",
setcolors = c("#D35A45", "#3B668B","#E0812E"))#brewer.pal(12, "Set3")[4:6])#c("red","gray", "black"))


#AUROC
auroc_plot_PSCPSCUC_UC <- plot_pROC_rocs(proclist = list(PSC_PSCUC_Diagnose_results$pROC_object,
                                                         PSCUCPSC_UC_Diagnose_results$pROC_object), 
                                         procnames = c("(iv) PSC alone vs PSC/UC",
                                                       "(iii) PSC alone + PSC/UC vs UC"),
                                         plot_tag = element_blank(), 
                                         plot_title="PSC alone vs PSC/UC and PSC alone + PSC/UC vs UC",
                                         setcolors = c("#5E912A","#D69FB3","#A6A6A6")#brewer.pal(12, "Set3")[7:9]#c("green","dimgray","black")
)


arranged_figure4 <- ggarrange(ncol = 1, nrow=2, align = "v",
                              ggarrange(ncol=2,nrow=1,align="v",auroc_plot_PSC_UC_COHORT_controls+ theme(legend.text=element_text(size=12))
                                        , UC_CON_plot+ theme(legend.text=element_text(size=12))
                                        ,labels= c("a)", "b)"),font.label = list(size = 24, color = "black", face = "plain", 
                                                                                 family = NULL),
                                        label.y = .97),
                              ggarrange(ncol=2,nrow=1,align="v",
                                        PSC_CON_plot+ theme(legend.text=element_text(size=12)),
                                        auroc_plot_PSCPSCUC_UC+ theme(legend.text=element_text(size=12)), 
                                        labels= c("c)", "d)"),font.label = list(size = 24, color = "black", face = "plain",
                                                                                family = NULL),
                                        label.y = .97),
                              # ggarrange(ncol=1,nrow=1,umap_grid_plot, labels= c("e)"),font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                              #           label.y = .97),
                              #umap_plot_marginal,
                              labels= c("","c)"),
                              label.y = .97,
                              font.label = list(size = 24, color = "black", face = "plain",
                                                family = NULL)
)+
  theme(panel.background = element_rect(fill="white", color = NA))

arranged_figure4


fig2_plot_width=12
fig2_plot_height=12
fig2_plotname <- "output/figure4"
ggsave(arranged_figure4, filename = paste0(fig2_plotname, ".pdf"),device = "pdf",width = fig2_plot_width, height = fig2_plot_height)
ggsave(arranged_figure4, filename = paste0(fig2_plotname, ".png"),device = "png",width = fig2_plot_width, height = fig2_plot_height)
ggsave(arranged_figure4, filename = paste0(fig2_plotname, ".svg"),device = "svg",width = fig2_plot_width, height = fig2_plot_height)



#### Supplementary FIGURE 5: Overrepresentation factors in DEG and coexpression modules ####
##### RF and DESeq enrichment in cemitool modules barplot #####
gene_background <- cem_merged_RF@selected_genes
DT_module<-as.data.table(cem_merged_RF@module)
UC_CON_RF_set <- mostimportantfeatures(variable_importance_table = CON_UC_Diagnose_results$variable_importance, limit = 50)$feature
PSCPSCUC_CON_RF_set <- mostimportantfeatures(variable_importance_table = PSCUCPSC_CON_Diagnose_results$variable_importance, limit = 50)$feature
PSCPSCUC_UC_RF_set <- mostimportantfeatures(variable_importance_table = PSCUCPSC_UC_Diagnose_results$variable_importance, limit = 50)$feature
UC_CON_DE_set <- sign_UC
PSCPSCUC_CON_DE_set <- sign_PSCUCPSC_CON
PSCPSCUC_UC_DE_set <- sign_UC_PSCUCPSC
model_names <- c("total","UC_CON_DE","PSC_CON_DE","PSC_UC_DE","UC_CON_RF","PSC_CON_RF","PSC_UC_RF")
model_gene_sets <- list(gene_background,UC_CON_DE_set,PSCPSCUC_CON_DE_set,PSCPSCUC_UC_DE_set,UC_CON_RF_set,PSCPSCUC_CON_RF_set,PSCPSCUC_UC_RF_set)
DT_module_sizes <- DT_module[,.N,by=modules]
DT_genesets_modules <- data.table(model_name=rep(model_names,length(unique(DT_module$modules))),
                                  module_name=rep(unique(DT_module$modules),each=length(model_names)),
                                  model_name_index=rep(1:length(model_names),length(unique(DT_module$modules))))
DT_genesets_modules[,geneset_size:=unlist(lapply(model_gene_sets[model_name_index],length))]
for (i in 1:nrow(DT_genesets_modules)) {
  set(DT_genesets_modules,i=i,j="module_size",value = DT_module_sizes[modules== DT_genesets_modules[i,module_name] ,N])
}
for (i in 1:nrow(DT_genesets_modules)) {
  set(DT_genesets_modules,i=i,j="intersect_size",value = length(intersect(x=DT_module[modules == DT_genesets_modules[i,module_name] ,genes], # genes in module
                                                                          y=model_gene_sets[[DT_genesets_modules[i,model_name_index]]])))
}
DT_genesets_modules[,percentage_of_geneset:=intersect_size/geneset_size]
DT_genesets_modules$model_name <- factor(DT_genesets_modules$model_name,levels = rev(model_names), labels = c("PSC alone + PSC/UC vs UC\ntop 50% importance","PSC alone + PSC/UC vs CON\ntop 50% importance","UC vs CON\ntop 50% importance",
                                                                                                              "PSC alone +PSC/UC vs UC\nDEG","PSC alone + PSC/UC vs CON\nDEG","UC vs CON\nDEG",
                                                                                                              "total"))
DT_genesets_modules$module_name <- factor(DT_genesets_modules$module_name,levels = rev(c(paste0("M",1:length(unique(DT_genesets_modules$module_name))-1),"Not.Correlated")))
DT_genesets_modules[,module_enrichment:= percentage_of_geneset/ .SD[model_name=="total",percentage_of_geneset]
                    ,by="module_name"]
DT_genesets_modules[,module_annotation:=ifelse(module_enrichment>2.5,as.character(module_name),NA_character_)]
DT_genesets_modules[,module_annotation_small:=ifelse(module_enrichment>2,as.character(module_name),NA_character_)]


ggplot(data = DT_genesets_modules) + 
  geom_col(mapping = aes(x=model_name,y=percentage_of_geneset,fill=module_name),position=position_dodge2()) + 
  coord_flip()
#
ggplot(data = DT_genesets_modules[module_name%in%c(paste0("M",1:9))]) + 
  geom_col(mapping = aes(x=module_name,y=module_enrichment,fill=model_name),position=position_dodge2()) + 
  coord_flip()
ggplot(data = DT_genesets_modules[module_name%in%c(paste0("M",10:18))]) + 
  geom_col(mapping = aes(x=module_name,y=module_enrichment,fill=model_name),position=position_dodge2()) + 
  coord_flip()


DT_genesets_modules$model_category <- ifelse(grepl("DEG",DT_genesets_modules$model_name), "DEG","RF")

gg_genesets_modules <- ggarrange(nrow=2,
  ggplot(data = DT_genesets_modules[model_category=="DEG"]) + 
    geom_col(mapping = aes(x=model_name,y=module_enrichment,fill=module_name),position=position_dodge2()) + 
    geom_label_repel(mapping = aes(
      x=model_name,y=module_enrichment  ,label=module_annotation,group=module_name),position = position_dodge(0.9),
      min.segment.length = 0.1,   size = 5,direction = 'x', ylim = c(30, 60) #direction = c("both")#label.size = 5.5,
      #nudge_y = 1
    ) +
    scale_fill_manual(values = rainbow(18)[c(18,1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12)]) +
    labs(x="gene set", y="overrepresentation factor")+
    theme_linedraw()+
    
    theme(legend.position = "none", 
          plot.margin = margin(5.5,20,5.5,27, "pt"),
          #legend.position="none",
          #axis.title.x = element_blank(),        
          legend.title = element_blank(),
          plot.tag.position = c(-.01, .97),
          #plot.tag = element_text(size = rel(1.5 * size.rel)),
          legend.text=element_text(size=rel(1.5 * size.rel)),
          axis.text = element_text(size=unit(17,"points")),
          axis.title = element_text(size=unit(17,"points")),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.tag =element_text(size=unit(22,"points")),
          plot.subtitle = element_text(size=unit(22,"points"))
          #axis.title.y = element_blank()
    )+coord_flip()+
    labs(tag="a)"),
  ggplot(data = DT_genesets_modules[model_category=="RF"][model_name!="total"]) + 
    geom_col(mapping = aes(x=model_name,y=module_enrichment,fill=module_name),position=position_dodge2()) + 
    geom_label_repel(mapping = aes(
      x=model_name,y=module_enrichment  ,label=module_annotation,group=module_name),position = position_dodge(0.9),
      min.segment.length = 0.1,   size = 5, hjust = 0,direction = 'x', ylim = c(5.5, 7.5) #direction = c("both")#label.size = 5.5,
      #nudge_y = 1
    ) + ylim(0,7)+# coord_cartesian(ylim = c(0, 60), expand = F) +
    #  scale_fill_manual(values = rainbow(22)[c(22,1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,17,21)]) +
    scale_fill_manual(values = rainbow(18)[c(18,1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12)]) +
    labs(x="gene set", y="Overrepresentation factor")+
    theme_linedraw()+
    
    theme(legend.position = "none", 
          plot.margin = margin(5.5,20,5.5,27, "pt"),
          #legend.position="none",
          #axis.title.x = element_blank(),        
          legend.title = element_blank(),
          plot.tag.position = c(-.01, .97),
          #plot.tag = element_text(size = rel(1.5 * size.rel)),
          legend.text=element_text(size=rel(1.5 * size.rel)),
          axis.text = element_text(size=unit(17,"points")),
          axis.title = element_text(size=unit(17,"points")),
          axis.title.y = element_blank(),
          plot.tag =element_text(size=unit(22,"points")),
          plot.subtitle = element_text(size=unit(22,"points"))
          #axis.title.y = element_blank()
    )+coord_flip()+
    labs(tag="b)")
)

# fwrite(DT_genesets_modules, file = "")

ggsave(gg_genesets_modules, filename = paste0("output/Supfig5_genesets", ".svg"),device = "svg",width = 10, height = 12)
ggsave(gg_genesets_modules, filename = paste0("output/Supfig5_genesets", ".pdf"),device = "pdf",width = 10, height = 12)
ggsave(gg_genesets_modules, filename = paste0("output/Supfig5_genesets", ".png"),device = "png",width = 10, height = 12)


svg(filename = "output/Supfig5_genesets.svg",height = 8,width=6)
gg_genesets_modules
dev.off()
png(filename = "output/Supfig5_genesets.png",height = 8,width=6,units = "in",res = 400)
gg_genesets_modules
dev.off()

gg_genesets_modules_zoom <- ggplot(data = DT_genesets_modules[model_name!="total"]) + 
  geom_col(mapping = aes(x=model_name,y=module_enrichment,fill=module_name),position=position_dodge2()) + 
  geom_text_repel(mapping = aes(
    x=model_name,y=ifelse(module_enrichment < 4.8,module_enrichment,4.8) ,label=module_annotation_small,group=module_name),position = position_dodge(0.9),min.segment.length = 0.1) +
  #  scale_fill_manual(values = rainbow(22)[c(22,1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,17,21)]) +
  scale_fill_manual(values = rainbow(18)[c(18,1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12)]) +
  theme_minimal() +
  #theme(legend.position="none") +
  labs(x="gene set", y="overrepresentation factor")+
  coord_flip(ylim = c(0,5)) 
svg(filename = "output/genesets_modules_zoom_signed.svg",height = 8,width=8)
gg_genesets_modules_zoom
dev.off()
png(filename = "output/genesets_modules_zoom_signed.png",height = 8,width=8,units = "in",res = 400)
gg_genesets_modules_zoom
dev.off()


##### RF and DESeq enrichment in cemitool modules barplot #####
gene_background <- cem_merged_RF@selected_genes
DT_module<-as.data.table(cem_merged_RF@module)
UC_CON_RF_set <- mostimportantfeatures(variable_importance_table = CON_UC_Diagnose_results$variable_importance, limit = 50)$feature
PSCPSCUC_CON_RF_set <- mostimportantfeatures(variable_importance_table = PSCUCPSC_CON_Diagnose_results$variable_importance, limit = 50)$feature
PSCPSCUC_UC_RF_set <- mostimportantfeatures(variable_importance_table = PSCUCPSC_UC_Diagnose_results$variable_importance, limit = 50)$feature
UC_CON_DE_set <- sign_UC
PSCPSCUC_CON_DE_set <- sign_PSCUCPSC_CON
PSCPSCUC_UC_DE_set <- sign_UC_PSCUCPSC
model_names <- c("total","UC_CON_DE","PSC_CON_DE","PSC_UC_DE","UC_CON_RF","PSC_CON_RF","PSC_UC_RF")
model_gene_sets <- list(gene_background,UC_CON_DE_set,PSCPSCUC_CON_DE_set,PSCPSCUC_UC_DE_set,UC_CON_RF_set,PSCPSCUC_CON_RF_set,PSCPSCUC_UC_RF_set)
DT_module_sizes <- DT_module[,.N,by=modules]
DT_genesets_modules <- data.table(model_name=rep(model_names,length(unique(DT_module$modules))),
                                  module_name=rep(unique(DT_module$modules),each=length(model_names)),
                                  model_name_index=rep(1:length(model_names),length(unique(DT_module$modules))))
DT_genesets_modules[,geneset_size:=unlist(lapply(model_gene_sets[model_name_index],length))]
for (i in 1:nrow(DT_genesets_modules)) {
  set(DT_genesets_modules,i=i,j="module_size",value = DT_module_sizes[modules== DT_genesets_modules[i,module_name] ,N])
}
for (i in 1:nrow(DT_genesets_modules)) {
  set(DT_genesets_modules,i=i,j="intersect_size",value = length(intersect(x=DT_module[modules == DT_genesets_modules[i,module_name] ,genes], # genes in module
                                                                          y=model_gene_sets[[DT_genesets_modules[i,model_name_index]]])))
}
DT_genesets_modules[,percentage_of_geneset:=intersect_size/geneset_size]
DT_genesets_modules$model_name <- factor(DT_genesets_modules$model_name,levels = rev(model_names), labels = c("PSC+PSCUC vs UC\ntop 50% importance","PSC+PSCUC vs CON\ntop 50% importance","UC vs CON\ntop 50% importance",
                                                                                                              "PSC+PSCUC vs UC\nDEG","PSC+PSCUC vs CON\nDEG","UC vs CON\nDEG",
                                                                                                              "total"))
DT_genesets_modules$module_name <- factor(DT_genesets_modules$module_name,levels = rev(c(paste0("M",1:length(unique(DT_genesets_modules$module_name))-1),"Not.Correlated")))
DT_genesets_modules[,module_enrichment:= percentage_of_geneset/ .SD[model_name=="total",percentage_of_geneset]
                    ,by="module_name"]
DT_genesets_modules[,module_annotation:=ifelse(module_enrichment>2.5,as.character(module_name),NA_character_)]
DT_genesets_modules[,module_annotation_small:=ifelse(module_enrichment>2,as.character(module_name),NA_character_)]

ggplot(data = DT_genesets_modules) + 
  geom_col(mapping = aes(x=model_name,y=percentage_of_geneset,fill=module_name),position=position_dodge2()) + 
  coord_flip()
#
ggplot(data = DT_genesets_modules[module_name%in%c(paste0("M",1:9))]) + 
  geom_col(mapping = aes(x=module_name,y=module_enrichment,fill=model_name),position=position_dodge2()) + 
  coord_flip()
ggplot(data = DT_genesets_modules[module_name%in%c(paste0("M",10:18))]) + 
  geom_col(mapping = aes(x=module_name,y=module_enrichment,fill=model_name),position=position_dodge2()) + 
  coord_flip()


gg_genesets_modules <- ggplot(data = DT_genesets_modules[model_name!="total"]) + 
  geom_col(mapping = aes(x=model_name,y=module_enrichment,fill=module_name),position=position_dodge2()) + 
  geom_text_repel(mapping = aes(
    x=model_name,y=module_enrichment  ,label=module_annotation,group=module_name),position = position_dodge(0.9),min.segment.length = 0.1) +
  #  scale_fill_manual(values = rainbow(22)[c(22,1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,17,21)]) +
  scale_fill_manual(values = rainbow(18)[c(18,1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12)]) +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x="gene set", y="overrepresentation factor")+
  coord_flip()
svg(filename = "output/Supfig5_genesets.svg",height = 8,width=6)
gg_genesets_modules
dev.off()
png(filename = "output/Supfig5_genesets.png",height = 8,width=6,units = "in",res = 400)
gg_genesets_modules
dev.off()

gg_genesets_modules_zoom <- ggplot(data = DT_genesets_modules[model_name!="total"]) + 
  geom_col(mapping = aes(x=model_name,y=module_enrichment,fill=module_name),position=position_dodge2()) + 
  geom_text_repel(mapping = aes(
    x=model_name,y=ifelse(module_enrichment < 4.8,module_enrichment,4.8) ,label=module_annotation_small,group=module_name),position = position_dodge(0.9),min.segment.length = 0.1) +
  #  scale_fill_manual(values = rainbow(22)[c(22,1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,17,21)]) +
  scale_fill_manual(values = rainbow(18)[c(18,1,7,13,2,8,14,3,9,15,4,10,16,5,11,17,6,12)]) +
  theme_minimal() +
  #theme(legend.position="none") +
  labs(x="gene set", y="overrepresentation factor")+
  coord_flip(ylim = c(0,5)) 
svg(filename = "output/genesets_modules_zoom_signed.svg",height = 8,width=8)
gg_genesets_modules_zoom
dev.off()
png(filename = "output/genesets_modules_zoom_signed.png",height = 8,width=8,units = "in",res = 400)
gg_genesets_modules_zoom
dev.off()
#heatmap nach module skaliert



# GRAPHICAL ABSTRACT ####

graph_abstract_figure2 <- ggarrange(ncol = 2, nrow=2,
                                   EnhancedVolcano::EnhancedVolcano(res_PSCandPSCUC_CON_untidy,
                                                                    title=element_blank(),
                                                                    subtitle="PSC alone + PSC/UC \nvs CON",
                                                                    caption="",
                                                                    selectLab = c(""),
                                                                    lab = paste0("italic('", rownames(res_PSCandPSCUC_CON_untidy), "')"),
                                                                    x = 'log2FoldChange',
                                                                    y = 'padj',
                                                                    xlim = c(min(res[['log2FoldChange']], na.rm = TRUE) - 0.1, max(res[['log2FoldChange']], na.rm = TRUE) +
                                                                               0.1),
                                                                    ylab = bquote(~- ~ Log[10] ~ italic(P- Adjusted)),
                                                                    #                  selectLab = selectLab_italics,
                                                                    xlab = bquote(~Log[2]~ 'fold change'),
                                                                    max.overlaps = 30,
                                                                    pCutoff = 1e-02,
                                                                    #                  FCcutoff = 1.0,
                                                                    pointSize = 3.0,
                                                                    labSize = 3,
                                                                    labCol = 'black',
                                                                    labFace = 'bold',
                                                                    legendLabels = c("NS", expression(Log[2] ~ FC), "p - Adjusted", expression(p - ~ Adjusted ~ and
                                                                                                                                               ~ log[2] ~ FC)), #c("NS", expression(Log[2] ~ FC), "p-value", paste0("P-Adjusted and ",expression(Log[2] ~ FC))),#expression(p - value ~ and ~ log[2] ~ FC)),
                                                                    boxedLabels = FALSE,
                                                                    parseLabels = TRUE,
                                                                    col = c('black', 'pink', 'purple', 'red3'),
                                                                    colAlpha = 4/5,
                                                                    legendPosition = 'right',
                                                                    legendLabSize = 14,
                                                                    legendIconSize = 4.0,
                                                                    drawConnectors = FALSE,
                                                                    widthConnectors = .50,
                                                                    lengthConnectors = unit(0.01, "npc"),
                                                                    colConnectors = 'black') +
                                     coord_flip() +theme_linedraw()+
                                     theme(plot.title = element_text(hjust = 0.0),
                                           legend.position= "none",
                                           plot.margin = margin(5.5,20,5.5,27, "pt"),
                                           #axis.title.x = element_blank(),        
                                           legend.title = element_blank(),
                                           plot.tag.position = c(-.01, .97),
                                           #plot.tag = element_text(size = rel(1.5 * size.rel)),
                                           legend.text=element_text(size=rel(1.5 * size.rel)),
                                           #axis.text = element_text(size=unit(17,"points")),
                                           axis.text = element_blank(),
                                           # axis.title = element_text(size=unit(17,"points")),
                                           axis.title = element_text(size=unit(18,"points")),
                                           plot.tag =element_text(size=unit(22,"points")),
                                           plot.subtitle = element_text(size=unit(22,"points"))
                                           #axis.title.y = element_blank()),
                                     ),
                                   EnhancedVolcano::EnhancedVolcano(res_PSCUC_PSC_untidy[!(rownames(res_PSCUC_PSC_untidy) %in% "MAN2A1"),],
                                                                    title=element_blank(),
                                                                    subtitle="PSC/UC vs \nPSC alone",
                                                                    caption="",
                                                                    selectLab = c(""),
                                                                    lab = paste0("italic('", rownames(res_PSCUC_PSC_untidy[!(rownames(res_PSCUC_PSC_untidy) %in% "MAN2A1"),]), "')"),
                                                                    x = 'log2FoldChange',
                                                                    y = 'padj',
                                                                    xlim = c(min(res[['log2FoldChange']], na.rm = TRUE) - 0.1, max(res[['log2FoldChange']], na.rm = TRUE) +
                                                                               0.1),
                                                                    ylab = bquote(~- ~ Log[10] ~ italic(P- Adjusted)),
                                                                    #                  selectLab = selectLab_italics,
                                                                    xlab = bquote(~Log[2]~ 'fold change'),
                                                                    max.overlaps = 30,
                                                                    pCutoff = 1e-02,
                                                                    #                  FCcutoff = 1.0,
                                                                    pointSize = 3.0,
                                                                    labSize = 3,
                                                                    labCol = 'black',
                                                                    labFace = 'bold',
                                                                    legendLabels = c("NS", expression(Log[2] ~ FC), "p - Adjusted", expression(p - ~ Adjusted ~ and
                                                                                                                                               ~ log[2] ~ FC)), #c("NS", expression(Log[2] ~ FC), "p-value", paste0("P-Adjusted and ",expression(Log[2] ~ FC))),#expression(p - value ~ and ~ log[2] ~ FC)),
                                                                    boxedLabels = FALSE,
                                                                    parseLabels = TRUE,
                                                                    col = c('black', 'pink', 'purple', 'red3'),
                                                                    colAlpha = 4/5,
                                                                    legendPosition = 'right',
                                                                    legendLabSize = 14,
                                                                    legendIconSize = 4.0,
                                                                    drawConnectors = FALSE,
                                                                    widthConnectors = .50,
                                                                    lengthConnectors = unit(0.01, "npc"),
                                                                    colConnectors = 'black') +
                                     coord_flip() +
                                     theme_linedraw()+
                                     theme(plot.title = element_text(hjust = 0.0),
                                           legend.position= "none",
                                           plot.margin = margin(5.5,20,5.5,27, "pt"),
                                           #axis.title.x = element_blank(),        
                                           legend.title = element_blank(),
                                           plot.tag.position = c(-.01, .97),
                                           #plot.tag = element_text(size = rel(1.5 * size.rel)),
                                           legend.text=element_text(size=rel(1.5 * size.rel)),
                                           #axis.text = element_text(size=unit(17,"points")),
                                           axis.text = element_blank(),
                                           # axis.title = element_text(size=unit(17,"points")),
                                           axis.title = element_text(size=unit(18,"points")),
                                           plot.tag =element_text(size=unit(22,"points")),
                                           plot.subtitle = element_text(size=unit(22,"points"))
                                           #axis.title.y = element_blank()),
                                     ),
                                   plot_pROC_rocs(
                                     proclist = list(#PSC_CONTROL_Diagnose_results$pROC_object
                                                     #PSCUC_CON_Diagnose_results$pROC_object,
                                                     PSCUCPSC_CON_Diagnose_results$pROC_object
                                     ), 
                                     procnames = c("PSC alone + PSC/UC\nvs Control"
                                                   #"PSC/UC vs Control",
                                                   #"(ii)  PSC + PSC/UC vs Control"
                                     ),
                                     plot_tag = element_blank(), 
                                     plot_title="",#PSC alone + PSC/UC \nvs Control",
                                     setcolors = c("#D35A45", "#3B668B","#E0812E")
                                   )+
                                     theme(plot.title = element_text(hjust = 0.0,
                                                                     size=unit(22,"points")),
                                           #legend.position= "none",
                                           plot.margin = margin(5.5,20,5.5,27, "pt"),
                                           #axis.title.x = element_blank(),        
                                           #legend.title = element_blank(),
                                           plot.tag.position = c(-.01, .97),
                                           #plot.tag = element_text(size = rel(1.5 * size.rel)),
                                           legend.text=element_text(size=rel(1.5 * size.rel)),
                                           #axis.text = element_text(size=unit(17,"points")),
                                           axis.text = element_blank(),
                                           # axis.title = element_text(size=unit(17,"points")),
                                           axis.title = element_text(size=unit(18,"points")),
                                           plot.tag =element_text(size=unit(22,"points")),
                                           plot.subtitle = element_text(size=unit(22,"points"))
                                           #axis.title.y = element_blank()),
                                     ),
                                   plot_pROC_rocs(
                                     proclist = list(PSC_PSCUC_Diagnose_results$pROC_object
                                                     #PSCUCPSC_UC_Diagnose_results$pROC_object
                                     ), 
                                     procnames = c("(iv) PSC alone vs PSC/UC"
                                                   #"(iii) PSC+PSC/UC vs UC"
                                     ),
                                     plot_tag = element_blank(), 
                                     plot_title="",#"PSC alone vs \nPSC/UC",
                                     setcolors = c("#5E912A","#D69FB3","#A6A6A6")
                                     #brewer.pal(12, "Set3")[7:9]#c("green","dimgray","black")
                                   )+#theme(plot.title = element_text(size=unit(22,"points")))
                                     theme(plot.title = element_text(hjust = 0.0,
                                                                     size=unit(22,"points")),
                                           #legend.position= "none",
                                           plot.margin = margin(5.5,20,5.5,27, "pt"),
                                           #axis.title.x = element_blank(),        
                                           #legend.title = element_blank(),
                                           plot.tag.position = c(-.01, .97),
                                           #plot.tag = element_text(size = rel(1.5 * size.rel)),
                                           legend.text=element_text(size=rel(1.5 * size.rel)),
                                           #axis.text = element_text(size=unit(17,"points")),
                                           axis.text = element_blank(),
                                           # axis.title = element_text(size=unit(17,"points")),
                                           axis.title = element_text(size=unit(18,"points")),
                                           plot.tag =element_text(size=unit(22,"points")),
                                           plot.subtitle = element_text(size=unit(22,"points"))
                                           #axis.title.y = element_blank()),
                                     )
                                   
)
graph_abstract_figure2
ggsave(graph_abstract_figure2, filename = paste0("output/graphical_abstract_figure_alt", ".svg"),device = "svg",width = 9, height = 9,bg="white")
ggsave(graph_abstract_figure2, filename = paste0("output/graphical_abstract_figure_alt", ".pdf"),device = "pdf",width = 9, height = 9,bg="white")




# Ostrowski_PSC Validation Preparation####
validation_ostrowski_PSC_cohortname<- "Ostrowski_PSC"
validation_cohorts <- list(Ostrowski ="UCAI Based\nDisease Severity",
                           Mo_UC="Diagnose",
                           Planell="Endoscopic\nMayo Score",
                           Ostrowski_PSC="Diagnose")
vsd_validation_ostrowski_PSC <- getCohortsVSD(cohortname=validation_ostrowski_PSC_cohortname)

# ostrowski_res <- getCohortsRES(cohortname = "Ostrowski_PSC")

severityscale <- validation_cohorts[[validation_ostrowski_PSC_cohortname]]
dds_validation_ostrowski_PSC <- getCohortsDDS(cohortname=validation_ostrowski_PSC_cohortname)
vsd_validation_ostrowski_PSC$severity <- vsd_validation_ostrowski_PSC$Diagnose
vsd_validation_ostrowski_PSC$severity[is.na(vsd_validation_ostrowski_PSC$severity)] <- "Control"  

# ggplot(as.data.table(t(assay(vsd_validation_ostrowski_PSC)[,])),aes(x=colData(vsd_validation_ostrowski_PSC)$severity,y=S100A6))+geom_boxplot()+geom_jitter(width=0.05)

# res_validation_ostrowski_PSC <- results(dds_validation_ostrowski_PSC, c("Diagnose","Control","UC"),tidy=TRUE)
# res_validation_ostrowski_PSC[res_validation_ostrowski_PSC$row %in% top25$feature,]

validation_ostrowski_PSC_controls <- ifelse(dds_validation_ostrowski_PSC$Diagnose == "Control",TRUE,FALSE)
validation_ostrowski_PSC_vst_counts <- data.table((assay(vsd_validation_ostrowski_PSC)[,]))

validation_ostrowski_PSC_z_con_stabilised <- data.table(rn=c(rownames(assay(vsd_validation_ostrowski_PSC))),
                                                        outlierfiltered_control_transformation(validation_ostrowski_PSC_vst_counts, controlsvector = validation_ostrowski_PSC_controls)
)
validation_ostrowski_PSC_z_con_stabilised <- transpose_datatable(validation_ostrowski_PSC_z_con_stabilised)

# Ostrowski PSC RF for validation_ostrowski_PSC Model####
Ostrowski_genefeatures <- c("rn",rownames(assay(vsd_validation_ostrowski_PSC)))
sum(colnames(filtered_merged_RF) %in% Ostrowski_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo

Ostrowski_PSC_dataset <- data.table(validation_ostrowski_PSC_z_con_stabilised, Diagnose=vsd_validation_ostrowski_PSC$Diagnose)


Ostrowski_PSC_dataset <- Ostrowski_PSC_dataset[Diagnose %in% c("PSC","Control"),]
Ostrowski_PSC_dataset$Diagnose <- droplevels(Ostrowski_PSC_dataset$Diagnose)

filtered_merged_PSCCON_RF <- merged_RF[merged_RF$Cohort == 0,]
#Diagnose factor: 0=Control, 1=PSC, 2=PSCUC, 3=UC

filtered_merged_PSCCON_RF$Diagnose <- droplevels(factor(filtered_merged_PSCCON_RF$Diagnose))
#0 = Control, 1= UC
levels(filtered_merged_PSCCON_RF$Diagnose) <- c(  "0","1","1")#0 = Control, 1= PSC or PSCUC


reduced4_merged_RF <- filtered_merged_PSCCON_RF[,intersect(colnames(Ostrowski_PSC_dataset), colnames(filtered_merged_PSCCON_RF)),with=F]
reduced3_Ostrowski_PSC_RF <- Ostrowski_PSC_dataset[,intersect(colnames(Ostrowski_PSC_dataset), colnames(filtered_merged_PSCCON_RF)),with=F]



levels(reduced3_Ostrowski_PSC_RF[["Diagnose"]]) <- c(0L,1L)

Ostrowski_CON_PSCPSCUC_Diagnose_results <- performML(dataset = reduced4_merged_RF[,-c("rn")],testing_dataset = reduced3_Ostrowski_PSC_RF,  splitfactor = 0.5, outcome_column = "Diagnose",seed = seednr)
reducedforOstrowski_CON_PSCPSCUC_Diagnose_results <- predictwithatunedmodel(tuned_model = Ostrowski_CON_PSCPSCUC_Diagnose_results$tuned_model, testing_dataset = Ostrowski_CON_PSCPSCUC_Diagnose_results$test_df, outcome_column = "Diagnose",seed = seednr)
#Auroc for external validation dataset
pROC::auc(Ostrowski_CON_PSCPSCUC_Diagnose_results$pROC_object)
pROC::ci.auc(Ostrowski_CON_PSCPSCUC_Diagnose_results$pROC_object)
# pROC::ci.auc(Ostrowski_CON_PSCPSCUC_Diagnose_results$pROC_object)
#validation with own test dataset
pROC::auc(reducedforOstrowski_CON_PSCPSCUC_Diagnose_results$pROC_object)


# Ostrowski UC RF for validation Model####
Ostrowski_UC_genefeatures <- c("rn",rownames(assay(vsd_validation)))
#sum(colnames(filtered_merged_RF) %in% Ostrowski_UC_genefeatures)
#we do not find all gene features from our dataset in the mo dataset, so we intersect, rerun RF on our merged dataset and then try to predict outcome in Mo
Ostrowski_UC_dataset <- data.table(validation_ostrowski_PSC_z_con_stabilised, Diagnose=vsd_validation_ostrowski_PSC$Diagnose)

Ostrowski_UC_dataset <- Ostrowski_UC_dataset[Diagnose %in% c("UC","Control"),]
Ostrowski_UC_dataset$Diagnose <- droplevels(Ostrowski_UC_dataset$Diagnose)


reduced_merged_RF <- filtered_merged_RF[,intersect(colnames(Ostrowski_UC_dataset), colnames(filtered_merged_RF)),with=F]
reduced_Ostrowski_UC_RF <- Ostrowski_UC_dataset[,intersect(colnames(Ostrowski_UC_dataset), colnames(filtered_merged_RF)),with=F]


Ostrowski_UC_CON_UC_Diagnose_results <- performML(dataset = reduced_merged_RF[,-c("rn")], testing_dataset = reduced_Ostrowski_UC_RF,  splitfactor = 0.5, outcome_column = "Diagnose", seed = seednr)

# pROC::auc(Ostrowski_UC_CON_UC_Diagnose_results$pROC_object)
# pROC::ci.auc(Ostrowski_UC_CON_UC_Diagnose_results$pROC_object)

reducedforOstrowski_CON_UC_Diagnose_results <- predictwithatunedmodel(tuned_model = Ostrowski_UC_CON_UC_Diagnose_results$tuned_model, testing_dataset = Ostrowski_UC_CON_UC_Diagnose_results$test_df, outcome_column = "Diagnose",seed = seednr)
#Auroc for external validation dataset
pROC::auc(Ostrowski_UC_CON_UC_Diagnose_results$pROC_object)
pROC::ci.auc(Ostrowski_UC_CON_UC_Diagnose_results$pROC_object)
# pROC::ci.auc(Ostrowski_CON_PSCPSCUC_Diagnose_results$pROC_object)
#validation with own test dataset
pROC::auc(reducedforOstrowski_CON_UC_Diagnose_results$pROC_object)

pROC::ci.auc(reducedforOstrowski_CON_UC_Diagnose_results$pROC_object)


# Supplement Figure 3 ####


# nes_matrix2 <- as.matrix(cem_merged_RF@enrichment$nes[,-1])
# rownames(nes_matrix2) <- cem_merged_RF@enrichment$nes$pathway 

nes_heatmap2 <- Heatmap(nes_matrix[,c("PSC","PSCUC","UC","Control")], 
                       name="NES",
                       heatmap_legend_param = list(
                         legend_direction = "horizontal", 
                         legend_width = unit(6, "cm")),
                       height = unit(10, "cm"),
                       cluster_columns = FALSE,
                       row_km = 5,
                       row_km_repeats = 10,
                       row_names_side = "left",
                       width = unit(14, "cm"),
                       column_names_rot = 0, 
                       row_title = NULL,
                       # column_names_max_height = unit(2, "cm"),
                       column_names_side = "bottom",
                       column_names_gp = gpar(fontsize = 10)
)

# 

nes_hmp <- draw(nes_heatmap2, heatmap_legend_side="top")
NES_grob <- grid.grabExpr(draw(nes_hmp))
#grob = grid.grabExpr(draw(Heatmap(...)))
df_nes_split <- stack(row_order(nes_hmp))

df_nes_split$modules <- rownames(nes_matrix[,c("PSC","PSCUC","UC","Control")])[df_nes_split$values]
rownames(df_nes_split) <- df_nes_split$modules
#df_nes_split[rownames(nes_matrix),"ind"]

DT_ora <- as.data.table(cem_merged_RF@ora)[p.adjust<0.05]
catchphrases <- c("neutrophil","T cell", "B cell", "killer cell", "platelet", "eosinophil", "macrophage", "dendritic", #cell types
                  "cell cycle","chemotaxis","cell death","defense","immune response","stress","signaling","splicing","metabolic","differentiation","migration", #cell states and functions
                  "RNA","mRNA","tRNA","rRNA","ncRNA","expression","ubiquitin","mitochond","ER","golgi","transcript","ribosom", #cell biology
                  "interleukin","interferon","NF-kappaB","wnt" #signaling componds, pathways
)
#modules <- c(paste0("M",1:21),"Not.Correlated")
modules <- c(paste0("M",1:18))#,"Not.Correlated")
DT_ora_catchphrases <- data.table(module=rep(modules,length(catchphrases)),catchphrase=rep(catchphrases,each=length(modules)))

for (i in 1:length(catchphrases)) {
  DT_ora_catchphrases[catchphrase==catchphrases[i],found:=case_when(catchphrase==catchphrases[i] & module %in% DT_ora[grepl(catchphrases[i],DT_ora$ID,ignore.case = T),Module] ~ 1,
                                                                    TRUE ~ 0)]
}
DT_ora_catchphrases[catchphrase=="ER",found := case_when(catchphrase=="ER" & module %in% DT_ora[grepl("ER",DT_ora$ID,ignore.case = F),Module] ~ 1,
                                                         TRUE ~ 0)]
DT_ora_catchphrases_wide <- pivot_wider(data=DT_ora_catchphrases,names_from = "catchphrase",values_from = "found") %>%as.data.table()
ora_matrix <- as.matrix(DT_ora_catchphrases_wide[,-1],dimnames=list(DT_ora_catchphrases_wide$module,colnames(DT_ora_catchphrases_wide)[-1]))
dimnames(ora_matrix) <- list(DT_ora_catchphrases_wide$module,colnames(DT_ora_catchphrases_wide)[-1])
ora_heatmap <- Heatmap(ora_matrix[df_nes_split$modules,], 
                       col=c("white","turquoise"),
                       name="ORA",
                       show_heatmap_legend = FALSE,
                       #heatmap_legend_param = list(
                       # legend_direction = "horizontal", 
                       #legend_width = unit(6, "cm")),
                       height = unit(10, "cm"),
                       cluster_rows = F,
                       clustering_distance_columns = "manhattan",
                       row_names_side = "left",
                       width = unit(14, "cm"),
                       column_names_rot = 90, 
                       border="grey70",
                       row_split = df_nes_split[#rownames(nes_matrix),
                         "ind"],
                       rect_gp = gpar(col = "grey70"),
                       # column_names_max_height = unit(2, "cm"),
                       #column_title = "GO term catchphrase",
                       column_title = element_blank(),
                       column_title_gp = gpar(fontsize=12),
                       column_title_side = "bottom",
                       column_names_side = "bottom",
                       column_names_gp = gpar(fontsize = 10),
                       row_title = NULL
                       # column_ti
)


# png("output/ORA_heatmap_signed.png",width=7.5,height=7.5,units="in",res=1200)
# draw(ora_heatmap)

GOcatchterm_grob <- grid.grabExpr(draw(ora_heatmap))#, heatmap_legend_side="top"))
# dev.off()
# DT_ora_catchphrases_wide[1,-1] %>% unlist()
# DT_ora_catchphrases_wide_pasted <-data.table(module=DT_ora_catchphrases_wide$module,catchphrase_concat=NA_character_)
# for (i in 1:nrow(DT_ora_catchphrases_wide_pasted)) {
#   set(DT_ora_catchphrases_wide_pasted,i=i,j="catchphrase_concat",value=paste(names(unlist(DT_ora_catchphrases_wide[i,-1]))[unlist(DT_ora_catchphrases_wide[i,-1])==1],collapse = ", "))
# }
# fwrite(DT_ora_catchphrases_wide_pasted,file="output/DT_ora_catchphrases_wide_pasted_signed.tsv",sep = "\t")

#c

# Cell Blueprints ####
#Cell Blueprints LM22 for every module created from cemitools:
enrichmentslist <- list()
for (x in unique(cem_merged_RF@module$modules)){
  print(x)
  genesinput <- data.table(cem_merged_RF@module)[modules == x,]$genes
  result <- blueprintenrichments(genes=genesinput)[[1]]
  # result <- result[result$enrichment >= 0.1,]
  enrichmentslist[[x]] <- result
}

enrichment_df <- data.table(bind_rows(enrichmentslist, .id = "column_label")) %>% pivot_wider(names_from = celltype, values_from = enrichment) %>% as.data.table() %>% as.matrix(., rownames="column_label")

htmp2 <- Heatmap(enrichment_df[df_nes_split$modules,], 
                name="Z-score",
                heatmap_legend_param = list(
                  legend_direction = "horizontal",
                  legend_width = unit(6, "cm")),
                # height = unit(10, "cm"),
                # # row_km = 1,
                # # row_km_repeats = 10,
                # row_names_side = "left",
                # width = unit(14, "cm"),
                # column_names_rot = 90, 
                # row_split = df_nes_split[#rownames(nes_matrix),
                #   "ind"],
                # cluster_rows = F,
                # # column_names_max_height = unit(2, "cm"),
                # column_names_side = "bottom",
                # column_names_gp = gpar(fontsize = 10),
                # row_title = NULL
                show_heatmap_legend = T,
                #heatmap_legend_param = list(
                # legend_direction = "horizontal", 
                #legend_width = unit(6, "cm")),
                height = unit(10, "cm"),
                cluster_rows = F,
                clustering_distance_columns = "manhattan",
                row_names_side = "left",
                width = unit(14, "cm"),
                column_names_rot = 90, 
                border="grey70",
                row_split = df_nes_split[#rownames(nes_matrix),
                  "ind"],
                rect_gp = gpar(col = "grey70"),
                # column_names_max_height = unit(2, "cm"),
                #column_title = "GO term catchphrase",
                column_title = element_blank(),
                column_title_gp = gpar(fontsize=12),
                column_title_side = "bottom",
                column_names_side = "bottom",
                column_names_gp = gpar(fontsize = 10),
                row_title = NULL
)
draw(htmp2, heatmap_legend_side="right")
CelltypeEnrichment_grob <- grid.grabExpr(draw(htmp2, heatmap_legend_side="top"))#, heatmap_legend_side="top"))


fig3_plot_width=20
fig3_plot_height=7
fig3_plotname <- "output/supplement_figure3_heatmaps"


figure3_supplement <- ggarrange(ncol = 2, nrow=1, align = "hv",
                                GOcatchterm_grob,
                                CelltypeEnrichment_grob,
                                labels= c("a)","b)"),
                                label.y = .97,
                                font.label = list(size = 24, color = "black", face = "bold", family = NULL))
ggsave(figure3_supplement, filename = paste0(fig3_plotname, ".pdf"),device = "pdf",width = fig3_plot_width-6.6, height = fig3_plot_height,bg="white")
ggsave(figure3_supplement, filename = paste0(fig3_plotname, ".png"),device = "png",width = fig3_plot_width-6.6, height = fig3_plot_height,bg="white")
ggsave(figure3_supplement, filename = paste0(fig3_plotname, ".svg"),device = "svg",width = fig3_plot_width-6.6, height = fig3_plot_height,bg="white")


##### Cibersort: rawcounts: run all samples of a cohort in one run for the batch effect correction to apply correctly #####

cibersortx_results_fractions_psc_cohort <- fread("../00_RawData/cibersort_results/psc_cohort_fractions_only/CIBERSORTx_Job32_Adjusted.txt")
cibersortx_results_fractions_psc_cohort$Cohort <- 0
cibersortx_results_fractions_psc_cohort[,Diagnose:=case_when(Mixture %in% merged_metadata[Diagnose=="PSCUC",SampleID] ~ "PSCUC",
                                                             Mixture %in% merged_metadata[Diagnose=="PSC",SampleID] ~ "PSC",
                                                             Mixture %in% merged_metadata[Diagnose=="Control",SampleID] ~ "Control")]
cibersortx_results_fractions_uc_cohort <- fread("../00_RawData/cibersort_results/uc_cohort_fractions_only/CIBERSORTx_Job33_Adjusted.txt")
cibersortx_results_fractions_uc_cohort$Cohort <- 1
cibersortx_results_fractions_uc_cohort[,Diagnose:=case_when(Mixture %in% merged_metadata[Diagnose=="UC",SampleID] ~ "UC",
                                                            Mixture %in% merged_metadata[Diagnose=="Control",SampleID] ~ "Control")]
cibersortx_results_cell_fractions_DT<-rbindlist(list(cibersortx_results_fractions_psc_cohort,cibersortx_results_fractions_uc_cohort))

cibersortx_results_cell_fractions_DT_long <- pivot_longer(cibersortx_results_cell_fractions_DT, cols = 2:23, values_to = "cell_fraction", names_to = "cell_type") %>% as.data.table()
cibersortx_results_cell_fractions_DT_long$Cohort <- cibersortx_results_cell_fractions_DT_long$Cohort %>% factor()
cibersortx_results_cell_fractions_DT_long[,CohortDiagnose:=factor(paste0(Cohort,Diagnose))]

# plotting
gg_cibersort_cell_fractions_PSC <- ggplot(data = cibersortx_results_cell_fractions_DT_long[Cohort==0],
                                          mapping = aes(x=Diagnose,y=cell_fraction,fill=Diagnose)) +
  geom_boxplot(position = position_dodge(0.8),coef=5) +
  scale_fill_manual(values = c("Control"="white","PSC"="#66C2A5","PSCUC"="#A6D854"))+
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank() ,legend.position=c(0.9,0.6), legend.text=element_text(size=12),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
        plot.margin = margin(0,0,0,0,"cm")) +
  stat_compare_means(comparisons = list(c("Control","PSC"),c("Control","PSCUC"),c("PSC","PSCUC")) ,label =  "p.signif", label.x = 1.5,  hide.ns = TRUE, method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"),na=FALSE),
                     step.increase = 0.05,tip.length = 0.01,vjust=0.65,size=3) +
  facet_grid(cols = vars(cell_type),scales = "free_x",switch="x") +
  labs(tag = "a") + ylab("cell fraction") + xlab("Diagnose")
gg_cibersort_cell_fractions_UC <- ggplot(data = cibersortx_results_cell_fractions_DT_long[Cohort==1],
                                         mapping = aes(x=Diagnose,y=cell_fraction,fill=Diagnose)) +
  geom_boxplot(position = position_dodge(0.8),coef=5) +
  scale_fill_manual(values = c("Control"="white","UC"="#8DA0CB"))+
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank() ,legend.position=c(0.9,0.6), legend.text=element_text(size=12),
        strip.text = element_text(angle = 90, hjust = 1, size=9,vjust = 0.5,margin = margin(0.1,0.5,0.0,0.5,"cm")),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
        plot.margin = margin(0,0,0,0,"cm")) +
  stat_compare_means(comparisons = list(c("Control","UC")) ,label =  "p.signif", label.x = 1.5,  hide.ns = TRUE, method = "wilcox.test",
                     step.increase = 0.05,tip.length = 0.01,vjust=0.65,size=3) +
  facet_grid(cols = vars(cell_type),scales = "free_x",switch = "x") +
  labs(tag = "b") + ylab("cell fraction") + xlab("Diagnose") 
ggarrange(plotlist = list(gg_cibersort_cell_fractions_PSC,gg_cibersort_cell_fractions_UC),nrow=2,align="none",heights = c(0.4,0.6))
svg(filename = "output/SupplementaryFigure6_cibersortx_fractions_by_cohort.svg",width = 8,height=8)
ggarrange(plotlist = list(gg_cibersort_cell_fractions_PSC,gg_cibersort_cell_fractions_UC),nrow=2,align="none",heights = c(0.4,0.6))
dev.off()
# manually remove the "NS." string from the plots, since the statistical test sometimes gives NA due to all-zero vectors compared against each other 




# sessionInfo()
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.2 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=de_DE.UTF-8
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#   [1] org.Hs.eg.db_3.14.0         geneplotter_1.76.0          annotate_1.76.0             XML_3.99-0.14               lattice_0.20-45             genefilter_1.80.3           topGO_2.50.0
# [8] SparseM_1.81                GO.db_3.14.0                AnnotationDbi_1.56.2        Rgraphviz_2.42.0            graph_1.76.0                EnhancedVolcano_1.16.0      ggrepel_0.9.3
# [15] rmarkdown_2.20              CEMiTool_1.22.0             ggpubr_0.6.0                ComplexHeatmap_2.14.0       InformationValue_1.2.3      Information_0.0.9           ranger_0.14.1
# [22] RColorBrewer_1.1-3          DESeq2_1.38.3               SummarizedExperiment_1.28.0 Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_0.63.0          GenomicRanges_1.50.2
# [29] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0         yardstick_1.1.0             workflowsets_1.0.0          workflows_1.1.3
# [36] tune_1.0.1                  rsample_1.1.1               recipes_1.0.5               parsnip_1.0.4               modeldata_1.1.0             infer_1.0.4                 dials_1.1.0
# [43] scales_1.2.1                broom_1.0.4                 tidymodels_1.0.0            lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0               dplyr_1.1.1
# [50] purrr_1.0.1                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1                ggplot2_3.4.1               tidyverse_2.0.0             plyr_1.8.8
# [57] data.table_1.14.8
# 
# loaded via a namespace (and not attached):
#   [1] ggthemes_4.2.4         coda_0.19-4            ragg_1.2.5             ggpmisc_0.5.2          bit64_4.0.5            knitr_1.42             DelayedArray_0.24.0    rpart_4.1.19
# [9] KEGGREST_1.38.0        hardhat_1.2.0          RCurl_1.98-1.10        doParallel_1.0.17      generics_0.1.3         GPfit_1.0-8            preprocessCore_1.60.2  cowplot_1.1.1
# [17] RSQLite_2.3.0          shadowtext_0.1.2       future_1.32.0          bit_4.0.5              tzdb_0.3.0             enrichplot_1.18.3      ggpp_0.5.1             viridis_0.6.2
# [25] gower_1.0.1            xfun_0.38              hms_1.1.3              evaluate_0.20          fansi_1.0.4            igraph_1.4.1           DBI_1.1.3              htmlwidgets_1.6.2
# [33] backports_1.4.1        vctrs_0.6.1            quantreg_5.94          abind_1.4-5            cachem_1.0.7           withr_2.5.0            ggforce_0.4.1          HDO.db_0.99.1
# [41] checkmate_2.1.0        sna_2.7-1              treeio_1.22.0          svglite_2.1.1          cluster_2.1.4          DOSE_3.24.2            ape_5.7-1              lazyeval_0.2.2
# [49] crayon_1.5.2           labeling_0.4.2         pkgconfig_2.0.3        tweenr_2.0.2           nlme_3.1-162           nnet_7.3-18            rlang_1.1.0            globals_0.16.2
# [57] lifecycle_1.0.3        MatrixModels_0.5-1     downloader_0.4         polyclip_1.10-4        Matrix_1.5-1           aplot_0.1.10           carData_3.0-5          base64enc_0.1-3
# [65] GlobalOptions_0.1.2    png_0.1-8              viridisLite_0.4.1      rjson_0.2.21           bitops_1.0-7           gson_0.1.0             pROC_1.18.0            Biostrings_2.66.0
# [73] blob_1.2.4             shape_1.4.6            qvalue_2.30.0          parallelly_1.35.0      rstatix_0.7.2          gridGraphics_0.5-1     ggsignif_0.6.4         memoise_2.0.1
# [81] magrittr_2.0.3         zlibbioc_1.44.0        compiler_4.2.3         scatterpie_0.1.8       intergraph_2.0-2       clue_0.3-64            cli_3.6.1              XVector_0.38.0
# [89] DiceDesign_1.9         listenv_0.9.0          patchwork_1.1.2        htmlTable_2.4.1        Formula_1.2-5          MASS_7.3-58.3          WGCNA_1.72-1           tidyselect_1.2.0
# [97] stringi_1.7.12         textshaping_0.3.6      GOSemSim_2.24.0        locfit_1.5-9.7         fastmatch_1.1-3        tools_4.2.3            timechange_0.2.0       future.apply_1.10.0
# [105] parallel_4.2.3         circlize_0.4.15        rstudioapi_0.14        foreach_1.5.2          foreign_0.8-82         gridExtra_2.3          prodlim_2019.11.13     farver_2.1.1
# [113] ggraph_2.1.0           digest_0.6.31          pracma_2.4.2           lava_1.7.2.1           Rcpp_1.0.10            car_3.1-1              httr_1.4.5             ggdendro_0.1.23
# [121] colorspace_2.1-0       splines_4.2.3          yulab.utils_0.0.6      tidytree_0.4.2         graphlayouts_0.8.4     ggplotify_0.1.0        systemfonts_1.0.4      xtable_1.8-4
# [129] jsonlite_1.8.4         ggtree_3.6.2           dynamicTreeCut_1.63-1  tidygraph_1.2.3        timeDate_4022.108      ggfun_0.0.9            ipred_0.9-14           R6_2.5.1
# [137] Hmisc_5.0-1            lhs_1.1.6              pillar_1.9.0           htmltools_0.5.5        glue_1.6.2             fastmap_1.1.1          clusterProfiler_4.6.2  DT_0.27
# [145] BiocParallel_1.32.6    class_7.3-21           codetools_0.2-19       fgsea_1.24.0           furrr_0.3.1            utf8_1.2.3             network_1.18.1         survival_3.5-3
# [153] statnet.common_4.8.0   munsell_0.5.0          GetoptLong_1.0.5       fastcluster_1.2.3      GenomeInfoDbData_1.2.9 iterators_1.0.14       impute_1.72.3          reshape2_1.4.4
# [161] gtable_0.3.3