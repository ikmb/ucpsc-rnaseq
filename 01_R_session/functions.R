require(data.table)
require(tidyverse)
require(DESeq2)
require(pheatmap)
require(rvcheck)
require(CEMiTool)
require(M3C)
require(rmarkdown)
require(RColorBrewer)
require(gridExtra)
require(ggpubr)
require(ggplotify)
require(VennDiagram)
require(extrafont)
require(SCORPIUS)
#########FUNCTIONS##########
getnormalizedDESeq2object <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname=NULL){
  require(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
 
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file=paste0(projectdir,"/00_RawData/",counttable))
  
  colnames(UC) <- gsub("_1Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("\\-L1_S[0-9]+_L00[0-9]_R1_001$","",colnames(UC))
  
  #countmatrix:
  if(colnames(UC)[1] %in% "Geneid"){UC <- UC[,-1]}
  countmatrix <- data.matrix(UC[,-1])
  mode(countmatrix) <- "double"
  rownames(countmatrix) <- unlist(UC[,1])
  rownames(countmatrix) <- make.names(rownames(countmatrix), unique = TRUE)
  rownames(countmatrix) <- gsub("[\\-]",".",rownames(countmatrix), perl=T)
  rm(UC)
  
  #readin metadata  
  metadata <- fread(paste0(projectdir,"/00_RawData/",metatable), header=metatableheader)
  #Diagnose must be "UC" or "Control for every sample
  #colnames(metadata) <- c("SampleID", "Diagnose", colnames(metadata)[-1:-2])
  
  #globin genes removal:
  countmatrix <- countmatrix[!(row.names(countmatrix) %in%  c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")),]
  
  #all samples are in metadata and countmatrix and same order
  metadata <- metadata[metadata$SampleID %in% colnames(countmatrix),]
  
  countmatrix <- countmatrix[,colnames(countmatrix) %in% metadata$SampleID]
  countmatrix <- countmatrix[,metadata$SampleID]
  #filtering to 10 counts per gene in at least 10% of samples
  countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
  
  metadata$PlateNr <- factor(metadata$PlateNr)
  metadata$Diagnose <- factor(metadata$Diagnose)
  
  require(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~PlateNr + Diagnose)
  dds <- DESeq(dds)
  #DESeqperformed <- TRUE
  #vsde <- vst(dds, blind=FALSE)
  #res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
  vsd <- vst(dds, blind=FALSE)
  return(vsd)
}

getrawcounts <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname="Our"){
  require(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file=paste0(projectdir,"/00_RawData/",counttable))
  
  colnames(UC) <- gsub("_1Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("\\-L1_S[0-9]+_L00[0-9]_R1_001$","",colnames(UC))
  #countmatrix:
  if(colnames(UC)[1] %in% "Geneid"){UC <- UC[,-1]}
  countmatrix <- data.matrix(UC[,-1])
  mode(countmatrix) <- "double"
  rownames(countmatrix) <- unlist(UC[,1])
  rownames(countmatrix) <- make.names(rownames(countmatrix), unique = TRUE)
  rownames(countmatrix) <- gsub("[\\-]",".",rownames(countmatrix), perl=T)
  rm(UC)
  
  #readin metadata  
  metadata <- fread(paste0(projectdir,"/00_RawData/",metatable), header=metatableheader)
  #Diagnose must be "UC" or "Control for every sample
  colnames(metadata) <- c("SampleID", "Diagnose", colnames(metadata)[-1:-2])
  
  #globin genes removal:
  countmatrix <- countmatrix[!(row.names(countmatrix) %in%  c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")),]
  
  #all samples are in metadata and countmatrix and same order
  metadata <- metadata[metadata$SampleID %in% colnames(countmatrix),]
  
  countmatrix <- countmatrix[,colnames(countmatrix) %in% metadata$SampleID]
  countmatrix <- countmatrix[,metadata$SampleID]
  #filtering to 10 counts per gene in at least 10% of samples
  countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]

  return(countmatrix)
}

getmetadata <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname="Our"){
  require(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)

  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file=paste0(projectdir,"/00_RawData/",counttable))
  
  colnames(UC) <- gsub("_1Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("\\-L1_S[0-9]+_L00[0-9]_R1_001$","",colnames(UC))
  #countmatrix:
  if(colnames(UC)[1] %in% "Geneid"){UC <- UC[,-1]}
  countmatrix <- data.matrix(UC[,-1])
  mode(countmatrix) <- "double"
  rownames(countmatrix) <- unlist(UC[,1])
  rownames(countmatrix) <- make.names(rownames(countmatrix), unique = TRUE)
  rownames(countmatrix) <- gsub("[\\-]",".",rownames(countmatrix), perl=T)
  rm(UC)
  
  #readin metadata  
  metadata <- fread(paste0(projectdir,"/00_RawData/",metatable), header=metatableheader)
  #Diagnose must be "UC" or "Control for every sample
  #colnames(metadata) <- c("SampleID", "Diagnose", colnames(metadata)[-1:-2])
  
  #globin genes removal:
  # countmatrix <- countmatrix[!(row.names(countmatrix) %in%  c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")),]
  
  #all samples are in metadata and countmatrix and same order
  metadata <- metadata[metadata$SampleID %in% colnames(countmatrix),]
  metadata$PlateNr <- as.factor(metadata$PlateNr)
  return(metadata)
}

vsdofallknownseverity <- function(meta=metadata, countmatrix=rawcounts, batch="PlateNr"){
  #throwing out all Unknown severity patients:
  meta <- meta[!(meta$severity %in% "Unknown"),]
  countmatrix <- countmatrix[,meta$SampleID]
  cohortname <- "Discovery"
  colnames(meta)[colnames(meta)==batch] <- "batch"
  #filtering to 10 counts per gene in at least 10% of samples
  countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
  #deseq2 of mild/con, severe/con:
  meta$batch <- as.factor(meta$batch)
  dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=meta, design= ~batch + severity)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind=FALSE)
  return(vsd)
  }

ddsofallknownseverity <- function(meta=metadata, countmatrix=rawcounts, batch="PlateNr",...){
  #throwing out all Unknown severity patients:
  meta <- meta[!(meta$severity %in% "Unknown"),]
  countmatrix <- countmatrix[,meta$SampleID]
  cohortname <- "Discovery"
  colnames(meta)[colnames(meta)==batch] <- "batch"
  meta$batch <- as.factor(meta$batch)
  #filtering to 10 counts per gene in at least 10% of samples
  countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
  #deseq2 of mild/con, severe/con:
  dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=meta, design= ~batch + severity)
  dds <- DESeq(dds)
  return(dds)
  }


getCohortsVSD <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname=NULL){
  require(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  #cohorttable <- fread("/home/sukmb465/Documents/Eike/Nextflow/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
  #cohortname <- "Our"
  ##cohortname <- dlg_list(cohorttable[,Cohort], multiple = FALSE)$res
  
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file= paste0(projectdir,"/00_RawData/",counttable))
 
  colnames(UC) <- gsub("_1Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(UC))
  
  #countmatrix:
  if(colnames(UC)[1] %in% "Geneid"){UC <- UC[,-1]}
  countmatrix <- data.matrix(UC[,-1])
  mode(countmatrix) <- "double"
  rownames(countmatrix) <- unlist(UC[,1])
  rownames(countmatrix) <- make.names(rownames(countmatrix), unique = TRUE)
  rownames(countmatrix) <- gsub("[\\-]",".",rownames(countmatrix), perl=T)
  rm(UC)
  
  #readin metadata  
  metadata <- fread(paste0(projectdir,"/00_RawData/",metatable), header=metatableheader)
  #Diagnose must be "UC" or "Control for every sample
  colnames(metadata) <- c("SampleID", "Diagnose", colnames(metadata)[-1:-2])
  
  #globin genes removal:
  countmatrix <- countmatrix[!(row.names(countmatrix) %in%  c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")),]
  
  #all samples are in metadata and countmatrix and same order
  metadata <- metadata[metadata$SampleID %in% colnames(countmatrix),]
  
  countmatrix <- countmatrix[,colnames(countmatrix) %in% metadata$SampleID]
  countmatrix <- countmatrix[,metadata$SampleID]
  #filtering to 10 counts per gene in at least 10% of samples
  countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
  require(DESeq2)
  #deseq2
  #check for countmatrix with only integer values, then check if countmatrix is already log-normalised:
  if (all(countmatrix%%1==0) == TRUE) {
    print("Values are integer")
    #filtering for at least 10 counts in at least 10% of all samples:######
    countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
    dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~Diagnose)
    dds <- DESeq(dds)
    vsd <- vst(dds, blind=FALSE)
    res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
    DESeqperformed<- TRUE
  }else{
    print("Values are double values, checking for normalisation...")
    if(max(countmatrix)>100){
      countmatrix <- round(countmatrix)
      #filtering for at least 10 counts in at least 10% of all samples:######
      countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
      dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~Diagnose)
      dds <- DESeq(dds)
      vsd <- vst(dds, blind=FALSE)
      res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
      DESeqperformed<- TRUE
    }else{
      dds <- DESeqDataSetFromMatrix(countData=round(countmatrix), colData=metadata, design= ~Diagnose)
      #dds <- DESeq(dds)
      print("DESeq cannot be performed with normalised Data")
      #vsd <- vst(dds, blind=FALSE)
      vsd <- dds
      assay(vsd) <- countmatrix
      #res <- results(dds,  contrast=c("Diagnose","UC","Control"))
      DESeqperformed<- FALSE
    }
  }

  return(vsd)
}

getCohortsDDS <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname=NULL){
  require(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file= paste0(projectdir,"/00_RawData/",counttable))
  
  colnames(UC) <- gsub("_1Aligned.sortedByCoord.out.bam","",colnames(UC))
  colnames(UC) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(UC))
  
  #countmatrix:
  if(colnames(UC)[1] %in% "Geneid"){UC <- UC[,-1]}
  countmatrix <- data.matrix(UC[,-1])
  mode(countmatrix) <- "double"
  rownames(countmatrix) <- unlist(UC[,1])
  rownames(countmatrix) <- make.names(rownames(countmatrix), unique = TRUE)
  rownames(countmatrix) <- gsub("[\\-]",".",rownames(countmatrix), perl=T)
  rm(UC)
  
  #readin metadata  
  metadata <- fread(paste0(projectdir,"/00_RawData/",metatable), header=metatableheader)
  #Diagnose must be "UC" or "Control for every sample
  colnames(metadata) <- c("SampleID", "Diagnose", colnames(metadata)[-1:-2])
  
  #globin genes removal:
  countmatrix <- countmatrix[!(row.names(countmatrix) %in%  c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")),]
  
  #all samples are in metadata and countmatrix and same order
  metadata <- metadata[metadata$SampleID %in% colnames(countmatrix),]
  
  countmatrix <- countmatrix[,colnames(countmatrix) %in% metadata$SampleID]
  countmatrix <- countmatrix[,metadata$SampleID]
  #filtering to 10 counts per gene in at least 10% of samples
  countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
  require(DESeq2)
  #deseq2
  #check for countmatrix with only integer values, then check if countmatrix is already log-normalised:
  if (all(countmatrix%%1==0) == TRUE) {
    print("Values are integer")
    #filtering for at least 10 counts in at least 10% of all samples:######
    countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
    dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~Diagnose)
    dds <- DESeq(dds)
    vsd <- vst(dds, blind=FALSE)
    res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
    DESeqperformed<- TRUE
  }else{
    print("Values are double values, checking for normalisation...")
    if(max(countmatrix)>100){
      countmatrix <- round(countmatrix)
      #filtering for at least 10 counts in at least 10% of all samples:######
      countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
      dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~Diagnose)
      dds <- DESeq(dds)
      vsd <- vst(dds, blind=FALSE)
      res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
      DESeqperformed<- TRUE
    }else{
      dds <- DESeqDataSetFromMatrix(countData=round(countmatrix), colData=metadata, design= ~Diagnose)
      #dds <- DESeq(dds)
      print("DESeq cannot be performed with normalised Data")
      #vsd <- vst(dds, blind=FALSE)
      vsd <- dds
      assay(vsd) <- countmatrix
      #res <- results(dds,  contrast=c("Diagnose","UC","Control"))
      DESeqperformed<- FALSE
    }
  }
  return(dds)
}

remove_outlier_filter <- function(x=NULL){
  calcedSD <- sd(x)
  calcedmean <- mean(x)
  booleanvector <- (x > (calcedmean - 3*calcedSD)) & (x < (calcedmean + 3*calcedSD))
  return(booleanvector)
}

# Substract the mean of all controls from any row and then divide by the SD of outlierfiltered controls. Then return a datatable
outlierfiltered_control_transformation <- function(dataset=NULL,controlsvector=NULL){
  cols <-  colnames(dataset)
  transformed_dataset <- apply(dataset,1,function(x){
                              y <- x[controlsvector];
                              z <- y[remove_outlier_filter(y)];
                              (x-mean(y))/sd(z)}
                              )
  transformed_dataset <- t(transformed_dataset)
  transformed_dataset <- data.table(transformed_dataset, keep.rownames=TRUE)
  colnames(transformed_dataset) <- cols
  return(transformed_dataset)
  }

#see https://stackoverflow.com/questions/28653867/best-way-to-transpose-data-table
transpose_datatable<- function(datatable){
  tdatatable <- t(datatable[,-1,with=F])
  colnames(tdatatable) <- datatable[[1]]
  tdatatable <- data.table(tdatatable, keep.rownames=T)
  setnames(tdatatable, 1, names(datatable)[1])
  return(tdatatable)
}

LEGACYRandomForestFeatureSelection <- function(table=NULL,set_seed=123, splittestratio=0.5,featuretablesize=50){
  require(data.table)
  require(rsample)      # data splitting 
  require(randomForest) # basic implementation
  require(tidyverse)
  require(ranger)
  require(Information)
  require(InformationValue)
  
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