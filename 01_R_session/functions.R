require(data.table)
require(tidyverse)
require(DESeq2)
# require(rvcheck)
require(CEMiTool)
# require(M3C)
require(rmarkdown)
# require(RColorBrewer)
# require(gridExtra)
# require(ggpubr)
# require(ggplotify)
# require(VennDiagram)
# require(extrafont)

# FUNCTIONS ####
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

#Wrapping CEMiTool Analysis with Report-Output:
setGeneric('generate_report', function(cem, ...) {
  standardGeneric('generate_report')
})

#' @rdname generate_report
setMethod('generate_report', signature('CEMiTool'),
          function(cem, max_rows_ora=50, title="Report", directory=output_location, force=FALSE, nameofcohort=cohortname, ...) {
            if(is.null(unique(cem@module$modules))){
              stop("No modules in CEMiTool object! Did you run find_modules()?")
            }
            if(dir.exists(directory)){
              if(!force){
                stop("Stopping analysis: ", directory, " already exists! Use force=TRUE to overwrite.")
              }
            }else{
              dir.create(directory, recursive=TRUE)
            }
            rmd <- system.file("report", "report.Rmd", package = "CEMiTool")
            rmarkdown::render(rmd,output_file = nameofcohort, output_dir=directory, intermediates_dir=directory, quiet=TRUE, ...)
          })

CEMiwrapper <- function(expressionmatrix=assay(vsd), ID=metadata$SampleID, Groups=metadata$selectedgenesm3c, gmt_location=gmt_file, reportname="report", beta20=FALSE, applyfiltering=TRUE, ...){
  if(isFALSE(beta20)){
    foreach::registerDoSEQ()
    sample_annot <- data.frame(ID, Groups)
    colnames(sample_annot)<-  c("SampleName", "Class")
    CountsforCemi <- data.frame(expressionmatrix[,ID], row.names = rownames(expressionmatrix))
    #Geneset enrichment analysis of CEMI:
    gmt_in <- CEMiTool::read_gmt(gmt_location)
    gmt_in <- gmt_in[!gmt_in$gene=="",]
    #Adding Interactions: 
    #int_df <- fread(interaction_location)
    #class(int_df) <- "data.frame"    
    #All together:
    cem <- CEMiTool::cemitool(CountsforCemi, sample_annot, gmt_in, 
                              filter=applyfiltering, plot=TRUE, verbose=TRUE,gsea_max_size = 10000)}
  else{
    foreach::registerDoSEQ()
    sample_annot <- data.frame(ID, Groups)
    colnames(sample_annot)<-  c("SampleName", "Class")
    CountsforCemi <- data.frame(expressionmatrix[,ID], row.names = rownames(expressionmatrix))
    #Geneset enrichment analysis of CEMI:
    gmt_in <- CEMiTool::read_gmt(gmt_location)
    gmt_in <- gmt_in[!gmt_in$gene=="",]
    
    #All together:
    cem <- CEMiTool::cemitool(CountsforCemi, sample_annot, gmt_in, 
                              filter=TRUE, plot=TRUE, verbose=TRUE,set_beta = 20,gsea_max_size = 10000)
  }
  generate_report(cem, force=T, title = cohortname, nameofcohort=reportname,...)
  return(cem)
}

# return a heatmap for the expression of the input-genes in defined cell types 
blueprintplot <- function(genes=NULL){
  require(pheatmap)
  require(ggplotify)
  bprint_expression <- fread("../data_rescources/E-MTAB-3827-query-results.tpms.tsv")
  bprint_expression_l<-pivot_longer(bprint_expression, 3:29, names_to = "cell_type", values_to = "TPM")%>%as.data.table()
  setnames(bprint_expression_l,c("gene_ID","gene_name","cell_type","TPM"))
  celltypes_to_keep<-c("mature neutrophil", "erythroblast", "conventional dendritic cell", "macrophage", "CD4-positive, alpha-beta T cell",
                       "CD8-positive, alpha-beta T cell", "CD38-negative naive B cell","CD8-positive", "CD14-positive, CD16-negative classical monocyte")
  
  length_genes_blueprint_intersect<-length(intersect(genes,bprint_expression$`Gene Name`))
  genes_blueprint_intersect<-intersect(genes,bprint_expression$`Gene Name`)
  bprint_expression_mat<-as.matrix(bprint_expression[`Gene Name`%in%genes,-2:-1],rownames = bprint_expression[`Gene Name`%in%genes,`Gene Name`])
  bprint_expression_mat[(bprint_expression_mat)==0]<-min(bprint_expression_mat,na.rm = T)
  #make it a z score of the log values
  z_bprint_expression_mat<-scale(t(log10(bprint_expression_mat)))
  #kick NAs
  z_bprint_expression_mat[is.na(z_bprint_expression_mat)]<-min(z_bprint_expression_mat,na.rm = T)
  z_bprint_expression_mat_reduced<-z_bprint_expression_mat#[rownames(z_bprint_expression_mat)%in%celltypes_to_keep,]
  rownames(z_bprint_expression_mat_reduced) <- gsub("-negative","-",gsub("-positive","+",rownames(z_bprint_expression_mat_reduced)))
  
  blueprint <- as.ggplot(pheatmap(t(z_bprint_expression_mat_reduced),cluster_rows = F,cluster_cols = F,fontsize=7,border_color=NA,legend=F,show_rownames=F,main="genes",breaks=seq(from=-3,to=3,length.out = 101),silent=F))#,filename="marker_expression_BLUEPRINT_celltypes.png", width = 7, height =7)
  return(blueprint)
  }

# #return a data table with z-scores for enrichment of specific cell type of given input-genes
# blueprintenrichments <- function(genes=NULL){
#   bprint_expression <- fread("../data_rescources/E-MTAB-3827-query-results.tpms.tsv")
#   bprint_expression_l<-pivot_longer(bprint_expression, 3:29, names_to = "cell_type", values_to = "TPM")%>%as.data.table()
#   setnames(bprint_expression_l,c("gene_ID","gene_name","cell_type","TPM"))
#   celltypes_to_keep<-c("mature neutrophil", "erythroblast", "conventional dendritic cell", "macrophage", "CD4-positive, alpha-beta T cell",
#                        "CD8-positive, alpha-beta T cell", "CD38-negative naive B cell","CD8-positive", "CD14-positive, CD16-negative classical monocyte")
#   
#   length_genes_blueprint_intersect<-length(intersect(genes,bprint_expression$`Gene Name`))
#   genes_blueprint_intersect<-intersect(genes,bprint_expression$`Gene Name`)
#   bprint_expression_mat<-as.matrix(bprint_expression[`Gene Name`%in%genes,-2:-1],rownames = bprint_expression[`Gene Name`%in%genes,`Gene Name`])
#   bprint_expression_mat[(bprint_expression_mat)==0]<-min(bprint_expression_mat,na.rm = T)
#   #make it a z score of the log values
#   z_bprint_expression_mat<-scale(t(log10(bprint_expression_mat)))
#   #kick NAs
#   z_bprint_expression_mat[is.na(z_bprint_expression_mat)]<-min(z_bprint_expression_mat,na.rm = T)
#   z_bprint_expression_mat_reduced<-z_bprint_expression_mat#[rownames(z_bprint_expression_mat)%in%celltypes_to_keep,]
#   rownames(z_bprint_expression_mat_reduced) <- gsub("-negative","-",gsub("-positive","+",rownames(z_bprint_expression_mat_reduced)))
#   
#   enrichmentstats <- apply(z_bprint_expression_mat_reduced,1,function(x) sum(x)/length(x))
#   enrichment_df <- data.frame(celltype=names(enrichmentstats),enrichment=enrichmentstats)
#   # blueprint <- as.ggplot(pheatmap(t(z_bprint_expression_mat_reduced),cluster_rows = F,cluster_cols = F,fontsize=7,border_color=NA,legend=F,show_rownames=F,main="genes",breaks=seq(from=-3,to=3,length.out = 101),silent=F))#,filename="marker_expression_BLUEPRINT_celltypes.png", width = 7, height =7)
#   return(enrichment_df)
#   }
#return a list: First element: data frame with average z-scores for enrichment of specific cell type of given input-genes
# Second element: matrix with z-scores for enrichment of specific cell type of given input-genes
blueprintenrichments <- function(genes=NULL){
  bprint_expression <- fread("../data_rescources/LM22_source_GEPs.txt")
  #bprint_expression <- fread("../data_rescources/E-MTAB-3827-query-results.tpms.tsv")
  bprint_expression$genesinput <- gsub("-",".",bprint_expression$genesinput,fixed = TRUE)
  bprint_expression_l<-pivot_longer(bprint_expression, 2:23, names_to = "cell_type", values_to = "expression")%>%as.data.table()
  #bprint_expression_l<-pivot_longer(bprint_expression, 3:29, names_to = "cell_type", values_to = "TPM")%>%as.data.table()
  setnames(bprint_expression_l,c("gene_name","cell_type","expression"))
  #setnames(bprint_expression_l,c("gene_ID","gene_name","cell_type","expression"))
  celltypes_to_keep<-bprint_expression_l$cell_type %>% unique()
  #celltypes_to_keep<-c("mature neutrophil", "erythroblast", "conventional dendritic cell", "macrophage", "CD4-positive, alpha-beta T cell",
  #                     "CD8-positive, alpha-beta T cell", "CD38-negative naive B cell","CD8-positive", "CD14-positive, CD16-negative classical monocyte")
  length_genes_blueprint_intersect<-length(intersect(genes,bprint_expression$genesinput))
  #length_genes_blueprint_intersect<-length(intersect(genes,bprint_expression$`Gene Name`))
  genes_blueprint_intersect<-intersect(genes,bprint_expression$genesinput)
  #genes_blueprint_intersect<-intersect(genes,bprint_expression$`Gene Name`)
  bprint_expression_mat<-as.matrix(bprint_expression[genesinput%in%genes,-1],rownames = bprint_expression[genesinput%in%genes,genesinput])
  #bprint_expression_mat<-as.matrix(bprint_expression[`Gene Name`%in%genes,-2:-1],rownames = bprint_expression[`Gene Name`%in%genes,`Gene Name`])
  bprint_expression_mat[(bprint_expression_mat)==0]<-min(bprint_expression_mat,na.rm = T)
  #make it a z score of the log values
  z_bprint_expression_mat<-scale(t(log10(bprint_expression_mat)))
  #kick NAs
  z_bprint_expression_mat[is.na(z_bprint_expression_mat)]<-min(z_bprint_expression_mat,na.rm = T)
  z_bprint_expression_mat_reduced<-z_bprint_expression_mat#[rownames(z_bprint_expression_mat)%in%celltypes_to_keep,]
  rownames(z_bprint_expression_mat_reduced) <- gsub("-negative","-",gsub("-positive","+",rownames(z_bprint_expression_mat_reduced)))
  
  enrichmentstats <- apply(z_bprint_expression_mat_reduced,1,function(x) sum(x)/length(x))
  enrichment_df <- data.frame(celltype=names(enrichmentstats),enrichment=enrichmentstats)
  # blueprint <- as.ggplot(pheatmap(t(z_bprint_expression_mat_reduced),cluster_rows = F,cluster_cols = F,fontsize=7,border_color=NA,legend=F,show_rownames=F,main="genes",breaks=seq(from=-3,to=3,length.out = 101),silent=F))#,filename="marker_expression_BLUEPRINT_celltypes.png", width = 7, height =7)
  return(list(enrichment_df,z_bprint_expression_mat_reduced))
}

#perform Random Forest in a tidymodels workflow with tuning of parameters. Takes some time. Added multithreading for some performance improvements. 
#Needs the training result of training() from the rsample library and the column name of the outcome variable.
return_tuned_RF_model <- function(training_df=NULL,outcome_col=NULL){
  require(doParallel)
  require(tidyverse)
  require(tidymodels)
  require(tune)
  
  #Registering all cores: 
  all_cores <- parallel::detectCores(logical = FALSE)
  
  cl <- makePSOCKcluster(all_cores)
  registerDoParallel(cl)
  
  # Models
  
  model_rf <- 
    rand_forest(mtry = tune::tune(), trees = tune::tune(), min_n = tune::tune()) %>% 
    set_engine("ranger", importance = "impurity") %>% 
    set_mode("classification")
  
  # model_xgboost <- 
  #   boost_tree(mtry = tune::tune(), trees = tune::tune(), min_n = tune::tune()) %>% 
  #   set_engine("xgboost", importance = "impurity") %>% 
  #   set_mode("classification")
  
  # Grid of hyperparameters
  
  grid_rf <- 
    grid_max_entropy(
      mtry(range = c(1, 20)), 
      trees(range = c(500, 1000)),
      min_n(range = c(2, 10)),
      size = 30)
  
  # Formula
  formula_rf <- paste0(outcome_col, " ~ .") %>% as.formula()
  
  # Workflow
  
  wkfl_rf <- 
    workflow() %>% 
    add_formula(formula_rf) %>%
    add_model(model_rf)
  
  # wkfl_wgboost <-
  #   workflow() %>%
  #   add_formula(Diagnose ~ .) %>%
  #   add_model(model_xgboost)
  
  # Cross validation method
  
  cv_folds <- vfold_cv(training_df, v = 5)
  cv_folds
  
  # Choose metrics
  
  my_metrics <- metric_set(yardstick::roc_auc, yardstick::accuracy, yardstick::sens, yardstick::spec, yardstick::precision, yardstick::recall)
  
  # Tuning
  
  rf_fit <- tune::tune_grid(
    wkfl_rf,
    resamples = cv_folds,
    grid = grid_rf,
    metrics = my_metrics,
    control = control_grid(verbose = TRUE) # don't save prediction (imho)
  )
  
  # Fit best model 
  
  tuned_model <-
    wkfl_rf %>% 
    finalize_workflow(select_best(rf_fit, metric = "accuracy")) %>% 
    fit(data = training_df)
  
  #unregister all cores
  stopCluster(cl)
  return(tuned_model)
}

# 
# dataset = reduced2_merged_RF[,-c("rn")]
# testing_dataset = reduced2_Planell_RF
# splitfactor = 0.5
# outcome_column = "Diagnose"
# seed = seednr
#make sure outcome_col is a binary factor, there is no sanity test yet in place for this.
performML <- function(dataset=NULL,testing_dataset=NULL,splitfactor=NULL,outcome_column=NULL,seed=123){
  require(rsample)
  require(pROC)
  require(InformationValue)
    
  dataset[[outcome_column]] <- base::droplevels(dataset[[outcome_column]])
  levels(dataset[[outcome_column]]) <- c(0,1)
  set.seed(seednr)
  splitted <- rsample::initial_split(dataset, splitfactor)
  rftrain <- training(splitted)
  rftest <- testing(splitted)
  if(is.null(testing_dataset)){testing_dataset <- rftest}
  Resultitem <- outcome_column
  
  

  tuned_model <- return_tuned_RF_model(training_df=rftrain,outcome_col=Resultitem)
  
  
  # predicted_train <- predict(tuned_model, rftrain, type = "prob")
  predicted_test <- predict(tuned_model, testing_dataset, type = "prob")
  predicted <- as.data.table(predicted_test)[,2]%>%unlist()%>%as.vector()%>%as.numeric()
  
  realvalue <- as.integer(testing_dataset[[Resultitem]])-1
  
  optCutOff <- InformationValue::optimalCutoff(realvalue, predicted, optimiseFor = "Both")[1]
  
  misClassErrorResult <- misClassError(testing_dataset[[Resultitem]], predicted)#, threshold = optCutOff)
  
  # ROCplotResult <- plotROC(rftest[[Resultitem]], predicted, returnSensitivityMat=T)
  proc_predicted_test <- predict(tuned_model, testing_dataset, type= "prob")[,2] %>% unlist() %>% as.vector()
  proc_predicted_test_integer <- as.integer(ifelse(proc_predicted_test>=optCutOff,1,0))
  
  realvalue <- as.integer(testing_dataset[[Resultitem]])-1
  rocobject <- pROC::roc(realvalue ~ proc_predicted_test)#proc_predicted_test_integer
  
  concordanceResult <- Concordance(testing_dataset[[Resultitem]], predicted)
  sensitivityResult <- InformationValue::sensitivity(testing_dataset[[Resultitem]], predicted, threshold = optCutOff)
  specificityResult <- InformationValue::specificity(testing_dataset[[Resultitem]], predicted, threshold = optCutOff)
  confusionMatrixResult <- confusionMatrix(testing_dataset[[Resultitem]], predicted, threshold = optCutOff)
  # confusionMatrixResult
  variable_importance <- data.table(feature=names(extract_fit_parsnip(tuned_model)$fit$variable.importance), importance=extract_fit_parsnip(tuned_model)$fit$variable.importance)
  variable_importance <- variable_importance[order(variable_importance$importance,decreasing = T),]
  # top25 <- head(variable_importance[order(variable_importance$importance, decreasing=T),],25)
  return_list <- list(rocobject,tuned_model,concordanceResult,sensitivityResult,specificityResult,confusionMatrixResult,variable_importance, rftrain, rftest)
  names(return_list) <- c("pROC_object","tuned_model","concordanceResult","sensitivityResult","specificityResult","confusionMatrixResult","variable_importance","train_df","test_df")
  return(return_list)
}

predictwithatunedmodel <- function(tuned_model=NULL,testing_dataset=NULL,outcome_column = NULL,seed = 123){
  require(pROC)
  Resultitem <- outcome_column
  predicted_test <- predict(tuned_model, testing_dataset, type = "prob")
  predicted <- as.data.table(predicted_test)[,2]%>%unlist()%>%as.vector()%>%as.numeric()
  
  optCutOff <- optimalCutoff(testing_dataset[[Resultitem]], predicted, optimiseFor = "Both")[1]
  
  misClassErrorResult <- misClassError(testing_dataset[[Resultitem]], predicted)#, threshold = optCutOff)
  
  # ROCplotResult <- plotROC(rftest[[Resultitem]], predicted, returnSensitivityMat=T)
  proc_predicted_test <- predict(tuned_model, testing_dataset, type=c("prob"))[,1] %>% unlist() %>% as.vector() #%>% as.integer()
  realvalue <- as.integer(testing_dataset[[Resultitem]])-1
  rocobject <- pROC::roc(realvalue ~ proc_predicted_test)
  
  concordanceResult <- Concordance(testing_dataset[[Resultitem]], predicted)
  sensitivityResult <- InformationValue::sensitivity(testing_dataset[[Resultitem]], predicted, threshold = optCutOff)
  specificityResult <- InformationValue::specificity(testing_dataset[[Resultitem]], predicted, threshold = optCutOff)
  confusionMatrixResult <- confusionMatrix(testing_dataset[[Resultitem]], predicted, threshold = optCutOff)
  # confusionMatrixResult
  variable_importance <- data.table(feature=names(extract_fit_parsnip(tuned_model)$fit$variable.importance), importance=extract_fit_parsnip(tuned_model)$fit$variable.importance)
  variable_importance <- variable_importance[order(variable_importance$importance,decreasing = T),]
  # top25 <- head(variable_importance[order(variable_importance$importance, decreasing=T),],25)
  return_list <- list(rocobject,tuned_model,concordanceResult,sensitivityResult,specificityResult,confusionMatrixResult,variable_importance)
  names(return_list) <- c("pROC_object","tuned_model","concordanceResult","sensitivityResult","specificityResult","confusionMatrixResult","variable_importance")
  return(return_list)
}

topGO_enrichment <- function(genelist = NULL, 
                             weights = NULL,
                             weight_threshold = 0,
                             statistical_test="Fisher", 
                             algorithm_topGO = "classic",
                             DESeq_result = NULL,
                             match_by_expression = FALSE,
                             gene_background = NULL,
                             ontology_type = "BP",
                             draw_plot = FALSE,
                             output_dir = "output/", 
                             debug_mode = FALSE) {
  # this function is purposed to take a list of genes, which is then analysed by the TOPGO algorithm using GO libraries
  # desired output is a tabular enrichment result
  # also, if the user wishes, a pdf containing a graph of the enrichment hierarchy is written
  # genelist  character vector coontaining gene names. If gene_background is given, genes not in genelist will be ignored.
  # If DESeq_result is given, genes not occuring in the DESeq-assay table will be ignored.
  # weigths are assumed to increasing with importance: Gene A: 5,  Gene B: 2.3: Gene A is more important and will be ranked higher
  # by the Kolmogorov-Smirnov test.
  require(Rgraphviz)
  if (debug_mode == TRUE) {
    #this is to intialize for a test run
    genelist <- c("S100A12","ETS1","ACTB")
    weights <- c(2,1,0.8)
    gene_background <- c("VSTM1","TARM1","OSCAR","NDUFA3","TFPT","PRPF31","AC012314.8","CNOT3","LENG1","TMC4","MBOAT7","TSEN34","S100A12","ETS1","ACTB")
  }
  # checks
  stopifnot( "no genelist provided" = !is.null(genelist) )
  stopifnot( "weights and genelist length differ" = is.null(weights) | length(weights) == length(genelist) )
  if (is.null(DESeq_result) & is.null(gene_background)) {
    stop("TopGO needs a gene universe (background). Provide either DESeq_result or gene_background.")
  }
  if ( !is.null(DESeq_result) & !is.null(gene_background) ) {
    message("DESeq_result and gene_background are both input. Using gene_background...")
  }
  stopifnot("library missing" = require(DESeq2) & require(topGO) & require(genefilter) & require(geneplotter) & require(tidyverse) & require(data.table) & require(org.Hs.eg.db))
  
  
  # to generate the gene background, prefer gene_background if provided
  # else if no expression matching is wished, use any expressed gene from DESeq2 result
  # else, calculate a background based on similarly expressed genes
  if ( !is.null(gene_background)) {
    gene_background <- gene_background
    backG <- NULL
  } else {
    
    if (match_by_expression == FALSE) {
      DESeq_result <- subset(DESeq_result, DESeq_result$baseMean > 0) 
      gene_background <- DESeq_result$row
      backG <- NULL
    } else {
      DESeq_result <- subset(DESeq_result, DESeq_result$baseMean > 0) 
      DESeq_result <- subset(DESeq_result, !is.na(DESeq_result$row) & !DESeq_result$row == "" )
      rownames(DESeq_result) <- DESeq_result$row
      gene_background <- rownames(DESeq_result)
      overallBaseMean <- as.matrix(DESeq_result[, "baseMean", drop = F])
      rownames(overallBaseMean) <- rownames(DESeq_result)
      sig_idx <- match(genelist, rownames(overallBaseMean))
      
      backG <- c()
      
      for(i in sig_idx){
        ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
        backG <- c(backG, ind)
        
      }
      
      backG <- unique(backG)
      backG <- rownames(overallBaseMean)[backG]
      backG <- setdiff(backG, genelist)
      
      multidensity( list( 
        all= log2(overallBaseMean[,"baseMean"]) ,
        foreground =log2(overallBaseMean[genelist, "baseMean"]), 
        background =log2(overallBaseMean[backG, "baseMean"])), 
        xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")
      
    }
    
  }
  
  
  if (!is.null(backG)) {
    inUniverse <- gene_background %in% c(genelist, backG)
    inSelection <- gene_background %in% genelist 
  } else {
    inUniverse <- gene_background %in% gene_background
    inSelection <- gene_background %in% genelist
  }
  # create a named vector with either weights or factor values
  
  if ( !is.null(weights) ) {
    alg <- as.numeric( inSelection[inUniverse] )
    names(alg) <- gene_background[inUniverse]
    alg[genelist] <- weights
  } else {
    alg <- factor( as.integer( inSelection[inUniverse] ) )
    names(alg) <- gene_background[inUniverse]
  }
  
  # function to pick genes which are used when weights are given. 
  
  geneSelFunction <- function(allScore) {
    return(allScore > weight_threshold)
  }
  ## prepare topGO object
  tgd <- new( "topGOdata", ontology = ontology_type, allGenes = alg, 
              geneSel = ifelse( !is.null(weights), geneSelFunction, NULL),
              nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "symbol" )
  
  ## run tests
  #resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO <- runTest(tgd, algorithm = algorithm_topGO, statistic = statistical_test , scoreOrder = "decreasing") #scoreOrder is handdown for the ks test
  
  ## look at results
  if(length(nodes(graph(tgd))) < 200){
    # tab <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
    #                  Fisher.classic = resultTopGO.classic,
    #                  orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
    tab <- GenTable( tgd, 
                     test.pvalue = resultTopGO,
                     orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
  }else{
    # tab <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
    #                  Fisher.classic = resultTopGO.classic,
    #                  orderBy = "Fisher.elim" , topNodes = 200)
    tab <- GenTable( tgd, 
                     test.pvalue = resultTopGO,
                     orderBy = "Fisher.elim" , topNodes = 200)
  }
  
  
  #printGraph(tgd, resultTopGO.elim, firstSigNodes = 5, fn.prefix = output_dir, useInfo = "all", pdfSW = TRUE)
  if(draw_plot == TRUE) 
  { printGraph(tgd, resultTopGO, firstSigNodes = 15, fn.prefix = output_dir, useInfo = "all", pdfSW = TRUE) }
  
  
  return(tab)
}


EnhancedVolcano_function <- function(res=NULL,plottitle=element_blank(),plotnumber=element_blank()){
  require(EnhancedVolcano)
  lab_italics <- paste0("italic('", rownames(res), "')")
  #selectLab_italics = paste0(
  #  "italic('",
  #  c('VCAM1','KCTD12','ADAM12', 'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA'),
  #  "')")
  EV <- EnhancedVolcano::EnhancedVolcano(res,
                                         title=plotnumber,
                                         subtitle=plottitle,
                                         caption="",
                                         lab = lab_italics,
                                         x = 'log2FoldChange',
                                         y = 'pvalue',
                                         #                  selectLab = selectLab_italics,
                                         xlab = bquote(~Log[2]~ 'fold change'),
                                         max.overlaps = 30,
                                         #                  pCutoff = 10e-14,
                                         #                  FCcutoff = 1.0,
                                         pointSize = 3.0,
                                         labSize = 3.0,
                                         labCol = 'black',
                                         labFace = 'bold',
                                         #legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", paste0( expression(~ p - value ~ and),"\n",expression(~ log[2] ~ FC))),
                                         boxedLabels = FALSE,
                                         parseLabels = TRUE,
                                         col = c('black', 'pink', 'purple', 'red3'),
                                         colAlpha = 4/5,
                                         legendPosition = 'right',
                                         legendLabSize = 14,
                                         legendIconSize = 4.0,
                                         drawConnectors = TRUE,
                                         widthConnectors = .50,
                                         lengthConnectors = unit(0.01, "npc"),
                                         colConnectors = 'black') +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.0))
  return(EV)
}

plot_pROC_rocs <- function(proclist=list(), procnames=NULL, size.rel=1, plot_tag="a)", plot_title=element_blank()){
  require(RColorBrewer)
  require(ggplot2)
  require(pROC)
  # Define the number of colors we need
  nb.cols <- length(proclist)
  mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
  
  # Plot rocs via ggroc
  rocs_plotted <- ggroc(proclist, size=1.5)+
    labs(title=plot_title, tag=plot_tag)+
    scale_colour_manual(values = mycolors, label= procnames) +
    theme(#legend.position = c(0.70, 0.19),
      # legend.position = c(1., 0.),
      # legend.justification = c(0, 1),
      legend.position = c(0.98,0.02), # bottom right position with 0.02 margin
      legend.justification = c(1, 0), # bottom right justification
      legend.box.margin = margin(1,1,1,1, unit = "mm"), # small margin
      legend.box.background = element_rect(colour = "black"),
      panel.background = element_rect(fill = "white", 
                                      colour = NA), panel.border = element_rect(fill = NA, 
                                                                                colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
      panel.grid.minor = element_line(size = rel(0.5)), 
      #panel.grid.minor = element_blank(),
      #panel.grid.major = element_blank(),
      #legend.title = element_text(size = rel(1.8 * size.rel)),
      #legend.spacing = unit(0.1, "cm"),
      #legend.box.margin = margin(0,0,0,0, "pt"),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.8 * size.rel)), 
      axis.title = element_text(size = rel(1.5 * size.rel)), 
      axis.text = element_text(size = rel(1.2 * size.rel)), 
      strip.text.x = element_text(size = rel(1.8 * size.rel)), 
      strip.text.y = element_text(size = rel(1.8 * size.rel)), 
      legend.key.height = unit(1.3 * size.rel,"line"), 
      legend.key.width = unit(1.3 * size.rel, "line"),
      plot.tag.position = c(-.01, 1.05),
      plot.tag = element_text(size = rel(3 * size.rel)),
      plot.margin = margin(40,20,5.5,27, "pt"),
      strip.background = element_rect(fill = "grey85", 
                                      colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                    colour = NA), complete = TRUE) +
    labs(x = "False positive rate (1-specificity)",
         y = "Detection rate (sensitivity)") +
    scale_x_reverse(#limits = c(0, 100),
      expand = expansion(mult = c(.004, .004)),
      labels = rev(c( "0%", "25%", "50%", "75%", "100%"))) +
    scale_y_continuous(#limits = c(0, 100),
      expand = expansion(mult = c(0.004, .004)),
      labels = (c( "0%", "25%", "50%", "75%", "100%")))+
    guides(#colour=guide_legend(ncol=2), 
      linetype = "none")+
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1, linetype="NULL_line"),color="grey85")+
    scale_linetype_manual(name = "NULL_line_name", values = c("NULL_line" = 1))+
    coord_fixed()
  return(rocs_plotted)
}

mostimportantfeatures <- function(variable_importance_table = NULL, limit= 50){
  i=as.double(0)
  j=1L
  feature_df <- data.table(feature=c(),importance=c())
  normed_importance <- (100/sum(variable_importance_table$importance))*variable_importance_table$importance
  while (i < limit){
    feature_df <- rbind(feature_df, variable_importance_table[j,])
    i = i + as.double(normed_importance[j])
    j = j + 1
        #print(i)
  }
  return(feature_df)
}
