require(data.table)
require(tidyverse)
require(DESeq2)
# require(pheatmap)
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
# require(SCORPIUS)
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
                              filter=applyfiltering, plot=TRUE, verbose=TRUE)}
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
                              filter=TRUE, plot=TRUE, verbose=TRUE,set_beta = 20)
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

#return a data table with z-scores for enrichment of specific cell type of given input-genes
blueprintenrichments <- function(genes=NULL){
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
  
  enrichmentstats <- apply(z_bprint_expression_mat_reduced,1,function(x) sum(x)/length(x))
  enrichment_df <- data.frame(celltype=names(enrichmentstats),enrichment=enrichmentstats)
  # blueprint <- as.ggplot(pheatmap(t(z_bprint_expression_mat_reduced),cluster_rows = F,cluster_cols = F,fontsize=7,border_color=NA,legend=F,show_rownames=F,main="genes",breaks=seq(from=-3,to=3,length.out = 101),silent=F))#,filename="marker_expression_BLUEPRINT_celltypes.png", width = 7, height =7)
  return(enrichment_df)
  }

#perform Random Forest in a tidymodels workflow with tuning of parameters. Takes some time. Added multithreading for some performance improvements. 
#Needs the training result of training() from the rsample library and the column name of the outcome variable.
return_tuned_RF_model <- function(training_df=NULL,outcome_col=NULL){
  require(doParallel)
  require(tidyverse)
  require(tidymodels)
  
  #Registering all cores: 
  all_cores <- parallel::detectCores(logical = FALSE)
  
  cl <- makePSOCKcluster(all_cores)
  registerDoParallel(cl)
  
  # Models
  
  model_rf <- 
    rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>% 
    set_engine("ranger", importance = "impurity") %>% 
    set_mode("classification")
  
  # model_xgboost <- 
  #   boost_tree(mtry = tune(), trees = tune(), min_n = tune()) %>% 
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
  
  my_metrics <- metric_set(roc_auc, accuracy, sens, spec, precision, recall)
  
  # Tuning
  
  rf_fit <- tune_grid(
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


# dataset = merged_RF[merged_RF$Diagnose == 0,-c("rn","Diagnose")]
# splitfactor = 0.5
# outcome_column = "Cohort"
# seed = seednr
#make sure outcome_col is a binary factor, there is no sanity test yet in place for this.
performML <- function(dataset=NULL,testing_dataset=NULL,splitfactor=NULL,outcome_column=NULL,seed=123){
  require(rsample)
  require(pROC)
  
  dataset[[outcome_column]] <- droplevels(dataset[[outcome_column]])
  levels(dataset[[outcome_column]]) <- c(0,1
  )
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
  
  optCutOff <- optimalCutoff(testing_dataset[[Resultitem]], predicted, optimiseFor = "Both")[1]
  
  misClassErrorResult <- misClassError(testing_dataset[[Resultitem]], predicted)#, threshold = optCutOff)
  
  # ROCplotResult <- plotROC(rftest[[Resultitem]], predicted, returnSensitivityMat=T)
  proc_predicted_test <- predict(tuned_model, testing_dataset) %>% unlist() %>% as.vector() %>% as.integer()
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
