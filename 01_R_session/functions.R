library(data.table)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(CEMiTool)
library(M3C)
library(rmarkdown)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggplotify)
library(VennDiagram)
library(extrafont)
library(SCORPIUS)
#########FUNCTIONS##########
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

CEMiwrapper <- function(expressionmatrix=assay(vsd), ID=metadata$SampleID, Groups=metadata$selectedgenesm3c, gmt_location=gmt_file, reportname="report", beta20=FALSE, ...){
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
                              filter=TRUE, plot=TRUE, verbose=TRUE)}
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
#interaction_location=interaction_file,
# #Adding Interactions: 
# int_df <- fread(interaction_location)
# class(int_df) <- "data.frame"    

#interactions=int_df, 
# set_beta=20, force_beta = TRUE


#NESclustering from CEMiTool file:
NESheatmap <- function(CEMiToolfile=cem){
  require(pheatmap)
  require(viridis)
  NES <- t(as.matrix(CEMiToolfile@enrichment$nes[,-1]))
  NES[is.na(NES)] <- 0
  colnames(NES) <- CEMiToolfile@enrichment$nes[,1]
  #ordering:
  NES <- NES[,order(nchar(colnames(NES)))]
  heatmap <- pheatmap::pheatmap(
    mat=NES,
    color= viridis::viridis(length(seq(min(NES, na.rm = T), max(NES,na.rm = T), length.out = 10)) -1),
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    cluster_cols = F,
    fontsize = 16,
    main="Coexpression Modules",
    legend_breaks = c(-4, -2, 0, 2, 4, max(NES)),
    legend_labels= c("-4","-2","0","2","4", "NES\n")                    
  )
  
  return(heatmap)
}
#running M3C clustering with only a subset of genes:
m3cassignments<- function (mydata,genes=rownames(mydata),... ){
  res2 <- M3C::M3C(mydata[rownames(mydata) %in% genes,], method = 1, cores=detectCores(), seed = 135, repsref=100, repsreal = 100)
  return(res2)
}

selectedgeneumap <- function(mydata = assay(vsd)[,metadata[,]$SampleID],
                             labels = FALSE,
                             seed = FALSE,
                             shape=NULL,
                             axistextsize = 6,
                             legendtextsize = 6,
                             dotsize = 1,
                             textlabelsize = 4,
                             legendtitle = "Group",
                             genes = rownames(mydata),
                             n_neighbors=15,
                             ellipse = F,
                             shapelegendtitle=FALSE,
                             metric="euclid",
                             colvec = c("skyblue", "gold","violet", "darkorchid", "slateblue", "forestgreen", "violetred", "orange", "midnightblue", "grey31", "black"),
                             text = FALSE){
  if (seed != FALSE) {set.seed(seed)}
  labels <- as.vector(labels)
  labels[is.na(labels)] <- "Unknown"
  labels <- as.factor(labels)
  scores <- data.frame(uwot::umap(t(as.matrix(mydata[rownames(mydata[rownames(mydata) %in% genes,]), ])), n_neighbors=n_neighbors, metric=metric))
  if(is.null(shape)){
    p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
      geom_point(aes(colour = labels), size = dotsize)+
      theme_bw() + theme(legend.position="bottom",legend.spacing.x = unit(0.1,"mm"),legend.margin=ggplot2::margin(0,0,0,0),legend.box.margin=ggplot2::margin(0,0,0,0),text=element_text(),panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize), 
                         axis.title.y = element_text(size = axistextsize), 
                         legend.title = element_text(size = legendtextsize), 
                         legend.text = element_text(size = legendtextsize))+
      labs(colour = legendtitle) + scale_colour_manual(values = colvec)+guides(colour=guide_legend(label.position = "bottom"))}
  else{
    p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
      geom_point(aes(colour = labels, shape=shape), size = dotsize)+ scale_colour_manual(values = colvec) + scale_shape_manual(values=c(16,15,17,18,4,5,6,3,8,9,10,11,12,13,14,21,22,23,24,25),labels=paste0("C",unique(shape)))+
      theme_bw() + theme(legend.position="bottom",legend.spacing.x = unit(0.1,"mm"),legend.margin=ggplot2::margin(0,0,0,0),legend.box.margin=ggplot2::margin(0,0,0,0),legend.direction = "horizontal", panel.grid.major = element_blank(), text=element_text(), 
                         panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                        colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                      colour = "black"), axis.title.x = element_text(size = axistextsize), 
                         axis.title.y = element_text(size = axistextsize), 
                         legend.title = element_text(size = legendtextsize), 
                         legend.text = element_text(size = legendtextsize)) + 
      labs(colour = legendtitle, shape = shapelegendtitle)+guides(color=guide_legend(label.position = "bottom",nrow=1),shape=guide_legend(label.position = "bottom",nrow=1))  }
  
  if(ellipse==TRUE){p <- p+stat_ellipse(aes(colour = labels), level=0.30,size=1, type="euclid")}
  if(text[1]!=FALSE){p <- p + geom_text(aes(label=text))}
  if(isFALSE(labels[1])){p <- p +guides(color=FALSE)}
  #p <- p+coord_fixed()
  p <- p+xlab("UMAP X1") + ylab("UMAP X2")
  return(p)
}
#values = c(16,17,18,19,15,7,8,9,10,11,as.numeric(paste0(0:6))),

#scale_fill_manual(values=c(RColorBrewer::brewer.pal(length(unique(metadata$severity)),"Set2")),name="",labels=c("Control","mild","severe"))

#ORA plot for all modules in one:
ORAplot <- function(cemi=NULL, outputfile="ORAplot",filetype="png",... ){
  #gsub("(GO:","\n(GO:",as.character(cemDiagnose@ora$ID), fixed=T)
  ora <- data.frame("Module"=cemi@ora$Module, "ID"=gsub("([^ ]+ [^ ]+ [^ ]+ [^ ]+) ", "\\1\n", as.character(cemi@ora$ID)), "GeneRatio"=cemi@ora$GeneRatio, "p.adjust"=cemi@ora$p.adjust)
  ora <- ora[  order( ora[,1], ora[,4] ),  ]
  oratop5 <-
    ora %>% 
    group_by(Module) %>% 
    slice_min(p.adjust, n = 5, with_ties=F)
  oratop5 <- as.data.frame(oratop5)
  oratop5[,"log10pvalue"] <- -log10(oratop5$p.adjust)
  oratop5[,"alpha"] <- 1
  if(length(oratop5[oratop5$p.adjust > 0.05,]$alpha)>0){oratop5[oratop5$p.adjust > 0.05,]$alpha <- 0}
  oratop5$ID <- ave(oratop5$ID, FUN = function(i) paste(strrep(" ", rowid(oratop5$ID)-1),i, sep=""))
  oratop5$ID <- factor(oratop5$ID, levels = unique(oratop5[order(oratop5[,1], oratop5[,4], decreasing=T),]$ID))
  #ordering to sort M10 after M9:
  oratop5$Module <- ordered(oratop5$Module , levels = c(paste0("M",1:(nmodules(cemi)-1)),"Not.Correlated"))
  
  ORAplot <-ggplot(oratop5)+geom_bar(stat="identity",aes(y=ID, x=log10pvalue,alpha=alpha, fill=interaction(Module)))+theme_bw()+
    theme(text=element_text( size = 12),legend.box.margin=ggplot2::margin(0,1,1,1), legend.text = element_text(size =12),
          axis.text=element_text(size=12, colour = "black"),  legend.justification="left", legend.title=element_blank(),axis.title.y = element_blank(),axis.text.y=element_text(margin = margin(r = 10)),legend.position="bottom", legend.direction = "horizontal", legend.key.height = unit(0.1, "mm"), legend.margin=margin(0,0,0,0))+
    scale_alpha(range=c(0.4, 1), guide="none") +
    geom_vline(xintercept=-log10(0.05), colour="grey", linetype="dotted")+ 
    scale_x_continuous(limits=c(0,5),oob = scales::rescale_none)+
    guides(fill=guide_legend(label.position = "bottom"))+
    labs(x="-log10(adjusted p-value)", title=element_blank(), y=element_blank())
  ORAplot
  # ggsave(plot=ORAplot, filename=paste0(output_location,"/",cohortname,"_",outputfile,".",filetype),dpi="retina",width = 78, height =120, units = "mm", device=filetype)
  return(ORAplot)
}
#+scale_fill_manual(labels=c(paste0("M",1:(nmodules(cem)-1)),"Not Cor."))
#combine NES heatmap of two CEMiTool runs if they share identical co-expression modules:
combinedNESheatmap <- function(cem1=cemDiagnose, cem2=cemim3c, groupnameforcem2="C",outputfile=paste0(cohortname,"_","combinedNESheatmap"), filetype = "png", ...){ 
  if(all(module_genes(cem1) == module_genes(cem2))){
    #building a heatmap plot of NES of severity and M3C in one:
    NES <- t(as.matrix(cem2@enrichment$nes[,-1]))
    colnames(NES) <- cem2@enrichment$nes[,1]
    rownames(NES) <-paste0(groupnameforcem2,unique(cem2@sample_annotation$Class))
    NES<- rbind(NES,t(as.matrix(cem1@enrichment$nes[,c(-1)])) )
    NES[is.na(NES)] <- 0
    #ordering:
    NES <- NES[,order(nchar(colnames(NES)))]
    heatmap <- pheatmap::pheatmap(
      mat=NES,
      cluster_cols = F,
      color= viridis::viridis(length(seq(min(NES, na.rm = T), max(NES, na.rm = T), length.out = 10)) -1),
      show_colnames     = TRUE,
      show_rownames     = TRUE,
      treeheight_row = 10,
      fontsize = 12,
      legend_breaks = c(-4, -2, 0, 2, 4, max(NES)),
      legend_labels= c("-4","-2","0","2","4", "NES\n"),
      width = 10,
      height = 10)
    heatmap <- as.ggplot(heatmap)+theme(plot.margin = unit(c(1,0,0,0), "cm"))
    ggsave(plot=heatmap, filename=paste0(output_location,"/",outputfile,".",filetype),dpi=400,width = 120, height = 120, units = "mm", device=filetype)
    return(heatmap)
  }else{print("Error: Modules are not identical")}}
getnormalizedDESeq2object <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname="Our"){
  require(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  #cohorttable <- fread("/home/sukmb465/Documents/Eike/Nextflow/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
  #cohortname <- "Our"
  ##cohortname <- dlg_list(cohorttable[,Cohort], multiple = FALSE)$res
  
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file=paste0(projectdir,"/00_RawData/",counttable))
  
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
  
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~PlateNr + Diagnose)
  dds <- DESeq(dds)
  #DESeqperformed <- TRUE
  #vsde <- vst(dds, blind=FALSE)
  res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
  vsd <- vst(dds, blind=FALSE)
  return(vsd)
}
getrawcounts <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname="Our"){
  library(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  #cohorttable <- fread("/home/sukmb465/Documents/Eike/Nextflow/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
  #cohortname <- "Our"
  ##cohortname <- dlg_list(cohorttable[,Cohort], multiple = FALSE)$res
  
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file=paste0(projectdir,"/00_RawData/",counttable))
  
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

  return(countmatrix)
}
getmetadata <- function(corhorttablelocation=paste0(projectdir,"/00_RawData/cohorts.csv"),cohortname="Our"){
  library(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  #cohorttable <- fread("/home/sukmb465/Documents/Eike/Nextflow/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
  #cohortname <- "Our"
  ##cohortname <- dlg_list(cohorttable[,Cohort], multiple = FALSE)$res
  
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file=paste0(projectdir,"/00_RawData/",counttable))
  
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
  # 
  # library(DESeq2)
  # dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~PlateNr + Diagnose)
  # dds <- DESeq(dds)
  # #DESeqperformed <- TRUE
  # #vsde <- vst(dds, blind=FALSE)
  # res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
  # vsd <- vst(dds, blind=FALSE)
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
return(vsd)}
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
  return(dds)}

#Taken from github.com/federicomarini/pcaExplorer/blob/master/R/pcaplot.R:
pcaplot <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE,title = NULL,
                     pcX = 1, pcY = 2, text_labels = TRUE, point_size = 3,
                     ellipse = TRUE, ellipse.prob = 0.95) # customized principal components
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(x)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, pcX], PC2 = pca$x[, pcY], group = group,
                  intgroup.df, names = colnames(x))
  colnames(d)[1] <- paste0("PC", pcX)
  colnames(d)[2] <- paste0("PC", pcY)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  # clever way of positioning the labels - worked good, then no need with ggrepel
  d$hjust <- ifelse((sign(d[, paste0("PC", pcX)]) == 1), 0.9, 0.1)# (1 + varname.adjust * sign(PC1))/2)
  
  g <- ggplot(data = d, aes_string(x = paste0("PC", pcX), y = paste0("PC", pcY), color = "group")) +
    geom_point(size = point_size) +
    xlab(paste0("PC", pcX, ": ", round(percentVar[pcX] * 100, digits = 2), "% variance")) +
    ylab(paste0("PC", pcY ,": ", round(percentVar[pcY] * 100, digits = 2), "% variance"))
  
  ## plot confidence ellipse
  # credit to vince vu, author of ggbiplot
  if (ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(d, "group", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x[[paste0("PC", pcX)]], x[[paste0("PC", pcY)]]))
      mu <- c(mean(x[[paste0("PC", pcX)]]), mean(x[[paste0("PC", pcY)]]))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$group[1])
    })
    # names(ell)[1:2] <- c('xvar', 'yvar')
    if (nrow(ell) > 0) {
      g <- g + geom_path(data = ell, aes_string(x = "X1", y = "X2", color = "groups", group = "groups"))
    }
  }
  
  if (text_labels)
    g <- g + geom_label_repel(mapping = aes_string(label = "names", fill = "group"),
                              color = "white", show.legend = TRUE) 
  if (!is.null(title)) g <- g + ggtitle(title)
  g <- g + theme_bw()
  # as in http://www.huber.embl.de/msmb/Chap-Graphics.html
  # "well-made PCA plots usually have a width that’s larger than the height"
  g <- g + coord_fixed()
  g
}

severityspineplot <- function(){
  metadata$severity <- droplevels(metadata$severity)
  library(tidyverse)
  library(dplyr)
  mosaicplotdata <- metadata %>%
    group_by(m3c, severity, .drop=F) %>%
    dplyr::summarise(count = n()) %>%
    mutate(cut.count = sum(count),
           prop = count/sum(count)) %>%
    ungroup()
  
  # supp.labs <- c("Cluster 1 (n=301)", "Cluster 2 (n=311)", "Cluster 3 (n=126)")
  # names(supp.labs) <- c("1", "2","3")
  
  supp.labs <- paste0("C",unique(metadata$m3c),"\n n=", table(metadata$m3c))
  names(supp.labs) <- unique(metadata$m3c)
  
  mosaicplotdata$m3c <- factor(mosaicplotdata$m3c, ordered = T, levels=as.data.table(mosaicplotdata)[order(-prop),][severity %in% "Control",]$m3c)
  
  mosaicplot <- ggplot(mosaicplotdata,
                       aes(x = m3c, y = prop, width = cut.count, fill = as.factor(severity))) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +
    # geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5)) + # if labels are desired
    facet_grid(~m3c, scales = "free_x", space = "free_x", labeller = labeller(m3c=supp.labs)) +
    scale_fill_manual(values=c(RColorBrewer::brewer.pal(length(unique(metadata$severity)),"Set2")),name="",labels=c("Control","mild","severe", "Unknown"))+
    # theme(panel.spacing.x = unit(0, "npc")) + # if no spacing preferred between bars
    theme_bw() +ylab("Frequency")+labs(fill="Observed\nSeverity")+theme(legend.margin=ggplot2::margin(0,0,0,0),
                                                                        legend.box.margin=ggplot2::margin(0,-0,0,0),
                                                                        legend.position="bottom",legend.key.height = unit(0.1, "mm"),panel.grid.major = element_blank(), text=element_text( size=6),legend.title=element_text( size=7), legend.text = element_text(size=7), 
                                                                        panel.grid.minor = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank())+ylab("Proportion")+scale_x_discrete(breaks=NULL)+scale_y_continuous(limits=c(0,1), expand=c(0,0))+guides(fill=guide_legend(label.position = "bottom",override.aes = list(size=.5)))
  #, "#F0F0EB"
  #RColorBrewer::brewer.pal(length(unique(metadata$supervised)),"YlOrRd")
  return(mosaicplot)
}

getCohortsVSD <- function(corhorttablelocation="/home/sukmb465/Documents/Eike/Nextflow/cohorts-nextflow2/cohorts/cohorts.csv",cohortname="Our"){
  library(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  #cohorttable <- fread("/home/sukmb465/Documents/Eike/Nextflow/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
  #cohortname <- "Our"
  ##cohortname <- dlg_list(cohorttable[,Cohort], multiple = FALSE)$res
  
  counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
  metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
  metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
  Disease <- as.vector(unlist(cohorttable[Cohort %in% cohortname, 5]))
  UC <- fread(file=counttable)
  
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
  metadata <- fread(metatable, header=metatableheader)
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
  library(DESeq2)
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
  
  # 
  # 
  # library(DESeq2)
  # dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~Diagnose)
  # dds <- DESeq(dds)
  # #DESeqperformed <- TRUE
  # #vsde <- vst(dds, blind=FALSE)
  # res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
  # vsd <- vst(dds, blind=FALSE)
  return(vsd)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  if(coef(m)[2]>=0){
  eq <- substitute(italic(y) == a ~+~~  b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(abs(coef(m)[2])), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  }else{
    eq <- substitute(italic(y) == a ~-~~  b %.% italic(x)*","~~italic(r)^2~"="~r2,
                         list(a = format(unname(coef(m)[1]), digits = 2),
                              b = format(unname(abs(coef(m)[2])), digits = 2),
                              r2 = format(summary(m)$r.squared, digits = 3)))
    
  }
  as.character(as.expression(eq));
}

