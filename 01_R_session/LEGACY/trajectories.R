library(data.table)
#library(svDialogs)

getnormalizedDESeq2object <- function(corhorttablelocation="/home/sukmb465/Documents/cohorts-nextflow2/cohorts/cohorts.csv",cohortname="Our"){
library(data.table)
cohorttable <- fread(corhorttablelocation, header=T)
#cohorttable <- fread("/home/sukmb465/Documents/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
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
dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~PlateNr + Diagnose)
dds <- DESeq(dds)
#DESeqperformed <- TRUE
#vsde <- vst(dds, blind=FALSE)
res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
vsd <- vst(dds, blind=FALSE)
return(vsd)
}
getrawcounts <- function(corhorttablelocation="/home/sukmb465/Documents/cohorts-nextflow2/cohorts/cohorts.csv",cohortname="Our"){
  library(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  #cohorttable <- fread("/home/sukmb465/Documents/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
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
  # 
  # library(DESeq2)
  # dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~PlateNr + Diagnose)
  # dds <- DESeq(dds)
  # #DESeqperformed <- TRUE
  # #vsde <- vst(dds, blind=FALSE)
  # res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
  # vsd <- vst(dds, blind=FALSE)
  return(countmatrix)
}
getmetadata <- function(corhorttablelocation="/home/sukmb465/Documents/cohorts-nextflow2/cohorts/cohorts.csv",cohortname="Our"){
  library(data.table)
  cohorttable <- fread(corhorttablelocation, header=T)
  #cohorttable <- fread("/home/sukmb465/Documents/cohorts-nextflow2/cohorts/cohorts.csv", header=T)
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
  # 
  # library(DESeq2)
  # dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=metadata, design= ~PlateNr + Diagnose)
  # dds <- DESeq(dds)
  # #DESeqperformed <- TRUE
  # #vsde <- vst(dds, blind=FALSE)
  # res <- results(dds,  contrast=c("Diagnose",Disease,"Control"))
  # vsd <- vst(dds, blind=FALSE)
  return(metadata)
}


vsd <- getnormalizedDESeq2object()
# geneset <- as.vector(unlist(fread("//home/sukmb465/Documents/Eike/Rprojects/UC/Pivot/a_geneset_rf15_optimal.txt", header=F)))
geneset <- as.vector(unlist(fread("//home/sukmb465/Documents/Eike/Rprojects/UC/Pivot/markers.txt", header=F)))
geneset <- gsub("-",".",geneset, perl=T)
counts <-assay(vsd)[,]
rawcounts <- getrawcounts()
metadata <- getmetadata()
# rawcounts <- rawcounts[,]

# trawcounts <-t(rawcounts) 
# tcounts <- t(counts)
# 
# library(dyno)
library(tidyverse)
# library(dynwrap)
# library(dynmethods)
# dataset <- wrap_expression(
#   counts = rawcounts,
#   expression = counts)

# dataset <- add_grouping(
#   dataset,
#   example_dataset$grouping
# )

# #Check for best method to use on dataset:
# library(dynguidelines)
# guidelines_shiny(dataset = dataset)


# test <- infer_trajectory(dataset, ti_scorpius(), verbose = TRUE)

# 
# model <- infer_trajectory(dataset, ti_slingshot(), verbose=TRUE)
# plot_dimred(model)
library(data.table)
# library(slingshot)
# library(M3C)
library(uwot)
# library(SingleCellExperiment)
# rawcounts[geneset,]
# counts[geneset,]


# rd2 <- uwot::umap(t(counts[geneset,]))
# colnames(rd2) <- c('UMAP1', 'UMAP2')
# 
# plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
# 
# m3cassignments<- function (mydata,genes=rownames(mydata),... ){
#   res2 <- M3C::M3C(mydata[rownames(mydata) %in% genes,], method = 1, cores=detectCores(), seed = 135, repsref=150, repsreal = 150)
#   return(res2)
# }
# sim <- SingleCellExperiment::SingleCellExperiment(assays = List(counts = counts))
# 
# 
# pca <- prcomp(t(counts[geneset,]), scale. = FALSE)
# rd1 <- pca$x[,1:2]
# 
# plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
# cluster <- m3cassignments(counts, geneset)
# 
# 
# reducedDims(sim) <- SimpleList(PCA = rd1, UMAP = rd2)
# colData(sim)$M3C <- cluster$assignments
# 
# 
# library(RColorBrewer)
# plot(rd1, col = brewer.pal(9,"Set1")[cluster$assignments], pch=16, asp = 1)
# 
# sim <- slingshot(sim, clusterLabels = 'M3C', reducedDim = 'UMAP')
# 
# summary(sim$slingPseudotime_1)
# 
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
# 
# plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sim), lwd=2, col='black')
# 
# 
# # p <- ggplot(data = scores, aes(x = X1, y = X2)) + 
# #   geom_point(aes(colour = labels, shape=shape), size = dotsize)+ scale_colour_manual(values = colvec) + scale_shape_manual(values=c(16,15,17,18,4,5,6,3,8,9,10,11,12,13,14,21,22,23,24,25),labels=paste0("C",unique(shape)))+
# #   theme_bw() + theme(legend.position="bottom",legend.spacing.x = unit(0.1,"mm"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),legend.direction = "horizontal", panel.grid.major = element_blank(), text=element_text(family="Arial"), 
# #                      panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
# #                                                                                     colour = "black"), axis.text.x = element_text(size = axistextsize, 
# #                                                                                                                                   colour = "black"), axis.title.x = element_text(size = axistextsize), 
# #                      axis.title.y = element_text(size = axistextsize), 
# #                      legend.title = element_text(size = legendtextsize), 
# #                      legend.text = element_text(size = legendtextsize)) + 
# #   labs(colour = legendtitle, shape = shapelegendtitle)+guides(color=guide_legend(label.position = "bottom",nrow=1),shape=guide_legend(label.position = "bottom",nrow=1))  }
# 
# 
# ggplot(data=data.frame(rd2), aes(x = UMAP1, y= UMAP2))+
#   geom_point(aes(shape=as.factor(metadata$severity), color=sim$slingPseudotime_1))+
#   geom_smooth(se=FALSE,method="gam", color='black')+
#   coord_fixed()
# 
# 
# df <- data.table(t(counts[geneset,]),Pseudotime=sim$slingPseudotime_1, severity=metadata$severity)

library(SCORPIUS)
rd2 <- uwot::umap(t(counts[geneset,]))
colnames(rd2) <- c('UMAP1', 'UMAP2')
traj <- SCORPIUS::infer_trajectory(rd2)
draw_trajectory_plot(rd2, metadata$severiy, traj$path, contour = TRUE)

#traj$time


ggplot(data=data.frame(rd2), aes(x = UMAP1, y= UMAP2))+
  geom_point(aes(shape=as.factor(metadata$severity), color=traj$time))+
  geom_line(data=data.frame(traj$path), aes(x=Comp1, y=Comp2))+
  coord_fixed()

ggplot(data=data.frame(rd2), aes(x = UMAP1, y= UMAP2))+
  geom_point(aes(shape=as.factor(cut(as.vector(traj$time), seq(0,1,0.2), include.lowest=T)), color=as.factor(metadata$severity)))+
  geom_line(data=data.frame(traj$path), aes(x=Comp1, y=Comp2))+
  coord_fixed()+ theme(legend.position="right", legend.box="horizontal")+labs(shape="Pseudotime", color="Disease severity")

df <- data.table(t(counts[geneset,]),Pseudotime=traj$time, severity=metadata$severity)
if(median(df[df$severity %in% "Control"]$Pseudotime)>=0.5){df$Pseudotime <- 1-df$Pseudotime}
library(reshape2)
# df <- merge(counts[geneset,], pseudotimes, by="SampleID")
melted <- melt(df[,colnames(df) %in% c("Pseudotime", geneset, "severity"), with=FALSE], id.vars=c('Pseudotime','severity'))
expressiontopseudotime <-  ggplot(melted,aes(x=Pseudotime,y=value,group='severity'))+
  geom_point(aes(color=severity))+geom_quantile(colour="black")+ 
  #geom_density(color="darkblue", fill="lightblue")+
  facet_wrap(~variable,ncol = 2,scales = 'free')+theme_minimal()+scale_color_manual(values=RColorBrewer::brewer.pal(length(unique(melted$severity)),"Accent"))
expressiontopseudotime
+ geom_density(color="darkblue", fill="lightblue")

gimp <- gene_importances(
  data.table(t(counts[geneset,])),
  traj$time,
  num_permutations = 10,
  num_threads = 8,
  ntree = 10000,
  ntree_perm = 1000
)
gimp
