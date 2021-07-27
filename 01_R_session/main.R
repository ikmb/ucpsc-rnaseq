library(data.table)
library(plyr)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(CEMiTool)
library(M3C)
library(rmarkdown)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggridges)
library(ggplotify)
library(VennDiagram)
library(extrafont)
library(SCORPIUS)
library(reshape2)

library(ggrepel)
library(ggvenn)
source("functions.R")

projectdir <- dirname(getwd())


#####Inputs#####
corhorttablelocation <- paste0(projectdir,"/00_RawData/cohorts.csv")
#corhorttablelocation <- "/home/sukmb465/Documents/cohorts-nextflow2/cohorts/cohorts.csv"
cohorttable <- fread(corhorttablelocation, header=T)
cohortname <- "Our"
# counttable <-unlist(cohorttable[Cohort %in% cohortname,2])
# metatable <- unlist(cohorttable[Cohort %in% cohortname,3])
# metatableheader <- as.logical(unlist(cohorttable[Cohort %in% cohortname,4]))
Disease <- unlist(cohorttable[Cohort %in% cohortname, 5])

gmt_file <- paste0(projectdir,"/data_rescources/GO_Biological_Process_2018.gmt")
output_location <- paste0("output")

rawcounts <- getrawcounts(corhorttablelocation=corhorttablelocation,cohortname=cohortname)
# vsd <- getnormalizedDESeq2object(corhorttablelocation=corhorttablelocation,cohortname=cohortname)
metadata <- getmetadata(corhorttablelocation=corhorttablelocation,cohortname=cohortname)

########Severity Defintion########
metadata$severity <- "Control"
metadata[metadata$supervised %in% "Unknown",]$severity <- "Unknown"
metadata[metadata$supervised %in% c(2,3,4,5),]$severity <- rep("mild", length(metadata[metadata$supervised %in% c(2,3,4,5),]$SampleID))
metadata[metadata$supervised %in% c(6,7,8),]$severity <- rep("severe", length(metadata[metadata$supervised %in% c(6,7,8),]$SampleID))
metadata$severity <- ordered(metadata$severity, c("Control", "mild", "severe", "Unknown"))
metadata$severity <- factor(metadata$severity, ordered = F)

#Discovery dataset contains only patients with known severity:
discoveryvsd <- vsdofallknownseverity(meta=metadata, countmatrix=rawcounts, batch="PlateNr")
vsd <- discoveryvsd

####Overview PCA####
OverviewPCA <- pcaplot(vsd, intgroup = "severity", text_labels = F, ntop=length(vsd@rowRanges))+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#999999"),name = "Phenotype", labels = c("Control", "Mild UC", "Severe UC"))+
  guides(colour=guide_legend(label.position = "bottom",keywidth=0.5,
                             keyheight=0.1,
                             default.unit="inch"))+theme_bw()+
  theme(text = element_text(size=16),legend.position="bottom",legend.spacing.x = unit(0.1,"mm"),legend.box.margin=ggplot2::margin(0,0,0,0),legend.direction = "horizontal")
ggsave(plot=OverviewPCA, filename=paste0(output_location,"/",cohortname,"_","OverviewPCA",".","pdf"),dpi="retina",width = 120, height =120, units = "mm", device="pdf")

#SCREE Plot of PCA#
library(pcaExplorer)
screeplot <- pcaExplorer::pcascree(prcomp(t(assay(vsd)[, ])), type="pev",pc_nr = 20)+theme_bw()
ggsave(plot=screeplot, filename=paste0(output_location,"/",cohortname,"_","OverviewPCAscree",".","pdf"),dpi="retina",width = 120, height =120, units = "mm", device="pdf")

####Overview UMAP####
# OverviewUMAP <- selectedgeneumap(assay(vsd),seed=135, n_neighbors=15, labels=as.factor(vsd$severity), colvec=c("#00AFBB", "#E7B800", "#FC4E07","#999999") ,text=FALSE,ellipse=T,legendtitle = "Phenotype",axistextsize = 20, legendtextsize = 20, dotsize = 5 )+ 
#   theme_bw()+scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#999999"),name = "Phenotype", labels = c("Control", "Mild UC", "Severe UC"))+theme(legend.position="bottom")+coord_fixed()

OverviewUMAP <- selectedgeneumap(assay(vsd),seed=135, n_neighbors=15, labels=as.factor(vsd$severity), colvec=c("#00AFBB", "#E7B800", "#FC4E07","#999999") ,text=FALSE,ellipse=T,legendtitle = "Phenotype",axistextsize = 20, legendtextsize = 20, dotsize = 3 )+
guides(colour=guide_legend(label.position = "bottom",keywidth=0.5,
                           keyheight=0.1,
                           default.unit="inch"))+theme_bw()+
  theme(text = element_text(size=16),legend.position="bottom",legend.spacing.x = unit(0.1,"mm"),legend.box.margin=ggplot2::margin(0,0,0,0),legend.direction = "horizontal")
ggsave(plot=OverviewUMAP, filename=paste0(output_location,"/",cohortname,"_","OverviewUMAP",".","pdf"),dpi="retina",width = 120, height =120, units = "mm", device="pdf")

####CEMiTool of Diagnose, all Samples####
cemDiagnose <- CEMiwrapper(expressionmatrix=assay(vsd), ID=vsd$SampleID, Groups=vsd$severity, reportname="Our_severity")
#NES heatmap of Diagnose CemiTool:
Diagnoseheatmap <- NESheatmap(cemDiagnose)
ggsave(plot=Diagnoseheatmap, filename=paste0(output_location,"/",cohortname,"_","Diagnoseheatmap",".","pdf"),dpi="retina",width = 120, height =120, units = "mm", device="pdf")

#ORA plot of Diagnose CemiTool:
ORAp <- ORAplot(cemi = cemDiagnose, outputfile ="ORAplot", filetype = "pdf")
ggsave(plot=ORAp, filename=paste0(output_location,"/",cohortname,"_","ORA",".","pdf"),dpi="retina",width = 160, height =260, units = "mm", device="pdf")

#####Markergeneset Identification######
discoverydds <- ddsofallknownseverity(meta=metadata, countmatrix=rawcounts, batch="PlateNr")
ressevcon <- results(discoverydds,  contrast=c("severity","severe","Control"),tidy = T)
resmilcon <- results(discoverydds,  contrast=c("severity","mild","Control"),tidy = T)
ressevmil <- results(discoverydds,  contrast=c("severity","severe","mild"),tidy = T)
allsigngenes <- list(CEMiTool::module_genes(cemDiagnose)[!(CEMiTool::module_genes(cemDiagnose)$modules %in% "Not.Correlated"),]$genes,
                     ressevcon[(abs(ressevcon$log2FoldChange)>0.25)&(ressevcon$padj < 0.01), ]$row,
                     ressevmil[(abs(ressevmil$log2FoldChange)>0.25)&(ressevmil$padj < 0.01), ]$row,
                     resmilcon[(abs(resmilcon$log2FoldChange)>0.25)&(resmilcon$padj < 0.01), ]$row,
                     as.vector(unlist(fread("rf500_with_seed.txt", header=F)   )))
allsigngenes <- lapply(allsigngenes, function(x) gsub("-",".",x,perl=T))
library("ggvenn")
DEVENN <- ggvenn(list("Severe UC vs.\nControl"=ressevcon[(abs(ressevcon$log2FoldChange)>0.25)&(ressevcon$padj < 0.01), ]$row,
            "Severe UC vs.\nMild UC"=ressevmil[(abs(ressevmil$log2FoldChange)>0.25)&(ressevmil$padj < 0.01), ]$row,
            "Mild UC vs.\nControl"=resmilcon[(abs(resmilcon$log2FoldChange)>0.25)&(resmilcon$padj < 0.01), ]$row), stroke_size=0.5, text_size = 5, set_name_size = 5,fill_color = rev(c("#00AFBB", "#E7B800", "#FC4E07")), fill_alpha = .5)

MMVENN <- ggvenn(list("Filtered by\nCEMiTool"=CEMiTool::module_genes(cemDiagnose)[!(CEMiTool::module_genes(cemDiagnose)$modules %in% "Not.Correlated"),]$genes,
            "Differential Expressed between \n severe UC, mild UC and Controls"=Reduce(intersect,lapply(list("Severe UC vs.\nControl"=ressevcon[(abs(ressevcon$log2FoldChange)>0.25)&(ressevcon$padj < 0.01), ]$row,
                                                                                                                       "Severe UC vs.\nMild UC"=ressevmil[(abs(ressevmil$log2FoldChange)>0.25)&(ressevmil$padj < 0.01), ]$row,
                                                                                                                       "Mild UC vs.\nControl"=resmilcon[(abs(resmilcon$log2FoldChange)>0.25)&(resmilcon$padj < 0.01), ]$row), function(x) gsub("-",".",x,perl=T))),
            "Feature Selected by\nRandom Forest"=as.vector(unlist(fread("rf500_with_seed.txt", header=F)   ))),stroke_size=0.5, text_size = 5, set_name_size = 5) 
ggsave(plot=DEVENN, filename=paste0(output_location,"/",cohortname,"_","DEVENN",".","svg"),dpi="retina",width = 150, height = 150, units = "mm", device="svg")
ggsave(plot=MMVENN, filename=paste0(output_location,"/",cohortname,"_","MMVENN",".","svg"),dpi="retina",width = 150, height = 150, units = "mm", device="svg")

MMGS <- Reduce(intersect,allsigngenes)
MMGS
fwrite(data.table(MMGS),file=paste0(output_location,"/","markergenes.txt"),col.names = FALSE)

RFGS <- as.vector(unlist(fread("rf50_with_seed.txt", header=F)))
RFGS

#DEGS:
# ressevcon[order(ressevcon$padj),][1:23,]$row

#Discovery CEMiTool:
cemDiscovery <- CEMiwrapper(expressionmatrix=assay(discoveryvsd), ID=discoveryvsd$SampleID, Groups=discoveryvsd$severity, reportname=paste0(cohortname,"_Severity"))
Discoveryheatmap <- NESheatmap(cemDiscovery)
ggsave(plot=Discoveryheatmap, filename=paste0(output_location,"/",cohortname,"_","Discoveryheatmap",".","svg"),dpi="retina",width = 120, height =120, units = "mm", device="svg")
DiscoveryORA <-  ORAplot(cemi = cemDiscovery, outputfile ="Discovery_ORA", filetype = "pdf")
ggsave(plot=DiscoveryORA, filename=paste0(output_location,"/",cohortname,"_","DiscoveryORA",".","svg"),dpi="retina",width = 160, height =230, units = "mm", device="svg")



#What sets to continue with?
marker <- list(DEGS=as.vector(unlist(fread("/home/sukmb465/Desktop/UCRNA/DE23.txt", header=F))),RFGS=RFGS,MMGS=MMGS)
# markergenes <- marker[["MMGS"]]
# geneset <- "MMGS"
m3clusters <- list()
pseudocorrelation <- list()
for (geneset in names(marker)) {
  print(geneset)
  markergenes <- marker[[geneset]]
  
  print(marker[[geneset]])
  





####Clustering based on markergenes####
  
# assignments <- m3cassignments(assay(vsd)[,], genes=markergenes, seed=135)
# #selectedgeneumap(assay(vsd)[,metadata[,]$SampleID],seed=135, n_neighbors=15, genes=markergenes, labels=as.factor(assignments$assignments), colvec=c(RColorBrewer::brewer.pal(ifelse(3<length(levels(as.factor(assignments$assignments)))-1, length(levels(as.factor(assignments$assignments)))-1,3),"YlOrRd"), "#F0F0EB") , dotsize = 2,ellipse = T )
# # metadata$m3c <- assignments$assignments
# vsd$m3c <- assignments$assignments
# 
# m3clusters[[geneset]] <- assignments$assignments
# 
# ggsave(plot=grid.arrange(as.ggplot(assignments$plots[[1]]),as.ggplot(assignments$plots[[2]]),as.ggplot(assignments$plots[[3]]),as.ggplot(assignments$plots[[4]])), filename=paste0(output_location,"/",cohortname,"_",geneset,"_M3C",".","png"),dpi="retina",width = 2*178, height =2*178, units = "mm", device="png")

MarkergeneUMAP <- selectedgeneumap(assay(vsd)[,],seed=135, n_neighbors=15, genes=markergenes, labels=as.factor(vsd$severity), colvec=c("#00AFBB", "#E7B800", "#FC4E07","#999999") ,text=FALSE,ellipse=T, shape=as.factor(vsd$m3c), shapelegendtitle = element_blank(), legendtitle = element_blank())+theme_bw()+coord_fixed()
ggsave(plot=MarkergeneUMAP, filename=paste0(output_location,"/",cohortname,"_",geneset,"_UMAP.svg"),dpi="retina",width = 160, height =160, units = "mm", device='svg')

#Plot of Disease-Severity Distribution in found Clusters: 
mosaicplot <- severityspineplot()+theme_bw()+ 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#999999"))
ggsave(plot=mosaicplot, filename=paste0(output_location,"/",cohortname,"_",geneset,"_mosaicplot.svg"),dpi="retina",width = 60, height =60, units = "mm", device='svg')



#####SCORPIUS Trajectory#######
rd2 <- uwot::umap(t(assay(vsd)[rownames(assay(vsd)) %in% markergenes,]))
colnames(rd2) <- c('UMAP1', 'UMAP2')
traj <- SCORPIUS::infer_trajectory(rd2)
draw_trajectory_plot(rd2, vsd$severiy, traj$path, contour = TRUE)




df <- data.table(data.table(t(assay(vsd)[rownames(assay(vsd)) %in% markergenes,])),Pseudotime=traj$time, severity=vsd$severity)
if(median(df[df$severity %in% "Control"]$Pseudotime)>=0.5){df$Pseudotime <- 1-df$Pseudotime}

# df <- merge(counts[geneset,], pseudotimes, by="SampleID")
melted <- melt(df[,colnames(df) %in% c("Pseudotime", markergenes, "severity"), with=FALSE], id.vars=c('Pseudotime','severity'))
expressiontopseudotime <-  ggplot(melted,aes(x=Pseudotime,group='severity'))+
  geom_point(aes(color=severity,y=value))+geom_quantile(aes(y=value),colour="black")+theme_bw()+
  facet_wrap(~variable,ncol = ifelse(length(unique(melted$variable))>20,3,2),scales = 'free')+theme_minimal()+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
expressiontopseudotime
ggsave(plot=expressiontopseudotime, filename=paste0(output_location,"/",cohortname,"_",geneset,"_expressiontopseudotimeUMAP.svg"),dpi="retina",width = ifelse(length(unique(expressiontopseudotime$data$variable))<20,160*2,160*3), height = 55*round(length(markergenes)/ifelse(length(unique(expressiontopseudotime$data$variable))<20,2,3)), units = "mm", device='svg')


# marginalPseudotimeUMAP <- ggExtra::ggMarginal(
#   ggplot(data=data.frame(rd2), aes(x = UMAP1, y= UMAP2))+
#     geom_point(aes( color=as.factor(vsd$severity)), size=4)+
#     geom_line(data=data.frame(traj$path), aes(x=Comp1, y=Comp2))+theme_bw()+
#     coord_fixed()+ theme(legend.position="bottom", legend.box="horizontal", text=element_text(family="Arial", size=20))+labs(shape="Pseudotime", color="Disease severity")+ 
#     scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")), type="density", groupColour = TRUE)
# ggsave(plot=marginalPseudotimeUMAP, filename=paste0(output_location,"/",cohortname,"_",geneset,"_marginalPseudotimeUMAP.svg"),dpi=300,width = 150, height = 150, units = "mm", device='svg')





marginalPseudotimeUMAP <- ggExtra::ggMarginal(
  ggplot(data=data.frame(rd2), aes(x = UMAP1, y= UMAP2))+
    geom_point(aes( color=as.factor(vsd$severity)),size=4)+
    geom_path(aes_string("Comp1", "Comp2"), data.frame(traj$path), 
              size = 0.5, alpha = 1)+ 
    theme_bw()+
    theme(legend.position="bottom",text=element_text(family="Arial", size=18), legend.box="horizontal")+
    labs(shape="Pseudotime", color="UC Severity")+ 
    scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(vsd$severity)),"Set2")),type="density", groupColour = TRUE)
ggsave(plot=marginalPseudotimeUMAP, filename=paste0(output_location,"/",cohortname,"_",geneset,"_marginalPseudotimeUMAP.svg"),dpi=300,width = 150, height = 150, units = "mm", device='svg')


dg <- data.table(Pseudotime=ifelse(
                              rep(
                                median(
                                data.table(Pseudotime=traj$time, Severity=vsd$severity)[Severity %in% "Control",]$Pseudotime
                                )>=0.5, times=length(traj$time)
                              ),
                              as.vector(1-traj$time),
                              as.vector(traj$time)),
                 Supervised=ordered(vsd$supervised))
dg$Supervised <- revalue(dg$Supervised, c("0"="Control"))

PseudotimeSupervisedCorrelation <- ggplot(data = dg, aes(x=Pseudotime, y = Supervised)) + 
  geom_density_ridges(aes(fill= Supervised,y=Supervised,point_fill=Supervised,point_color=Supervised),jittered_points=TRUE,alpha=0.5) +#,stat="binline", bins=15
  theme_bw() + 
  labs(x = "Pseudotime", y = "Supervised\nSeverity Score") +
  geom_smooth(aes(Pseudotime,as.numeric(Supervised) ), method = "lm",formula="y~x")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(dg$Supervised )),"Set2"))+
                      #c("#00AFBB", rep("#E7B800",times=4), rep("#FC4E07", times=3))
  theme(legend.position="none", text=element_text(family="Arial", size=20))+
  xlim(-0.1,1.1)+
  scale_discrete_manual(aesthetics = "point_color", values =RColorBrewer::brewer.pal(length(unique(dg$Supervised )),"Set2")
                          #c("#00AFBB", rep("#E7B800",times=4), rep("#FC4E07", times=3))
                        )

  
ggsave(plot=PseudotimeSupervisedCorrelation, filename=paste0(output_location,"/",cohortname,"_",geneset,"_PseudotimeSupervisedCorrelation.svg"),dpi=300,width = 120, height =120, units = "mm", device='svg')
pseudocorrelation[[geneset]]<- cor.test(as.numeric(dg$Supervised), dg$Pseudotime, method="spearman")
}





















# lm_eqn <- function(df){
#   m <- lm(y ~ x, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(a = format(unname(coef(m)[1]), digits = 2),
#                         b = format(unname(coef(m)[2]), digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));
# }
# p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)
# ggplot(data.frame(Supervised=(as.numeric(replace(vsd$supervised, vsd$supervised==0,2))-2)/length(unique((as.numeric(replace(vsd$supervised, vsd$supervised==0,2))-2))),Pseudotime=1-traj$time), aes(x=Pseudotime, y=Supervised, color=Supervised))+geom_jitter()+geom_smooth(method = "glm", method.args=list(family="quasibinomial"))
# fit <- glm(formula= Supervised ~ Pseudotime, data=data.frame(Supervised=(as.numeric(replace(vsd$supervised, vsd$supervised==0,2))-2)/length(unique((as.numeric(replace(vsd$supervised, vsd$supervised==0,2))-2))),Pseudotime=1-traj$time), family = "quasibinomial")
# fit <-lm(formula= Supervised ~ Pseudotime, data=data.frame(Supervised=(as.numeric(replace(vsd$supervised, vsd$supervised==0,2))-2)/length(unique((as.numeric(replace(vsd$supervised, vsd$supervised==0,2))-2))),Pseudotime=1-traj$time))
# paste(coef(fit), names(coef(fit)), sep = ' * ', collapse = ' + ')
# 
# 
#     
# # Reduce(intersect,marker[-2])
# 
MarkerVenn <- ggvenn(marker,fill_color = rev(c("#00AFBB", "#E7B800", "#FC4E07")), fill_alpha = .5)
ggsave(plot=MarkerVenn, filename=paste0(output_location,"/",cohortname,"_","MarkerVenn",".","svg"),dpi="retina",width = 150, height = 150, units = "mm", device="svg")




#####celltypeblueprints#######
bprint_expression<-fread("E-MTAB-3827-query-results.tpms.tsv")
bprint_expression_l<-pivot_longer(bprint_expression, 3:29, names_to = "cell_type", values_to = "TPM")%>%as.data.table()
setnames(bprint_expression_l,c("gene_ID","gene_name","cell_type","TPM"))
celltypes_to_keep<-c("mature neutrophil", "erythroblast", "conventional dendritic cell", "macrophage", "CD4-positive, alpha-beta T cell",
                     "CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta thymocyte","CD8-positive, alpha-beta thymocyte")
expression_plots_of_module_genes<-function(expression_matrix,gene_module_table,output,celltypes_to_keep){
  n_modules<-length(unique(gene_module_table$modules))
  module_names<-unique(gene_module_table$modules)
  blueprints <- list()
  for (i in 1:n_modules) {
    #pick the Genes of the module from the gene-module list
    markergenes<-gene_module_table[modules==module_names[[i]],genes]
    #next step (matrix conversion of a data.table) bugs with modules length 1 because of a bug in as.matrix.data.table(). make a manual solution for this special case
    length_genes_blueprint_intersect<-length(intersect(markergenes,expression_matrix$`Gene Name`))
    genes_blueprint_intersect<-intersect(markergenes,expression_matrix$`Gene Name`)
    if(length_genes_blueprint_intersect==1){bprint_expression_mat<-matrix(data=as.vector(as.matrix(expression_matrix[`Gene Name`%in%markergenes,-2:-1],mode="numeric",rownames.value=genes_blueprint_intersect)),
                                                                          nrow = 1,
                                                                          ncol=dim(expression_matrix)[[2]]-2,
                                                                          dimnames = list(genes_blueprint_intersect,colnames(expression_matrix)[3:29]))}
    #the standard solution
    else{bprint_expression_mat<-as.matrix(expression_matrix[`Gene Name`%in%markergenes,-2:-1],rownames = expression_matrix[`Gene Name`%in%markergenes,`Gene Name`])}
    #kick zero values to avoid division by zero and stuff
    bprint_expression_mat[(bprint_expression_mat)==0]<-min(bprint_expression_mat,na.rm = T)
    #make it a z score of the log values
    z_bprint_expression_mat<-scale(t(log10(bprint_expression_mat)))
    #kick NAs
    z_bprint_expression_mat[is.na(z_bprint_expression_mat)]<-min(z_bprint_expression_mat,na.rm = T)
    #blueprints[[i]] <- (pheatmap(t(z_bprint_expression_mat),cluster_rows = ifelse(length_genes_blueprint_intersect==1,F,T),cluster_cols = F,fontsize=5,silent=F,cellheight = 5))
    #start a pdf device, which is very large, but possibly not large enough: couple the height to the number of blueprint intersect genes
    #pdf(file=paste0(output,"/",cohortname,"_",module_names[[i]],".pdf"), width=10,height=length_genes_blueprint_intersect*(5/72)+3)
    #pheatmap(t(z_bprint_expression_mat),cluster_rows = ifelse(length_genes_blueprint_intersect==1,F,T),cluster_cols = F,fontsize=5,silent=F,cellheight = 5)#,filename="marker_expression_BLUEPRINT_celltypes.png", width = 7, height =7)
    #dev.off()
    #make mini-heatmaps to show in the paper
    z_bprint_expression_mat_reduced<-z_bprint_expression_mat[rownames(z_bprint_expression_mat)%in%celltypes_to_keep,]
    
    #Reduce annotation length for bigger panels: 
    rownames(z_bprint_expression_mat_reduced)<- gsub("-positive","+",rownames(z_bprint_expression_mat_reduced))
    rownames(z_bprint_expression_mat_reduced)<- gsub("-negative","-",rownames(z_bprint_expression_mat_reduced))
    #pdf(file=paste0(output,"/",cohortname,"_",module_names[[i]],"_small.pdf"), width=1.5,height=3)
    #pheatmap(t(z_bprint_expression_mat_reduced),cluster_rows = F,cluster_cols = F,fontsize=8,border_color=NA,legend=F,show_rownames=F,main=module_names[[i]],breaks=seq(from=-3,to=3,length.out = 101),silent=F)#,filename="marker_expression_BLUEPRINT_celltypes.png", width = 7, height =7)
    #dev.off()
    blueprints[[i]] <- as.ggplot(as.ggplot(pheatmap((z_bprint_expression_mat_reduced),cluster_rows = F,cluster_cols = F,fontsize=8,border_color=NA,legend=F,show_rownames=T,show_colnames = F,main=module_names[[i]],breaks=seq(from=-3,to=3,length.out = 101),silent=F))+theme(text=element_text(family="Arial")))#,filename="marker_expression_BLUEPRINT_celltypes.png", width = 7, height =7)
    
  }
  names(blueprints) <- paste0(module_names)
  return(blueprints)}
blueprinter <- expression_plots_of_module_genes(expression_matrix=bprint_expression,
                                                gene_module_table=data.table(CEMiTool::module_genes(cemDiagnose)[order(CEMiTool::module_genes(cemDiagnose)$modules),]),
                                                output=output_location,
                                                celltypes_to_keep=c("mature neutrophil", "erythroblast", "conventional dendritic cell", "macrophage", "CD4-positive, alpha-beta T cell",
                                                                    "CD8-positive, alpha-beta T cell", "CD38-negative naive B cell","CD8-positive", "CD14-positive, CD16-negative classical monocyte"))

ggsave(ggpubr::ggarrange(plotlist=blueprinter,nrow=1),file=paste0(output_location,"/blueprints.pdf"),units="mm", height=100, width=300, dpi=300)

ggsave(ggpubr::ggarrange(plotlist=blueprinter,norw=5,ncol=1),file=paste0(output_location,"/blueprints.pdf"),units="mm", height=230, width=100, dpi=300)


blueprinter[1]$M1+theme(text=element_text(family="Arial", size = 12))


#######Pseudotime Age Relationship#######
cordf <- data.table(Age=cut(vsd$Alter_bei_Probennahme,breaks=seq(10,80,by=10)),
                    Severity=vsd$severity,
                    Pseudotime=ifelse(
                      rep(
                        median(
                          data.table(Pseudotime=traj$time, Severity=vsd$severity)[Severity %in% "Control",]$Pseudotime
                        )>=0.5, times=length(traj$time)
                      ),
                      as.vector(1-traj$time),
                      as.vector(traj$time)))[!is.na(Age),]
cor.test(method="spearman",x=cordf$Pseudotime,y=as.numeric(cordf$Age))


LMagepseudo <-ggplot(data=cordf,aes(x=Age,y=Pseudotime))+
  geom_jitter(aes(color=Severity))+ 
  geom_smooth(aes(x=as.numeric(Age)),method="lm",formula = "y~x")+
  geom_text(x = 4, y = 1.05, label = lm_eqn(data.table(x=as.numeric(cordf$Age),y=cordf$Pseudotime)), parse = T)+
  ylim(0,1.1)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  xlab("Age [years]")+
  theme_bw()+
  theme(legend.position = "bottom",text=element_text(family="Arial", size=18))
ggsave(plot=LMagepseudo, filename=paste0(output_location,"/",cohortname,"_",geneset,"_LMagepseudo.svg"),dpi=300,width = 160, height =120, units = "mm", device='svg')


#####UNUSED######
# ######Identifying wrongly labled gender in samples ########
# scoring <- fread("/home/sukmb465/Documents/Eike/Data/UC/scoringmetadata.csv")
# mscoring <- merge(colData(vsd), scoring, by="SampleID")
# coords <- selectedgeneumap(assay(vsd),seed=135, n_neighbors=15, labels=as.factor(mscoring$gender), colvec=c("#00AFBB", "#E7B800", "#FC4E07","#999999") ,text=FALSE,ellipse=T,legendtitle = "Phenotype",axistextsize = 20, legendtextsize = 20, dotsize = 5, genes = c("XIST", "USP9Y", "UTY") )
# mscoring <- data.frame(mscoring, X1=coords$data$X1, X2=coords$data$X2)
# mscoring[mscoring$X1 >= 0 & mscoring$gender %in% "Female",]
# c("DE09NGSUKBR201383","DE13NGSUKBR200729","DE13NGSUKBR200729")
# selectedgeneumap(assay(vsd),seed=135, n_neighbors=15, labels=as.factor(mscoring$gender), colvec=c("#00AFBB", "#E7B800", "#FC4E07","#999999") ,text=FALSE,ellipse=T,legendtitle = "Phenotype",axistextsize = 20, legendtextsize = 20, dotsize = 5, genes = c("XIST", "USP9Y", "UTY") )
  
# #####Geneset Foldchange and Module check######
# markergenes <- RFGS
# # data.table(genes=markergenes,
#                         FCsevcon=data.frame(row.names=ressevcon$row,ressevcon)[markergenes,]$log2FoldChange,
#                         PADJsevcon=data.frame(row.names=ressevcon$row,ressevcon)[markergenes,]$padj,
#                         FCsevmild=data.frame(row.names=ressevmil$row,ressevmil)[markergenes,]$log2FoldChange,
#                         PADJsevmild=data.frame(row.names=ressevmil$row,ressevmil)[markergenes,]$padj,
#                         FCmildcon=data.frame(row.names=resmilcon$row,resmilcon)[markergenes,]$log2FoldChange,
#                         PADJmildcon=data.frame(row.names=resmilcon$row,resmilcon)[markergenes,]$padj,
#                         CEMiModule=data.frame(row.names=module_genes(cemDiagnose)$genes,module_genes(cemDiagnose))[markergenes,]$modules)



sc <-results(discoverydds,  contrast=c("severity","severe","Control"),alpha=0.01)
mc <- results(discoverydds,  contrast=c("severity","mild","Control"),alpha=0.01)
sm <- results(discoverydds,  contrast=c("severity","severe","mild"),alpha=0.01)
length(Reduce(intersect,list(
  sc[sc$padj<0.01& abs(sc$log2FoldChange)>= 0.25,]$row,
  mc[mc$padj<0.01& abs(mc$log2FoldChange)>= 0.25,]$row,
  sm[sm$padj<0.01& abs(sm$log2FoldChange)>= 0.25,]$row
))
)
results(discoverydds,  contrast=c("Diagnose","UC","Control"),alpha=0.01,tidy=T)
summary(sm)
length(summary(sc[sc$padj<0.01,]$row))

dfg <- DESeqDataSetFromMatrix(countData=getrawcounts(), colData=metadata, design= ~PlateNr + Diagnose)
dfg <- DESeq(dfg)


ts <- vsdofallknownseverity(meta=metadata, countmatrix=getrawcounts(), batch="PlateNr")
ddsofallknownseverity3 <- function(meta=metadata, countmatrix=rawcounts, batch="PlateNr",...){
  #throwing out all Unknown severity patients:
  meta <- meta[!(meta$severity %in% "Unknown"),]
  countmatrix <- countmatrix[,meta$SampleID]
  cohortname <- "Discovery"
  colnames(meta)[colnames(meta)==batch] <- "batch"
  meta$batch <- as.factor(meta$batch)
  #filtering to 10 counts per gene in at least 10% of samples
  countmatrix <- countmatrix[as.vector(rowSums(countmatrix>10)>(dim(countmatrix)[2]/10)),]
  #deseq2 of mild/con, severe/con:
  dds <- DESeqDataSetFromMatrix(countData=countmatrix, colData=meta, design= ~batch + Diagnose)
  dds <- DESeq(dds)
  return(dds)}
tv <- ddsofallknownseverity3(meta=metadata, countmatrix=getrawcounts(), batch="PlateNr")


summary(results(tv,  contrast=c("Diagnose","UC","Control"),alpha=0.01))



###remission umap####
remissionfactor <- vsd$supervised
remissionfactor[(remissionfactor >= 3)] <- 3
remissionfactor <- as.factor(remissionfactor)
marginalerplot <-ggExtra::ggMarginal(
  ggplot(data=data.frame(rd2), aes(x = UMAP1, y= UMAP2))+
    geom_point(aes( color=remissionfactor),size=4)+
    geom_path(aes_string("Comp1", "Comp2"), data.frame(traj$path), 
              size = 0.5, alpha = 1)+ 
    theme_bw()+
    theme(legend.position="bottom",text=element_text(family="Arial", size=18), legend.box="horizontal")+
    labs(shape="Pseudotime", color="UC Supervised")+ 
    scale_color_manual(labels= c("Control", "Remission", "Active UC"),values = RColorBrewer::brewer.pal(length(unique(remissionfactor)),"Set2")),type="density", groupColour = TRUE)
ggsave(plot=marginalerplot, filename=paste0(output_location,"/",cohortname,"_",geneset,"_marginalPseudotimeUMAPREMISSION.svg"),dpi=300,width = 150, height = 150, units = "mm", device='svg')
###supplement figure 1####
UMAPscreen(mydata=assay(vsd),meta=fread("/home/sukmb465/Documents/Eike/Rprojects/UC/Projectmetadata/Projectmetadata/metadata.csv"))
