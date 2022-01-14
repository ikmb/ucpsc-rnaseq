#Needed libraries:
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(uwot)
library(gridExtra)
library(svglite)
library(DESeq2)

#At startup please run:
#####Function#####
UMAPscreen <- function(mydata=assay(vsd), input=rownames(mydata), n_neighbors=15, meta=metadata[fread("/home/sukmb465/Documents/Eike/Rprojects/UC/Cohorts/metadata.csv")[,!c("Diagnose","PlateNr")], on="SampleID"], seed = 135, 
                       axistextsize = 7, legendtextsize = 9, dotsize = 1, textlabelsize = 3, 
                       legendtitle = "Group", controlscale = FALSE){
  if (seed != FALSE) {
    set.seed(seed)
  }
  

  mydata <- as.matrix(mydata, rownames="genes")
  
  #ordering factors and combining duplicates
  meta$CRPcategory <- as.factor(meta$CRPcategory)
  meta$finding <- as.factor(meta$finding)
  meta$Stuhlfrqz <- as.factor(meta$Stuhlfrqz)
  meta$BlutimStuhl <- as.factor(meta$BlutimStuhl)
  

  meta$CRPcategory  <- ordered(meta$CRPcategory , levels = c("Control","normal", "high", "missing"))
  meta$finding   <- ordered(meta$finding  , levels = c("Control", "remission", "chronically active" ,"acute flare",  "other"    ,"missing"))
  meta$Stuhlfrqz   <- ordered(meta$Stuhlfrqz  , levels = c("Control","1: normal", "2: 1-2 more than normal" ,"3: 3-4 more than normal" ,"4: >5 more than normal" ,"missing" ))
  meta$BlutimStuhl<- ordered(meta$BlutimStuhl, levels = c("Control", "none"   , "little(<50%)","much(>50%)" , "heavy(onlyblood)" ))
  
  
  
  scores <- data.frame(uwot::umap(t(as.matrix(mydata[rownames(mydata[rownames(mydata) %in% input,]), ])), n_neighbors=n_neighbors))
  
  diagnoseplot <- ggplot(data = scores, aes(x = X1, y = X2)) + 
    geom_point(aes(colour = meta$Diagnose), size = dotsize) + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                      colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                    colour = "black"), axis.title.x = element_text(size = axistextsize), 
                       axis.title.y = element_text(size = axistextsize), 
                       legend.title = element_text(size = legendtextsize), 
                       legend.text = element_text(size = legendtextsize)) + 
    labs(colour = "Diagnose") + scale_colour_manual(values = as.vector(RColorBrewer::brewer.pal(8,"Set2")))+stat_ellipse(aes(colour = meta$Diagnose), level=0.95,size=dotsize, type="euclid")
  #c("lightblue", heat.colors(length(levels(as.factor(meta$Diagnose)))-1)
  genderplot <- ggplot(data = scores, aes(x = X1, y = X2)) + 
    geom_point(aes(colour = meta$gender), size = dotsize) + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                      colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                    colour = "black"), axis.title.x = element_text(size = axistextsize), 
                       axis.title.y = element_text(size = axistextsize), 
                       legend.title = element_text(size = legendtextsize), 
                       legend.text = element_text(size = legendtextsize)) + 
    labs(colour = "Gender") + scale_colour_manual(values = as.vector(RColorBrewer::brewer.pal(8,"Set2")))+stat_ellipse(aes(colour = meta$gender), level=0.95,size=dotsize, type="euclid")
  
  findingplot <- ggplot(data = scores[!(meta$finding %in% c("missing", "other")),], aes(x = X1, y = X2)) + 
    geom_point(aes(colour = meta[!(meta$finding %in% c("missing", "other")),]$finding), size = dotsize) + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                      colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                    colour = "black"), axis.title.x = element_text(size = axistextsize), 
                       axis.title.y = element_text(size = axistextsize), 
                       legend.title = element_text(size = legendtextsize), 
                       legend.text = element_text(size = legendtextsize)) + 
    labs(colour = "General Medical Finding") + scale_colour_manual(values = as.vector(RColorBrewer::brewer.pal(8,"Set2")))+stat_ellipse(aes(colour = meta[!(meta$finding %in% c("missing", "other")),]$finding), level=0.95,size=dotsize, type="euclid")
  
  Stuhlfrqzplot <- ggplot(data = scores[!(meta$Stuhlfrqz %in% "missing"),], aes(x = X1, y = X2)) + 
    geom_point(aes(colour = meta[!(meta$Stuhlfrqz %in% "missing"),]$Stuhlfrqz), size = dotsize) + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                      colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                    colour = "black"), axis.title.x = element_text(size = axistextsize), 
                       axis.title.y = element_text(size = axistextsize), 
                       legend.title = element_text(size = legendtextsize), 
                       legend.text = element_text(size = legendtextsize)) + 
    labs(colour = "Bowel Movements") + scale_colour_manual(values = as.vector(RColorBrewer::brewer.pal(8,"Set2")))+stat_ellipse(aes(colour = meta[!(meta$Stuhlfrqz %in% "missing"),]$Stuhlfrqz), level=0.95,size=dotsize, type="euclid")
  BlutimStuhlplot <- ggplot(data = scores, aes(x = X1, y = X2)) + 
    geom_point(aes(colour = meta$BlutimStuhl), size = dotsize) + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                      colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                    colour = "black"), axis.title.x = element_text(size = axistextsize), 
                       axis.title.y = element_text(size = axistextsize), 
                       legend.title = element_text(size = legendtextsize), 
                       legend.text = element_text(size = legendtextsize)) + 
    labs(colour = "Bloody Stool") + scale_colour_manual(values = as.vector(RColorBrewer::brewer.pal(8,"Set2")))+stat_ellipse(aes(colour = meta[!(meta$BlutimStuhl %in% "missing"),]$BlutimStuhl), level=0.95,size=dotsize, type="euclid")
  
  CRPcategoryplot <- ggplot(data = scores[!(meta$CRPcategory %in% "missing"),], aes(x = X1, y = X2)) + 
    geom_point(aes(colour = meta[!(meta$CRPcategory %in% "missing"),]$CRPcategory), size = dotsize) + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize, 
                                                                                      colour = "black"), axis.text.x = element_text(size = axistextsize, 
                                                                                                                                    colour = "black"), axis.title.x = element_text(size = axistextsize), 
                       axis.title.y = element_text(size = axistextsize), 
                       legend.title = element_text(size = legendtextsize), 
                       legend.text = element_text(size = legendtextsize)) + 
    labs(colour = "CRP level") + scale_colour_manual(values = as.vector(RColorBrewer::brewer.pal(8,"Set2")))+stat_ellipse(aes(colour = meta[!(meta$CRPcategory %in% "missing"),]$CRPcategory), level=0.95,size=dotsize, type="euclid")
  
  
  umapplot <- grid.arrange(diagnoseplot, genderplot, findingplot, Stuhlfrqzplot, BlutimStuhlplot, CRPcategoryplot)
  ggsave(file="umapplots.svg", plot=umapplot, width=178, height = 100,unit="mm")
  return(umapplot)
}


#####Script#####
#set working directory, please adjust:
# setwd("/home/ewacker/Desktop/UMAPselection/")

##read in a gene list, with a single gene per row, no other separations, no header:
#300 most significant DE genes:
#inputgenes <- fread("inputgenes.csv", header=F)
inputgenes <- markergenes
#your gene list:
#inputgenes <- fread(file.choose(), header=F)
#or manually into this script:
#inputgenes <- c("CD177", "S100A12","S100A9")

#Creates vector of read in data.frame structured genes:
inputgenes <- as.vector(unlist(inputgenes))


dds1 <- DESeqDataSetFromMatrix(countData=assay(discoverydds), colData=colData(discoverydds), design= ~batch + Diagnose)
dds1 <- DESeq(dds1)
resUCCON<- results(dds1,  contrast=c("Diagnose","UC","Control"),tidy = T)
vsd1 <- vst(dds1, blind=FALSE)

meta1<- metadata[fread("/home/sukmb465/Documents/Eike/Rprojects/UC/Cohorts/metadata.csv")[,!c("Diagnose","PlateNr")], on="SampleID"]
meta1$Stuhlfrqz<-factor(meta1$Stuhlfrqz)
meta1$Stuhlfrqz<- revalue(meta1$Stuhlfrqz,c("0: normal"="1: normal","1: 1-2 more than normal"="2: 1-2 more than normal","2: 3-4 more than normal"="3: 3-4 more than normal","3: >5 more than normal"="4: >5 more than normal"))
#Function that will create UMAP plots, automatically saved in directoy as umapplots.svg:
UMAPscreen(input=data.table(resUCCON[order(resUCCON$padj),])[1:300,row],meta=meta1, n_neighbors=15, dotsize= .8, legendtextsize = 7)


#metadata[fread("/home/sukmb465/Documents/Eike/Rprojects/UC/Cohorts/metadata.csv")[,!c("Diagnose","PlateNr")], on="SampleID"]

