setwd(paste0("/ssdpool/ewacker/RNAseq","/01_R_session"))
# setwd(paste0("/home/sukmb465/Documents/Eike/Projects/uc-rnaseq/01_R_session"))

library(data.table)
library(plyr)
library(tidyverse)
library(tidymodels)
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
library(ggpubr)
library(EnhancedVolcano)
source("functions.R")

##Lookup single gene expressions:
#ggplot(test_vst_counts)+geom_jitter(aes(x=merged_metadata$Diagnose,y=ZFP36L2))+facet_wrap(~merged_metadata$Cohort)
size.rel = 1L
# Volcano Plots ####
#PSC vs CON
# res_PSC_untidy <- results(dds,  contrast=c("Diagnose","PSC","Control"))
#resLFC_PSC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

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

# EnhancedVolcano_function(res_PSC_untidy)



volcano_UCvsCON <- EnhancedVolcano_function(res_UC_untidy, plottitle="UC vs CON",plotnumber=element_blank())+ 
  theme(legend.position = "none", 
        plot.margin = margin(5.5,20,5.5,27, "pt"),
        axis.title.x = element_blank(),        
        plot.tag.position = c(-.01, .97),
        plot.tag = element_text(size = rel(1.5 * size.rel))) +
  labs(tag="a)")# axis.title=element_text(size=12,face="bold"),
volcano_PSCvsCON <- EnhancedVolcano_function(res_PSC_CON_untidy, plottitle="PSC vs CON",plotnumber=element_blank())+  
  theme(legend.position = "none", 
                                                                                                                            plot.margin = margin(5.5,20,5.5,27, "pt"),
                                                                                                                            axis.title.x = element_blank(),        
                                                                                                                            plot.tag.position = c(-.01, .97),
                                                                                                                            plot.tag = element_text(size = rel(1.5 * size.rel)))+
  labs(tag="b)")
volcano_UCvsPSC <- EnhancedVolcano_function(res_UC_PSC_untidy, plottitle="UC vs PSC",plotnumber=element_blank())+ 
  theme(legend.position = "none", 
        plot.margin = margin(5.5,20,5.5,27, "pt"),
        axis.title.x = element_blank(),        
        plot.tag.position = c(-.01, .97),
        plot.tag = element_text(size = rel(1.5 * size.rel)))+
  labs(tag="c)")
volcano_PSCUCvsPSC <- EnhancedVolcano_function(res_PSCUC_PSC_untidy, plottitle="PSCUC vs PSC_noUC",plotnumber=element_blank())+
  labs(tag="d)")+  
  theme(legend.position = "right", 
        plot.margin = margin(5.5,20,5.5,27, "pt"),
        #axis.title.x = element_blank(),
        plot.tag.position = c(-.01, .97),
        plot.tag = element_text(size = rel(1.5 * size.rel)))

arranged_volcano <- ggarrange( ncol = 1,nrow=4,
  volcano_UCvsCON,volcano_PSCvsCON,volcano_UCvsPSC,volcano_PSCUCvsPSC)
arranged_volcano
plot_width=8
plot_height=18
plotname <- "output/VolcanoplotRNAseq"
# paste0(plotname, ".pdf")
ggsave(arranged_volcano, filename = paste0(plotname, ".pdf"),device = "pdf",width = plot_width, height = plot_height)
ggsave(arranged_volcano, filename = paste0(plotname, ".png"),device = "png",width = plot_width, height = plot_height)
ggsave(arranged_volcano, filename = paste0(plotname, ".svg"),device = "svg",width = plot_width, height = plot_height)



#TopGO enrichments for sign genesets from volcano pvalue and log2 FC:
sign_UC <- res_UC_tidy[res_UC_tidy$padj < 0.0001 & abs(res_UC_tidy$log2FoldChange)>1,]$row
sign_UC_PSC <-res_UC_PSC_tidy[res_UC_PSC_tidy$padj < 0.0001 & abs(res_UC_PSC_tidy$log2FoldChange)>1,]$row
sign_PSC_CON <-res_PSC_CON_tidy[res_PSC_CON_tidy$padj < 0.0001 & abs(res_PSC_CON_tidy$log2FoldChange)>1,]$row
sign_PSCUC_PSC <-res_PSCUC_PSC_tidy[res_PSCUC_PSC_tidy$padj < 0.0001 & abs(res_PSCUC_PSC_tidy$log2FoldChange)>1,]$row #is empty

sign_UC <- sign_UC[!is.na(sign_UC)]
sign_UC_PSC <- sign_UC_PSC[!is.na(sign_UC_PSC)]
sign_PSC_CON <- sign_PSC_CON[!is.na(sign_PSC_CON)]


topGO_enrichment(genelist = sign_UC, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                 DESeq_result = res, match_by_expression = TRUE, ontology_type = "MF")
topGO_enrichment(genelist = sign_UC_PSC, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                 DESeq_result = res, match_by_expression = TRUE, ontology_type = "MF")
topGO_enrichment(genelist = sign_PSC_CON, statistical_test = "Fisher", #gene_background = PSCUC_UC_Diagnose_results$variable_importance$feature
                 DESeq_result = res, match_by_expression = TRUE, ontology_type = "MF")
# test1 <- assays(dds)[["cooks"]]
# mcols(dds)$maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
# plot(mcols(dds)$baseMean, mcols(dds)$maxCooks)
# # this requires you not filter or order the 'res' object
# stopifnot(all.equal(rownames(dds), rownames(res_PSC_untidy)))
# plot(res_PSCUC_PSC_untidy$log2FoldChange, mcols(dds)$maxCooks)


Planell_CON_UC_Diagnose_results
Mo_CON_UC_Diagnose_results
Ostrowski_CON_UC_Diagnose_results

Ostrowski_CON_PSC_Diagnose_results
PSCUCPSC_CON_Diagnose_results
PSC_CONTROL_Diagnose_results

#Our Dataset:
CON_UC_Diagnose_results
PSC_UC_COHORT_results
PSC_UC_Diagnose_results
PSCUC_UC_Diagnose_results
PSC_PSCUC_Diagnose_results

PSCUCPSC_UC_Diagnose_results


pROC::auc(Ostrowski_CON_PSC_Diagnose_results$pROC_object)
pROC::ggroc(Ostrowski_CON_PSC_Diagnose_results$pROC_object)
pROC::ci.auc(Ostrowski_CON_PSC_Diagnose_results$pROC_object)



