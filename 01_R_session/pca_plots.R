## Author: Florian Uellendahl-Werth & Eike Matthias Wacker
##
## IKMB Kiel, 2023
## Email: e.wacker at ikmb.uni-kiel.de



library(PCAtools)
library(tidyverse)
library(DESeq2)
library(ggExtra)
vsd_PCA <- pca(assay(vsd), metadata = data.frame(row.names=colData(vsd)$SampleID, as.data.table(colData(vsd)), Cohort_Diagnose=paste0(colData(vsd)$Cohort,"_",colData(vsd)$Diagnose)), removeVar = 0.1)

screeplot(vsd_PCA, axisLabSize = 18, titleLabSize = 22)
pre_z_plot <- biplot(vsd_PCA, 
                     showLoadings = F, 
                     colkey = c('0_Control' = 'red', '1_Control' = 'darkred', '0_PSC' = 'lightgreen', '0_PSCUC' = 'forestgreen', '1_UC' = 'purple'),
                     colLegendTitle = 'Cohort_Diagnose',
                     legendPosition = 'top', 
                     legendLabSize = 16, 
                     legendIconSize = 8.0,
                     lab = NULL,
                     colby = "Cohort_Diagnose",
                     alpha= 0.5,
                     labSize = 5, 
                     #shape = "Cohort",
                     pointSize = 5, 
                     sizeLoadingsNames = 5)
legend_saved <- get_legend(pre_z_plot)

pre_z_plot <- ggMarginal(pre_z_plot+theme(legend.position = "none"), type="boxplot",groupFill=T,groupColour = T)


zscore_matrix <- t(as.matrix(merged_RF[,-c("rn","Diagnose","Cohort")]))
colnames(zscore_matrix) <- merged_RF$rn
rownames(zscore_matrix) <- colnames(merged_RF[,-c("rn","Diagnose","Cohort")])  

# zs_metadata <- data.table(colData(vsd), Cohort_Diagnose=paste0(colData(vsd)$Cohort,"_",colData(vsd)$Diagnose)
zs_PCA <- pca(zscore_matrix, metadata = data.frame(row.names=colData(vsd)$SampleID, as.data.table(colData(vsd)), Cohort_Diagnose=paste0(colData(vsd)$Cohort,"_",colData(vsd)$Diagnose)), removeVar = 0.1)

screeplot(zs_PCA, axisLabSize = 18, titleLabSize = 22)

pairsplot(zs_PCA)

after_z_plot <- biplot(zs_PCA, 
                       showLoadings = F, 
                       colkey = c('0_Control' = 'red', '1_Control' = 'darkred', '0_PSC' = 'lightgreen', '0_PSCUC' = 'forestgreen', '1_UC' = 'purple'),
                       colLegendTitle = 'Cohort_Diagnose',
                       legendPosition = 'top', 
                       legendLabSize = 16, 
                       legendIconSize = 8.0,
                       lab = NULL,
                       colby = "Cohort_Diagnose",
                       
                       #set_alpha= 0.5,
                       #shape = "Cohort",
                       labSize = 5, 
                       pointSize = 5, 
                       sizeLoadingsNames = 5)#+geom_point(aes(fill="Diagnose"))#+geom_point()



after_z_plot <- ggMarginal(after_z_plot+theme(legend.position = "none"), type="boxplot",groupFill=T,groupColour = T)





# 
# biplot(zs_PCA, showLoadings = TRUE,
# #       legendPosition = "top",
# #       lab = merged_RF$Diagnose, 
#        pointSize = 2, 
#        sizeLoadingsNames = 5)
# # ?biplot
# # 
# pairsplot(zs_PCA,
#           components = getComponents(zs_PCA, c(1:10)),
#           triangle = TRUE, trianglelabSize = 12,
#           hline = 0, vline = 0,
#           pointSize = 0.4,
#           gridlines.major = FALSE, gridlines.minor = FALSE,
#           colby = 'Grade',
#           title = 'Pairs plot', plotaxes = FALSE,
#           margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))




pca_pre_and_after_z <- ggpubr::ggarrange(cowplot::plot_grid(pre_z_plot,after_z_plot, 
                                                            nrow=1,
                                                            ncol=2,
                                                            labels = c("a)","b)"),
                                                            label_y = .97,
                                                            label_size=24), 
                                         legend.grob = legend_saved, 
                                         #labels = c("a)","b)"),
                                         #label.y = .97,
                                         legend = "bottom",
                                         font.label = list(size = 24, color = "black", face = "bold", family = NULL))

pca_pre_and_after_z

ggsave(pca_pre_and_after_z, filename = "output/pca_plots.svg", height = 7.5, width = 16, device = "svg",bg="white")
ggsave(pca_pre_and_after_z, filename = "output/pca_plots.png", height = 7.5, width = 16, device = "png",bg="white")




after_z_plot2 <- biplot(zs_PCA,
                        showLoadings = F,
                        colkey = c('Control' = 'red', 'PSC' = 'lightgreen', 'PSCUC' = 'forestgreen', 'UC' = 'purple'),
                        colLegendTitle = 'Diagnose',
                        legendPosition = 'top',
                        legendLabSize = 16,
                        legendIconSize = 8.0,
                        lab = NULL,
                        colby = "Diagnose",
                        shape = "Cohort",
                        labSize = 3,
                        pointSize = 3,
                        sizeLoadingsNames = 5)
# 
# ggsave(after_z_plot2, filename = "output/pca_plot_afterbatchcorrection.svg", height = 7.5, width = 7.5, device = "svg",bg="white")



highvariance_genes <- keep_high_variance_genes(matrix_input = as.matrix(merged_RF[,-c("rn","Diagnose","Cohort")]),
                                               remove_var = 0.9)

genes_modules <- cem_merged_RF@module[order(cem_merged_RF@module$genes),]
genes_modules[genes_modules$genes %in% highvariance_genes,]
# cem_merged_RF@module[order(cem_merged_RF@module$genes)%in% highvariance_genes,]#[cem_merged_RF@module$genes %in% highvariance_genes,]
# ha = HeatmapAnnotation(Modules = genes_modules[genes_modules$genes %in% highvariance_genes,]$modules)#anno_block(gp = genes_modules[genes_modules$genes %in% highvariance_genes,]$modules)
#)#genes_modules[genes_modules$genes %in% highvariance_genes,]$modules)

heat <- Heatmap(as.matrix(merged_RF[,-c("rn","Diagnose","Cohort")])[,highvariance_genes[order(highvariance_genes)]],
                # top_annotation = ha, 
                column_split = genes_modules[genes_modules$genes %in% highvariance_genes,]$modules,
                row_split = merged_RF$Diagnose,
                name = "Genes with top 1% variance",
                row_labels = NULL,
                show_column_names = FALSE,
                show_heatmap_legend = FALSE
                # width = unit(8, "cm"),
                # height = unit(8, "cm")
)
draw(heat)
