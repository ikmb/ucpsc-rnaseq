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
library(PCAtools)
library(cowplot)

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

merged_metadata_yearsbinned <- rbind(PSC_metadata[,c("SampleID","Diagnose","PlateNr", "yearssincediagnose")],UC_metadata[,c("SampleID","Diagnose","PlateNr", "yearssincediagnose")])
merged_metadata_yearsbinned$Cohort <- factor(c(rep(0,length(PSC_metadata$SampleID)),rep(1,length(UC_metadata$SampleID))))

merged_metadata_yearsbinned$yearssincebinned <- case_when(merged_metadata_yearsbinned$yearssincediagnose == 0 ~ "newly",
          merged_metadata_yearsbinned$yearssincediagnose > 0 & merged_metadata_yearsbinned$yearssincediagnose < 4 ~ "newly",
          merged_metadata_yearsbinned$yearssincediagnose >= 4 & merged_metadata_yearsbinned$yearssincediagnose < 10 ~ "medium",
          merged_metadata_yearsbinned$yearssincediagnose >= 10 & merged_metadata_yearsbinned$yearssincediagnose < 30 ~ "long",
          merged_metadata_yearsbinned$yearssincediagnose >= 30 & merged_metadata_yearsbinned$yearssincediagnose < 40 ~ "long",
          merged_metadata_yearsbinned$yearssincediagnose >= 50 & merged_metadata_yearsbinned$yearssincediagnose < 60 ~ "SixthDecade",
          is.na(merged_metadata_yearsbinned$yearssincediagnose) ~ "Missing")

dds_yearssincediagnose_binned <- DESeqDataSetFromMatrix(countData=merged_rawcounts, colData=merged_metadata_yearsbinned, design= ~Diagnose + PlateNr + yearssincebinned)
dds_yearssincediagnose_binned <- DESeq(dds_yearssincediagnose_binned)



# results(dds_yearssincediagnose_binned, contrast = c("yearssincebinned", "newly", "long"), tidy = T) %>% arrange(padj)
# results(dds_yearssincediagnose_binned, contrast = c("yearssincebinned", "newly", "medium"), tidy = T) %>% arrange(padj)
# results(dds_yearssincediagnose_binned, contrast = c("yearssincebinned", "long", "medium"), tidy = T) %>% arrange(padj)
results(dds_yearssincediagnose_binned, contrast = c("yearssincebinned", "newly", "long"), tidy = T) %>% arrange(padj)

# #lookup significant genes:
# ggplot(as.data.table(t(assay(dds_yearssincediagnose_binned))))+geom_boxplot(aes(x=merged_metadata_yearsbinned$Diagnose, y= TMEM216))+geom_jitter(aes(x=merged_metadata_yearsbinned$Diagnose,y=TMEM216))+facet_wrap(~merged_metadata_yearsbinned$Cohort)





# zscore_matrix <- t(as.matrix(merged_RF[,-c("rn","Diagnose","Cohort")]))
# colnames(zscore_matrix) <- merged_RF$rn
# rownames(zscore_matrix) <- colnames(merged_RF[,-c("rn","Diagnose","Cohort")])  

zscore_matrix <- t(as.matrix(merged_RF[Diagnose %in% c(1,2),-c("rn","Diagnose","Cohort")]))
colnames(zscore_matrix) <- merged_RF[Diagnose %in% c(1,2),]$rn
rownames(zscore_matrix) <- colnames(merged_RF[,-c("rn","Diagnose","Cohort")])  

zs_PCA <- pca(zscore_matrix, metadata = colData(vsd)[colData(vsd)$Diagnose %in% c("PSC","PSCUC"),], removeVar = 0.1)

pca_timesincediagnose_plot <- biplot(zs_PCA, 
       showLoadings = F, 
       colkey = c('newly' = 'red', 'medium' = 'lightgreen', 'long' = 'forestgreen', 'Missing' = 'purple'),
       colLegendTitle = 'Binned time since PSC diagnose',
       legendPosition = 'top', 
       legendLabSize = 16, 
       legendIconSize = 8.0,
       lab = NULL,
       colby = "yearssincebinned",
       # shape = "yearssincebinned",
       labSize = 5, 
       pointSize = 5, 
       sizeLoadingsNames = 5,
       # encircle = T
       ellipse = T
       )+
  guides(color=guide_legend(nrow=2,byrow=T, title.position = "top"))

pca_timesincediagnose_plot+ theme(legend.position = "none")

pca_timesincediagnose_plot

# 
# ggsave(pca_timesincediagnose_plot, filename = "output/pca_plot_psc_timesincediagnose.svg", height = 7.5, width = 7.5, device = "svg",bg="white")
# ggsave(pca_timesincediagnose_plot, filename = "output/pca_plot_psc_timesincediagnose.png", height = 7.5, width = 7.5, device = "png",bg="white")



arranged_pca_psc_plot <- cowplot::plot_grid(
  pca_timesincediagnose_plot+theme(legend.position = "none")+
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))),
  biplot(zs_PCA, 
         showLoadings = F, 
         colkey = c('newly' = 'red', 'medium' = 'lightgreen', 'long' = 'forestgreen', 'Missing' = 'purple'),
         colLegendTitle = 'Binned time since PSC diagnose',
         legendPosition = 'top', 
         legendLabSize = 16, 
         legendIconSize = 8.0,
         lab = NULL,
         colby = "yearssincebinned",
         # shape = "yearssincebinned",
         labSize = 5, 
         pointSize = 5, 
         sizeLoadingsNames = 5,
         # encircle = T
         ellipse = T,
         x = "PC3",
         y= "PC4"
  )+theme(legend.position = "none")+
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))),
  
  biplot(zs_PCA, 
         showLoadings = F, 
         colkey = c('newly' = 'red', 'medium' = 'lightgreen', 'long' = 'forestgreen', 'Missing' = 'purple'),
         colLegendTitle = 'Binned time since PSC diagnose',
         legendPosition = 'top', 
         legendLabSize = 16, 
         legendIconSize = 8.0,
         lab = NULL,
         colby = "yearssincebinned",
         # shape = "yearssincebinned",
         labSize = 5, 
         pointSize = 5, 
         sizeLoadingsNames = 5,
         # encircle = T
         ellipse = T,
         x = "PC5",
         y= "PC6"
  )+theme(legend.position = "none")+
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))),
  
  biplot(zs_PCA, 
         showLoadings = F, 
         colkey = c('newly' = 'red', 'medium' = 'lightgreen', 'long' = 'forestgreen', 'Missing' = 'purple'),
         colLegendTitle = 'Binned time since PSC diagnose',
         legendPosition = 'top', 
         legendLabSize = 16, 
         legendIconSize = 8.0,
         lab = NULL,
         colby = "yearssincebinned",
         # shape = "yearssincebinned",
         labSize = 5, 
         pointSize = 5, 
         sizeLoadingsNames = 5,
         # encircle = T
         ellipse = T,
         x = "PC7",
         y= "PC8"
  )+theme(legend.position = "none")+
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))),
  
  biplot(zs_PCA, 
           showLoadings = F, 
           colkey = c('newly' = 'red', 'medium' = 'lightgreen', 'long' = 'forestgreen', 'Missing' = 'purple'),
           colLegendTitle = 'Binned time since PSC diagnose',
           legendPosition = 'top', 
           legendLabSize = 16, 
           legendIconSize = 8.0,
           lab = NULL,
           colby = "yearssincebinned",
           # shape = "yearssincebinned",
           labSize = 5, 
           pointSize = 5, 
           sizeLoadingsNames = 5,
           # encircle = T
           ellipse = T,
           x = "PC9",
           y= "PC10"
  )+theme(legend.position = "none")+
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))),
  
  get_legend(pca_timesincediagnose_plot)
)
arranged_pca_psc_plot

ggsave(arranged_pca_psc_plot, filename = "output/arrange_pca_plot_psc_timesincediagnose.svg", height = 14, width = 16, device = "svg",bg="white")
ggsave(arranged_pca_psc_plot, filename = "output/arrange_pca_plot_psc_timesincediagnose.png", height = 14, width = 16, device = "png",bg="white")
