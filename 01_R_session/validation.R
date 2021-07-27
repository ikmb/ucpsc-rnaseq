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
library(ggridges)
library(ggplotify)
library(VennDiagram)
library(extrafont)
library(SCORPIUS)
library(reshape2)
library(plyr)
library(ggrepel)
source("functions.R")

output_location <- paste0("output")
validationcohorts <- list(Ostrowski ="UCAI Based\nDisease Severity",
                          Mo_UC="Diagnose",
                          Planell="Endoscopic\nMayo Score")
# marker <- list(RFGS=as.vector(unlist(fread("/home/sukmb465/Desktop/UCRNA/RFGS.txt", header=F))),
#                DEGS=as.vector(unlist(fread("/home/sukmb465/Desktop/UCRNA/DE23.txt", header=F))),
#                MMGS=as.vector(unlist(fread("/home/sukmb465/Desktop/UCRNA/MMGS.txt", header=F))))

cortestlist <- list()
k <- 1
for (cohort in names(validationcohorts)){
  vsd <- getCohortsVSD(cohortname = cohort)
  cohortname <- cohort
  severityscale <- validationcohorts[[cohort]]
  
    if (cohortname=="Planell"){
      vsd$severity <- vsd$endoscopic_mayo
      vsd$supervised <- vsd$endoscopic_mayo
      vsd$severity <- ordered(replace(vsd$supervised, vsd$supervised%in%NA,"Control"), levels=c("Control","0","1","2","3"))
      #plotcolors <- c("#00AFBB", rep("#E7B800",times=2), rep("#FC4E07", times=2))
      plotcolors <- RColorBrewer::brewer.pal(length(levels(vsd$severity)),"Set2")
      MyPalette <-setNames(plotcolors,as.list(levels(vsd$severity)))


    }else if (cohortname=="Ostrowski"){
      vsd$severity[is.na(vsd$severity)] <- "Control"  
      vsd$severity <- ordered(vsd$severity, levels=c("Control","0","1","2","3"))
      #plotcolors <- c("#00AFBB", rep("#E7B800",times=2), rep("#FC4E07", times=3))
      plotcolors <- RColorBrewer::brewer.pal(length(levels(vsd$severity)),"Set2")
      MyPalette <-setNames(plotcolors,as.list(levels(vsd$severity)))
    }else if (cohortname=="Mo_UC"){
      vsd$severity <- vsd$Diagnose
      vsd$severity[is.na(vsd$severity)] <- "Control"  
      #plotcolors <- c("#00AFBB", rep("#E7B800",times=0), rep("#FC4E07", times=1))
      plotcolors <- RColorBrewer::brewer.pal(length(levels(vsd$severity)),"Set2")
      MyPalette <-setNames(plotcolors,as.list(levels(vsd$severity)))
      }
  for (geneset in names(marker)) {
    print(geneset)
    markergenes <- marker[[geneset]]
    print(marker[[geneset]])
    
    #####SCORPIUS Trajectory#######
    rd2 <- uwot::umap(t(assay(vsd)[rownames(assay(vsd)) %in% markergenes,]))
    colnames(rd2) <- c('UMAP1', 'UMAP2')
    traj <- SCORPIUS::infer_trajectory(rd2)
    #draw_trajectory_plot(rd2, vsd$Diagnose, traj$path, contour = TRUE)

    
    #####Plots based on Trajectory#######
    marginalPseudotimeUMAP <- ggExtra::ggMarginal(
      ggplot(data=data.frame(rd2), aes(x = UMAP1, y= UMAP2))+
        geom_point(aes( color=as.factor(vsd$severity)),size=4)+
        geom_path(aes_string("Comp1", "Comp2"), data.frame(traj$path), 
                  size = 0.5, alpha = 1)+ theme_bw()+
        theme(legend.position="bottom",text=element_text(family="Arial", size=18), legend.box="horizontal")+labs(shape="Pseudotime", color=severityscale)+ 
        scale_color_manual(values = MyPalette),type="density", groupColour = TRUE)
    ggsave(plot=marginalPseudotimeUMAP, filename=paste0(output_location,"/",cohortname,"_",geneset,"_marginalPseudotimeUMAP.svg"),dpi=300,width = 150, height = 150, units = "mm", device='svg')
    
    dg <- data.table(Pseudotime=ifelse(
      rep(
        median(
          data.table(Pseudotime=traj$time, Severity=vsd$severity)[Severity %in% "Control",]$Pseudotime
        )>=0.5, times=length(traj$time)
      ),
      as.vector(1-traj$time),
      as.vector(traj$time)),
      severity=vsd$severity)
    # try(dg$severity <- revalue(dg$severity, c("0"="Control")))
    
    PseudotimeSupervisedCorrelation <- ggplot(data = dg, aes(x=Pseudotime, y = severity)) + 
      geom_density_ridges(aes(fill= severity,y=severity,point_fill=severity,point_color=severity),jittered_points=TRUE,alpha=0.5) +#,stat="binline", bins=15
      theme_bw() + 
      labs(x = "Pseudotime", y = severityscale) +
      geom_smooth(aes(Pseudotime,as.numeric(severity) ), method = "lm",formula="y~x")+
      scale_fill_manual(values=MyPalette)+
      theme(legend.position="none", text=element_text(family="Arial", size=20))+
      xlim(-0.1,1.1)+
      scale_discrete_manual(aesthetics = "point_color", values = MyPalette)+ scale_y_discrete(drop=FALSE)
    ggsave(plot=PseudotimeSupervisedCorrelation, filename=paste0(output_location,"/",cohortname,"_",geneset,"_PseudotimeSupervisedCorrelation.svg"),dpi=300,width = 120, height =120, units = "mm", device='svg')    
    
    df <- data.table(data.table(t(assay(vsd)[rownames(assay(vsd)) %in% markergenes,])),Pseudotime=traj$time, Diagnose=vsd$Diagnose)
    if(median(df[df$Diagnose %in% "Control",]$Pseudotime)>=0.5){df$Pseudotime <- 1-df$Pseudotime}
    melted <- melt(df[,colnames(df) %in% c("Pseudotime", markergenes, "Diagnose"), with=FALSE], id.vars=c('Pseudotime','Diagnose'))
    expressiontopseudotime <-  ggplot(melted,aes(x=Pseudotime,group='Diagnose'))+
      geom_point(aes(color=Diagnose,y=value))+
      geom_quantile(aes(y=value),colour="black", formula="y~x")+
      facet_wrap(~variable,ncol = ifelse(length(unique(melted$variable))>20,3,2),scales = 'free')+
      theme_bw()+ 
      scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
    expressiontopseudotime
    ggsave(plot=expressiontopseudotime, filename=paste0(output_location,"/",cohortname,"_",geneset,"_expressiontopseudotimeUMAP.svg"),dpi="retina",width = ifelse(length(unique(expressiontopseudotime$data$variable))<20,160*2,160*3), height = 55*round(length(markergenes)/ifelse(length(unique(expressiontopseudotime$data$variable))<20,2,3)), units = "mm", device='svg')
    
    cortestlist[[k]] <- cor.test(x=as.numeric(dg$severity),
                                 y=as.numeric(dg$Pseudotime),method=c("spearman"))
    names(cortestlist)[[k]] <- paste0(cohortname,"_",geneset )
    k<-k+1
    }
}

# cor.test(x=as.numeric(ordered(replace(vsd$severity, vsd$severity%in%NA,"Control"), levels=c("Control","0","1","2","3"))),
#          y=1-traj$time,method=c("spearman"))
