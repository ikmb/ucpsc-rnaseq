roclist <- list(
Planell_CON_UC_Diagnose_results$pROC_object,
Mo_CON_UC_Diagnose_results$pROC_object,
Ostrowski_CON_UC_Diagnose_results$pROC_object,
Ostrowski_CON_PSC_Diagnose_results$pROC_object,


#Our Dataset:
CON_UC_Diagnose_results$pROC_object,

PSC_UC_COHORT_results$pROC_object,
PSC_UC_Diagnose_results$pROC_object,
PSCUC_UC_Diagnose_results$pROC_object,
PSC_PSCUC_Diagnose_results$pROC_object,
PSC_CONTROL_Diagnose_results$pROC_object,
PSCUCPSC_UC_Diagnose_results$pROC_object,
PSCUCPSC_CON_Diagnose_results$pROC_object,
PSCUC_CON_Diagnose_results$pROC_object)

roclist_names <- c(
  #"NULL",
  "Planell_CON_UC",
  "Mo_CON_UC",
  "Ostrowski_CON_UC",
  "Ostrowski_CON_PSC",
    #Our Dataset:
  "CON_UC",
  "PSC_UC_COHORT",
  "PSC_UC",
  "PSCUC_UC",
  "PSC_PSCUC",
  "PSC_CONTROL",
  "PSCUCPSC_UC",
  "PSCUCPSC_CON",
  "PSCUC_CON")


plot_pROC_rocs <- function(proclist=list(), procnames=NULL, size.rel=1,plot_tag="a)", plot_title=element_blank()){
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

roclist_UC_CON <- list(
  CON_UC_Diagnose_results$pROC_object,
  Planell_CON_UC_Diagnose_results$pROC_object,
  Mo_CON_UC_Diagnose_results$pROC_object,
  Ostrowski_CON_UC_Diagnose_results$pROC_object
 # Ostrowski_CON_PSC_Diagnose_results$pROC_object,
 )

roclist_names_UC_CON <- c(
  #"NULL",
  "This study: UC vs Control ",
  "Planell: UC vs Control",
  "Mo: UC vs Control",
  "Ostrowski: UC vs Control"
  )

UC_CON_plot <- plot_pROC_rocs(proclist = roclist_UC_CON, procnames = roclist_names_UC_CON,plot_tag = "a)", plot_title="UC vs Control")
PSC_CON_plot <- plot_pROC_rocs(proclist = list(PSC_CONTROL_Diagnose_results$pROC_object,
                                               PSCUC_CON_Diagnose_results$pROC_object,
                                               PSCUCPSC_CON_Diagnose_results$pROC_object
                                               ), 
                               procnames = c("PSC vs Control",
                                             "PSC/UC vs Control",
                                             "PSC+PSC/UC vs Control"),
                               plot_tag = "b)", 
                               plot_title="PSC and PSC/UC vs Control")


rocs_a_b <- ggarrange(UC_CON_plot,PSC_CON_plot)
rocs_a_b
output <- "output/ROCS_arranged"
ggsave(rocs_a_b, filename = paste0(output,".svg"), height = 7.5, width = 15, device = "svg")
ggsave(rocs_a_b, filename = paste0(output,".png"), height = 7.5, width = 15, device = "png")


ggsave(plot_pROC_rocs(proclist = list(PSC_UC_COHORT_results$pROC_object), procnames = c("UC Control vs PSC Control"), plot_tag = "a)", plot_title=c("UC Control vs PSC Control")), filename = paste0("output/ROC_UC_PSC_COHORT",".png"), height = 7.5, width = 7.5, device = "png")
ggsave(plot_pROC_rocs(proclist = list(PSC_UC_COHORT_results$pROC_object), procnames = c("UC Control vs PSC Control"), plot_tag = "a)", plot_title=c("UC Control vs PSC Control")), filename = paste0("output/ROC_UC_PSC_COHORT",".svg"), height = 7.5, width = 7.5, device = "svg")

ggsave(plot_pROC_rocs(proclist = list(PSC_PSCUC_Diagnose_results$pROC_object), procnames = c("PSC vs PSC/UC"), plot_tag = "c)", plot_title=c("PSC vs PSC/UC")), filename = paste0("output/ROC_PSC_PSCUC_DIAGNOSE",".png"), height = 7.5, width = 7.5, device = "png")
ggsave(plot_pROC_rocs(proclist = list(PSC_PSCUC_Diagnose_results$pROC_object), procnames = c("PSC vs PSC/UC"), plot_tag = "c)", plot_title=c("PSC vs PSC/UC")), filename = paste0("output/ROC_PSC_PSCUC_DIAGNOSE",".svg"), height = 7.5, width = 7.5, device = "svg")


