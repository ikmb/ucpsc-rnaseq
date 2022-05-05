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

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

rocs_plotted <- ggroc(roclist, size=1.5)+
  # annotate("text", x = .30, y=.39,
  #          label = paste0("AUC with PRS = ",formatC(Test_LDpred2__aucwith[[1]],format="g",digits=2),
  #                         "\n",
  #                         "CI95%: ",
  #                         formatC(Test_LDpred2__aucwithCI[[1]],format="g",digits=2),
  #                         "-",
  #                         formatC(Test_LDpred2__aucwithCI[[3]],format="g",digits=2),
  #                         "\n",
  #                         "AUC without PRS = ",formatC(Test_LDpred2__aucwithout[[1]],format="g",digits = 2),
  #                         "\n","CI95%: ",formatC(Test_LDpred2__aucwithoutCI[[1]],format="g",digits=2),
  #                         "-",
  #                         formatC(Test_LDpred2__aucwithoutCI[[3]],format="g",digits=2))
  #          ,size = 6)+
  # annotate("text", x = .70, y=.98, 
  #          label = paste0("GLM PRS p-value = ",
  #                         formatC(summary(data_table__auto.model2)$coefficients["SCORE",4],format="e",digits=2)),
  #          size = 6)+
  labs(title=element_blank(), tag="(a)")+
  #geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1, colour="NULL"))+
  scale_colour_manual(values = mycolors, label= roclist_names) +
  theme(legend.position = c(0.70, 0.19),
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
        plot.margin = margin(5.5,20,5.5,27, "pt"),
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
  guides(colour=guide_legend(ncol=2))+
  coord_fixed()

output <- "output/ROCS_plotted_all"
ggsave(rocs_plotted, filename = paste0(output,".svg"), height = 7.7, width = 7.5, device = "svg")

