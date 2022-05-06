ggplot(aes(x=PSCUCPSC_UC_Diagnose_results$variable_importance[order(PSCUCPSC_UC_Diagnose_results$variable_importance$importance,decreasing = T),]$feature,y=PSCUCPSC_UC_Diagnose_results$variable_importance[order(PSCUCPSC_UC_Diagnose_results$variable_importance$importance,decreasing = T),]$importance))+geom_point()


feature_importance_df <- PSCUCPSC_UC_Diagnose_results$variable_importance
feature_importance_df <- feature_importance_df[order(feature_importance_df$importance, decreasing = T),]
feature_importance_df$feature <- ordered(feature_importance_df$feature, levels=feature_importance_df[order(feature_importance_df$importance, decreasing = T),]$feature)
ggplot(data=feature_importance_df, aes(x=feature,y=importance))+geom_point()





results_list <- list(
  CON_UC_Diagnose_results$variable_importance,
  PSC_UC_COHORT_results$variable_importance,
  PSC_UC_Diagnose_results$variable_importance,
  PSCUC_UC_Diagnose_results$variable_importance,
  PSC_PSCUC_Diagnose_results$variable_importance,
  PSC_CONTROL_Diagnose_results$variable_importance,
  PSCUCPSC_UC_Diagnose_results$variable_importance,
  PSCUCPSC_CON_Diagnose_results$variable_importance,
  Planell_CON_UC_Diagnose_results$variable_importance,
  Mo_CON_UC_Diagnose_results$variable_importance,
  Ostrowski_CON_UC_Diagnose_results$variable_importance
)
results_names<- c(
  "CON_UC_Diagnose",
  "PSC_UC_COHORT",
  "PSC_UC_Diagnose",
  "PSCUC_UC_Diagnose",
  "PSC_PSCUC_Diagnose",
  "PSC_CONTROL_Diagnose",
  "PSCUCPSC_UC_Diagnose",
  "PSCUCPSC_CON_Diagnose",
  "Planell_CON_UC_Diagnose_results",
  "Mo_CON_UC_Diagnose_results",
  "Ostrowski_CON_UC_Diagnose_results"
)

importance_table <- data.table(features=CON_UC_Diagnose_results$variable_importance$feature)

for(i in seq(1,length(results_list))){
  j <- data.table(results_list[[i]])
  colnames(j) <- c("features", results_names[i])
  importance_table <- merge(importance_table,j,by="features", all.x = T)
}
fwrite(importance_table, file="RF_importance_table.csv")

##### AUROC plots
results_list <- list(
  CON_UC_Diagnose_results$variable_importance,
  PSC_UC_COHORT_results$variable_importance,
  PSC_UC_Diagnose_results$variable_importance,
  PSCUC_UC_Diagnose_results$variable_importance,
  PSC_PSCUC_Diagnose_results$variable_importance,
  PSC_CONTROL_Diagnose_results$variable_importance,
  PSCUCPSC_UC_Diagnose_results$variable_importance,
  PSCUCPSC_CON_Diagnose_results$variable_importance,
  Planell_CON_UC_Diagnose_results$variable_importance,
  Mo_CON_UC_Diagnose_results$variable_importance,
  Ostrowski_CON_UC_Diagnose_results$variable_importance
)

pROC::ggroc(PSC_UC_COHORT_results$pROC_object)

library(pROC)
library(RColorBrewer)
Test_ROCplot <- ggroc(CON_UC_Diagnose_results$pROC_object, size=1.5)+
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
  labs(title=paste("ROC of PSC and UC Cohort, only Controls"), tag="(a)")+
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1, colour="NULL"))+
  scale_color_brewer(palette="Set2", label= c("NULL", "PSCvsUC Cohort, only CON")) +
  theme(legend.position = c(0.83, 0.13),
        panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.border = element_rect(fill = NA, 
                                                                                  colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        #legend.title = element_text(size = rel(1.8 * size.rel)),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(1.3 * size.rel)), 
        axis.title = element_text(size = rel(1.5 * size.rel)), 
        axis.text = element_text(size = rel(1.2 * size.rel)), 
        strip.text.x = element_text(size = rel(1.8 * size.rel)), 
        strip.text.y = element_text(size = rel(1.8 * size.rel)), 
        legend.key.height = unit(1.3 * size.rel,"line"), 
        legend.key.width = unit(1.3 * size.rel, "line"),
        plot.tag.position = c(-.01, .97),
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
  coord_fixed()
