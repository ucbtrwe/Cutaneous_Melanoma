
#### This code creates a Kaplan_Meier Survival Curve 


###Survival Curves

#install.packages('survival')
library('survival')

#install.packages('survminer')
library('survminer')

Clinical_data_groups$DSS_months <- Clinical_data_groups$DSS.time/31
Clinical_data_groups$DSS_years <- Clinical_data_groups$DSS.time/365


fit = survfit(Surv(DSS_years, DSS) ~ Sample_Group_Split, data = Clinical_data_groups  )
s<-ggsurv<-ggsurvplot(
  fit,
  xlab = "Years following diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.03),
  legend.title = "",
  legend.labs = c("Metastatic A1" , "Metastatic A2", "Metastatic B1", "Metastatic B2"),
  legend = "right",
  pval.size = 8,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.27,
  palette = c(   "#bdeaea" , "#5acbca", '#3e1453', "#a47ab9" ),
  font.legend = c(20),
  font.main = c(20),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

ggsurv$table <- ggrisktable(fit,
                            data = Clinical_data_groups,
                            
                            y.text = T,   ylab = "",  xlab = "Years",
                            x.text = T,
                            #   legend.labs = c("1", "2"),
                            tables.theme = theme_survminer(base_size = (10),
                                                           font.tickslab = c(10),
                                                           font.x = c(15),
                                                           font.y = c(15),
                                                           font.legend= c(15),
                                                           risk.table.fontsize = (18),
                                                           risk.table.y.text.col	= TRUE,
                                                           legend.title = "Sample Group: "
                                                           
                                                           
                            )
                            
)


pdf(("SKCM_Metastatic_Disease_Specific_survival_by_reg_cluster_August.pdf"), width = 11, height = 8.5)

print(s)
dev.off()


survdiff(Surv(DSS_years, DSS) ~ Sample_Group_Split, data = Clinical_data_groups )

