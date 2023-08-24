### This code Forest plot showing cox proportional analysis of relative effects of different 
####factors affecting differences of survival in 

library(survminer)

Clinical_data_groups$Sample_Group <- Clinical_data_groups$Final_Sample_Group_Paper
model <- coxph( Surv(DSS.time, DSS) ~ Age + Stage +
                  
                  Gender + Sample_Group ,data = Clinical_data_groups)
ggforest(model)

pdf(("Disease_Specific_survival_by_cluster_cox_forrest_regressed_Met_August_big_no_na.pdf"), width = 27.5, height = 18)
ggforest(model, fontsize = 3)
ggforest


dev.off()
