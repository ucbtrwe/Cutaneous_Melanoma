
### this code identifies the optimal cutpoint to divide samples by expression to create two 
#### groups with m greatest differences in survival and calculates significance of differences with 
### p values and adjusted p values calculated with different methods

library (survminer)
### Read in expression data #### 
PN_gene_exp <- read_excel(  'SKCM_Metastatic_clinical_exp_with_new_clusters_DS_status.xlsx')

### remove rows with NA in DSS column

PN_gene_exp <- PN_gene_exp [which(!is.na( PN_gene_exp$DSS)),]
PN_gene_exp <- PN_gene_exp [which(!is.na( PN_gene_exp$DSS.time )),]

##For loop
names <-names(PN_gene_exp)
# 1. Determine the optimal cutpoint of variables

for (i in 42:441 ){ try({
 
  results = surv_cutpoint(PN_gene_exp, time = "DSS.time", event = "DSS"
                          ,variables = names[i] )
  print(results)
  
  res.cut <- surv_cutpoint(PN_gene_exp, time = "DSS.time", event = "DSS",
                           variables = names[i])
  
  summary(res.cut)
  print(res.cut)
  res.cat <- surv_categorize(res.cut)
  print(head(res.cat))
  
  # print curve
  pdf(paste0(names[i],"_DSS_survival_curve_met.pdf"), width = 10, height = 7)
  
  fit <- survfit(Surv(DSS.time, DSS) ~ res.cat[,3], data = res.cat)
  
  
  s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE, pval = TRUE,
                      ggtheme = theme_classic2(base_size=20),
                      legend.title = names[i],
                      legend.labs = c("high", "low"))
  
  
  print(s_plot)
  dev.off()
  print(s_plot)
  print(i)
})
}


#### calculate adjusgted p values ####
#### survival cutpoint p_values manually entered into spreadsheet due to difficulities writing loop
#### load .csv with p values 

gene_surv_data <-  read.csv('Primary_and_Metastatic_Gene_Clusters_with_survival_cutpoints_p_val_March.csv')

### calculate adjusted p values

gene_surv_data$adj_p_values_met <- p.adjust(gene_surv_data$met_P_value)

gene_surv_data$adj_p_values_hoch_met <- p.adjust(gene_surv_data$met_P_value,method="hochberg" )

gene_surv_data$adj_p_values_BH_met <- p.adjust(gene_surv_data$met_P_value,method="BH" )

gene_surv_data$adj_p_values_holm_met <- p.adjust(gene_surv_data$met_P_value,method="holm" )
gene_surv_data$adj_p_values_hommel_met <- p.adjust(gene_surv_data$met_P_value,method="hommel" )
gene_surv_data$adj_p_values_BY_met <- p.adjust(gene_surv_data$met_P_value,method="BY" )

###add column to indicate if both adjusted p values are significant

gene_surv_data$significance<- ifelse (gene_surv_data$adj_p_values_BH_.prim.<0.1 &
                                        gene_surv_data$adj_p_values_BH_met<0.1 ,"sig", "non_sig")

###add column to indicate which survival is longer
gene_surv_data$survival <- paste(gene_surv_data$better.survival_.prim.,'_',gene_surv_data$higher..survival)



# save as Rdata and CSV
save(gene_surv_data, file ='Primary_and_Metastatic_Gene_Clusters_with_survival_cutpoints_p_val_Met_and_Prim.Rdata')
write.xlsx(gene_surv_data, file ='Primary_and_Metastatic_Gene_Clusters_with_survival_cutpoints_p_val_Met_and_Prim.xlsx')
