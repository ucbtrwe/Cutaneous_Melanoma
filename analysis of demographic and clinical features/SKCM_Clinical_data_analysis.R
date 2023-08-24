## this code merges clinical data from TCGA and sample groups identified based on analysis of PN gene expression
### and creates survival curves showing survival of groups of patients based on different demographic features and
### clinical features and box plots and barcharts showing distribution of demographic and clincial features between groups 
#### and tests for significance of differences in survival or distribution


library('survival')
library('survminer')
library('ggplot2')

##  read csv of clinical data
clinical_data <- read.csv('SKCM_exp_info_clinical_November.csv')


###merge sample clusters with clinical data
sample_group_clinical<- merge(new_sample_clusters,clinical_data, by = 'Patient_Barcode')

sample_group_clinical$Gender <-  ifelse ((sample_group_clinical$gender == "MALE") ,"Male", "Female")
sample_group_clinical$Age <-as.numeric(sample_group_clinical$age_at_initial_pathologic_diagnosis)

sample_group_clinical$Subsequent_Metastasis_Category <-  
  ifelse ((sample_group_clinical$new_tumor_event_type == "None"), "None",
          ifelse ((sample_group_clinical$new_tumor_event_type %in% c("Distant Metastasis",
                                                                               "Distant Metastasis|Regional lymph node",
                                                                               "Regional lymph node|Distant Metastasis|Distant Metastasis",
                                                                               "Regional lymph node|Distant Metastasis" ,
                                                                               "Distant Metastasis|Locoregional Recurrence",
                                                                               "Regional lymph node|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Locoregional Recurrence|Distant Metastasis" ,
                                                                               "Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis",
                                                                               "Regional lymph node|Regional lymph node|Distant Metastasis|Regional lymph node|Distant Metastasis|Distant Metastasis"  ,
                                                                               "Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis",
                                                                               "Locoregional Recurrence|Locoregional Recurrence|Locoregional Recurrence|Distant Metastasis|Distant Metastasis",
                                                                               "Locoregional Recurrence|Distant Metastasis", "Distant Metastasis|Distant Metastasis|Distant Metastasis|Distant Metastasis", "Distant Metastasis|Distant Metastasis|Regional lymph node")), "Distant Metastasis",    
                  #ifelse ((sample_group_clinical$new_tumor_event_type == "NA"), "NA",
                  ifelse ((sample_group_clinical$new_tumor_event_type == "New primary melanoma"), 
                          "New Primary Melanoma", 
                          ifelse ((sample_group_clinical$new_tumor_event_type %in% c("Locoregional Recurrence" ,"Locoregional Recurrence", "Locoregional Recurrence","Regional lymph node" ,"Locoregional Recurrence|Regional lymph node" )), "Regional Metastasis", "Other"))))
sample_group_clinical$Subsequent_Metastasis_Category<- factor (sample_group_clinical$Subsequent_Metastasis_Category,
                                                                         levels<- rev(c("None", "New Primary Melanoma","Regional Metastasis", "Distant Metastasis","Other")))             



#Make Sample Group a factor
sample_group_clinical$Sample_Group_F<- factor (sample_group_clinical$Sample_Group_F,
                                               levels<- (c("Low PN Prim" , "High PN Prim" )))             
## add DSS by month
sample_group_clinical$DSS_months <- sample_group_clinical$DSS.time/31




#####Primary Survival ####
## select primary samples 
primary_regressed <- subset(sample_group_regressed_expr,Sample_Group %in% c("Low PN Prim" , "High PN Prim"  ))


# create table of expression

primary_regressed_exp <- data.frame(primary_regressed, row.names = 1)
# select columms with expression data and transpose dataframe so genes are rows and samples are columns
primary_regressed_exp <- primary_regressed_exp [,-c(1,2)]
primary_regressed_exp_t <- as.data.frame(t(primary_regressed_exp))




colours2 <- list('Sample_Group' = c('Low_PN' = "#296a8c", 'High_PN' = '#fde725'))
                 
  

###Survival Curves

### Create list of Low PN samples with subsequent mets

metastasised_low <- subset(sample_group_clinical$Patient_Barcode, 
                           ((sample_group_clinical$Subsequent_Metastasis_Category =="Distant Metastasis") &
                              (sample_group_clinical$Sample_Group =="Low PN")))

# ##### Save regressed sample groups and clinical data
save(sample_group_clinical,  file = 'Sample_groups_clinical.RData')
write.xlsx( sample_group_clinical, file = "Sample_groups_clinical.xlsx")

### read in clinical and sample table

sample_group_clinical <- read.csv (file = "Sample_groups_clinical.csv")

#Create subset of clinical + supplementary  data for primary samples
sample_group_clinical_and_supplementary_prim <- subset (clinical_and_supplementary_data_no_dups , clinical_and_supplementary_data_no_dups$Sample_Group %in% c("Low PN Prim","High PN Prim") )
##remove duplicates
sample_group_clinical_and_supplementary_prim <-  sample_group_clinical_and_supplementary_prim[!duplicated(sample_group_clinical_and_supplementary_prim$Sample.x), ]

#Create subset of clinical  data for primary samples
sample_group_clinical_prim <- subset (sample_group_clinical , sample_group_clinical$Sample_Group_Name.x %in% c("Low PN Prim","High PN Prim") )
##remove duplicates
sample_group_clinical_prim <-  sample_group_clinical_prim[!duplicated(sample_group_clinical_prim$Patient_Barcode), ]


#Create subset of clinical  data for Metastatic samples
sample_group_clinical_met <- subset (sample_group_clinical , sample_group_clinical$Sample_Group_Name.x %in% c("Low PN Met","High PN Met") )
##remove duplicates
sample_group_clinical_met <-  sample_group_clinical_met[!duplicated(sample_group_clinical_met$Patient_Barcode), ]



fit = survfit(Surv(DSS_months_3y, DSS_3y) ~ Sample_Group, data = sample_group_clinical_prim  )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months Following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  legend.title = "",
  legend.labs = c("High PN Prim", "Low PN Prim"),
  legend = "right",
  pval.size = 8,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.2,
  palette = c(  "#fde725",  '#296a8c'  ) ,
  font.legend = c(20),
  font.main = c(20),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  axes.offset = TRUE,
  tables.theme = clean_theme())

ggsurv$table <- ggrisktable(fit,
                            data = sample_group_clinical_prim,
                            fontsize = 30,
                            y.text = T,   ylab = "",  xlab = "Months",
                            x.text = T,
                            
                            #legend.labs = c("1", "2"),
                            tables.theme = theme_survminer(base_size = (30),
                                                           font.tickslab = c(30),
                                                           font.x = c(30),
                                                           font.y = c(30),
                                                           font.legend= c(30),
                                                           risk.table.fontsize = 30,
                                                           risk.table.y.text.col	= TRUE,
                                                           fontsize = 30,
                                                           legend.title = "Sample Group: "
                                                           
                                                           
                            )
                            
)




pdf(("SKCM_Primary_Disease_Specific_survival_by_sample_group_dec_3y.pdf"), width = 9 ,height = 7)

print(s)
dev.off()


survdiff(Surv(DSS_months_3y, DSS_3y) ~ Sample_Group, data = sample_group_clinical_prim )



####Drug Response primary
### read in edited response data

drug_response <- read.xlsx('drug_treatment_and_sample_group_edited.xlsx', sheetName = "Response")

drug_response_with_sample_group <- merge(new_sample_clusters,drug_response, by.x = 'Sample',
                                         by.y = "Name")


drug_response_with_sample_group$Sample_Group <- factor(drug_response_with_sample_group$Sample_Group_split,
                                           levels = c( "Low PN Prim" , "High PN Prim"))

####  bar chart drug response
barchart_df <- subset(drug_response_with_sample_group,!is.na(measure_of_response))

barchart_df$measure_of_response<- factor(barchart_df$measure_of_response,
                                         levels = rev(c("Complete Response", "Partial Response" ,
                                                    #"stable disease" ,
                                                    "Clinical Progressive Disease")))


t <- ggplot(barchart_df, aes(Sample_Group, fill = measure_of_response))+
  geom_bar(position="fill") +
  theme_bw() +
  
  scale_fill_manual(name = 'Drug Response',
                    values=rev(c(  
                      "Complete Response"= "#d5e1df",
                      "Partial Response" = "#b5e7a0",
                     # "Stable Disease" = "#e3eaa7",
                      "Clinical Progressive Disease" = "#86af49"))) +

  
  # ggtitle ("SKCM Primary Samples Tumour Stage by Cluster") 
  
  labs(
    #title="Age at Diagnosis",
    x ="", y = "Proportion") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20) ,
    # panel.background = element_rect(fill = "white"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16,face="bold"),
    panel.background = element_rect(fill = 'white', colour = 'black'))

pdf ("Drug_response_by_sample_group_prim_dec.pdf", width = 8, height = 4)
print(t)
dev.off()

###make list of samples with non-complete response

non_responding_samples<- subset()


### Breslow Thickness #####

##create df with only samples with breslow thickness

TCGA_data_sample_clusters_bres<- subset(TCGA_data_sample_clusters, !(TCGA_data_sample_clusters$CURATED_BRESLOW %in% c("[Not Available]","-")))

TCGA_data_sample_clusters_bres$breslow_thickness <-as.numeric(TCGA_data_sample_clusters_bres$CURATED_BRESLOW)

TCGA_data_sample_clusters_bres$Breslow_Category <- ifelse(TCGA_data_sample_clusters_bres$breslow_thickness< 5, 'Thin', 'Thick' ) 
TCGA_data_sample_clusters_bres$Breslow_Category<- factor(TCGA_data_sample_clusters_bres$Breslow_Category,
                                                           levels = c("Thick","Thin"))
TCGA_data_sample_clusters_bres$Sample_Group <- factor(TCGA_data_sample_clusters_bres$Sample_Group,
                                                      levels = c('Low_PN', 'High_PN'))
 ##Remove duplicates by patient barcode    
##Save
save(TCGA_data_sample_clusters_bres, file ='Breslow_Thickness_Clinical_Data.Rdata')
write.xlsx(TCGA_data_sample_clusters_bres, file ='Breslow_Thickness_Clinical_Data.xlsx')



####Breslow Thickness Bar Chart#####
t <-  ggplot(TCGA_data_sample_clusters_bres , aes(Sample_Group,
                                             fill = Breslow_Category))+
  geom_bar(position="fill") +
  theme_bw() +
  
 # scale_fill_manual(values=c("None" = "darkblue",
 # 
 #                            "New Primary Melanoma"  = "chocolate4",
 #                            "Regional Metastasis" = "yellow",
 #                            "Distant Metastasis" =  "red3"
 # )) +
  labs(fill = "Breslow Thickness", 
       x ="Sample Group", y = "Proportion of Patients")

# ggtitle ("New Tumour Event Type by Cluster")

pdf (("SKCM_Breslow_Thickness_By_Final_Sample_Group.pdf"),width=4, height=2)

print(t)
dev.off()
# 
# Fisher test for progression

fisher.test(TCGA_data_sample_clusters_bres$Sample_Group,TCGA_data_sample_clusters_bres$Breslow_Category)


####Survival Curve Breslow ###
fit = survfit(Surv(DSS_months_3y.y, DSS_3y.y) ~ Breslow_Category, data = TCGA_data_sample_clusters_bres  )
s<-ggsurv <- ggsurvplot(
  fit,
  xlab = "Months",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.03),
  legend.title = "",
 # legend.labs = c("C (Metastatic)","D (Metastatic)","E (Metastatic)" , "F (Metastatic)" ),
  legend = "right",
  pval.size = 7,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.27,
#  palette = rev(c("#005a00", "#99c199",'#a47ab9',  "#3e1453" ) ),
  font.legend = c(20),
  font.main = c(20),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

ggsurv$table <- ggrisktable(fit,
                            data = TCGA_data_sample_clusters_bres ,
                            
                            y.text = T,   ylab = "",  xlab = "Months",
                            x.text = T,
                            #   legend.labs = c("1", "2"),
                            tables.theme = theme_survminer(base_size = (10),
                                                           font.tickslab = c(10),
                                                           font.x = c(10),
                                                           font.y = c(10),
                                                           font.legend= c(10),
                                                           risk.table.fontsize = (10),
                                                           risk.table.y.text.col	= TRUE,
                                                           legend.title = "Sample Group: "
                                                           
                                                           
                            )
                            
)
ggsurv

pdf(("SKCM_Primary_Disease_Specific_survival_by_Breslow_Category_Nov.pdf"), width = 14, height = 6)

print(s)
dev.off()


survdiff(Surv(DSS_months_3y.y, DSS_3y.y) ~ Breslow_Category, data = TCGA_data_sample_clusters_bres  )

#####Box Plot Breslow Thickness
j <- ggplot(TCGA_data_sample_clusters_bres , aes(x= Sample_Group, y= breslow_thickness, fill= Sample_Group) )+ 
  geom_boxplot() +
  
  scale_fill_manual( name = "Sample Group", values=c(  "#296a8c",
                                                       "#fde725" )) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  #scale_fill_manual(values=c( '', '')) +
  # scale_x_discrete(labels=c('1'=' Cluster 1' ,  '2' ='Cluster 2'  )) +
  stat_compare_means(aes(group = Sample_Group), label = "p.format")  +
  # geom_jitter() +
  
  labs(
    #title="Age at Diagnosis",
    x ="Sample Group", y = "Breslow Thickness") +
  theme(
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20) ,
    # panel.background = element_rect(fill = "white"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16,face="bold"),
    panel.background = element_rect(fill = 'white', colour = 'black')
    
  ) 
pdf(("Breslow_THickness_boxplot_Nov.pdf"), width = 8, height = 5)

print(j)
dev.off()

#### Pigmentation ####
# ## create Bar chart of pigmentation scores 

###Primary 
clinical_and_supplementary_data_prim <- subset(clinical_and_supplementary_data_no_dups,Sample_Group
                                               %in% c("Low PN Prim",  "High PN Prim"))

clinical_and_supplementary_data_prim$Sample_Group_Final <- factor(clinical_and_supplementary_data_prim$Sample_Group_Final,
                                                            levels = c( "Low PN Prim","High PN Prim" ))


###Bar chart by pigment score
barchart_df <- clinical_and_supplementary_data_prim 

barchart_df$PIGMENT.SCORE<- factor(barchart_df$PIGMENT.SCORE,
                                   levels = c("3", "2","1","0" ))

###

t <- ggplot(barchart_df, aes(Sample_Group_Final, fill = PIGMENT.SCORE))+
  geom_bar(position="fill") +
  
  
  scale_fill_manual(name = 'Pigmentation Score', 
                    values=c(  "0" = "#f3ece9",
                               "1"= "#c39f91",
                               "2" = "#873e23",  
                               "3"= "#512515" )) +
  
  
  # ggtitle ("SKCM Primary Samples Tumour Stage by Cluster") 
  
  labs(
    #title="Age at Diagnosis",
    x ="", y = "Proportion") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20) ,
        # panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        panel.background = element_rect(fill = 'white', colour = 'black'))

pdf ("SKCM Primary Samples Pigmentation by Sample Group.pdf", width = 8, height = 5)
print(t)
dev.off()


#### Pigmentation ####
# ## create Bar chart of pigmentation scores 

###Primary 
clinical_and_supplementary_data_prim <- subset(clinical_and_supplementary_data_no_dups,Sample_Group
                                               %in% c("Low PN Prim",  "High PN Prim"))
###Bar chart by pigment score
barchart_df <- clinical_and_supplementary_data_prim 

barchart_df$pigmentation_level<- factor(barchart_df$pigmentation_level,
                                   levels = c("Dark", "Light" ))

###

t <- ggplot(barchart_df, aes(Sample_Group_Final, fill = pigmentation_level))+
  geom_bar(position="fill") +
  
  
  scale_fill_manual(name = 'Pigmentation', 
                    values=c(  "Light" = "#f3ece9",
                               
                               "Dark" = "#873e23"  
                                )) +
  
  
  # ggtitle ("SKCM Primary Samples Tumour Stage by Cluster") 
  
  labs(
    #title="Age at Diagnosis",
    x ="", y = "Proportion") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20) ,
        # panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        panel.background = element_rect(fill = 'white', colour = 'black'))

pdf ("SKCM Primary Samples Pigmentation category by Sample Group.pdf", width = 8, height = 5)
print(t)
dev.off()




##fisher
fisher.test(barchart_df$Sample_Group_Final,barchart_df$pigmentation_level)
####Survival analysis Primary Pigmentation

# Survival curves
fit = survfit(Surv(DSS_months, DSS) ~ PIGMENT.SCORE , data = clinical_and_supplementary_data_prim )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation Score:",
  legend.labs = c("0", "1", "2","3" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in Primary Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c( "#f3ece9",
               "#c39f91",
               "#873e23",  
               "#512515") ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                            data = clinical_and_supplementary_data_prim,
                            
                            y.text = T,   ylab = "",  xlab = "Months",
                            x.text = T, 
                            legend.labs = c("0", "1", "2" , "3" ),
                            tables.theme = theme_survminer(base_size = (16),
                                                           font.tickslab = c(16),
                                                           font.x = c(14),
                                                           font.y = c(14),
                                                           font.legend= c(16),
                                                           risk.table.fontsize = 20, 
                                                           risk.table.y.text.col	= TRUE,
                                                          legend.title = "Pigmentation Score: "
                                                           
                                                           
                            )
                            
)
pdf(("survival_by_pigmentation_score_prim.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS.time, DSS) ~ PIGMENT.SCORE, data = clinical_and_supplementary_data_prim )

# Survival curves Pigmentation
fit = survfit(Surv(DSS_months, DSS) ~ Pigmentation , data = clinical_and_supplementary_data_prim )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation",
  legend.labs = c("Not Pigmented", "Pigmented" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in Primary Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c( "#f3ece9",
               
               "#873e23"  
               ) ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                       data = clinical_and_supplementary_data_prim,
                       
                       y.text = T,   ylab = "",  xlab = "Months",
                       x.text = T, 
                       legend.labs = c("Not Pigmented", "Pigmented" ),
                       tables.theme = theme_survminer(base_size = (16),
                                                      font.tickslab = c(16),
                                                      font.x = c(14),
                                                      font.y = c(14),
                                                      font.legend= c(16),
                                                      risk.table.fontsize = 20, 
                                                      risk.table.y.text.col	= TRUE,
                                                      # legend.title = "Pigmentation Score: "
                                                      
                                                      
                       )
                       
)
pdf(("survival_by_pigmentation_presence_prim.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS.time, DSS) ~ Pigmentation , data = clinical_and_supplementary_data_prim )

# Survival curves light dark primary 3 years
fit = survfit(Surv(DSS_months_3y, DSS_3y) ~ pigmentation_level , data = clinical_and_supplementary_data_prim )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation Level:",
  legend.labs = c( "Dark","Light" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in primary Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c("#873e23", "#f3ece9") ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                       data = clinical_and_supplementary_data_prim,
                       
                       y.text = T,   ylab = "",  xlab = "Months",
                       x.text = T, 
                       legend.labs = c( "Dark" , "Light"),
                       tables.theme = theme_survminer(base_size = (16),
                                                      font.tickslab = c(16),
                                                      font.x = c(14),
                                                      font.y = c(14),
                                                      font.legend= c(16),
                                                      risk.table.fontsize = 20, 
                                                      risk.table.y.text.col	= TRUE,
                                                      # legend.title = "Pigmentation Score: "
                                                      
                                                      
                       )
                       
)
pdf(("survival_by_pigmentation_level_2_prim_3y.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS_months_3y, DSS_3y) ~ pigmentation_level , data = clinical_and_supplementary_data_prim )






# Survival curves light dark primary 2 years
fit = survfit(Surv(DSS_months_2y, DSS_2y) ~ pigmentation_level , data = clinical_and_supplementary_data_prim )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation Level:",
  legend.labs = c( "Dark","Light" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in primary Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c("#873e23", "#f3ece9") ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                       data = clinical_and_supplementary_data_prim,
                       y.text = T,   ylab = "",  xlab = "Months",
                       x.text = T, 
                       legend.labs = c( "Dark" , "Light"),
                       tables.theme = theme_survminer(base_size = (16),
                                                      font.tickslab = c(16),
                                                      font.x = c(14),
                                                      font.y = c(14),
                                                      font.legend= c(16),
                                                      risk.table.fontsize = 20, 
                                                      risk.table.y.text.col	= TRUE,
                                                      # legend.title = "Pigmentation Score: "
                                                      
                                                      
                       )
                       
)
pdf(("survival_by_pigmentation_level_2_prim.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS_months_2y, DSS_2y) ~ pigmentation_level , data = clinical_and_supplementary_data_prim )


# Survival curves light dark primary 1 years
fit = survfit(Surv(DSS_months_1y, DSS_1y) ~ pigmentation_level , data = clinical_and_supplementary_data_prim )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation Level:",
  legend.labs = c( "Dark","Light" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in primary Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c("#873e23", "#f3ece9") ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                       data = clinical_and_supplementary_data_prim,
                       
                       y.text = T,   ylab = "",  xlab = "Months",
                       x.text = T, 
                       legend.labs = c( "Dark" , "Light"),
                       tables.theme = theme_survminer(base_size = (16),
                                                      font.tickslab = c(16),
                                                      font.x = c(14),
                                                      font.y = c(14),
                                                      font.legend= c(16),
                                                      risk.table.fontsize = 20, 
                                                      risk.table.y.text.col	= TRUE,
                                                      # legend.title = "Pigmentation Score: "
                                                      
                                                      
                       )
                       
)
pdf(("survival_by_pigmentation_level_2_prim_1y.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS_months_1y, DSS_1y) ~ pigmentation_level , data = clinical_and_supplementary_data_prim )








#####Metastatic Pigmentation ####

clinical_and_supplementary_data_met <- subset(clinical_and_supplementary_data_no_dups,Sample_Group
                                               %in% c("Low PN Met" ,  "High PN Met" ))
clinical_and_supplementary_data_met <- subset(clinical_and_supplementary_data_met,clinical_and_supplementary_data_met$PIGMENT.SCORE !="-")

###Bar chart by pigment score
barchart_df <- subset(clinical_and_supplementary_data_met,clinical_and_supplementary_data_met$PIGMENT.SCORE !="-")
barchart_df$PIGMENT.SCORE<- factor(barchart_df$PIGMENT.SCORE,
                                   levels = c("3", "2","1","0" ))

###

t <- ggplot(barchart_df, aes(Sample_Group_Final, fill = PIGMENT.SCORE))+
  geom_bar(position="fill") +
  
  
  scale_fill_manual(name = 'Pigmentation Score', 
                    values=c(  "0" = "#f3ece9",
                               "1"= "#c39f91",
                               "2" = "#873e23",  
                               "3"= "#512515" )) +
  
  
  # ggtitle ("SKCM metastatic Samples Tumour Stage by Cluster") 
  
  labs(
    #title="Age at Diagnosis",
    x ="", y = "Proportion") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20) ,
        # panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        panel.background = element_rect(fill = 'white', colour = 'black'))

pdf ("SKCM metastatic Samples Pigmentation by Sample Group.pdf", width = 8, height = 5)
print(t)
dev.off()

##fisher
fisher.test(barchart_df$Sample_Group_Final,barchart_df$PIGMENT.SCORE)


#bookmarkc####

barchart_df$Sample_Group_Final<- factor(barchart_df$Sample_Group_Final,
                                        levels = c("Low PN Met" , "High PN Met" ))
barchart_df$pigmentation_level<- factor(barchart_df$pigmentation_level,
                                        levels = c("Dark", "Light" ))

###

t <- ggplot(barchart_df, aes(Sample_Group_Final, fill = pigmentation_level))+
  geom_bar(position="fill") +
  
  
  scale_fill_manual(name = 'Pigmentation', 
                    values=c(  "Light" = "#f3ece9",
                               
                               "Dark" = "#873e23"  
                    )) +
  
  
  # ggtitle ("SKCM Primary Samples Tumour Stage by Cluster") 
  
  labs(
    #title="Age at Diagnosis",
    x ="", y = "Proportion") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20) ,
        # panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        panel.background = element_rect(fill = 'white', colour = 'black'))

pdf ("SKCM Metastatic Samples Pigmentation category by Sample Group.pdf", width = 8, height = 5)
print(t)
dev.off()


##fisher
fisher.test(barchart_df$Sample_Group_Final,barchart_df$pigmentation_level)






####Survival analysis metastatic Pigmentation

# Survival curves
fit = survfit(Surv(DSS_months, DSS) ~ PIGMENT.SCORE , data = clinical_and_supplementary_data_met )

s<-ggsurvplot(
  fit,
  break.time.by = 48,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation Score:",
  legend.labs = c("0", "1", "2","3" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in metastatic Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c( "#f3ece9",
               "#c39f91",
               "#873e23",  
               "#512515") ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                       data = clinical_and_supplementary_data_met,
                       
                       y.text = T,   ylab = "",  xlab = "Months",
                       x.text = T, 
                       legend.labs = c("0", "1", "2" , "3" ),
                       tables.theme = theme_survminer(base_size = (16),
                                                      font.tickslab = c(16),
                                                      font.x = c(14),
                                                      font.y = c(14),
                                                      font.legend= c(16),
                                                      risk.table.fontsize = 20, 
                                                      risk.table.y.text.col	= TRUE,
                                                      # legend.title = "Pigmentation Score: "
                                                      
                                                      
                       )
                       
)
pdf(("survival_by_pigmentation_score_met.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS_months, DSS) ~ PIGMENT.SCORE, data = clinical_and_supplementary_data_met )

# Survival curves light dark metastatic
fit = survfit(Surv(DSS_months, DSS) ~ pigmentation_level , data = clinical_and_supplementary_data_met )

s<-ggsurvplot(
  fit,
  break.time.by = 48,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation Level:",
  legend.labs = c( "Dark","Light" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in Primary Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c("#873e23", "#f3ece9") ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                       data = clinical_and_supplementary_data_prim,
                       
                       y.text = T,   ylab = "",  xlab = "Months",
                       x.text = T, 
                       legend.labs = c( "Dark" , "Light"),
                       tables.theme = theme_survminer(base_size = (16),
                                                      font.tickslab = c(16),
                                                      font.x = c(14),
                                                      font.y = c(14),
                                                      font.legend= c(16),
                                                      risk.table.fontsize = 20, 
                                                      risk.table.y.text.col	= TRUE,
                                                      # legend.title = "Pigmentation Score: "
                                                      
                                                      
                       )
                       
)
pdf(("survival_by_pigmentation_level_2_met.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS.time, DSS) ~ pigmentation_level , data = clinical_and_supplementary_data_met )

# Survival curves Pigmentation
fit = survfit(Surv(DSS_months, DSS) ~ Pigmentation , data = clinical_and_supplementary_data_met )

s<-ggsurvplot(
  fit,
  break.time.by = 48,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  
  legend.title = "Pigmentation",
  legend.labs = c("Not Pigmented", "Pigmented" ),
  legend = "right",
  pval.size = 7,
  #title = "Disease Specific Survival Probability in metary Samples by Pigmentation Score",
  risk.table = TRUE,
  risk.table.fontsize = 25,
  risk.table.height =0.27,
  palette = c( "#f3ece9",
               
               "#873e23"  
  ) ,
  font.legend = c(18),
  font.main = c(18),
  font.x = c(16),
  font.y = c(16),
  font.tickslab = c(18),
  fontsize = 25,
  risk.table.y.text.col	=TRUE,
  
  tables.theme = clean_theme())

s$table <- ggrisktable(fit,
                       data = clinical_and_supplementary_data_met,
                       
                       y.text = T,   ylab = "",  xlab = "Months",
                       x.text = T, 
                       legend.labs = c("Not Pigmented", "Pigmented" ),
                       tables.theme = theme_survminer(base_size = (16),
                                                      font.tickslab = c(16),
                                                      font.x = c(14),
                                                      font.y = c(14),
                                                      font.legend= c(16),
                                                      risk.table.fontsize = 20, 
                                                      risk.table.y.text.col	= TRUE,
                                                      # legend.title = "Pigmentation Score: "
                                                      
                                                      
                       )
                       
)
pdf(("survival_by_pigmentation_presence_met.pdf"), width = 8.5, height = 7)
print(s)

dev.off()


survdiff(Surv(DSS.time, DSS) ~ Pigmentation , data = clinical_and_supplementary_data_met )

#### survival by metastasis ### metastatic samples ####
clinical_and_supplementary_data_met <- subset(clinical_and_supplementary_data_no_dups,Sample_Group
                                              %in% c("Low PN Met" ,  "High PN Met" ))
#remove samples without subsequent metastasis data
clinical_and_supplementary_data_met <- subset(clinical_and_supplementary_data_met,!(Subsequent_Metastasis_Category
                                              %in% c("Other" )))

clinical_and_supplementary_data_met$Subsequent_Metastasis_Category <- factor (clinical_and_supplementary_data_met$Subsequent_Metastasis_Category,
           levels<- rev(c("None", "New Primary Melanoma","Regional Metastasis", "Distant Metastasis")))             

clinical_and_supplementary_data_met$Sample_Group<- factor(clinical_and_supplementary_data_met$Sample_Group,
                                                    levels = c("Low PN Met" , "High PN Met"))





fit = survfit(Surv(DSS_months, DSS) ~ Subsequent_Metastasis_Category, data = clinical_and_supplementary_data_met)

s<-ggsurvplot(
  fit,
  break.time.by = 48,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  legend.title = "",
  legend.labs = c("Distant\nMetastasis", "Regional\nMetastasis", "New Primary\nMelanoma", "None"),
  legend = "right",
  pval.size = 10,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.25,
  palette = rev (c(  "darkblue","chocolate4","yellow", "red3"  )) ,
  font.legend = c(20),
  font.main = c(20),
  font.x = c(19),
  font.y = c(19),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  axes.offset = TRUE,
  tables.theme = clean_theme())

ggrisktable(fit,
            data = clinical_and_supplementary_data_met,
            fontsize = 30,
            y.text = T,   ylab = "",  xlab = "Months",
            x.text = T,
            
            legend.labs = c("Distant Metastasis", "Regional Metastasis", "New Primary Melanoma", "None"),
            tables.theme = theme_survminer(base_size = (30),
                                           font.tickslab = c(30),
                                           font.x = c(30),
                                           font.y = c(30),
                                           font.legend= c(30),
                                           risk.table.fontsize = 30,
                                           risk.table.y.text.col	= TRUE,
                                           fontsize = 30,
                                           legend.title = "Subsequent Metastasis "
                                           
                                           
            )
            
)

pdf(("SKCM_Metastatic_Disease_Specific_survival_by_metastasis_nov.pdf"), width = 10 ,height = 8)

print(s)
dev.off()


survdiff(Surv(DSS_months, DSS) ~ Subsequent_Metastasis_Category, data = clinical_and_supplementary_data_met )

#### Survival by sample group in samples with no metastasis

sample_group_clinical_surv_metastatic_no_subs_met <- subset (clinical_and_supplementary_data_met ,clinical_and_supplementary_data_met$Subsequent_Metastasis_Category == "None")

fit = survfit(Surv(DSS_months, DSS) ~ Sample_Group, data = sample_group_clinical_surv_metastatic_no_subs_met  )

s<-ggsurv<-ggsurvplot(
  fit,
  break.time.by = 48,
  size.est = 1,
  xlab = "Months Following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  legend.title = "",
  legend.labs = c("High PN Met", "Low PN Met"),
  legend = "right",
  pval.size = 8,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.2,
  palette = c(  "#fde725",  '#296a8c'  ) ,
  font.legend = c(20),
  font.main = c(20),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  axes.offset = TRUE,
  tables.theme = clean_theme())

ggsurv$table <- ggrisktable(fit,
                            data = sample_group_clinical_surv_no_met,
                            fontsize = 30,
                            y.text = T,   ylab = "",  xlab = "Months",
                            x.text = T,
                            legend.labs = c("High PN Met", "Low PN Met"),
                            tables.theme = theme_survminer(base_size = (30),
                                                           font.tickslab = c(30),
                                                           font.x = c(30),
                                                           font.y = c(30),
                                                           font.legend= c(30),
                                                           risk.table.fontsize = 30,
                                                           risk.table.y.text.col	= TRUE,
                                                           fontsize = 30,
                                                           legend.title = "Sample Group: "
                                                           
                                                           
                            )
                            
)




pdf(("SKCM_metastatic_Disease_Specific_survival_by_Sample_Group_No_Metastasis_nov.pdf"), width = 8.5 ,height = 7)

print(s)
dev.off()


survdiff(Surv(DSS_months, DSS) ~ Sample_Group, data = sample_group_clinical_surv_metastatic_no_subs_met )


#### Survival by sample group in samples with metastasis

sample_group_clinical_surv_metastatic_subs_met <- subset (clinical_and_supplementary_data_met ,clinical_and_supplementary_data_met$Subsequent_Metastasis_Category != "None")

fit = survfit(Surv(DSS_months, DSS) ~ Sample_Group, data = sample_group_clinical_surv_metastatic_subs_met  )

s<-ggsurv<-ggsurvplot(
  fit,
  break.time.by = 48,
  size.est = 1,
  xlab = "Months Following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  legend.title = "",
  legend.labs = c("High PN Met", "Low PN Met"),
  legend = "right",
  pval.size = 8,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.2,
  palette = c(  "#fde725",  '#296a8c'  ) ,
  font.legend = c(20),
  font.main = c(20),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  axes.offset = TRUE,
  tables.theme = clean_theme())

ggsurv$table <- ggrisktable(fit,
                            data = sample_group_clinical_surv_metastatic_subs_met,
                            fontsize = 30,
                            y.text = T,   ylab = "",  xlab = "Months",
                            x.text = T,
                            legend.labs = c("High PN Met", "Low PN Met"),
                            tables.theme = theme_survminer(base_size = (30),
                                                           font.tickslab = c(30),
                                                           font.x = c(30),
                                                           font.y = c(30),
                                                           font.legend= c(30),
                                                           risk.table.fontsize = 30,
                                                           risk.table.y.text.col	= TRUE,
                                                           fontsize = 30,
                                                           legend.title = "Sample Group: "
                                                           
                                                           
                            )
                            
)




pdf(("SKCM_metastatic_Disease_Specific_survival_by_Sample_Group_Metastasis_nov.pdf"), width = 8.5 ,height = 7)

print(s)
dev.off()


survdiff(Surv(DSS_months, DSS) ~ Sample_Group, data = sample_group_clinical_surv_metastatic_subs_met )





####Subsequent Metastasis Bar Chart#####
p <-  ggplot(clinical_and_supplementary_data_met , aes(Sample_Group,
                                             fill = Subsequent_Metastasis_Category))+
  geom_bar(position="fill") +
  theme_bw() +
  
  scale_fill_manual(values=c("None" = "darkblue",
                             "New Primary Melanoma"  = "chocolate4",
                             "Regional Metastasis" = "yellow",
                             "Distant Metastasis" =  "red3"
  )) +
  labs(fill = "Subsequent Metastasis", 
       x ="Sample Group", y = "Proportion of Patients")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggtitle ("New Tumour Event Type by Cluster")

pdf (("SKCM_Subsequent_Metastasis_By_Final_Sample_Group_Metastatic.pdf"),width=4.5, height=2)

print(p)
dev.off()
# 
# Fisher test for progression

fisher.test(clinical_and_supplementary_data_met$Sample_Group,clinical_and_supplementary_data_met$Subsequent_Metastasis_Category)


#### survival by metastasis primary ####
fit = survfit(Surv(DSS_months_3y, DSS_3y) ~ Subsequent_Metastasis_Category, data = clinical_and_supplementary_data_prim  )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  legend.title = "",
  legend.labs = c("Distant\nMetastasis", "Regional\nMetastasis", "New Primary\nMelanoma", "None"),
  legend = "right",
  pval.size = 10,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.25,
  palette = rev (c(  "darkblue","chocolate4","yellow", "red3"  )) ,
  font.legend = c(20),
  font.main = c(20),
  font.x = c(19),
  font.y = c(19),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  axes.offset = TRUE,
  tables.theme = clean_theme())

ggrisktable(fit,
            data = clinical_data,
            fontsize = 30,
            y.text = T,   ylab = "",  xlab = "Months",
            x.text = T,
            
            legend.labs = c("Distant Metastasis", "Regional Metastasis", "New Primary Melanoma", "None"),
            tables.theme = theme_survminer(base_size = (30),
                                           font.tickslab = c(30),
                                           font.x = c(30),
                                           font.y = c(30),
                                           font.legend= c(30),
                                           risk.table.fontsize = 30,
                                           risk.table.y.text.col	= TRUE,
                                           fontsize = 30,
                                           legend.title = "Subsequent Metastasis "
                                           
                                           
            )
            
)




pdf(("SKCM_Primary_Disease_Specific_survival_by_metastasis_3y_nov.pdf"), width = 10 ,height = 8)

print(s)
dev.off()


survdiff(Surv(DSS_months, DSS) ~ Subsequent_Metastasis_Category, data = sample_group_clinical_prim )

#### Survival by sample group in samples with no metastasis

sample_group_clinical_surv_no_met <- subset (sample_group_clinical_prim ,sample_group_clinical_prim$Subsequent_Metastasis_Category == "None")

fit = survfit(Surv(DSS_months_3y, DSS_3y) ~ Sample_Group, data = sample_group_clinical_surv_no_met  )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  legend.title = "",
  legend.labs = c("High PN Prim", "Low PN Prim"),
  legend = "right",
  pval.size = 8,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.2,
  palette = c(  "#fde725",  '#296a8c'  ) ,
  font.legend = c(20),
  font.main = c(20),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  axes.offset = TRUE,
  tables.theme = clean_theme())

ggsurv$table <- ggrisktable(fit,
                            data = sample_group_clinical_surv_no_met,
                            fontsize = 30,
                            y.text = T,   ylab = "",  xlab = "Months",
                            x.text = T,
                            
                            #legend.labs = c("1", "2"),
                            tables.theme = theme_survminer(base_size = (30),
                                                           font.tickslab = c(30),
                                                           font.x = c(30),
                                                           font.y = c(30),
                                                           font.legend= c(30),
                                                           risk.table.fontsize = 30,
                                                           risk.table.y.text.col	= TRUE,
                                                           fontsize = 30,
                                                           legend.title = "Sample Group: "
                                                           
                                                           
                            )
                            
)




pdf(("SKCM_Primary_Disease_Specific_survival_by_Sample_Group_No_Metastasis_nov.pdf"), width = 8.5 ,height = 7)

print(s)
dev.off()


survdiff(Surv(DSS_months_3y, DSS_3y) ~ Sample_Group, data = sample_group_clinical_surv_no_met )


#### Survival by sample group in samples with metastasis

sample_group_clinical_surv_met <- subset (sample_group_clinical_prim ,sample_group_clinical_prim$Subsequent_Metastasis_Category != "None")

fit = survfit(Surv(DSS_months_3y, DSS_3y) ~ Sample_Group, data = sample_group_clinical_surv_met  )

s<-ggsurvplot(
  fit,
  break.time.by = 6,
  size.est = 1,
  xlab = "Months following Diagnosis",
  ylab = "Disease Specific Survival Probability",
  pval = TRUE,
  pval.coord = c(0, 0.11),
  legend.title = "",
  legend.labs = c("High PN", "Low PN"),
  # legend = "right",
  pval.size = 8,
  # title = "Disease Specific Survival Probability by Cluster",
  risk.table = TRUE,
  risk.table.fontsize = 8,
  risk.table.height =0.2,
  palette = c(  "#fde725",  '#296a8c'  ) ,
  font.legend = c(20),
  font.main = c(20),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(20),
  fontsize = 20,
  risk.table.y.text.col	=TRUE,
  axes.offset = TRUE,
  tables.theme = clean_theme())

ggsurv$table <- ggrisktable(fit,
                            data = sample_group_clinical_surv_met,
                            fontsize = 30,
                            y.text = T,   ylab = "",  xlab = "Months",
                            x.text = T,
                            
                            #legend.labs = c("1", "2"),
                            tables.theme = theme_survminer(base_size = (30),
                                                           font.tickslab = c(30),
                                                           font.x = c(30),
                                                           font.y = c(30),
                                                           font.legend= c(30),
                                                           risk.table.fontsize = 30,
                                                           risk.table.y.text.col	= TRUE,
                                                           fontsize = 30,
                                                           legend.title = "Sample Group: "
                                                           
                                                           
                            )
                            
)




pdf(("SKCM_Primary_Disease_Specific_survival_by_Sample_Group_Metastasis_nov.pdf"), width = 8 ,height = 7)

print(s)
dev.off()


survdiff(Surv(DSS_months_3y, DSS_3y) ~ Sample_Group, data = sample_group_clinical_surv_met )
###bookmarkA####

####Subsequent Metastasis Bar Chart_Primary#####

###Primary 
sample_group_clinical_prim$Sample_Group_F<-factor(sample_group_clinical_prim$Sample_Group_F,
                                                                   levels =c("Low PN Prim" , "High PN Prim"))
#remove duplicated samples
p <-  ggplot(sample_group_clinical_prim , aes(Sample_Group_F,
                                             fill = Subsequent_Metastasis_Category))+
  geom_bar(position="fill") +
  theme_bw() +
  
  scale_fill_manual(values=c("None" = "darkblue",
                             
                             "New Primary Melanoma"  = "chocolate4",
                             "Regional Metastasis" = "yellow",
                             "Distant Metastasis" =  "red3"
  )) +
  labs(fill = "Subsequent Metastasis", 
       x ="Sample Group", y = "Proportion of Patients")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggtitle ("New Tumour Event Type by Cluster")

pdf (("SKCM_Subsequent_Metastasis_By_Final_Sample_Group_primary_1.pdf"),width=4.5, height=2)

print(p)
dev.off()
# 
# Fisher test for progression

fisher.test(sample_group_clinical_prim$Sample_Group_F,sample_group_clinical_prim$Subsequent_Metastasis_Category)



####Compare period from diagnosis to last contact prim ####
sample_group_clinical_prim$max_time <-ifelse(sample_group_clinical_prim$last_contact_days_to !="None", 
                                             sample_group_clinical_prim$last_contact_days_to,
                                             sample_group_clinical_prim$death_days_to
                                             )
sample_group_clinical_prim$max_time <- as.numeric(sample_group_clinical_prim$max_time)

sample_group_clinical_prim$prim_last_contact_years <-  sample_group_clinical_prim$max_time/365

p<- ggplot(data = sample_group_clinical_prim, aes(x=prim_last_contact_years, fill = factor(sample_group_clinical_prim$DSS))) + 
  geom_histogram()+ 
 # geom_histogram(binwidth=1)
# Change colors
# p<-ggplot(df, aes(x=weight)) + 
#   geom_histogram(color="black", fill="white")
  scale_fill_manual(name = "Status at last contact",
                    labels = c( "Alive", "Dead"),
                    values=c(c("0" = "darkgreen",
                               "1"  = "#838B83")))+
  theme_bw()+
  labs( x ="Years From Diagnosis to Death or Last Contact" , y = 'Count'   )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        )+
geom_vline(aes(xintercept=mean(prim_last_contact_years)),
             color="blue", linetype="dashed", linewidth=1)
  #####bookmark ####
  
pdf(("primary_last_contact_days_to_August.pdf"), width = 7 ,height = 4)

print(p)
dev.off()

####Compare period from diagnosis to last contact met ####

sample_group_clinical_met$max_time <-ifelse(sample_group_clinical_met$last_contact_days_to !="None", 
                                             sample_group_clinical_met$last_contact_days_to,
                                             sample_group_clinical_met$death_days_to
)
sample_group_clinical_met$max_time <- as.numeric(sample_group_clinical_met$max_time)

sample_group_clinical_met$met_last_contact_years <-  sample_group_clinical_met$max_time/365


sample_group_clinical_met<- subset(sample_group_clinical_met, DSS %in% c("0","1")) 
  
p<- ggplot(data = sample_group_clinical_met, aes(x=met_last_contact_years, fill = factor(sample_group_clinical_met$DSS))) + 
  geom_histogram()+ 
  # geom_histogram(binwidth=1)
  # Change colors
  # p<-ggplot(df, aes(x=weight)) + 
  #   geom_histogram(color="black", fill="white")
  scale_fill_manual(name = "Status at last contact",
                    labels = c( "Alive", "Dead"),
                    values=c(c("0" = "darkgreen",
                               "1"  = "#838B83")))+
  scale_x_continuous(breaks= c(5,10,15,20,25,30))+
theme_bw()+
  labs( x ="Years From Diagnosis to Death or Last Contact"  , y = 'Count'  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept=mean(met_last_contact_years, na.rm = TRUE)),
             color="blue", linetype="dashed", linewidth=1)
#####bookmark ####

pdf(("metastatic_last_contact_days_to_August.pdf"), width = 7 ,height = 4)

print(p)
dev.off()


###Meatastatic Breslow

##create clinical table with metastatic samples

clinical_and_supplementary_data_no_dups_met<- subset(clinical_and_supplementary_data_no_dups, 
                                                     (clinical_and_supplementary_data_no_dups$Final_Sample_Group.x %in% c("D (Metastatic)", "F (Metastatic)" ,"E (Metastatic)" ,"C (Metastatic)")))
###Add column with sample group name for paper

clinical_and_supplementary_data_no_dups_met$Paper_Sample_Group <-
  ifelse(clinical_and_supplementary_data_no_dups_met$Sample_Group_Name.x == "Low PN Met", "Metastatic A","Metastatic B" )

### Breslow Thickness #####

##create df with only samples with breslow thickness

clinical_and_supplementary_data_no_dups_met_bres<- subset(clinical_and_supplementary_data_no_dups_met, !(clinical_and_supplementary_data_no_dups_met$CURATED_BRESLOW %in% c("[Not Available]","-")))

clinical_and_supplementary_data_no_dups_met_bres$breslow_thickness <-as.numeric(clinical_and_supplementary_data_no_dups_met_bres$CURATED_BRESLOW)

clinical_and_supplementary_data_no_dups_met_bres$Breslow_Category <- ifelse(clinical_and_supplementary_data_no_dups_met_bres$breslow_thickness< 5, 'Thin', 'Thick' ) 
clinical_and_supplementary_data_no_dups_met_bres$Breslow_Category<- factor(clinical_and_supplementary_data_no_dups_met_bres$Breslow_Category,
                                                         levels = c("Thick","Thin"))
clinical_and_supplementary_data_no_dups_met_bres$Sample_Group <- factor(clinical_and_supplementary_data_no_dups_met_bres$Sample_Group,
                                                      levels = c('Low_PN', 'High_PN'))
##Remove duplicates by patient barcode    
##Save
save(clinical_and_supplementary_data_no_dups_met_bres, file ='Breslow_Thickness_Clinical_Data_met.Rdata')
write.xlsx(clinical_and_supplementary_data_no_dups_met_bres, file ='Breslow_Thickness_Clinical_Data_met.xlsx')

#####Box Plot Breslow Thickness
j <- ggplot(clinical_and_supplementary_data_no_dups_met_bres , aes(x= Paper_Sample_Group, y= breslow_thickness, fill= Paper_Sample_Group) )+ 
  geom_boxplot() +
  
  scale_fill_manual( name = "Sample Group", values=c(  "#7bd5d5",
                                                       "#7f4e97" )) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  #scale_fill_manual(values=c( '', '')) +
  # scale_x_discrete(labels=c('1'=' Cluster 1' ,  '2' ='Cluster 2'  )) +
  stat_compare_means(aes(group = Paper_Sample_Group), 
                     label = "p.format", size = 8,
                     label.y =45)  +
  # geom_jitter() +
  
  labs(
    #title="Age at Diagnosis",
    x ="Sample Group", y = "Breslow Thickness") +
  theme(
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20) ,
    # panel.background = element_rect(fill = "white"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16,face="bold"),
    panel.background = element_rect(fill = 'white', colour = 'black')
    
  ) 
pdf(("Breslow_THickness_boxplot_Met_Jan.pdf"), width = 8, height = 5)

print(j)
dev.off()



