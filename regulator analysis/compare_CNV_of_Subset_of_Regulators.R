### This code creates a faceted bar chart showing frequency of copy number variations in reguators, 
###based on data downloaded from cBioPortal


library('ggpubr')
library('ggplot2')
library(reshape2)
library (forcats)
library (xlsx)
setwd("~/iCloud_Drive_Archive/Desktop/R_Proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()

###load CNV data downloaded from cbioportal
 CNV_df<- read.xlsx('regenrich_CNV_data.xlsx', sheetIndex = 1)
 
 ####load group allocations
 
 ###read in list of sample group allocations  ###
 sample_allocations<- read.xlsx('Final_Expression_primary_and_metastatic_samples_with_groups.xlsx', sheetName = "Metastatic")
 
 ##merge regulator mutations with sample allocations
 
 CNV_df_with_groups <- merge(CNV_df,sample_allocations, by.x ='SAMPLE_ID', by.y = 'Sample')
 
 #### save mutations table
 
 write.xlsx(CNV_df_with_groups,file = "metastatic_CNV_with_sample_groups_July.xlsx")
 
 
# # convert first column to rowname
#  CNV_df<- data.frame(CNV_df, row.names = 1)
 
 ###melt into 3 columns
 
 CNV_table_melt <- as.data.frame(melt(CNV_df_with_groups, id= c(1,37)))
 

 
 ### Rename columns in melted dataframe
 colnames(CNV_table_melt)[colnames(CNV_table_melt)=="variable"] <- "Regulator"
 colnames(CNV_table_melt)[colnames(CNV_table_melt)=="value"] <- "CNV"
 
 
 ##Select subset of TFs
 CNV_table_melt <- subset( CNV_table_melt, Regulator %in% c('NR3C1','SP4', 'CREB1', 'GABPA','ATF2',  'MEF2A','ZFX', 'BCLAF1'))
 
 
 
 CNV_table_melt$Sample_Group <-factor(CNV_table_melt$Final_Sample_Group_Paper,
                                                  levels = c("Metastatic (A)", "Metastatic (B)") )
 
 CNV_table_melt$Sample_Group<- ifelse( CNV_table_melt$Sample_Group == "Metastatic (A)", "Metastatic A", "Metastatic B"  )
 
 CNV_table_melt$CNV <-factor(CNV_table_melt$CNV,
                                      levels = c( -2,-1, 1, 2, 0 ) )
 CNV_table_melt$CNV <-fct_rev( CNV_table_melt$CNV)
 
 CNV_table_melt<- as.data.frame( CNV_table_melt)
 ## create bar chart showing CNV

v <- ggplot(CNV_table_melt, aes(Sample_Group, fill = CNV))+
   geom_bar(position="fill") +
   scale_fill_manual(name = "", values=c("0" = "#d8d8d8",
                              "-2"  = "#3f93ce",
                              "-1"  = "#c5dff0",
                              "1"  = "#fca0ad",
                              "2"  = "#F8415A"))+
  facet_grid(~ Regulator)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20) ,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, face="bold"),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, 
                                   hjust=1,
                                   size =20,face="bold"),
        axis.text.y = element_text( size =20,face="bold"),
        panel.background = element_rect(fill = 'white', colour = 'black'),
                                        strip.text.x = element_text(size = 20),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
ylab("Proportion of samples")
  pdf ("CNV_in_subset_Metastatic_samples_margin.pdf", width=15, height=8)
print(v)
dev.off()




#investigate if significant difference between CNV in sample groups
install.packages('tables')
library(tables)
m <-tabular(Regulator+Sample_Group~CNV, data=CNV_table_melt)
 
###remove first colum
CNV_table <- CNV_table_melt[,-1]
m<-table(CNV_table) 
 
CNV_table_freq <- as.data.frame(m)
##convert from long to wide
library(tidyr)
CNV_table_wide <-spread(CNV_table_freq, CNV, Freq)

write.xlsx(CNV_table_wide, file = "Frequency of CNV in Regulators between_Sample_groups.xlsx")
 
## loop to compare frequencies of CNV in sample groups 

Regulators<- CNV_table_wide$Regulator[1:35]
Regulators<- as.factor(Regulators)


CNV_fisher_df <- data.frame (Regulator= character(0),   P.value=numeric(0))

n =0
 
while(n<=(length(Regulators)-1)) {
  n =n+1
  print(n)
  Regulator <- Regulators[n]
  print(Regulator)
m = matrix(c(CNV_table_wide$`0`[n],CNV_table_wide$`2` [n],CNV_table_wide$`1`[n],CNV_table_wide$`-1` [n],CNV_table_wide$`-2`[n],
             CNV_table_wide$`0`[n+35],CNV_table_wide$`2` [n+35],CNV_table_wide$`1`[n+35],CNV_table_wide$`-1` [n+35],CNV_table_wide$`-2`[n+35]), ncol =2)
  mt<-t(m)         
results<-fisher.test(mt, workspace = 2e8)    
P.value<-results$p.value
           
           temp_df <- data.frame (Regulator, P.value)
           CNV_fisher_df<- rbind(CNV_fisher_df, temp_df)
}     

CNV_fisher_df$p_adjust <- p.adjust(CNV_fisher_df$P.value, method =  "BH")
CNV_fisher_df$p_adjust <- formatC(CNV_fisher_df$p_adjust, format = "e", digits = 3)

CNV_fisher_df$p_adjust <- as.numeric(CNV_fisher_df$p_adjust)

### order fisher df by pvalue

CNV_fisher_df<-  CNV_fisher_df[order(CNV_fisher_df$p_adjust),]

write.xlsx(CNV_fisher_df,file = ("CNV_of_regulators_between_sample_groups_significance.xlsx"))


###m reorder barchart by significance

Sig_Regulators <- CNV_fisher_df$Regulator

## create bar chart showing CNV

w <- ggplot(CNV_table_melt, aes(Sample_Group, fill = CNV))+
  geom_bar(position="fill") +
  scale_fill_manual(values=c("0" = "#d8d8d8",
                             "-2"  = "#3f93ce",
                             "-1"  = "#c5dff0",
                             "1"  = "#fca0ad",
                             "2"  = "#F8415A"))+
  facet_grid(~factor(Regulator,levels= Sig_Regulators))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20) ,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, 
                                   hjust=1,
                                   size =13,face="bold"),
        axis.text.y = element_text( size =30,face="bold"),
        panel.background = element_rect(fill = 'white', colour = 'black') 
        
  )
pdf ("CNV_in_Metastatic_samples_by_significance.pdf", width=30, height=6)
print(w)
dev.off()
