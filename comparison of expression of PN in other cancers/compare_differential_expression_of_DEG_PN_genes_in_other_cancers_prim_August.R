##This script compares expression of key PN genes in different 
#cancers and produces a heatmap of expression for each cancer and a list of differentially expressed PN genes
### for each cancer

#install.packages("survminer")
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library("readxl")
library('xlsx')

#load Rdata ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##### SET WORKING DIRECTORY , use getwd to designate new variable'current_dir' to use in sub directory names

setwd("~/iCloud_Drive_Archive/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()

# #### Read in list of PN Genes
PN_Gene_Clusters <- read_xlsx('Difference_in_Exp_of_PN_Genes_in_primary_samples_corrected_expression_groups_August.xlsx', sheet = 1)
PN_Genes <- PN_Gene_Clusters$Gene

## Read in list of key SKCM PN genes
sig_SKCM_genes_df<- read_xlsx('TCGA_primary_PN_DSeq_Expression_Difference_edited.xlsx', sheet = 2)
sig_SKCM_genes <- subset(sig_SKCM_genes_df$Gene, sig_SKCM_genes_df$Expression %in% c("Lower","Higher"))
higher_sig_SKCM_genes <- subset(sig_SKCM_genes_df$Gene, sig_SKCM_genes_df$Expression == "Higher" )
lower_sig_SKCM_genes <- subset(sig_SKCM_genes_df$Gene, sig_SKCM_genes_df$Expression == "Lower")
#### Create list of studies
studies <- read.csv('TCGA_study_list.csv')
studies <- studies[,1]

## create completed studies vector
completed_studies<- vector()
####Other Cancers Loop ####
#### loop to produce heatmaps of expression of PN genes and table showing mean SKCM-up regulated 
##and down regulated genes for each sample then calculate if their is a significant difference in expression 
##for the genes in the allocated sample groups


Differentially_expressed_SKCM_Key_genes_df <- data.frame (Study= character(0), 
                                                          Lower_expressed_SKCM_Lower_genes=numeric(0), Higher_expressed_SKCM_Higher_genes=numeric(0))

Cancer_expression_df<- as.data.frame(sig_SKCM_genes)
Cancer_expression_df<- as.data.frame(Cancer_expression_df[-1,])



DEG_cancer_table_df <- data.frame (s= character(0), 
                                        Higher_genes_overlap=numeric(0),Lower_genes_overlap=numeric(0))


#set s to 0 then instruct programme to loop while s < number of columns of data
# for(s in studies){
#### start of loop #####
for (s in studies) {
print(s)
try({
  

#SKCM significant genes heatmap
#s<-'UVM'
  
 PN_exp<- read.csv(paste0('TCGA_data/',s,'/', s,'_PN_Gene_exp.csv'))  
#### Rename column X "Sample"
 names( PN_exp)[names( PN_exp) == 'X'] <- 'Sample'
 ### Select rows where sample begins with "TCGA" to get rid of extraneous rows
 
 PN_exp<-subset (PN_exp, Sample %like% "TCGA")
 
 ###Create list of sample type based on sample name
 PN_exp$SampleType<- sapply( PN_exp$Sample, function(x)
   strsplit(as.character(x),"\\.")[[1]][4])
 
 # remove non numeric elements of SampleType (letter refers to sample vial number)
 PN_exp$SampleType <- gsub("[^0-9]", "", PN_exp$SampleType)
 
 # Move SampleType after Patient ID 
 PN_exp <- PN_exp[moveme(names(PN_exp), "SampleType after Sample")]
 
 ####Select rows where sample type is 01 or 03 or 05   (primary samples)
 
 PN_exp<-subset (PN_exp, SampleType %in% c("01","03","05"))
 
 
 ###remove duplicates
 
 PN_exp <- PN_exp[!duplicated( PN_exp$Sample ), ]
 ###convert first column to rowname
 PN_exp <- data.frame( PN_exp, row.names = 1) 
 
##### select significant genes 
 Sig_PN_exp <- (subset(PN_exp, 
         select =c(names(PN_exp) %in% sig_SKCM_genes)) )
 
 
  #  order columns alphabetically
  Sig_PN_exp<-Sig_PN_exp[,order(colnames(Sig_PN_exp))]
  
  #transpose so samples are columns and genes are rows
  Sig_PN_exp_info_t <- t(Sig_PN_exp)
  
  Sig_PN_exp_info_t_num <- matrix(as.numeric(Sig_PN_exp_info_t),    # Convert to numeric matrix
                    ncol = ncol(Sig_PN_exp_info_t))
  ### add gene and sample names back in
  colnames(Sig_PN_exp_info_t_num) <- colnames(Sig_PN_exp_info_t )
  rownames(Sig_PN_exp_info_t_num) <- rownames(Sig_PN_exp_info_t )
  
  
  #log transform
  
  Sig_PN_Exp_log <- log2(Sig_PN_exp_info_t_num + 1)
  
  # Create new matrix of Z scores
  
  Sig_PN_Exp_Z_Score <- (apply(Sig_PN_Exp_log, 1, function(x) (x - mean(x)) / sd(x)))
  
  # convert to Dataframe
  
  Sig_PN_exp_df <- as.data.frame(Sig_PN_Exp_Z_Score)
  
  #remove columns with NA
  Sig_PN_exp_df <-Sig_PN_exp_df [ ,colSums(is.na(Sig_PN_exp_df)) == 0]
  
 # # ## remove genes not in primary exp df 
 Sig_PN_genes_table <- subset(sig_SKCM_genes_df, sig_SKCM_genes_df$Gene %in% names( Sig_PN_exp_df))
 
    # ### Sort so genes are in the same order as the genes in the heatmap
 Sig_PN_genes_table <- Sig_PN_genes_table[match(names(Sig_PN_exp), Sig_PN_genes_table$Gene),]   
  # 
  ######### HEATMAPS ######### 
  
  #create matrix for heatmap (gene expression columns transposed so genes are rows and samples are columns )
  heatmap_matrix <- 
    as.matrix(t(Sig_PN_exp_df))
  
  
  # ________________________________________
  #Create Colour function to distribute colours more evenly ####
  
  heatcolS <- circlize::colorRamp2 (breaks = c(-7,-0.1,0,0.1, 2.5,5, 20),
                                    c("dodgerblue4","#f1fbfe" ,"white", "#ffff55", "#ff4015", "#ff0000", '#5e0406'),
                                    transparency = .3)
  heatcolS(seq(-7, 20, 0.01))
  
    rowAnn <- data.frame (   Sig_PN_genes_table$Expression  , stringsAsFactors = FALSE)

  # # remove'Sample' row
  # create list of rownames in rowann not in matrix
  colnames(rowAnn) <- c('Expression_in_PN_Low_vs_PN')
  rowAnn$Expression_in_PN_Low_vs_PN <-factor(rowAnn$Expression_in_PN_Low_vs_PN, levels = c('Lower','Higher' ))


  colours3 <- list('Expression_in_PN_Low_vs_PN' = c('Lower' = '#3b3bcf','Higher' = '#cf3b3b'))

  rowAnn2 <- HeatmapAnnotation(df = rowAnn,
                               which = 'row',
                               col = colours3,
                               annotation_width = unit(c(2, 10), 'cm'),
                               show_annotation_name = FALSE,
                               annotation_legend_param = list(
                                 Expression_in_PN_Low_vs_PN = list(
                                   #   title_position = "lefttop-rot",
                                   title = "Difference in Expression in\nSKCM Low PN Primary, High PN Primary ",
                                   annotation_legend_side="right",
                                   annotation_legend_height = unit(160, "mm"),
                                   annotation_legend_width = unit(60, "mm"),
                                   grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
                                   labels_gp = gpar(fontsize = 20),
                                   title_gp = gpar(fontsize = 20, fontface = 'bold'))
                               ))
  #
  # ## Heatmap with ward clustering of samples and genes
  
     HM<-(Heatmap (heatmap_matrix, width = unit(260, "mm"), height = unit(150, "mm"),
                 # name ='',
                 rect_gp = gpar(col= FALSE),
                 
                 column_title_gp = gpar(fontsize = 30, fontface = "bold"),
                 # column_names_centered = TRUE,
                 #column_names_rot = 90,
                 clustering_method_rows = "ward.D2",
                 clustering_method_columns = "ward.D2",
                 # column_names_gp = gpar(fontsize = 25, fontface = "bold"),
                 row_names_gp = gpar(fontsize = 3, fontface = "bold"),
                 column_dend_height = unit(10,'mm'),
                 row_dend_width=unit(10,"mm"),
                 column_title = paste0(s),
                 
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 column_split = 2,
                #row_split = 2,
               row_gap = unit(5, "mm"),
                 column_gap = unit(5, "mm"),
                 col =  heatcolS,
                 # top_annotation=colAnn2,
              right_annotation = rowAnn2 ,
                 cluster_column_slices = FALSE,
                 cluster_row_slices = FALSE,
                 
                 heatmap_legend_param = list(
                   title = "Expression", 
                   # title_position = "lefttop-rot",
                   color_bar = "continuous",
                   legend_height = unit(120, "mm"),  
                   legend_width = unit(20, "mm"), 
                   fonts = 20, 
                   fontface = "bold",
                   title_gp = gpar(fontsize = 20, fontface = 'bold'),
                   labels_gp = gpar(col = "black", fontsize = 20),
                   heatmap_legend_side="right")
                 
  ))
  #ggsave(paste0("Sig_PN_Gene_Heatmaps/Sig_PN_heatmap_cancer_SKCM_gene_and_sample_clustered_",counter,".pdf", HM))
   
     pdf (paste0(s,"_PN_Gene_Expression_nov.pdf"),width=17.6, height=9)
     print(HM)
   dev.off()
  
  
   ### Extract column clusters
   
   ### Extract cluster allocation from heatmap
   c.dend <- column_dend(HM)  #If needed, extract row dendrogram
   ccl.list <- column_order(HM)  #Extract clusters (output is a list)
   
   lapply(ccl.list, function(x) length(x))  #check/confirm size gene clusters
   
   # loop to extract samples for each cluster.
   t_heatmap_matrix<- t(heatmap_matrix)
   for (i in 1:length(column_order(HM))){
     if (i == 1) {
       clu <- t(t(row.names(t_heatmap_matrix[column_order(HM)[[i]],])))
       out_col <- cbind(clu, paste("cluster", i, sep=""))
       colnames(out_col) <- c("Sample", "Cluster")
     } else {
       clu <- t(t(row.names(t_heatmap_matrix[column_order(HM)[[i]],])))
       clu <- cbind(clu, paste("cluster", i, sep=""))
       out_col <- rbind(out_col, clu)
     }
   }
   
   
   #check
   out_col
   
   new_sample_clusters  <- out_col
   
## convert out_col Sample Column to rowname
out_col <- data.frame(out_col, row.names = 1)   

###bind sample cluster allocation to PN Exp Table
Sig_PN_exp_cluster <- merge(Sig_PN_exp,out_col, by.x = 0, by.y = 0)

# ##### Save new sample clusters
save(Sig_PN_exp_cluster ,  file = paste( 'all_cancers/',s,'SKCM_sig_gene_expression_sample_groups_lower__oct.RData'))
write.csv( Sig_PN_exp_cluster , file = paste( 'all_cancers/',s,"SKCM_sig_gene_expression_sample_groups_lower__oct.csv"))

#### Identify how many genes are differentially expressed between clusters ####


Expression_difference_df <- data.frame (Gene= character(0), 
                                        P.value=numeric(0),Foldchange_A_over_B=numeric(0),mean_A =numeric(0), mean_B =numeric(0))

Genes<- (names(Sig_PN_exp))
n =0

while(n<=(length(Genes)-1)) {
  n =n+1
  print(n)
  Gene <- Genes[n]
  print(Gene)
  A <-subset (Sig_PN_exp_cluster[,Gene],(Sig_PN_exp_cluster$Cluster == "cluster2"))
  B<-subset (Sig_PN_exp_cluster[,Gene],(Sig_PN_exp_cluster$Cluster == "cluster1"))
  Foldchange_A_over_B <- (mean(as.numeric(A))/mean(as.numeric(B)))
  mean_A<-mean(as.numeric(A))
  mean_A<-format(mean_A, round(mean_A, 3),scientific = FALSE)
  mean_A<- as.numeric(mean_A)
  mean_B <-  mean(as.numeric(B))
  mean_B <-  format(mean_B, round(mean_A, 3),scientific = FALSE)
  mean_B<- as.numeric(mean_B)
  results = t.test((as.numeric(A)),(as.numeric(B)))
  P.value<-results$p.value
  formatC(P.value, format = "e", digits = 3)
  P.value<-as.numeric(P.value)
  temp_df <- data.frame (Gene, P.value,Foldchange_A_over_B,mean_A,mean_B)
  # temp_df$TF <- TF
  # temp_df$P.value <- results$p.value
  # temp_df$Foldchange_A_over_B <- Foldchange
  Expression_difference_df<- rbind(Expression_difference_df, temp_df)
}


##Add columnn listing higher mean and adjusted p value
Expression_difference_df$max_mean <- pmax(Expression_difference_df$mean_A,
                                          Expression_difference_df$mean_B)
Expression_difference_df$P.value <- formatC(Expression_difference_df$P.value, format = "e", digits = 3)

Expression_difference_df$P.value <- as.numeric(Expression_difference_df$P.value)

Expression_difference_df$p_adjust <- p.adjust(Expression_difference_df$P.value, method =  "hochberg")
Expression_difference_df$p_adjust <- formatC(Expression_difference_df$p_adjust, format = "e", digits = 3)

Expression_difference_df$p_adjust <- as.numeric(Expression_difference_df$p_adjust)


# ## sort expression difference table by p value and foldchange and alphabetically
# 
# Expression_difference_df <- Expression_difference_df[
#   order( desc(Expression_difference_df$max_mean) ,  Expression_difference_df$P.value, desc(Expression_difference_df$Foldchange_A_over_B) ),
#]

#Save expression data
save(Expression_difference_df, file = paste('all_cancers/',s,"_Difference_in_Expression_of_SKCM_Key_lower_Genes_oct.Rdata"))
write.xlsx(Expression_difference_df, file = paste('all_cancers/',s,"_Difference_in_Expression_of_SKCM_Key_lower_Genes_oct.xlsx"))

##Create list of Genes with significantly different expression and significant level FPKM >2
Sig_Gene_exp_X<- as.data.frame(subset (Expression_difference_df,
                       (Expression_difference_df$p_adjust <= 0.1) 
                      & (Foldchange_A_over_B <= 0.8)
                         # Expression_difference_df$max_mean >= 2
                        ))
Sig_Gene_exp_Y<- as.data.frame(subset (Expression_difference_df,
                         (Expression_difference_df$p_adjust <= 0.1) 
                         & (Foldchange_A_over_B >= 1.25))
                         # Expression_difference_df$max_mean >= 2
)

#### Save sig Gene list

Sig_Gene_exp_X$Expression <-if (nrow(Sig_Gene_exp_X)>nrow(Sig_Gene_exp_Y)){
  "Lower"
  }else{
    "Higher"}

Sig_Gene_exp_Y$Expression <-if (nrow(Sig_Gene_exp_X)<nrow(Sig_Gene_exp_Y)){
  "Lower"
}else{
  "Higher"}

Sig_Gene_exp<- rbind(Sig_Gene_exp_X,Sig_Gene_exp_Y)

Higher_genes <-if (nrow(Sig_Gene_exp_X)<nrow(Sig_Gene_exp_Y)){
  Sig_Gene_exp_X$Gene
}else{
  Sig_Gene_exp_Y$Gene}

Higher_genes_overlap <-intersect(Higher_genes,higher_sig_SKCM_genes)

Lower_genes <-if (nrow(Sig_Gene_exp_X)>nrow(Sig_Gene_exp_Y)){
  Sig_Gene_exp_X$`Gene`
}else{
  Sig_Gene_exp_Y$Gene}

Lower_genes_overlap <-intersect(Lower_genes,lower_sig_SKCM_genes)

save(Sig_Gene_exp, file = paste('all_cancers/',s,"significant_lower_PN SKCM_Key_Genes_Primary_SKCM_oct.Rdata"))
write.csv(Sig_Gene_exp, file = paste('all_cancers/',s,"significant_lower_PN SKCM_Key_Genes_Primary_SKCM_oct.csv"))

##Create temp data frame of study and number of significantly differentially expressed genes

temp_df_sig <- data.frame (s, length(Lower_genes_overlap),length(Higher_genes_overlap))

DEG_cancer_table_df <-rbind(DEG_cancer_table_df,temp_df_sig)

# Differentially_expressed_SKCM_Key_genes_df<- rbind(Differentially_expressed_SKCM_Key_genes_df, temp_df_sig)

#### Add study name to list of completed studies
completed_studies<- append(completed_studies,s) 

#### create new data frame of expression expression column to gene expression table

Expression_difference_column<- as.data.frame(ifelse(Expression_difference_df$Gene %in% Higher_genes_overlap,"Higher",
                                             ifelse(Expression_difference_df$Gene %in% Lower_genes_overlap,"Lower", "Other"   ) ) ) 
#rename column
colnames(Expression_difference_column)<- s

Differentially_expressed_SKCM_Key_genes_df$low_percent<- Differentially_expressed_SKCM_Key_genes_df[,2]/1.3
Differentially_expressed_SKCM_Key_genes_df$high_percent<- Differentially_expressed_SKCM_Key_genes_df[,3]/1.3



# ##### Save table of number of lower differentially expressed genes
save(Differentially_expressed_SKCM_Key_genes_df ,  file = paste( 'Number_of_SKCM_Key_lower_genes_Diff_exp_in_each_study_oct.RData'))
write.xlsx( Differentially_expressed_SKCM_Key_genes_df , file = paste( 'Number_of_SKCM_Key_lower_genes_Diff_exp_in_each_study_oct.xlsx'))

###bind column to cancer express

Cancer_expression_df<- cbind(Cancer_expression_df,Expression_difference_column )



})}


Cancer_expression_df[1] <- NULL

# ##### Save table of number of lower differentially expressed genes
save(Cancer_expression_df ,  file = paste( 'Number_of_DEGs_other_cancers.RData'))
write.xlsx( Cancer_expression_df , file = paste( 'Number_of_DEGs_other_cancers.xlsx'))

####save table of number of DEGs
###rename columns 
colnames(DEG_cancer_table_df) <-  c("Study","Lower", "Higher")
###sort table by lower then higher
DEG_cancer_table_df_order <- DEG_cancer_table_df[order(-DEG_cancer_table_df$Lower, -DEG_cancer_table_df$Higher),]


save(DEG_cancer_table_df ,  file = paste( 'Number_of_DEGs_other_cancers_table.RData'))
write.xlsx( DEG_cancer_table_df , file = paste( 'Number_of_DEGs_other_cancers_table.xlsx'))

####Read in table of DEGs ####


DEG_cancer_table_df <- read.xlsx(file = ( 'Number_of_DEGs_other_cancers_table.xlsx'), sheetIndex = 2)
DEG_cancer_table_df<- DEG_cancer_table_df[,1:3]

####convert from wide table to long
library(data.table)
DEG_cancer_table_df_long <- melt(setDT( DEG_cancer_table_df), id.vars = c("Study"))
###rename variable column "Expression"
names(DEG_cancer_table_df_long)[2] <- "Expression"
DEG_cancer_table_df_long$percent <- DEG_cancer_table_df_long$value/1.35

DEG_cancer_table_df_long$Expression<- factor(DEG_cancer_table_df_long$Expression,
                                             levels = c("Higher", "Lower"))

  save(DEG_cancer_table_df_long ,  file = paste( 'Number_of_DEGs_other_cancers_table_long.RData'))
write.xlsx( DEG_cancer_table_df_long , file = paste( 'Number_of_DEGs_other_cancers_table_long.xlsx'))


#####:::::::::::::::::::::::::::::::::::::
## plot number of DEGs ####
require(forcats)
DEG_barchart <- ggplot(DEG_cancer_table_df_long, aes(x = reorder(DEG_cancer_table_df_long$Study, -percent), y = percent, fill = DEG_cancer_table_df_long$Expression))+
  geom_bar(stat="identity") +
  theme_bw() +
  
  scale_fill_manual(name = "Expression in\nA vs B",
                    labels = c( "Higher", "Lower"),
                    values=c(c("Higher" = "#cf3b3b",
                               "Lower"  = "#3b3bcf")))+

  theme(text = element_text(size=50 ,face = "bold"))+
  xlab("")+
  ylab("Proportion of Genes")+
# labs(fill="Gene Set")
# coord_flip()
# ggtitle
# 
theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 20),
        panel.border = element_rect(color = "black", fill=NA ),
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 20),
        axis.text.y = element_text( size = 20),
        axis.title.y = element_text( size = 22)
  )+
  labs( 
    x ="", y = 
      "Percentage of CM DEGs with same expression pattern") + 
  # scale_y_discrete(
  #   labels = scales::number_format(accuracy = 1))
  
pdf ("DEGs_other_cancers_barchart_primary_August.pdf",width=14, height=9 )
print(DEG_barchart)
dev.off()






      
      