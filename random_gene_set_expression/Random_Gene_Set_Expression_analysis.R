## This script produces random lists of 427 genes and produces expression heatmaps to compare them


## 
# install.packages("pheatmap")
# install.packages("d3heatmap")
# install.packages("survival")
# install.packages("survminer")
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(sva)
library(matrixStats)
library(data.table)
library(dplyr)
library(reshape)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(factoextra)
library(MASS)
library(mclust)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms

library(dendextend) # for comparing two dendrograms
library(survival)
library(survminer)
library(tidyr)

### create moveme function to move columns
####  https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}



##### SET WORKING DIRECTORY , use getwd to designate new variable'current_dir' to use in sub directory names

setwd("~/iCloud Drive (Archive)/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()


#load Rdata ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# read csv of SKCM FKPM to sub folder
FKPM_data<- loadRData( file =  'SKCM_FKPM.Rdata')

allgenes <- FKPM_data$hgnc_symbol
# load SKCM clinical info ####

clinical_info_df <- loadRData('SKCM_clinical_info.Rdata')

# load primary sample cluster allocation table ####
# 
primary_clusters_df<- loadRData('Primary_SKCM_Samples_Cluster_Allocation.Rdata')


PN_cluster_list <- read.csv ('Primary_and_Metastatic_Gene_Clusters_June.csv')
# Create new subfolder to store heatmaps

#dir.create("Random_Gene_Heatmaps")
#set counter to 0 then instruct programme to loop until counter is 10
counter<- 16
#while(counter<=11){
  counter<- counter+1
  
  print (paste0('Counter = ', counter))

# select n random genes from columns list 

Random_Genes <- sample (allgenes, size=427, replace =F)

###Random ####

# reduce table to include only random genes 
Random_Gene_Expression <-  FKPM_data[FKPM_data$hgnc_symbol %in%Random_Genes ,]

# Remove Entrez id column 

Random_Gene_Expression <-Random_Gene_Expression[, -(2) ]

# remove duplicate rows by hgnc symbol
Random_Gene_Expression = Random_Gene_Expression[!duplicated(Random_Gene_Expression$hgnc_symbol),]


# convert hgnc_symbol column to row name

Random_Gene_Expression <- data.frame(Random_Gene_Expression, row.names = 1)

# transpose table so cases are rows and genes are columns
Random_Gene_exp <- t(Random_Gene_Expression)

# add rownames and column names back in

colnames(Random_Gene_exp) <- rownames(Random_Gene_Expression)
rownames(Random_Gene_exp) <- colnames(Random_Gene_Expression)


#  order columns alphabetically
Random_Gene_exp<-Random_Gene_exp[,order(colnames(Random_Gene_exp))]

#transpose so samples are columns and genes are rows
Random_exp_info_t <- t(Random_Gene_exp)

#log transform

Random_Exp_log <- log2(Random_exp_info_t + 1)

# Create new matrix of Z scores

Random_Exp_Z_Score <- (apply(Random_Exp_log, 1, function(x) (x - mean(x)) / sd(x)))

# convert to Dataframe

Random_exp_df <- as.data.frame(Random_Exp_Z_Score)

#convert rownames to column

Random_exp_df$Sample <- rownames(Random_exp_df)

## move Sample first
Random_exp_df <-Random_exp_df[moveme(names(Random_exp_df), "Sample first")]

##### Merge clinicial info df with expression ####

Random_clinical_exp_df <- merge(clinical_info_df,Random_exp_df,
                               by.x="Sample",by.y= "Sample" )


########### Create table of just primary sample expression ####
# # select just primary samples 

Random_exp_info_prim <-  subset(Random_clinical_exp_df,Random_clinical_exp_df$SampleTypeCateg == 'Cancer (Primary)')


#remove columns with NA
Random_exp_info_prim <-Random_exp_info_prim [ ,colSums(is.na(Random_exp_info_prim )) == 0]

## merge with primary cluster column

Random_exp_info_prim_cluster <- cbind(Random_exp_info_prim, primary_clusters_df ) 

## move Cluster after Sample 
Random_exp_info_prim_cluster <-Random_exp_info_prim_cluster[moveme(names(Random_exp_info_prim_cluster), "Cluster after Sample")]

## Create new dataframe of genes with column to say if they are in PN list
#create list of UV_Down genes
Random_genes <- names(Random_exp_df)
## remove genes not in primary exp df 
Random_genes <- intersect(Random_genes,names(Random_exp_info_prim_cluster))


# convert to dataframe

Random_genes_df <- data.frame(Random_genes )
# Merge with PN List gene cluster lists

Random_genes_PN_df <- merge (Random_genes_df,PN_cluster_list,by.x = 'Random_genes',by.y = 'Gene.Symbol',  all.x=TRUE)
# remove extra columns
Random_genes_PN_df <- subset(Random_genes_PN_df , select = c(1,12))


######### HEATMAPS ######### 

#create matrix for heatmap (gene expression columns transposed so genes are rows and samples are columns )
heatmap_matrix <- 
  as.matrix(t(Random_exp_info_prim_cluster[,57:(length(names(Random_exp_info_prim_cluster)))]))


# ________________________________________
#Create Colour function to distribute colours more evenly ####

heatcolS <- circlize::colorRamp2 (breaks = c(-7,-0.1,0,0.1, 2.5,5, 20),
                                  c("dodgerblue4","#f1fbfe" ,"white", "#ffff55", "#ff4015", "#ff0000", '#5e0406'),
                                  transparency = .3)
heatcolS(seq(-7, 20, 0.01))

# Create top annotation

ann2 <- data.frame (primary_clusters_df$Cluster)


colnames(ann2) <- c('Sample_Cluster_Primary')

ann2$Sample_Cluster_Primary <- factor(ann2$Sample_Cluster_Primary,
                                      levels = c( "SKCM(P)_A","SKCM(P)_B"))

colours2 <- list('Sample_Cluster_Primary' = c('SKCM(P)_A' = 'darkorchid4', 'SKCM(P)_B' = 'gold')
                 
)


colAnn2 <- HeatmapAnnotation(df = ann2,
                             which = 'col',
                             col = colours2,
                             annotation_width = unit(c(2, 10), 'cm'),
                             show_annotation_name = FALSE,
                             annotation_legend_param = list(
                               Sample_Cluster_Primary = list(
                                 title = "Sample Cluster",
                                 #  title_position = "lefttop-rot",
                                 annotation_legend_side="right",
                                 annotation_legend_height = unit(160, "mm"),  
                                 annotation_legend_width = unit(60, "mm"), 
                                 grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
                                 labels_gp = gpar(fontsize = 20),
                                 title_gp = gpar(fontsize = 20, fontface = 'bold'))
                             )
)


rowAnn <- data.frame (Random_genes_PN_df$Random_genes ,Random_genes_PN_df$Prim_Gene_Cluster, stringsAsFactors = FALSE)
rowAnn   <- data.frame(rowAnn, row.names = 1)
# remove'Sample' row
# create list of rownames in rowann not in matrix

remove<- setdiff(rownames(rowAnn), rownames(heatmap_matrix))
rowAnn   <- as.data.frame(rowAnn[!rownames(rowAnn) %in% remove, ])




colnames(rowAnn) <- c('PN_Cluster')
rowAnn$PN_Cluster <-factor(rowAnn$PN_Cluster, levels = c('PN_1(P)','PN_2(P)', 'PN_3(P)' ))


colours3 <- list('PN_Cluster' = c('PN_1(P)' = 'lightslateblue','PN_2(P)' = 'pink', 'PN_3(P)' = 'darkcyan'))

rowAnn2 <- HeatmapAnnotation(df = rowAnn,
                             which = 'row',
                             col = colours3,
                             annotation_width = unit(c(2, 10), 'cm'),
                             show_annotation_name = FALSE,
                             annotation_legend_param = list(
                               PN_Cluster = list(
                                 #   title_position = "lefttop-rot",
                                 title = "PN_Cluster",
                                 annotation_legend_side="right",
                                 annotation_legend_height = unit(160, "mm"),
                                 annotation_legend_width = unit(60, "mm"),
                                 grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
                                 labels_gp = gpar(fontsize = 20),
                                 title_gp = gpar(fontsize = 20, fontface = 'bold'))
                             ))

## Heatmap with ward clustering of samples

pdf (paste0("Random_Gene_Heatmaps/Random_heatmap_cancer_SKCM_gene_and_sample_clustered_",counter,".pdf"),width=17, height=7)
Heatmap (heatmap_matrix, width = unit(260, "mm"), height = unit(150, "mm"),
         # name ='',
         rect_gp = gpar(col= FALSE),
         
         #   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
         # column_names_centered = TRUE,
         #column_names_rot = 90,
         clustering_method_rows = "ward.D2",
         clustering_method_columns = "ward.D2",
         # column_names_gp = gpar(fontsize = 25, fontface = "bold"),
         row_names_gp = gpar(fontsize = 3, fontface = "bold"),
         column_dend_height = unit(10,'mm'),
         row_dend_width=unit(10,"mm"),
         #column_title = "Expression of Proteostasis Network Genes in Primary Skin Cutaneous Melanoma Samples",
         #column_title = "A",
         show_column_names = FALSE,
         show_row_names = FALSE,
         column_split = 2,
         row_split = 3,
         row_gap = unit(3, "mm"),
         column_gap = unit(3, "mm"),
         col =  heatcolS,
         top_annotation=colAnn2,
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
         
)
#ggsave(paste0("Random_Gene_Heatmaps/Random_heatmap_cancer_SKCM_gene_and_sample_clustered_",counter,".pdf", HM))
dev.off()

}

# #identify genes in each cluster
# 
# 
# # methods to assess
# m <- c( "average", "single", "complete", "ward")
# names(m) <- c( "average", "single", "complete", "ward")
# 
# # function to compute coefficient
# ac <- function(x) {
#   agnes(heatmap_matrix, method = x)$ac
# }
# 
# map_dbl(m, ac)
# 
# 
# ### ward method gives highest clustering coefficient (0.935)
# 
# hc3 <- agnes(heatmap_matrix, method = "ward")
# pltree(hc3, cex = 0.1, hang = -1, main = "Dendrogram of Genes")
# 
# # Cut tree into 3 groDowns
# sub_grp <- cutree(hc3, k = 3)
# 
# # Number of members in each cluster
# table(sub_grp)
# 
# #Add gene cluster number 
# heatmap_matrix_df <-heatmap_matrix
# heatmap_matrix_df$Cluster  <- sub_grp 
# 
# #Renumber clusters to match heatmap
# 
# 
# heatmap_matrix_df$Cluster <- sapply(heatmap_matrix_df$Cluster, function(x)
#   ifelse ( x %in% c('1'),"2",
#            ifelse (x %in% c('3' ),"1", '3')))
# 
# # convert cluster column to data frame
# gene_cluster_df <- as.data.frame(heatmap_matrix_df$Cluster)
# 
# #Convert Row Names to Column
# gene_cluster_df$Gene <- rownames(gene_cluster_df)
# 
# 
# #save Lindquist gene clusters 
# 
# save(gene_cluster_df, file = 'Random_only_gene_clusters.Rdata')
# 
# 
# write.csv(gene_cluster_df, file = 'Random_only_gene_clusters.csv')
