
## This script producesproduces expression heatmaps of targets of TFs and MiRNAs


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
library("readxl")

library('xlsx')

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

##Create matrix of ENrichr sources
## Read in Matrices of targets of TFS

Encode_Targets<- read_excel("ENCODE_TF_ChIP-seq_Sig_PN_Genes.xls", sheet = "Targets")

CHEA_Targets<- read_excel("ChEA_2016_TF_SIg_PN_Genes.xls", sheet = "Targets")


#### Add column with number of genes in overlap
Encode_Targets$Overlap_N <- as.numeric(gsub("\\/[.,0-9]*$", "", Encode_Targets$Overlap))

### delete eGFP- in TF
Encode_Targets$TF<- gsub("eGFP-", "", Encode_Targets$TF)


###Sort by Overlap_N
Encode_Targets <- Encode_Targets[
  order( desc(Encode_Targets[,('Overlap_N')]) ),
  ]
####Delete duplicate entries for TFS

Encode_Targets <- Encode_Targets[!duplicated(Encode_Targets$TF), ]


### Delete columnns 'Overlap', 'Overlap_N'

Encode_Targets <- subset(Encode_Targets , select = -c(Overlap, Overlap_N) )

# ## Replace NA with 0
# Encode_Targets[is.na(Encode_Targets)] <- 0

Encode_Targets <- as.data.frame(Encode_Targets)
Encode_Targets_Rank <- Encode_Targets 
Encode_Targets_Rank$EncodeRank <- (row_number(Encode_Targets_Rank)-1)
#
#CHEA
### delete eGFP- in TF
CHEA_Targets$TF<- gsub("eGFP-", "", CHEA_Targets$TF)

#### Add column with number of genes in overlap
CHEA_Targets$Overlap_N <- as.numeric(gsub("\\/[.,0-9]*$", "", CHEA_Targets$Overlap))
###Sort by Overlap_N
CHEA_Targets <- CHEA_Targets[
  order( desc(CHEA_Targets[,('Overlap_N')]) ),
  ]
####Delete duplicate entries for TFS

CHEA_Targets <- CHEA_Targets[!duplicated(CHEA_Targets$TF), ]


CHEA_Targets <- subset(CHEA_Targets , select = -c(Overlap, Overlap_N) )

# ## Replace NA with 0
# CHEA_Targets[is.na(CHEA_Targets)] <- 0


CHEA_Targets <- as.data.frame(CHEA_Targets)

####Merge two date frames

TF_Targets <- merge(Encode_Targets ,CHEA_Targets, by = 'TF', all = TRUE)

## First column to rownames

TF_Targets  <- data.frame(TF_Targets , row.names = 1)   


## Transpose
TF_Targets <- t(TF_Targets)

#Delete first column
TF_Targets <- TF_Targets[,-c(1)]  

#### Create matrix showing which genes are regulated by each TF
##Read in list of significant PN Genes
Sig_PN_Gene_Table <- read_excel('Significant_PN_Genes.xls')

Sig_PN_Genes <-Sig_PN_Gene_Table$Gene
x <- colnames(TF_Targets)
y <- Sig_PN_Gene_Table$Gene
  #subset(PN_cluster_list$Gene.Symbol, PN_cluster_list$Prim_Gene_Cluster %in% c('PN_1(P)','PN_2(P)'))
Sig_TF_target_df<- as.data.frame( matrix(ncol = length(x), nrow=length(y), dimnames = list(y,x)) )
Sig_TF_target_df$PN_Genes <- rownames(Sig_TF_target_df)
Sig_TF_target_df <- Sig_TF_target_df[moveme(names(Sig_TF_target_df), "PN_Genes first")]

#### create loop to write 'X'  if PN Gene is in the list of targets for the TF colname

## Read in csv of targets of C1C2TFs 

#C1C2_TF_Targets <- read.csv ('C1C2_TF_Crosstarget.csv')
t=0
for(t in 0:((length(x))-1))
  
{
 t <-t+1
  print(t)
  TF <- x[t]
  TF_target_genes <- TF_Targets [,TF]
  Sig_TF_target_df[,TF]<- sapply(Sig_TF_target_df$PN_Genes, function(x)
    ifelse (x %in% TF_target_genes,1,0))
}

###Delete PN Genes column
Sig_TF_target_df <- 
subset(Sig_TF_target_df, select = -c(PN_Genes) )

Sig_TF_target_df <- t(Sig_TF_target_df)
####
write.xlsx( Sig_TF_target_df, file = "Chipseq_PN_TF_targets.xlsx")

##########


## Read in Matrices of targets of TFS

Encode_Targets_Ranked<- read_excel("ENCODE_TF_ChIP-seq_Sig_PN_Genes.xls", sheet = "Targets_by_combined_score")

CHEA_Targets_Ranked<- read_excel("ChEA_2016_TF_SIg_PN_Genes.xls", sheet = "Targets_by_combined_score")



### delete eGFP- in TF
Encode_Targets_Ranked$TF<- gsub("eGFP-", "", Encode_Targets_Ranked$TF)

####Delete duplicate entries for TFS

Encode_Targets_Ranked <- Encode_Targets_Ranked[!duplicated(Encode_Targets_Ranked$TF), ]


# ## Replace NA with 0
# Encode_Targets_Ranked[is.na(Encode_Targets_Ranked)] <- 0

Encode_Targets_Ranked <- as.data.frame(Encode_Targets_Ranked)
Encode_Rank <- Encode_Targets_Ranked[,c(1:2)]
#CHEA
### delete eGFP- in TF
CHEA_Targets_Ranked$TF<- gsub("eGFP-", "", CHEA_Targets_Ranked$TF)

####Delete duplicate entries for TFS

CHEA_Targets_Ranked <- CHEA_Targets_Ranked[!duplicated(CHEA_Targets_Ranked$TF), ]

CHEA_Rank <- CHEA_Targets_Ranked[,c(1:2)]
# ## Replace NA with 0
# CHEA_Targets_Ranked[is.na(CHEA_Targets_Ranked)] <- 0


CHEA_Targets_Ranked <- as.data.frame(CHEA_Targets_Ranked)

####Merge two date frames

TF_Targets_Ranked <- merge(Encode_Targets ,CHEA_Targets_Ranked, by = 'TF', all = TRUE)

## First column to rownames

TF_Targets_Ranked  <- data.frame(TF_Targets_Ranked , row.names = 1)   


## Transpose
TF_Targets_Ranked <- t(TF_Targets_Ranked)


##Merge the dataframe of rank
TF_Ranks <- merge(Encode_Rank, CHEA_Rank, by = 'TF', all = TRUE)

#### Create matrix showing which genes are regulated by each TF
##Read in list of significant PN Genes
Sig_PN_Gene_Table <- read_excel('Significant_PN_Genes.xls')

Sig_PN_Genes <-Sig_PN_Gene_Table$Gene
x <- colnames(TF_Targets_Ranked)
y <- Sig_PN_Gene_Table$Gene
#subset(PN_cluster_list$Gene.Symbol, PN_cluster_list$Prim_Gene_Cluster %in% c('PN_1(P)','PN_2(P)'))
Sig_Ranked_TF_target_df<- as.data.frame( matrix(ncol = length(x), nrow=length(y), dimnames = list(y,x)) )
Sig_Ranked_TF_target_df$PN_Genes <- rownames(Sig_Ranked_TF_target_df)
Sig_Ranked_TF_target_df <- Sig_Ranked_TF_target_df[moveme(names(Sig_Ranked_TF_target_df), "PN_Genes first")]

#### create loop to write 'X'  if PN Gene is in the list of targets for the TF colname

## Read in csv of targets of C1C2TFs 

#C1C2_TF_Targets_Ranked <- read.csv ('C1C2_TF_Crosstarget.csv')
t=0
for(t in 0:((length(x))-1))
  
{
  t <-t+1
  print(t)
  TF <- x[t]
  TF_target_genes <- TF_Targets_Ranked [,TF]
  Sig_Ranked_TF_target_df[,TF]<- sapply(Sig_Ranked_TF_target_df$PN_Genes, function(x)
    ifelse (x %in% TF_target_genes,1,0))
}

###Delete PN Genes column
Sig_Ranked_TF_target_df <- 
  subset(Sig_Ranked_TF_target_df, select = -c(PN_Genes) )

Sig_Ranked_TF_target_df <- t(Sig_Ranked_TF_target_df)
Sig_Ranked_TF_target_df <-as.data.frame(Sig_Ranked_TF_target_df)
##Convert rowname to column 
Sig_Ranked_TF_target_df$TF <- rownames(Sig_Ranked_TF_target_df)

#merge with ranked
Sig_Ranked_TF_target_df_rank<- merge(TF_Ranks,Sig_Ranked_TF_target_df, by = 'TF', all = TRUE)
  
####
write.xlsx( Sig_Ranked_TF_target_df_rank, file = "Chipseq_PN_TF_Targets_Ranked.xlsx")


 



# ________________________________________
#Create Colour function to distribute colours more evenly ####

heatcolS <- circlize::colorRamp2 (breaks = c(0,1),
                                  c("white","blue"),
                                  transparency = .3)
heatcolS(seq(0, 1, 1))




## Heatmap with ward clustering of samples and genes

pdf (paste0("Encode_TF_target_matrix",".pdf"),width=27.5, height=15)
print(Heatmap (TF_heatmap_matrix, width = unit(650, "mm"), height = unit(300, "mm"),
               # name ='',
               rect_gp = gpar(col= FALSE),
               
               column_title_gp = gpar(fontsize = 25, fontface = "bold"),
               column_names_centered = TRUE,
               column_names_rot = 90,
               clustering_method_rows = "ward.D2",
               clustering_method_columns = "ward.D2",
               column_names_gp = gpar(fontsize = 12, fontface = "bold"),
               row_names_gp = gpar(fontsize = 12, fontface = "bold"),
               column_dend_height = unit(10,'mm'),
               row_dend_width=unit(10,"mm"),
               #column_title = "Expression of Proteostasis Network Genes in Primary Skin Cutaneous Melanoma Samples",
              # column_title = paste('Significant PN Gene Transcription Factor Targets (ENCODE Chip Seq)'),
               show_column_names = TRUE,
               show_row_names = TRUE,
              # column_split = 2,
               #  row_split = 1,
               #   row_gap = unit(3, "mm"),
               column_gap = unit(3, "mm"),
               col =  heatcolS,
              
               
               # cluster_column_slices = FALSE,
               # cluster_row_slices = FALSE,
               
               heatmap_legend_param = list(
                 show_legend = c( FALSE,FALSE,FALSE))
               
))
#ggsave(paste0("TF_Gene_Heatmaps/TF_heatmap_cancer_SKCM_gene_and_sample_clustered_",counter,".pdf", HM))
dev.off()


### Pheatmap matrix
install.packages(pheatmap)

# load package
library(pheatmap)
pdf (paste0("Encode_TF_target_matrix_pheatmap",".pdf"),width=27.5, height=15)
pheatmap(TF_heatmap_matrix, 
         color = colorRampPalette(c("white","lightslateblue"))(2),
         angle_col = "90",
        # legend =FALSE
        )
dev.off()
