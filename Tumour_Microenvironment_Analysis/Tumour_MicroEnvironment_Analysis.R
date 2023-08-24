options(repos = c(CRAN = 'https://cloud.r-project.org'))
##THis script will analyse the distribution of cells in the Tumour microenvironment based on 
## Gene expression data from TCGA and creates a table showing levels of significance of difference of prevalence
install.packages("devtools")
install.packages("readxl")

install.packages("ggpubr")
devtools::install_github("cansysbio/ConsensusTME")
BiocManager::install("edgeR")
options(install.packages.compile.from.source="interactive")

BiocManager::install("singscore")

library("ggpubr")

library("readxl")
library('xlsx')
library(dplyr) 
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

setwd("~/iCloud_Drive_Archive/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()


#load Rdata ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#### Read in expression ####
FKPM_data<- loadRData( file =  'SKCM_FKPM.Rdata')


#remove duplicates
FKPM_data<- FKPM_data[!duplicated(FKPM_data$hgnc_symbol), ]
#### Make gene column rowname
FKPM_data<- data.frame(FKPM_data, row.names = 1)
#remove first column
FKPM_data_exp<- FKPM_data[,-1]

FKPM_data_exp_t<- as.data.frame(t(FKPM_data_exp))

##create sample column 

FKPM_data_exp_t$Sample<- rownames(FKPM_data_exp_t)

FKPM_data_exp_t$Sample<- substr(FKPM_data_exp_t$Sample, start=1,stop=15)

## remove columns with NA

FKPM_data_exp_t_no_NA<-FKPM_data_exp_t[ , colSums(is.na(FKPM_data_exp_t)) ==0]





###Read in final group allocation

final_sample_groups <- read.xlsx('Final_primary_and_metastatic_sample_groups.xlsx', sheetIndex = 1)

# ##merge with sample group allocation
# 
# sample_group_FKPM_data_exp <- merge(final_sample_groups,FKPM_data_exp,by.x = 'Sample.' , by.y ='Sample')
# 


### Create a matrix for sample group Primary Low ####

Low_prim_samples <- subset(final_sample_groups$Sample., 
                           final_sample_groups$Sample_Group_F %in% "Low PN Prim" )

LP_Matrix<- (subset(FKPM_data_exp_t, 
                    FKPM_data_exp_t$Sample %in% Low_prim_samples ))

##remove Sample column
LP_Matrix<- LP_Matrix[,-(length(LP_Matrix))]


# ### transpose so genes are rows and samples are columns and convert to a matrix
# 
LP_Matrix<- as.matrix(LP_Matrix)
LP_Matrix<- t(LP_Matrix)
# Estimate proportion of cell types in Tumour Microenvironment
Cell_distribution_A <-as.data.frame(ConsensusTME::consensusTMEAnalysis(LP_Matrix, 
                                   cancer = "SKCM", statMethod = "ssgsea"))
### Transpose so samples are rows and cell types are columns

Cell_distribution_A_t <- as.data.frame(t(Cell_distribution_A))

###Add Column with Sample Group9
Cell_distribution_A_t$Sample_Group <- "Low PN Prim"
###Move Sample_Group column to beginning
Cell_distribution_A_t <- Cell_distribution_A_t[moveme(names(Cell_distribution_A_t), "Sample_Group first")]



### Create a matrix for sample group Primary High ####

High_prim_samples <- subset(final_sample_groups$Sample., 
                           final_sample_groups$Sample_Group_F %in% "High PN Prim" )

HP_Matrix<- (subset(FKPM_data_exp_t, 
                    FKPM_data_exp_t$Sample %in% High_prim_samples ))

##remove Sample column
HP_Matrix<- HP_Matrix[,-(length(HP_Matrix))]


# ### transpose so genes are rows and samples are columns and convert to a matrix
# 
HP_Matrix<- as.matrix(HP_Matrix)
HP_Matrix<- t(HP_Matrix)
# Estimate proportion of cell types in Tumour Microenvironment
Cell_distribution_B <-as.data.frame(ConsensusTME::consensusTMEAnalysis(HP_Matrix, 
                                                                       cancer = "SKCM", statMethod = "ssgsea"))
### Transpose so samples are rows and cell types are columns

Cell_distribution_B_t <- as.data.frame(t(Cell_distribution_B))

###Add Column with Sample Group9
Cell_distribution_B_t$Sample_Group <- "High PN Prim"
###Move Sample_Group column to beginning
Cell_distribution_B_t <- Cell_distribution_B_t[moveme(names(Cell_distribution_B_t), "Sample_Group first")]


### Create a matrix for sample group Metastatic Low ####

Low_Met_samples <- subset(final_sample_groups$Sample., 
                           final_sample_groups$Sample_Group_F %in% "Low PN Met" )

LP_Matrix<- (subset(FKPM_data_exp_t, 
                    FKPM_data_exp_t$Sample %in% Low_Met_samples ))

##remove Sample column
LP_Matrix<- LP_Matrix[,-(length(LP_Matrix))]


# ### transpose so genes are rows and samples are columns and convert to a matrix
# 
LP_Matrix<- as.matrix(LP_Matrix)
LP_Matrix<- t(LP_Matrix)
# Estimate proportion of cell types in Tumour Microenvironment
Cell_distribution_C <-as.data.frame(ConsensusTME::consensusTMEAnalysis(LP_Matrix, 
                                                                       cancer = "SKCM", statMethod = "ssgsea"))
### Transpose so samples are rows and cell types are columns

Cell_distribution_C_t <- as.data.frame(t(Cell_distribution_C))

###Add Column with Sample Group9
Cell_distribution_C_t$Sample_Group <- "Low PN Met"
###Move Sample_Group column to beginning
Cell_distribution_C_t <- Cell_distribution_C_t[moveme(names(Cell_distribution_C_t), "Sample_Group first")]

### Create a matrix for sample group Metastatic HIgh ####

High_Met_samples <- subset(final_sample_groups$Sample., 
                          final_sample_groups$Sample_Group_F %in% "High PN Met" )

LP_Matrix<- (subset(FKPM_data_exp_t, 
                    FKPM_data_exp_t$Sample %in% High_Met_samples ))

##remove Sample column
LP_Matrix<- LP_Matrix[,-(length(LP_Matrix))]


# ### transpose so genes are rows and samples are columns and convert to a matrix
# 
LP_Matrix<- as.matrix(LP_Matrix)
LP_Matrix<- t(LP_Matrix)
# Estimate proportion of cell types in Tumour Microenvironment
Cell_distribution_D <-as.data.frame(ConsensusTME::consensusTMEAnalysis(LP_Matrix, 
                                                                       cancer = "SKCM", statMethod = "ssgsea"))
### Transpose so samples are rows and cell types are columns

Cell_distribution_D_t <- as.data.frame(t(Cell_distribution_D))

###Add Column with Sample Group9
Cell_distribution_D_t$Sample_Group <- "High PN Met"
###Move Sample_Group column to beginning
Cell_distribution_D_t <- Cell_distribution_D_t[moveme(names(Cell_distribution_D_t), "Sample_Group first")]




####Amalgamate all Matrix into one

Cell_distribution_matrix<- rbind(Cell_distribution_A_t,Cell_distribution_B_t,Cell_distribution_C_t, Cell_distribution_D_t)

Cell_distribution_matrix$Sample_Group<- factor(Cell_distribution_matrix$Sample_Group,
                                               levels = c("Low PN Prim" , "High PN Prim" ,"Low PN Met" ,  "High PN Met" ))


Cell_distribution_matrix_prim <- subset(Cell_distribution_matrix, Cell_distribution_matrix$Sample_Group %in% c("Low PN Prim" , "High PN Prim"))

Cell_distribution_matrix_met <- subset(Cell_distribution_matrix, Cell_distribution_matrix$Sample_Group %in% c("Low PN Met" , "High PN Met"))



### save Cell_distributino_matrix
write.xlsx(Cell_distribution_matrix, file = "Immune_Cell_Distribution_by_corrected_Sample_Group_Nov.xlsx")

save(Cell_distribution_matrix, file = "Immune_Cell_Distribution_by__corrected_Sample_Group_Nov.RData")


write.xlsx(Cell_distribution_matrix_prim, file = "Immune_Cell_Distribution_by_corrected_Sample_Group_Nov_prim.xlsx")

save(Cell_distribution_matrix_prim, file = "Immune_Cell_Distribution_by__corrected_Sample_Group_Nov_prim.RData")

write.xlsx(Cell_distribution_matrix_met, file = "Immune_Cell_Distribution_by_corrected_Sample_Group_Nov_met.xlsx")

save(Cell_distribution_matrix_met, file = "Immune_Cell_Distribution_by__corrected_Sample_Group_Nov_met.RData")

Cell_distribution_matrix_met$Sample_Group<-factor(Cell_distribution_matrix_met$Sample_Group,
                                                  levels=c("Low PN Met" , "High PN Met"))

### Prepare box plots

cell_types <-  names(Cell_distribution_matrix_met)

cell_types<- cell_types[-c(1)]

##_________________________
for (c in 1:length(cell_types)){

type<- cell_types[c]
#i <- which(cell_types == type)
 
j <- ggplot(Cell_distribution_matrix_met, aes(x=Sample_Group, y = Cell_distribution_matrix_met[,(c+1)])) + 
 #scale_x_discrete(labels=c('chap_mean'='Chaperones' )) +
  geom_boxplot()+
  #scale_fill_manual(values=c("darkred","black" )) +
  stat_compare_means(aes(group = Sample_Group), label = "p.format", size = 10,
                    label.y = -0.1 ) + #+
 geom_jitter()+
  
 ylab(type)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # legend.title = element_text(size = 20),
        # legend.text = element_text(size = 20) ,
        # # panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"),
        panel.background = element_rect(fill = 'white', colour = 'black'))+
 labs( x = "", y = gsub("\\_", " ", type ))

 pdf(file = paste0(type, "_boxplot_met-.pdf"), width = 6, height = 5)
 print(j)
dev.off()}

##_________________________prim ####
for (c in 1:length(cell_types)){
  
  type<- cell_types[c]
  #i <- which(cell_types == type)
  
  j <- ggplot(Cell_distribution_matrix_prim, aes(x=Sample_Group, y = Cell_distribution_matrix_prim[,(c+1)])) + 
    #scale_x_discrete(labels=c('chap_mean'='Chaperones' )) +
    geom_boxplot()+
    #scale_fill_manual(values=c("darkred","black" )) +
    stat_compare_means(aes(group = Sample_Group), label = "p.format", size = 10,
                       label.y = 0.1 ) + #+
    geom_jitter()+
    
    ylab(type)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          # legend.title = element_text(size = 20),
          # legend.text = element_text(size = 20) ,
          # # panel.background = element_rect(fill = "white"),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25,face="bold"),
          panel.background = element_rect(fill = 'white', colour = 'black'))+
    labs( x = "")
  pdf(file = paste0(type, "_boxplot_prim.pdf"), width = 6, height = 5)
  print(j)
  dev.off()}



####import cell type significance table

cell_type_table <- read.xlsx('cell_type_table_nov.xlsx', sheetIndex = 1)
cell_type_table$Primary_asterisk <- sapply(cell_type_table$Primary , function(x)
  ifelse(x > 0.05,"ns" ,
         ifelse(x <= 0.05 & x > 0.01 ,"*",
                ifelse(x <= 0.01 & x > 0.001 ,"**", 
                       ifelse(x <= 0.001 & x > 0.0001 ,"***", "****"
                )))))
  
cell_type_table$Metastatic_asterisk <- sapply(cell_type_table$Metastatic , function(x)
  ifelse(x > 0.05,"ns" ,
         ifelse(x <= 0.05 & x > 0.01 ,"*",
                ifelse(x <= 0.01 & x > 0.001 ,"**", 
                       ifelse(x <= 0.001 & x > 0.0001 ,"***", "****"
                       )))))        

### save Cell_distributino_matrix
write.xlsx(cell_type_table, file = "cell_type_table_asterisk.xlsx")

save(cell_type_table, file = "cell_type_table_asterisk.RData")

###Generate Genesets
rawMethodSignatures <- ConsensusTME::methodSignatures

matchedSigs <- ConsensusTME::matchGeneSigs(rawMethodSignatures)

ConsensusTME::buildConsensusGenes(matchedSigs)


##Read PN Genes

PN_Genes<- read.xlsx('PN_Genes.xlsx', sheetIndex = 1)
PN_Genes <- PN_Genes$Genes

#PN_Endothelial <- intersect(PN_Genes, Endothelial_Cells)

###Read in Unfiltered Consensus Genes

Consensus_Genes<- read.xlsx('Unfiltered_Consensus_Genes.xlsx', sheetIndex = 1)

### Create matrix of Consensus genes that are PN genes

x <- colnames(Consensus_Genes)
y <- PN_Genes
#subset(PN_cluster_list$Gene.Symbol, PN_cluster_list$Prim_Gene_Cluster %in% c('PN_1(P)','PN_2(P)'))
Consensus_Genes_PN_Genes_df<- 
  as.data.frame( matrix(ncol = length(x), nrow=length(y), dimnames = list(y,x)) )
Consensus_Genes_PN_Genes_df$PN_Genes <- rownames(Consensus_Genes_PN_Genes_df)
Consensus_Genes_PN_Genes_df <- Consensus_Genes_PN_Genes_df[moveme(names(Consensus_Genes_PN_Genes_df), "PN_Genes first")]

#### create loop to write 'X'  if PN Gene is in the list of targets for the TF colname

## Read in csv of targets of C1C2TFs 

#C1C2_TF_Targets <- read.csv ('C1C2_TF_Crosstarget.csv')
t=0
for(t in 0:((length(x))-1))
  
{
  t <-t+1
  print(t)
  cell_type <- x[t]
  cell_type_target_genes <- Consensus_Genes [,cell_type]
  Consensus_Genes_PN_Genes_df[,cell_type]<- sapply(Consensus_Genes_PN_Genes_df$PN_Genes, function(x)
    ifelse (x %in% cell_type_target_genes,1,0))
}

### save Cell_distributino_matrix
write.xlsx(Consensus_Genes_PN_Genes_df, file = "Consensus_Genes_PN_Genes_df.xlsx")

save(Consensus_Genes_PN_Genes_df, file = "Consensus_Genes_PN_Genes_df.RData")
