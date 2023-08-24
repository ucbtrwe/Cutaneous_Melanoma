#https://github.com/hamidghaedi/RNA-seq-differential-expression
#cite Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
# sequence count data: removing the noise and preserving large differences.
# Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
#________Packgaes______________#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("IHW")
# BiocManager::install("apeglm")
# BiocManager::install("PCAtools")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
# 
# BiocManager::install("ExperimentHub")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
# BiocManager::install("TCGAbiolinks") # bioconductor package
# BiocManager::install("SummarizedExperiment") # bioconductor package
# BiocManager::install("DESeq2") # bioconductor package
# BiocManager::install("IHW") # bioconductor package
# BiocManager::install("biomaRt") # bioconductor package
# BiocManager::install("apeglm") # bioconductor package
# install.packages("pheatmap") # CRAN package
# install.packages("RColorBrewer") # CRAN package
# BiocManager::install("PCAtools") # bioconductor package
# install.packages('reshape2') # CRAN package
# BiocManager::install("TCGAbiolinks")
# install.packages('xlsx') # CRAN package
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("vsn")

#install.packages('survminer')
install.packages('survival')
library("survminer")
library("TCGAbiolinks") # bioconductor package
library("SummarizedExperiment") # bioconductor package
library("DESeq2") # bioconductor package
library("IHW") # bioconductor package
library("biomaRt") # bioconductor package
library("apeglm") # bioconductor package
library("pheatmap") # CRAN package
library("RColorBrewer") # CRAN package
library("PCAtools") # bioconductor package
library(reshape2) # CRAN package
library('xlsx')
library('readxl')
library('TCGAbiolinks')
library('DESeq2')
library('ggplot2')
library('vsn')
library('pheatmap')
library('RColorBrewer')
library('apeglm')
library('ggrepel')
library('EnsDb.Hsapiens.v79')
library ("SummarizedExperiment")
library ("vsn")
library ("ComplexHeatmap")
library ("survival")
#load Rdata ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}
  
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinks") 

setwd("~/iCloud_Drive_Archive/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()


# # #_______ Downloading_Data_______#
query_TCGA = GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  legacy = FALSE)

GDCdownload(query_TCGA)
rna <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
# exp matrix
rna <- as.data.frame(SummarizedExperiment::assay(dat))
# ###Save RNA
save(rna, file = 'TCGA_SKCM_raw_RNA_seq-reads_Dan.Rdata')
# write.csv(rna, file = 'TCGA_SKCM_raw_RNA_seq-reads_Dan.csv')

# ## read in raw RNA reads
# rna <-  read.csv('TCGA_SKCM_raw_RNA_seq-reads_Dan.csv')

## select only columns that are metastatic samples
rna<- loadRData('TCGA_SKCM_raw_RNA_seq-reads_Dan.Rdata')  
rna_met <- rna[ , grep("06A", colnames(rna))]  # Find matches
rna<- rna_met    
names(rna) <- substr(names(rna),start=1,stop=15)
## sort columns alphabetically
#sort data frame by name columns alphabetically

rna<- rna[,order(names(rna)) ]


# #### Read in list of PN Genes
# # #### Read in list of PN Genes
# PN_Genes_df <- read_xlsx('PN_Genes_FInal_LIst_428.xlsx', sheet = 1)
# PN_Genes <- PN_Genes_df$Gene
# 
# #### Select PN Genes from TCGA data
# rna <- subset(rna,
#                        rowData(rna)$gene_name %in% PN_Genes)
# ##save metastatic PN gene raw expression
# save(rna, file = 'TCGA_metastatic_SKCM_raw_PN_Gene_RNA_seq-reads_Dan.Rdata')
# 
# dim(rna)

### read in table of sample allocations
sample_allocation <- read.xlsx('Final_Corrected_Expression_primary_and_metastatic_sample_groups.xlsx', sheetIndex = 1)
sample_allocation$Sample_A <-paste0(sample_allocation$Sample,'A')
###create llist of samples that are in A_metastatic

A_metastatic <- sample_allocation$Sample_A[sample_allocation$Sample_Group_Name == "Low PN Met"]
B_metastatic <- sample_allocation$Sample_A[sample_allocation$Sample_Group_Name == "High PN Met"]


rna$Sample_Group <-ifelse(rna$sample %in% (A_metastatic),"Low PN Met","High PN Met"  )

rna$Sample_Group <-factor(rna$Sample_Group, levels = c("High PN Met", "Low PN Met"))

#_______Making_Expression_Object__________#

dds <- DESeqDataSet(rna, design = ~ Sample_Group)



####Check if any genes have 0 counts in all columns
table(rowSums(counts(dds)) > 0)

dds <- estimateSizeFactors(dds)

summary(dds$sizeFactor)

counts_normalized <- counts(dds, normalized = TRUE)

###save dds

save(dds, file = 'TCGA_metastatic_SKCM_raw_Gene_RNA_seq-reads_Dan_dds.Rdata')


# Compare normalized and non-normalized counts data
summary(colSums(counts(dds)))

summary(colSums(counts(dds, normalized = TRUE))) # == colSums(counts_normalized)


par(mfrow=c(1,2))
hist(colSums(counts(dds)), breaks = 20, xlab = 'Total Raw Counts', main = 'Raw Counts')
hist(colSums(counts_normalized), breaks = 20, xlab = 'Total Normalized Counts', main = 'Normalized Counts')

counts_log_normalized <- log2(counts_normalized + 1)

# Normalized counts: variance increases with mean.
SdPlot <- meanSdPlot(counts_normalized, ranks = FALSE, plot = FALSE)
SdPlot$gg + ggtitle('sequencing depth normalization') + ylab('standard deviation')


SdPlot_log <- meanSdPlot(counts_log_normalized, ranks = FALSE, plot = FALSE)
SdPlot_log$gg + ggtitle('sequencing depth normalized log2(read counts)') + ylab('standard deviation')

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# We see that the variance in expression is now essentially independent of mean expression, as displayed by a virtually horizontal fit line.
SdPlot_vsd <- meanSdPlot(assay(vsd), ranks = FALSE, plot = FALSE)
SdPlot_vsd$gg + ggtitle('variance stabilizing transformation') + ylab('standard deviation')
# 
par(mfrow = c(2,2))
# Raw data normalized by sequencing depth
plot(counts_normalized[,1:2], cex=.1, main = 'Normalized by sequencing depth')
# log2(x+1) transformation
plot(counts_log_normalized[,1:2], cex=.1, main = 'Normalized log2(read counts)')
# variance stabilizing transformation
plot(assay(vsd)[,1:2], cex=.1, main = 'Variance stabilizing transformation')


# Create the distance matrix
sampleDists_vsd <- dist(t(assay(vsd)))
sampleDistsMatrix_vsd <- as.matrix(sampleDists_vsd)
# For easier identification, let the rownames be a combination of patient and Sample_Type_Code
rownames(sampleDistsMatrix_vsd) <- paste(vsd$Sample_Type_Code, vsd$patient, sep = ' - ')
colnames(sampleDistsMatrix_vsd) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,'Blues')))(255)
par(mfrow=c(1,1))
heatmap <- pheatmap(sampleDistsMatrix_vsd,
                    clustering_distance_rows = sampleDists_vsd,
                    clustering_distance_cols = sampleDists_vsd,
                    col = colors)

# Plot PCA on Sample_Type
plotPCA(vsd, intgroup = 'Sample_Group')

# Run DeSeq2 DGE pipeline
dds_DGE <- DESeq(dds)

plotDispEsts(dds_DGE)

dds_DGE_results <- results(dds_DGE)
# write.table(dds_DGE_results,file='~/Documents/Rotation3/Scripts/DGE_Analysis/LUAD_full/DESeq2_results_LUAD.txt',
#             col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

dds_DGE_results <- dds_DGE_results[order(dds_DGE_results$padj),] # Re-orders results by significance
head(dds_DGE_results)



###Edit to add column to say if gene has lower or higher expresssion in Low PN Met
dds_DGE_results$Expression<-  ifelse(dds_DGE_results$padj > 0.1, "No Significant Difference",
                                     ifelse(dds_DGE_results$log2FoldChange <= -0.321928095, "Lower",
                                            ifelse(dds_DGE_results$log2FoldChange >= 0.321928095, "Higher","No Significant Difference")))

### remove "." from Ensemble Gene Name row name

rownames(dds_DGE_results) <-gsub(".", "", rownames(dds_DGE_results), fixed = TRUE)

###save dds

save(dds_DGE_results, file = 'TCGA_metastatic_SKCM_raw_Gene_RNA_seq-reads_Dan_dds_DGE_results.Rdata')
write.xlsx(dds_DGE_results, file = 'TCGA_metastatic_SKCM_raw_Gene_RNA_seq-reads_Dan_dds_DGE_results.xlsx')


####Merge eds results with gene name
## read in csv of Gene Decoder with ensmble and Hugo names
Gene_decoder <- read.csv("Gene_Decoder.csv")

dds_DGE_results_with_gene_name <- merge(Gene_decoder,as.data.frame(dds_DGE_results) , by.x = 'Ensembl', by.y = 0)

###save dds

save(dds_DGE_results_with_gene_name, file = 'TCGA_metastatic_SKCM_raw_Gene_RNA_seq-reads_Dan_dds_DGE_results_with_gene_name.Rdata')
write.csv(dds_DGE_results_with_gene_name, file = 'TCGA_metastatic_SKCM_raw_Gene_RNA_seq-reads_Dan_dds_DGE_results_with_gene_name.csv')





# MA plot
plotMA(dds_DGE_results, main = 'A (metastatic) vs B (metastatic)', ylim = c(-5,5))

##Volcano Plots
results_order <- as.data.frame(dds_DGE_results) # dds_DGE_results is a DESeqResults object. Need a simple data frame
results_order$label <- ifelse(results_order$padj < 0.1,
                            ifelse(abs(results_order$log2FoldChange) > 1, "log2(FC)>1","FDR<0.1"),"NotSig")
p <- ggplot(results_order, aes(x = log2FoldChange, y = -log10(padj))) +
 geom_point(aes(col = label)) +
  scale_color_manual(values = c('red','blue','black')) +
  ggtitle('PN Gene Expression A (metastatic) vs B (metastatic)')

# Add gene symbol names to most significant genes with log2FoldChange > 2
DEgenes_DESeq <- results_order[abs(results_order$log2FoldChange) > 1 & results_order$padj < 0.1 & !is.na(results_order$padj),]
DEgenes_DESeq <- DEgenes_DESeq[order(DEgenes_DESeq$padj),]

DEgenes_DESeq <- merge(Gene_decoder,DEgenes_DESeq , by.x = 'Ensembl', by.y = 0)

# DEgenes_DESeq$GeneSymbol <- ensembldb::select(EnsDb.Hsapiens.v79,
#                                               keys = rownames(DEgenes_DESeq),
#                                               keytype = 'GENENAME',
#                                               columns = 'SYMBOL')

#$SYMBOL
 q <- p + geom_text_repel(data=DEgenes_DESeq[1:100,] , aes(label = DEgenes_DESeq$hgnc_symbol[1:100]) )

 pdf ("PN_Genes_DSeq_Volcano_Plot.pdf", width=17.6, height=9)
 print(q)
 dev.off()
 
 #Heatmap
 topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE),428)
 mat <- assay(vsd)[topVarGenes,]
 mat <- mat - rowMeans(mat) # Normalize using mean expression of given gene
 anno <- as.data.frame(colData(vsd)[,'Sample_Group']) ; colnames(anno) = 'Sample_Group' ;
 rownames(anno) = colnames(mat)
 t <- pheatmap(mat, annotation_col = anno, show_colnames = FALSE, show_rownames = FALSE)
 pdf ("PN_Genes_DSeq_Heatmap.pdf", width=17.6, height=9)
 print(t)
 dev.off()
 
 #Heatmap top 150 variable genes
 topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 150)
 mat <- assay(vsd)[topVarGenes,]
 mat <- mat - rowMeans(mat) # Normalize using mean expression of given gene
 anno <- as.data.frame(colData(vsd)[,'Sample_Group']) ; colnames(anno) = 'Sample_Group' ;
 rownames(anno) = colnames(mat)
 t <- pheatmap(mat, annotation_col = anno, show_colnames = FALSE, show_rownames = FALSE)
 pdf ("PN_Genes_DSeq_top_150_DEG_Heatmap.pdf", width=17.6, height=9)
 print(t)
 dev.off()
 
 colnames(mat) <- str_sub(colnames(mat), 1, -3)
 rownames(df) <- colnames(mat)
 
 #Heatmap top 120 variable genes
 topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 120)
 mat <- assay(vsd)[topVarGenes,]
 mat <- mat - rowMeans(mat) # Normalize using mean expression of given gene
 anno <- as.data.frame(colData(vsd)[,'Sample_Group']) ; colnames(anno) = 'Sample_Group' ;
 rownames(anno) = colnames(mat)
 t <- pheatmap(mat, annotation_col = anno, show_colnames = FALSE, show_rownames = FALSE)
 pdf ("PN_Genes_DSeq_top_120_DEG_Heatmap.pdf", width=17.6, height=9)
 print(t)
 dev.off()
 
 
 
 
  
 
 
 