library("readxl")

library('xlsx')
library(ComplexHeatmap)
####  https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping


##### SET WORKING DIRECTORY , use getwd to designate new variable'current_dir' to use in sub directory names
setwd("~/iCloud_Drive_Archive/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()


#### Re-Create heatmap coloured by new clusters ####
sample_group_regressed_expr <- read.xlsx('SKCM_expression_regressed_by_purity_with_final_sample_groups.xlsx', sheetIndex = 1)

###Read in sample group allocations
sample_groups<- read.xlsx('Final_Corrected_Expression_primary_and_metastatic_sample_groups.xlsx',sheetIndex = 1)

#add patient barcode column to make merging with clinical tables easier
sample_groups$Patient_Barcode <- substr(sample_groups$Sample, start = 1, stop = 12 )
#save revised sample groups
write.xlsx(sample_groups, file ='Final_Corrected_Expression_primary_and_metastatic_sample_groups.xlsx')

##merge sample groups with expression
sample_group_regressed_expr <-  merge(sample_groups, sample_group_regressed_expr , by = "Sample"  )

#save expression and sample groups table
write.xlsx(sample_group_regressed_expr, file ='Final_Expression_primary_and_metastatic_samples_with_groups.xlsx')

###load table with no dups
sample_group_regressed_expr<- read.xlsx (file ='Final_Expression_primary_and_metastatic_samples_with_groups_no_dups.xlsx', sheetIndex = 1)
## select metastatic samples

Metastatic_regressed <- subset(sample_group_regressed_expr, Sample_Group_Split %like% c("Met"))

#### Heatmap
#Create table of group allocations
new_sample_clusters <- Metastatic_regressed[,c(2,4,5)]

# create table of expression

Metastatic_regressed_exp <- data.frame(Metastatic_regressed, row.names = 1)
# select columms with expression data and transpose dataframe so genes are rows and samples are columns
Metastatic_regressed_exp <- Metastatic_regressed_exp [,-c(1:7)]
Metastatic_regressed_exp <- as.data.frame(t(Metastatic_regressed_exp))
# Create new matrix of Z scores

Reg_Exp_Z_Score <- (apply(Metastatic_regressed_exp, 1, function(x) (x - mean(x)) / sd(x)))

# convert to Dataframe

Reg_exp_df <- as.data.frame(Reg_Exp_Z_Score)


######### HEATMAPS ######### to allocate samples to groups 

#create matrix for heatmap (gene expression columns transposed so genes are rows and samples are columns )
heatmap_matrix <- as.matrix(t(Reg_exp_df))
### Remove rows with NA
heatmap_matrix <-    heatmap_matrix[complete.cases( heatmap_matrix), ]
# ________________________________________
#Create Colour function to distribute colours more evenly ####

heatcolS <- circlize::colorRamp2 (breaks = c(-7,-0.1,0,0.1, 2.5,5, 20),
                                  c("dodgerblue4","#f1fbfe" ,"white", "#ffff55", "#ff4015", "#ff0000", '#5e0406'),
                                  transparency = .3)
heatcolS(seq(-7, 20, 0.01))


 


######### HEATMAPS ######### 

#create matrix for heatmap (gene expression columns transposed so genes are rows and samples are columns )
heatmap_matrix <- as.matrix(t(Reg_exp_df))
### Remove rows with NA
heatmap_matrix <-    heatmap_matrix[complete.cases( heatmap_matrix), ]
# ________________________________________
#Create Colour function to distribute colours more evenly ####

heatcolS <- circlize::colorRamp2 (breaks = c(-7,-0.1,0,0.1, 2.5,5, 20),
                                  c("dodgerblue4","#f1fbfe" ,"white", "#ffff55", "#ff4015", "#ff0000", '#5e0406'),
                                  transparency = .3)
heatcolS(seq(-7, 20, 0.01))

# Create top annotation

ann2 <- data.frame (new_sample_clusters$Sample_Group_split)
colnames(ann2) <- c('Sample_Group')

ann2$Sample_Group <- factor(ann2$Sample_Group,
                            levels = c( "Low PN Met (A)","Low PN Met (B)" , "High PN Met (A)", "High PN Met (B)"))

colours2 <- list('Sample_Group' = c("Low PN Met (A)" = '#bdeaea', 'Low PN Met (B)' ="#5acbca" , "High PN Met (A)" = "#3e1453" ,"High PN Met (B)" = "#a47ab9")
)

colAnn2 <- HeatmapAnnotation(df = ann2,
                             which = 'col',
                             col = colours2,
                             annotation_width = unit(c(2, 10), 'cm'),
                             show_annotation_name = FALSE,
                             annotation_legend_param = list(
                               Sample_Group = list(
                                 title = "Sample Group",
                                 #  title_position = "lefttop-rot",
                                 annotation_legend_side="right",
                                 annotation_legend_height = unit(160, "mm"),  
                                 annotation_legend_width = unit(60, "mm"), 
                                 grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
                                 labels_gp = gpar(fontsize = 20),
                                 title_gp = gpar(fontsize = 20, fonPIk3ace = 'bold'))
                             )
)


## Heatmap with ward clustering of samples and genes

x<- (Heatmap (heatmap_matrix, width = unit(260, "mm"), height = unit(200, "mm"),
              # name ='',
              rect_gp = gpar(col= FALSE),
              
              column_title_gp = gpar(fontsize = 15, fontface = "bold"),
              # column_names_centered = TRUE,
              #column_names_rot = 90,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              # column_names_gp = gpar(fontsize = 25, fonPIk3ace = "bold"),
              row_names_gp = gpar(fontsize = 6),
              column_dend_height = unit(10,'mm'),
              row_dend_width=unit(10,"mm"),
              #column_title = "Expression of Proteostasis Network Genes in Metastatic Skin Cutaneous Melanoma Samples",
              column_title = paste('Gene expression - Regressed for Tumour Purity'),
              show_column_names = FALSE,
              show_row_names = FALSE,
              column_split = 4,
              row_split = 3,
              row_gap = unit(3, "mm"),
              column_gap = unit(3, "mm"),
              col =  heatcolS,
              top_annotation=colAnn2,
              
              cluster_column_slices = FALSE,
              cluster_row_slices = FALSE,
              use_raster = FALSE,
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
#ggsave(paste0("Reg_Gene_Heatmaps/Reg_heatmap_cancer_SKCM_gene_and_sample_clustered_",counter,".pdf", HM))

pdf (paste0("Regressed_expression_heatmap_cancer_SKCM_Metastatic_Nov_23_regressed_groups",".pdf"),width=17, height=10)
print(x)
dev.off()









####Extract gene clusters
## Extract Row order
row_order(x)

### Extract row cluster allocation from heatmap
x = draw(x)
r.dend <- row_dend(x)  #If needed, extract row dendrogram
rcl.list <- row_order(x)  #Extract clusters (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters


# loop to extract genes for each cluster.
for (i in 1:length(row_order(x))){
  if (i == 1) {
    clu <- t(t(row.names(heatmap_matrix[row_order(x)[[i]],])))
    out <- cbind(clu, paste("reg_cluster_", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(heatmap_matrix[row_order(x)[[i]],])))
    clu <- cbind(clu, paste("reg_cluster_", i, sep=""))
    out <- rbind(out, clu)
  }
}


#check
out

# ##### Save new gene clusters
save(out,  file = 'Reg_Gene_Sets.RData')
write.xlsx( out, file = "Reg_Gene_Sets.xlsx")

### load list of previous genesets
previous_genesets <- read.xlsx('SKCM_Gene_Sets_May.xlsx', sheetIndex = 1)

gene_sets_may<- merge(previous_genesets, out, by.x = 'Gene.Symbol', by.y = 'GeneID')

save(gene_sets_may,  file = 'SKCM_Gene_Sets_regressed.RData')
write.xlsx( gene_sets_may, file = "SKCM_Gene_Sets_regressed.xlsx")



## read in list of HSP_Types

HSP_Types <- read.xlsx('HSP_Types.xlsx', sheetIndex = 1)

gene_sets_HSP_Types <- merge(HSP_Types, gene_sets_may, by.x = 'Gene.Symbol', by.y = 'Gene.Symbol')


save(gene_sets_HSP_Types,  file = 'SKCM_Gene_Sets_HSP_Types.RData')
write.xlsx( gene_sets_HSP_Types, file = "SKCM_Gene_HSP_Types.xlsx")
