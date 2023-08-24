library(data.table)
library('xlsx')

library('dplyr')
#load Rdata ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##### SET WORKING DIRECTORY , use getwd to designate new variable'current_dir' to use in sub directory names
setwd("~/iCloud_Drive_Archive/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()


## Read in UCSC normalised TCGA and GTEX FPKM data
TCGA_GTEX_FPKM <-fread("TcgaTargetGtex_rsem_gene_fpkm.gz")

# ###Read in non-sun_exposed_data_GTEX
# Not_Sun_GTEX_FPKM <-fread("gene_reads_2017-06-05_v8_skin_not_sun_exposed_suprapubic (1).gct.gz")
#
# #### Save
# save(Not_Sun_GTEX_FPKM, file ='GTEX_Not_Sun_Exposed_Gene_Expression.Rdata')
# write.csv(Not_Sun_GTEX_FPKM, file ='GTEX_Not_Sun_Exposed_Gene_Expression.csv')
#
#
#
# ###Read in non-sun_exposed_data_GTEX
# Sun_GTEX_FPKM <-fread("gene_reads_2017-06-05_v8_skin_sun_exposed_lower_leg.gct.gz")
#
# #### Save
# save(Sun_GTEX_FPKM, file ='GTEX_Sun_Exposed_Gene_Expression.Rdata')
# write.csv(Sun_GTEX_FPKM, file ='GTEX_Not_Sun_Exposed_Gene_Expression.csv')
#
# #### create list of sun gtex samples
# Sun_Gtex_samples <- names(Sun_GTEX_FPKM)
# Sun_Gtex_samples <-Sun_Gtex_samples [-c(1:3)]
# Sun_Sample_df <- as.data.frame(Sun_Gtex_samples)
# Sun_Sample_df$Sample_Group <- 'Sun_GTEX'
# colnames(Sun_Sample_df)[1]<- "Sample"
# #### create list of sun gtex samples
# Not_Sun_Gtex_samples <- names(Not_Sun_GTEX_FPKM)
# Not_Sun_Gtex_samples <-Not_Sun_Gtex_samples [-c(1:3)]
# Not_Sun_Sample_df <- as.data.frame(Not_Sun_Gtex_samples)
# Not_Sun_Sample_df$Sample_Group <- 'No_Sun_GTEX'
# colnames(Not_Sun_Sample_df)[1]<- "Sample"
#
# GTEX_Samples<- rbind(Sun_Sample_df,Not_Sun_Sample_df)
#
# ###Read in TCGA sample group list
#
# TCGA_Samples <- read.xlsx('Final_primary_and_metastatic_sample_groups.xlsx', sheetIndex = 1)
#
# ###Merge table of GTEX and TCGA
#
# Sample_table <- rbind(GTEX_Samples, TCGA_Samples )
#
# #### Save
# save(Sample_table, file ='GTEX_TCGA_Skin_Sample_table.Rdata')
# write.csv(Sample_table, file ='GTEX_TCGA_Skin_Sample_table.csv')


Sample_table<-read.xlsx('GTEX_TCGA_Skin_Sample_table_edited.xlsx',sheetIndex = 1)
###Select Skin samples
skin_samples <-Sample_table$Sample
skin_samples <- c(skin_samples,"sample")
names.use <- names(TCGA_GTEX_FPKM)[(names(TCGA_GTEX_FPKM) %in% skin_samples)]


TCGA_GTEX_FPKM_Skin <-TCGA_GTEX_FPKM[, ..names.use]


Gene_Decoder<- read.xlsx('Gene_Decoder_ENSG_Hugo.xlsx', sheetIndex = 1)
##Read in PN genes
PN_Genes_expression <-read.xlsx('TCGA_primary_PN_DSeq_Expression_Difference_edited.xlsx', sheetIndex = 2)
## add in ensembl ID
PN_Genes_expression<- merge(PN_Genes_expression,Gene_Decoder, by.x = "Gene", by.y ="Gene")

###translate PN Genes into ENSG



TCGA_GTEX_FPKM_Skin_with_gene_name  <-merge (TCGA_GTEX_FPKM_Skin, Gene_Decoder, by.x = "sample", by.y = "Gene_ENSG")

TCGA_GTEX_FPKM_Skin_with_gene_name <- as.data.frame(TCGA_GTEX_FPKM_Skin_with_gene_name )

##remove duplicated gene names

TCGA_GTEX_FPKM_Skin_with_gene_name<- TCGA_GTEX_FPKM_Skin_with_gene_name %>% distinct(Gene, .keep_all = TRUE)




rownames(TCGA_GTEX_FPKM_Skin_with_gene_name)<- TCGA_GTEX_FPKM_Skin_with_gene_name$Gene
## ##remove none-sample columns
TCGA_GTEX_FPKM_Skin_with_gene_name<- TCGA_GTEX_FPKM_Skin_with_gene_name[, -c(1,960)]


# #### Save
# save(TCGA_GTEX_FPKM_Skin_with_gene_name, file ='GTEX_TCGA_Skin_Sample_Expression.Rdata')
# write.csv(TCGA_GTEX_FPKM_Skin_with_gene_name, file ='GTEX_TCGA_Skin_Sample_Expression.csv')

##Read in TCGA GTEX expression data
TCGA_GTEX_FPKM_Skin_with_gene_name<- read.csv (file ='GTEX_TCGA_Skin_Sample_Expression.csv')
### create subset of data for PN genes
TCGA_GTEX_FPKM_Skin_PN_Genes<- subset(TCGA_GTEX_FPKM_Skin_with_gene_name, rownames(TCGA_GTEX_FPKM_Skin_with_gene_name)%in%PN_Genes_expression$Gene) 

# #### Save
# save(TCGA_GTEX_FPKM_Skin_PN_Genes, file ='GTEX_TCGA_Skin_Sample_PN_Gene_Expression.Rdata')
# write.csv(TCGA_GTEX_FPKM_Skin_PN_Genes, file ='GTEX_TCGA_Skin_Sample_PN_Gene_Expression.csv')

# ###merge with sample groups
# TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups<- as.data.frame(t(TCGA_GTEX_FPKM_Skin_PN_Genes))
# 
# ### load skin PN table
# TCGA_GTEX_FPKM_Skin_PN_Genes <- loadRData(file ='GTEX_TCGA_Skin_Sample_PN_Gene_Expression.Rdata')
# 
# ###turn first column into rownames
# TCGA_GTEX_FPKM_Skin_PN_Genes <- data.frame(TCGA_GTEX_FPKM_Skin_PN_Genes, row.names = 1)   
# 
# ##merge wiht sample table
# 
# TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups <- as.data.frame(t(TCGA_GTEX_FPKM_Skin_PN_Genes))
# TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample<- rownames(TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups)
# ### replace . with -
# TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample<- gsub(".", "-", TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample, fixed = TRUE)
# 
# TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample<- factor(TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample_Group,
#                                                                 levels =c("Low PN Prim","High PN Prim", "Low PN Met", "High PN Met" , "Normal" )  )
# 
# 
# TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups<- merge(TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups,Sample_table, by = "Sample")
# 
# # save(TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups, file ='GTEX_TCGA_with_Sample_Groups_and_Expression.Rdata')
# # write.csv(TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups, file ='GTEX_TCGA_with_Sample_Groups_and_Expression.csv')
# 
# 
# TCGA_GTEX_Sample_Groups_subset<- TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups[,c(1,392)]

# save(TCGA_GTEX_Sample_Groups_subset, file ='GTEX_TCGA_Sample_Groups_Subset.Rdata')
# write.csv(TCGA_GTEX_Sample_Groups_subset, file ='GTEX_TCGA_Sample_Groups_Subset.csv')

TCGA_GTEX_Sample_Groups_subset<- read.csv(file ='GTEX_TCGA_Sample_Groups_Subset.csv')


TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups_prim<- subset (TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups, 
                                                               TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample_Group 
                                                               %in% c("Normal",        
                                                                      "Cancer Low PN",           "Cancer High PN"           ))
### PCA https://towardsdatascience.com/principal-component-analysis-pca-101-using-r-361f4c53a9ff
wdbc.pr <- prcomp(TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups_prim[-c(1,392)], center = TRUE, scale = TRUE)
summary(wdbc.pr)

plot <- screeplot(wdbc.pr, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(wdbc.pr$sdev^2 / sum(wdbc.pr$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)

plot(wdbc.pr$x[,1],wdbc.pr$x[,2], xlab="PC1 (26.7%)", ylab = "PC2 (16.6%)", main = "PC1 / PC2 - plot")

library("factoextra")
plot <- fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 19, 
             pointsize = 1, 
             fill.ind = TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups_prim$Sample_Group, 
             col.ind = "black", 
             habillage= TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups_prim$Sample_Group,
            # palette = "jco", 
             addEllipses = FALSE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Sample Group") +
  scale_color_manual(values=c('#fde725',"#296a8c","#00e18b"))
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))






pdf ("PCA_Normal_cancer.pdf",width=3, height=2)
print(plot)
dev.off()

##met and prim
### PCA https://towardsdatascience.com/principal-component-analysis-pca-101-using-r-361f4c53a9ff
wdbc.pr <- prcomp(
  #TCGA_GTEX_Sample_Groups_subset,
TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups[-c(1,2,393,394)],
 center = TRUE, scale = TRUE)
summary(wdbc.pr)

plot <- screeplot(wdbc.pr, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(wdbc.pr$sdev^2 / sum(wdbc.pr$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)

plot(wdbc.pr$x[,1],wdbc.pr$x[,2], xlab="PC1 (28.8%)", ylab = "PC2 (14.6%)", main = "PC1 / PC2 - plot")

library("factoextra")
plot <- fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 19, 
                     pointsize = 0.6, 
                     fill.ind = TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample_Group, 
                     col.ind = "black", 
                     habillage= TCGA_GTEX_FPKM_Skin_PN_Genes_with_Sample_Groups$Sample_Group,
                     # palette = "jco", 
                     addEllipses = FALSE,
                     label = "var",
                     col.var = "black",
                     repel = TRUE,
                     legend.title = "Sample Group") +
scale_color_manual(values=c('#7f4e97','#fde725',"#7bd5d5", "#296a8c", "#00e18b"))+
  scale_fill_discrete(labels=c('Metastatic (B)', 'Primary (B)', 'Metastatic (A)',"Primary (A)", "Normal"))
ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))






pdf ("PCA_Normal_cancer_met_prim_April.pdf",width=3, height=2)
print(plot)
dev.off()




# ###PCA#####
# 
# library(Rtsne)
# library(ggplot2)
# library(RColorBrewer)
# library(sva)
# library(matrixStats)
# library(data.table)
# library(factoextra)
# library(matrixStats)
# library(devtools)
# install_github("vqv/ggbiplot")
# 
# library(ggbiplot)
# 
# library("FactoMineR")
# library("factoextra")
# install.packages(c("FactoMineR", "factoextra"))
# 
# # ## create moveme function to move columns
# # ####  https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping
# #
# moveme <- function (invec, movecommand) {
#   movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]],
#                                  ",|\\s+"), function(x) x[x != ""])
#   movelist <- lapply(movecommand, function(x) {
#     Where <- x[which(x %in% c("before", "after", "first",
#                               "last")):length(x)]
#     ToMove <- setdiff(x, Where)
#     list(ToMove, Where)
#   })
#   myVec <- invec
#   for (i in seq_along(movelist)) {
#     temp <- setdiff(myVec, movelist[[i]][[1]])
#     A <- movelist[[i]][[2]][1]
#     if (A %in% c("before", "after")) {
#       ba <- movelist[[i]][[2]][2]
#       if (A == "before") {
#         after <- match(ba, temp) - 1
#       }
#       else if (A == "after") {
#         after <- match(ba, temp)
#       }
#     }
#     else if (A == "first") {
#       after <- 0
#     }
#     else if (A == "last") {
#       after <- length(myVec)
#     }
#     myVec <- append(temp, values = movelist[[i]][[1]], after = after)
#   }
#   myVec
# }
# 
# 
# 
# 
# 
# # # # Set Seed to ensure consistent results when experiment is repeated
# set.seed(13)
# 
# # # specify range of p values to try
# # p_number <- seq(from = 200, to = 1500, by = 100)
# #transpose table with just expression values so samples are columns and genes are rows
# exp_table <- t(TCGA_GTEX_FPKM_Skin_with_gene_name)
# 
# # carry out principal component analysis
# gene_exp_cancer_PCA <-PCA(exp_table, graph = FALSE)
# 
# 
# 
# fviz_pca_ind(gene_exp_cancer_PCA,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )
# 
# #plot PCA
# ggbiplot(gene_exp_cancer_PCA)
# 
# 
# 
# 
# ####Scale data
# # calculate standard deviation of each gene and add as column to 
# # copy expression table
# scale_exp_table <- exp_table
# # # transpose expression table so genes are rows and samples are columns
# # 
# # scale_exp_table <- t(scale_exp_table)
# #convert back to data frame
# 
# scale_exp_table <- as.data.frame(scale_exp_table)
# 
# #Calculate Standard Deviation for each gene
# 
# scale_exp_table$SD <-  apply(scale_exp_table,1,sd)
# 
# 
# # move SD column to beginning of table
# 
# scale_exp_table <- scale_exp_table[moveme(names(scale_exp_table), "SD first")]
# 
# 
# # Divide values in other columns by SD
# 
# require(magrittr)
# scale_exp_table[,-(1)] %<>% sapply(`/`, scale_exp_table[,1])
# 
# 
# # remove standard deviation column
# 
# 
# scaled_exp <- scale_exp_table[,-(1)]
# 
# #Carry Out PCA of scaled data
# scale_gene_exp_cancer_PCA <-PCA(scaled_exp, graph = FALSE, scale=TRUE)
# 
# 
# #plot PCA of scaled data
# ggbiplot(scale_gene_exp_cancer_PCA)
# 
# 
# pca.data <- data.frame(Sample=rownames(scale_gene_exp_cancer_PCA$x),
#                        X=scale_gene_exp_cancer_PCA$x[,1],
#                        Y=scale_gene_exp_cancer_PCA$x[,2])
# pca.data
# 
# ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
#   geom_text() +
#   xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
#   ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
#   theme_bw() +
#   ggtitle("My PCA Graph")
# 

# ####### Using ComBat to remove batch effect of studies from expression data ####
# 
# 
# # create dataframe of samples and study
# 
# batch <-as.data.frame(gene_exp_cancer[,c(3,8)])
# 
# # name column of batch 'study'
# colnames(batch) <- c("sample","study")
# 
# # #create dataframe of expression data without non-expression columns e.g study 
# # 
# # combat_ready_exp <-   gene_exp_cancer[,15:ncol(gene_exp_cancer)]
# 
# # # # transpose expression dataframe so cases are columns and genes are rows (automatically converts it to matrix)
# # combat_ready_exp <- t(combat_ready_exp)
# 
# 
# # convert scaled data dataframe to matrix
# 
# scaled_exp<- as.matrix(scaled_exp)
# 
# # check if any rows have variance of 0 (as this stops ComBat working)
# 
# length(which ((apply(scaled_exp, 1, var)!=0) == "TRUE"))
# 
# length(which ((apply(scaled_exp, 1, var)!=0) == "FALSE"))


# # # # transpose expression dataframe so genes are columns and samples are rows 
# combat_ready_exp <- t(combat_ready_exp)
# 
# # remove gene / column(s) for which variance is 0
# 
# zeroVar <- function(combat_ready_exp, useNA = 'ifany') { 
#   
#   out <- apply(combat_ready_exp, 2, function(x) {length(table(x, useNA = useNA))}) 
#   
#   which(out==1) 
#   
# } 

# combat_ready_exp <- combat_ready_exp[,-zeroVar(combat_ready_exp)] 
# 
# #  transpose data frame again to apply ComBat so samples are columns and genes are rows 
# combat_ready_exp <- (t(combat_ready_exp)) 

#ComBat correction 
# 
# 
# # combat.expr <- ComBat(dat = scaled_exp, 
# #                       
# #                       batch = batch$study, 
# #                       
# #                       mod = NULL, 
# #                       
# #                       par.prior = TRUE, 
# #                       
# #                       prior.plots = FALSE) 
# 
# 
# # transpose matrix again to tSNE or PCA
# 
# combat.expr <-t(combat.expr)
# 
# # # write csv of combat expr
# # 
# # write.csv (combat.expr, file = "combat_cancer_only_exp.csv")
# 
# #Carry out PCA of combatted data 
# 
# combat.expr_PCA <-PCA(combat.expr, graph = FALSE)
# 
# #plot PCA
# ggbiplot(combat.expr_PCA)
# 

# # specify PN_gene categories
# gene_categories <- colnames(PN_gene_exp)
# gene_categories <- gene_categories [9:14]


# #carry out tSNE on expression data for several perplexities coloured for expression ####
# # for each group of PN genes 
# 
# p_number <-c((seq(from = 200, to = 100, by = -10)))
# 
# for (p in p_number) {
#   print(paste0('Running tSNE with perplexity = ',p))
#   set.seed(13)
#   PN_gene_exp_full <- Rtsne(scaled_exp, perplexity = p)
#   PN_gene_exp.tsne <- data.frame(TSNE1 = PN_gene_exp_full$Y[,1],
#                                  TSNE2 = PN_gene_exp_full$Y[,2])
#   Exp_tsne <- cbind(TCGA_GTEX_FPKM_Skin_with_gene_name, PN_gene_exp.tsne)
#   
#   
#   
#   pdf(file = paste0("GTEX_TCGA_PCA.pdf"), width = 8, height = 4)  
#   
#   g <- ggplot(data = Exp_tsne, aes(x = TSNE1, y = TSNE2
#                                    #, color = Exp_tsne$chap_mean, shape=Cancer_Status
#                                    )
#                                    ) + geom_point(size = 0.5) +
#    # scale_colour_gradient(low = 'white', high = 'navyblue') +
#     
#    # labs(title=paste0('Cancer tSNE (ComBatted) with perplexity = ',p, '. Mean Chaperone Gene Expression.')) +
#     theme(
#       legend.title = element_text(size = 14),
#       legend.text = element_text(size = 10) ,
#       panel.background = element_rect(fill = "gray86")) 
#   guides(colour = guide_legend(override.aes = list(size=4))) 
#   
#   print(g)
#   
#   dev.off()
#   
#   
#   
#   pdf(file = paste0("Cancer_CoChaperone_gene_exp_tsne_perp=",p,"_.pdf"), width = 8, height = 4)
#   
#   g <- ggplot(data = Exp_tsne, aes(x = TSNE1, y = TSNE2, color = Exp_tsne$cochap_mean, shape=Cancer_Status)) + geom_point(size = 0.5) +
#     scale_colour_gradient(low = 'white', high = 'dodgerblue4') +
#     
#     labs(title=paste0('Cancer tSNE (ComBatted) with perplexity = ',p, '. Mean CoChaperone Gene Expression.')) +
#     theme(
#       legend.title = element_text(size = 14),
#       legend.text = element_text(size = 10) ,
#       panel.background = element_rect(fill = "gray86")) 
#   guides(colour = guide_legend(override.aes = list(size=4))) 
#   print(g)
#   
#   dev.off()
#   
#   
# 
