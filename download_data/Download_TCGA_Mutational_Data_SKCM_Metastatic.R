if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maftools")


install.packages("DT")

install.packages("ensembldb")
  
  
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(DT)
library(ggplot2)
library(ensembldb)


##### SET WORKING DIRECTORY , use getwd to designate new variable'current_dir' to use in sub directory names

setwd("~/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()

# ## Specify gene set 
# 
# GeneSet <- 'DNAJ&HSP&CCT'
# 
# # Set genes of interest
GenesOfInterest_df <- read.xlsx("Regenrich_tables_met_and_prim_with_min.xlsx", sheetName = "top30")

GenesOfInterest <- GenesOfInterest_df$gene_name

## Read in Studies of Interest

# CT<-c("BRCA","CHOL","KIRC","LUSC","READ","PRAD","KICH","LIHC","BLCA","STAD","THCA","UCEC")

#for i in CT {
i <-'SKCM'

#current_dir<-getwd()
  
  print(i)
  
    cancer_type<-paste0("TCGA-",i,sep="")
##MAf = mutation annotation file
    query <- GDCquery(project = "TCGA-SKCM", 
                         data.category = "Simple Nucleotide Variation", # Simple nucleotide variation if legacy
                data.type = "Masked Somatic Mutation",
                         access = "open", 
                         legacy = F
                ,
                #   sample.type = "Metastatic"
                )
GDCdownload(query)
maf <- GDCprepare(query)

df.maf <- maf
df.maf.short <- df.maf[,c("Tumor_Sample_Barcode",
                          "Hugo_Symbol",
                           "Chromosome",
                           "Start_Position","End_Position",
                          "Variant_Classification",
                          "Variant_Type",
                          "Reference_Allele",
                         "Tumor_Seq_Allele1",
                          "Tumor_Seq_Allele2",
                          "Consequence",
                          "SIFT",
                          "PolyPhen"
                          #, 
                          #'Amino_acids'
                          )]


# snvs <- snvs[which((snvs$Tumor_Seq_Allele2 %in% 
#                       c("A","C","G","T"))&
#                      (snvs$Reference_Allele %in%
#                         c("A","C","G","T"))),]



df.maf.short$SampleType <- sapply(df.maf.short$Tumor_Sample_Barcode, function(x)
  paste(strsplit(x,"-")[[1]][4], collapse="-"))
df.maf.selected <- df.maf.short[df.maf.short$SampleType %like% c("06"), ] 
df.maf.final <- df.maf.selected
df.maf.selected <- df.maf.selected %>% read.maf
rm(maf)


#Save df.maf.seleced as RData

save(df.maf.final, file = 'SKCM_TCGA.metastatic.mutations_Simple somatic mutation.RData')

save(df.maf, file = 'SKCM_TCGA.metastatic.mutations_Simple somatic mutation_long.RData')


##Note the only variants that are predicted to be damaging by SIFT and PolyPhen are missense mutations 
# #and translation start site variants
# 
pdf (paste0(current_dir, "/mutation_data_summary_charts_met.pdf"))

plotsummary <- plotmafSummary(maf = df.maf.selected,
               rmOutlier = TRUE, addStat = 'median',
               dashboard = TRUE)

print (plotsummary )

dev.off()


# # assign name to gene list for each cluster
# assign(paste0(i,  "_mutations_plots"), plotsummary)

# Create Oncoplot of genes of interest
pdf ("SKCM_oncoplot_Met.pdf", width=60,height=50)
oncoplot <- oncoplot(maf = df.maf.selected, genes = GenesOfInterest, showTumorSampleBarcodes = TRUE, 
                     barcode_mar=30, gene_mar = 30, fontSize = 1.3, SampleNamefontSize = 1.8, legendFontSize = 4,titleFontSize = 5)
print (plotsummary )

dev.off()

#Download clinical data
query <- GDCquery(project = "TCGA-SKCM", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)








# Create oncoplot



skcm.maf <- system.file("extdata", "TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.maf", package = "maftools")
#skcm.clin = system.file('extdata', 'tcga_skcm_annot.tsv', package = 'maftools')
skcm <- read.maf(maf = skcm.maf
                 #, clinicalData = skcm.clin
                 )
#Basic oncoplot
oncoplot(maf = skcm, top = 3)
#Changing colors for variant classifications (You can use any colors, here in this example we will use a color palette from RColorBrewer)
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
#Color coding for FAB classification; try getAnnotations(x = skcm) to see available annotations.
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
fabcolors = list(FAB_classification = fabcolors)
oncoplot(maf = skcm, colors = col, clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, annotationColor = fabcolors)


