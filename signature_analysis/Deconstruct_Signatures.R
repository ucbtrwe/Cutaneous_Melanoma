###This code analysesth contribution of mutational signatures to each sample ###

# Enabling automated, flexible mutational signature analysis in all cancer types
setwd("~/iCloud_Drive_Archive/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()

# http://bioinformaticsfmrp.github.io/TCGAbiolinks/index.html

# Load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library('TCGAbiolinks')
install.packages("deconstructSigs")
library('deconstructSigs')
library('BSgenome.Hsapiens.UCSC.hg38')
library('xlsx')
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Initialize variables
# Studies list: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
studies = c(
  # 'LAML',
  # 'ACC',
  # 'BLCA',
  # 'LGG',
  # 'BRCA',
  # 'CESC',
  # 'CHOL',
  # 'LCML'
  # 'COAD',
  # 'CNTL',
  # 'ESCA',
  # 'FPPP',
  #'GBM',
  # 'HNSC',
  #'KICH',
  # 'KIRC',
  # 'KIRP',
  # 'LIHC',
  # 'LUAD',
  # 'LUSC',
  # 'DLBC',
  # 'MESO',
  # 'MISC',
  # 'OV',
  # 'PAAD',
  # 'PCPG',
  # 'PRAD',
  # 'READ',
  #'SARC',
  'SKCM'
  #,
  # 'STAD',
  # 'TGCT',
  # 'THYM',
  # 'THCA',
 # 'UCS'
  # 'UCEC',
  # 'UVM'
)
signature.set = 'Synapse2019'
min.mutation.count = 50

# Load signature data from Synapse2019
load('Synapse.SBS.Signatures.2019.Rdata')

###select signatures that are found in SKCM
signatures.synapse2019_SKCM <- subset(signatures.synapse2019, rownames(signatures.synapse2019) %in% c('SBS1','SBS2','SBS5','SBS7a','SBS7b','SBS7c',
                                                                                                      'SBS7d','SBS13','SBS38','SBS40','SBS43', 'SBS58'))

run.signature.analysis = function(study, signature.profile, mutation.minimum) {
  
  # Download TCGA mutation data for 'study' cancer study
#  mutations <- GDCquery_Maf(tumor = study, pipelines = 'mutect2')
  
#### new code for downloading mutations
  
  query <- GDCquery(
    project = "TCGA-SKCM", 
    data.category = "Simple Nucleotide Variation", 
    access = "open", 
    legacy = FALSE, 
    data.type = "Masked Somatic Mutation", 
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  GDCdownload(query)
  mutations <- GDCprepare(query)
  
  
  
  # Select columns of interest (from 120 variables)
  snvs <- data.frame(mutations[,c('Tumor_Sample_Barcode','Hugo_Symbol',
                                  'Chromosome','Start_Position','End_Position',
                                  'Variant_Classification','Variant_Type',
                                  'Reference_Allele','Tumor_Seq_Allele1',
                                  'Tumor_Seq_Allele2')])
  
  # Only interested in point mutations (i.e. VariantType == 'SNP' or Alleles %in% c('A','T','C','G'))
  snvs <- snvs[snvs$Reference_Allele %in% c('A','C','T','G') &
                 snvs$Tumor_Seq_Allele2 %in% c('A','C','T','G'),]
  
  # Use annotations to remove samples that are not primary tumours
  snvs$Sample_Type_Code <- sapply(snvs$Tumor_Sample_Barcode,
                                  function(tsb) substr(strsplit(tsb,'-')[[1]][4],1,2))
  snvs$Sample_Type <- sapply(snvs$Sample_Type_Code,
                             function(x) ifelse(x=='01','PrimaryTumour',
                                                ifelse(x =='03', 'PrimaryBloodDerivedCancer',
                                                       ifelse(x %in% c('10','11'),'Normal',
                                                              ifelse(x=='06','Metastatic','Other')))))
  snvs_prim <- snvs[snvs$Sample_Type == 'PrimaryTumour',]
  snvs_met <- snvs[snvs$Sample_Type == 'Metastatic',]
  
  # MUTATIONAL SIGNATURE ANALYSIS Primary
  
  # 1. Convert to deconstructSigs input:
  #   Assigns each mutation to one of 96 mutation trinucleotide profiles
  sigs.input.full <- mut.to.sigs.input(mut.ref = snvs_met, # mutation matrix
                                       sample.id = "Tumor_Sample_Barcode",
                                       chr = "Chromosome",
                                       pos = "Start_Position", # SNVs => Start_Position == End_Position)
                                       ref = "Reference_Allele",
                                       alt = "Tumor_Seq_Allele2",
                                       bsg = BSgenome.Hsapiens.UCSC.hg38 # specify reference genome
  )
  
  # 2. Remove samples with < mut.count.limit mutations
  sigs.input <- sigs.input.full[apply(sigs.input.full,1,sum) >= mutation.minimum,]
  
  # 3. Ensure correct signatures obtained
  signatures <- NULL
  if (signature.profile == 'Synapse2019') signatures = signatures.synapse2019_SKCM
  if (signature.profile == 'COSMIC') signatures = signatures.cosmic
  if (is.null(signatures)) print('signature.profile must be "COSMIC" or "Synapse2019"', quote=FALSE)
  
  # 4. Run signature identification
  sigs <- NULL
  for (sample in rownames(sigs.input)) {
    print(paste('Hello',': ',sample,': ',which(sample == rownames(sigs.input)),' of ',nrow(sigs.input)),quote=FALSE)
    sigs_1 = whichSignatures(tumor.ref = sigs.input, # converted mutation data: assignment to 96-element profiles
                             signatures.ref = signatures.synapse2019_SKCM, # known signatures
                             sample.id = sample, # name of sample - subject of looping
                             contexts.needed = TRUE, # normalizes mutation counts
                             signature.cutoff = 0, # discards sufficiently low weightings
                             tri.counts.method = 'default' # if further normalization is required
    )
    sigs <- rbind(sigs,sigs_1$weights)
  }
  
  return(data.frame(sigs))
}

# Run mutational signature analysis for all types listed in 'studies'
for (study.type in studies) {
  result = run.signature.analysis(study = study.type, signature.profile = signatures.synapse2019_SKCM, mutation.minimum = min.mutation.count)
  assign(paste('sigs',study.type,signature.set,sep='.'),result)
  save(list=paste('sigs',study.type,signature.set,sep='.'),file=('SKCM_metastatic_sigs_Filtered_2.Rdata'))
met_sigs_2 <- loadRData('SKCM_metastatic_sigs_Filtered_2.Rdata')
  
   write.xlsx(met_sigs_2,file=('SKCM_metastatic_sigs_Filtered_2.xlsx'))
  }

