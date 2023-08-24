# 
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
 library(ggpubr)
 library("readxl")
 library('xlsx')
##
## create moveme function to move columns
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

#load mutation data

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

SKCM_signatures<- loadRData('SKCM_TCGA.primary.mutations.RData')


#load clinical and expression data from 
signature_table <- loadRData('sigs_SKCM_allSamples.Rdata')

###load table of sample groups ####

sample_groups <- read.xlsx ('Final_Corrected_Expression_primary_and_metastatic_sample_groups.xlsx', sheetIndex = 1)

sample_groups_table <-sample_groups[,c("Sample","Final_Sample_Group_Paper")]
sample_groups_table<-data.frame(sample_groups_table, row.names = 1)

  
### merge sample groups and signatures tables
signatures_and_sample_groups <- merge(sample_groups_table,signature_table, by =0)

###create tables with just metastatic or primary samples

signatures_and_sample_groups_met <- subset(signatures_and_sample_groups, 
                                           Final_Sample_Group_Paper %in% c("Metastatic B", "Metastatic A"))

signatures_and_sample_groups_prim <- subset(signatures_and_sample_groups, 
                                           Final_Sample_Group_Paper %in% c("Primary B", "Primary A"))

####Primary ####


########Create balloon plot showing contribution of signatures to each patient ###

#Sort by cluster

signatures_and_sample_groups_prim <-  signatures_and_sample_groups_prim[order(signatures_and_sample_groups_prim$Final_Sample_Group_Paper ),]

##### Prepare table for balloon plot
library(dplyr)
#Create copy
SKCM_exp_info_clinical_primary_sig_sub <-signatures_and_sample_groups_prim

#Remove signature columns with no values above 0
SKCM_exp_info_clinical_primary_sig_sub <- SKCM_exp_info_clinical_primary_sig_sub[, colSums(SKCM_exp_info_clinical_primary_sig_sub != 0) > 0]

#calculate column means

mns <- (colMeans(SKCM_exp_info_clinical_primary_sig_sub [,3:14]))
#convert to dataframe
mns<- as.data.frame(mns)
#convert rownames to column
mns$Signature <- rownames(mns)

#sort by mean
mns  <-  mns [order(-mns$mns),]

mnst <- t(mns)
mnst <- as.data.frame(mnst)
#bring columns SBS 7b, 7c and 7d after SBS7a
mnst <-
  mnst[moveme(names(mnst), "SBS7a first ")]

mnst <-
  mnst[moveme(names(mnst), "SBS7b after SBS7a ")]
mnst <-
  mnst[moveme(names(mnst), "SBS7c after SBS7b ")]
mnst <-
  mnst[moveme(names(mnst), "SBS7d after SBS7c ")]

### rename Row name column Sample
names(SKCM_exp_info_clinical_primary_sig_sub)[1] <- "Sample"
names(SKCM_exp_info_clinical_primary_sig_sub)[2] <- "Sample_Group"

### specify column order


  col_order <- c('Sample', 'Sample_Group','SBS7a',
    'SBS7b',
    'SBS7c',
    'SBS7d',
    'SBS38',
    'SBS5',
    'SBS40',
    'SBS1',
    'SBS43',
    'SBS2',
    'SBS13',
    'SBS58')
  SKCM_exp_info_clinical_primary_sig_sub <- SKCM_exp_info_clinical_primary_sig_sub[, col_order]  
  
  
  ### add in total UV exposure coliumns
  
  SKCM_exp_info_clinical_primary_sig_sub_with_totals<-    SKCM_exp_info_clinical_primary_sig_sub
  
  SKCM_exp_info_clinical_primary_sig_sub_with_totals$SBS7_total <- rowSums(SKCM_exp_info_clinical_primary_sig_sub_with_totals[ ,c('SBS7a','SBS7b','SBS7c', 'SBS7d')])
  
 ###read in suplementary clinical data
  
  sup_clinical_df <- read.csv('TCGA_sup_clinical_data_with_sample_groups_and_main_clinical_data_no_dups.csv') 
 
  ###merge clinical and signature data
  
  SKCM_exp_info_supplementary_clinical_primary_sig_sub_with_totals_complete <- merge (SKCM_exp_info_clinical_primary_sig_sub_with_totals,sup_clinical_df, by = "Sample" )
 
   ###rename sample group column
  names(SKCM_exp_info_supplementary_clinical_primary_sig_sub_with_totals_complete)[2] <- "Sample_Group" 
                                                                        
    
  ###read in clinical data
  
clinical_df <- read.csv('SKCM_exp_info_clinical_November_22.csv') 
  
  ###merge clinical and signature data
  
  SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete <- merge (SKCM_exp_info_clinical_primary_sig_sub_with_totals,clinical_df, by = "Sample" )
  
  ###rename sample group column
  names(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete)[2] <- "Sample_Group" 
  
  ###remove duplicates by sample
  SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete <- (SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete[!duplicated(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete$Sample), ])
  
  ####box plots of proportions  of signatures ##€##
  sigs <- colnames (SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete[3:15])
  
  
  j <- ggplot(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete , aes(x= Sample_Group , y= SBS13, fill= Sample_Group) )+ 
    geom_boxplot() +
    
    scale_fill_manual( name = "Sample Group", values=c(  "#296a8c",
                                                         "#fde725" )) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    #scale_fill_manual(values=c( '', '')) +
    # scale_x_discrete(labels=c('1'=' Cluster 1' ,  '2' ='Cluster 2'  )) +
    stat_compare_means(aes(group = Sample_Group), label = "p.format",
                       size = 8,
                       label.y =0.95
                       )  +
    # geom_jitter() +
    
    labs(
      title= paste0 ("Proportion of ", 'SBS13', " Signature in Primary Samples"),
      x ="Sample Group", y = "Proportion") +
    theme(
      
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20) ,
      # panel.background = element_rect(fill = "white"),
      axis.text=element_text(size=14),
      axis.title=element_text(size=16,face="bold"),
      panel.background = element_rect(fill = 'white', colour = 'black')
      
    ) 
  pdf(paste0("Proportion_of_SBS13_in_primary_samples.pdf"), width = 8, height = 5)
  
  print(j)
  dev.off() 
  
  ### scatter plots ###
  
  # Scatter plot by group 7 vs sbs38

 t <-ggplot(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete, 
            aes(x = SBS7_total, y = SBS38,  color = Sample_Group))+
              geom_point()+
                scale_color_manual(values=c("#296a8c","#fde725"))+
   theme(
     
     legend.title = element_text(size = 20),
     legend.text = element_text(size = 20) ,
     # panel.background = element_rect(fill = "white"),
     axis.text=element_text(size=14),
     axis.title=element_text(size=16,face="bold"),
     panel.background = element_rect(fill = 'white', colour = 'black')
     
   ) 
    
  
 pdf("Scatter SBS7 vs SBS38.pdf", width = 8, height = 5)
 
 print(t)
 dev.off() 
 
 
 # Scatter plot by group 40 vs sbs38
 
 t <-ggplot(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete, 
            aes(x = SBS40, y = SBS38,  color = Sample_Group))+
   geom_point()+
   scale_color_manual(values=c("#296a8c","#fde725"))+
   theme(
     
     legend.title = element_text(size = 20),
     legend.text = element_text(size = 20) ,
     # panel.background = element_rect(fill = "white"),
     axis.text=element_text(size=14),
     axis.title=element_text(size=16,face="bold"),
     panel.background = element_rect(fill = 'white', colour = 'black')
     
   ) 
 
 
 pdf("Scatter SBS40 vs SBS38.pdf", width = 8, height = 5)
 
 print(t)
 dev.off() 
 
 # Scatter plot by group sbs38 vs age
 
 t <-ggplot(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete, 
            aes(x = age_at_initial_pathologic_diagnosis, y = SBS38,  color = Sample_Group))+
   geom_point()+
   scale_color_manual(values=c("#296a8c","#fde725"))+
   theme(
     
     legend.title = element_text(size = 20),
     legend.text = element_text(size = 20) ,
     # panel.background = element_rect(fill = "white"),
     axis.text=element_text(size=14),
     axis.title=element_text(size=16,face="bold"),
     panel.background = element_rect(fill = 'white', colour = 'black')
     
   ) 
 
 
 pdf("Scatter age vs SBS38.pdf", width = 20, height = 5)
 
 print(t)
 dev.off() 
 
 # ###histogram of ratio of SBS7/SBS38
 # SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete$SBS7_SBS38_ratio<-
 #   (100*SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete$SBS7_total) /
 #   (100*SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete$SBS38)
 # 
 # 
 # 
 # pdf("ratio of SBS7 to SBS38.pdf", width = 20, height = 5)
 # 
 # hist(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete$SBS7_SBS38_ratio)
 # 
 # dev.off() 
 # 
 # write.csv(SKCM_exp_info_clinical_primary_sig_sub_with_totals_complete, file = 'SKCM_mut_sigs_clinical.csv ')
 # 
  ###balloon plot#### 
#melt into 2 column dataframe

### melt dataframe so ID comprises Patient Barcode and SampleTypeCateg for all


SKCM_exp_info_clinical_primary_sig_sub_melt <- melt(SKCM_exp_info_clinical_primary_sig_sub, id=c(1,2))

### Rename columns in melted dataframe

colnames(SKCM_exp_info_clinical_primary_sig_sub_melt)[colnames(SKCM_exp_info_clinical_primary_sig_sub_melt)=="variable"] <- "Signature"
colnames(SKCM_exp_info_clinical_primary_sig_sub_melt)[colnames(SKCM_exp_info_clinical_primary_sig_sub_melt)=="value"] <- "Proportion"

# This script forces plot to maintain the row order of the source dataframe
SKCM_exp_info_clinical_primary_sig_sub_melt$Signature<-
  factor(SKCM_exp_info_clinical_primary_sig_sub_melt$Signature,
     levels= c('SBS7a',
               'SBS7b',
               'SBS7c',
               'SBS7d',
               'SBS38',
               'SBS5',
               'SBS40',
               'SBS1',
               'SBS43',
               'SBS2',
               'SBS13',
               'SBS58'))
SKCM_exp_info_clinical_primary_sig_sub_melt$Sample <- factor(SKCM_exp_info_clinical_primary_sig_sub_melt$Row.names,
                                                                levels= SKCM_exp_info_clinical_primary_sig_sub$Row.names)

##remove rownames column
SKCM_exp_info_clinical_primary_sig_sub_melt<- SKCM_exp_info_clinical_primary_sig_sub_melt[,-1]

# # Replace values == 0.000000000 with 0
#
# SKCM_exp_info_clinical_primary_sig_sub_melt$Proportion[SKCM_exp_info_clinical_primary_sig_sub_melt$Proportion == 0.000000000] <- NA
#

#create melted table with no patient Barcode


SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat <- melt(SKCM_exp_info_clinical_primary_sig_sub[,2:14], id=c(1))
### Rename columns in melted dataframe
#
colnames(SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat)[colnames(SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat)=="variable"] <- "Signature"
colnames(SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat)[colnames(SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat)=="value"] <- "Proportion"


# This script forces plot to maintain the row order of the source dataframe
SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat$Signature<-
  factor(SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat$Signature,
         levels= c(levels= c('SBS7a',
                             'SBS7b',
                             'SBS7c',
                             'SBS7d',
                             'SBS38',
                             'SBS5',
                             'SBS40',
                             'SBS1',
                             'SBS43',
                             'SBS2',
                             'SBS13',
                             'SBS58') ))




#Balloon Plot showing proportion of each signature with samples in cluster overlaid


hplot <- ggballoonplot(
  SKCM_exp_info_clinical_primary_sig_sub_melt_no_pat,
  x = 'Signature',
  y = 'Final_Sample_Group_Paper',
  size = "Proportion",
  facet.by = NULL,
  size.range = c(-1, 20),
  shape = 21,
  color = "black",
  fill ='Final_Sample_Group_Paper',
  palette = c('#296a8c', '#fde725'),
  main = "Proportion of each signature in each sample, coloured by sample group",
  show.label = FALSE,
  font.label = list(size = 20, color = "black"),
  rotate.x.text = TRUE,
  text.size = 20,
  zlab = "Fold Change",
  font("title", size = 40) +
    font("x.text", size = 22) +
    font("y.text", size = 30)


)
pdf("Signature_Representation_Balloon_Plot_primary_no_pat.pdf", width = 25, height = 5)

print(hplot)

dev.off()


#Balloon Plot showing proportion of each signature with samples in cluster separate
pdf("Signature_Representation_Balloon_Plot_sep_all_Samples_Jan_23.pdf", width = 40, height = 20)


hplot <- ggballoonplot(
  SKCM_exp_info_clinical_primary_sig_sub_melt,
  x = 'Sample',
  y = 'Signature',
  size = "Proportion",
  facet.by = NULL,
  size.range = c(-0.5, 13),
  shape = 21,
  color = "black",
  fill = "Final_Sample_Group_Paper",
  palette = c('#296a8c', '#fde725'),
  main = "Proportion of each signature in each sample, coloured by sample group",
  show.label = FALSE,
  font.label = list(size = 12, color = "black"),
  rotate.x.text = TRUE,

  text.size = 20,
  font("title", size = 40) +
    font("x.text", size = 22) +
    font("y.text", size = 30)

) +
(theme (panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "gray97"),
        legend.key.size  = unit(3,"cm"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 35, face = "bold")
        ))+
  guides(fill = guide_legend(override.aes = list(size=10)))+
labs(fill = "Sample Group", size ='Proportion')

print(hplot)

dev.off()

# ##Create histograms of proportion of each signature in samples coloured by cluster
# # Map smoke to fill, make the bars NOT stacked, and make them semitransparent
# ggplot(SKCM_exp_info_clinical_primary_sig_sub, aes(x = SBS7a
#                                                    #, fill = cluster4
#                                                    )) +
#   geom_histogram(position = "identity", alpha = 0.4, binwidth  = 0.02)
# 


####SIGNATURE SURVIVAL CUTPOINTS #####
###merge signatures with clinical data

##Determine survival cutpoints for SBS7b proportion

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "OS.time", event = "OS",
                         variables = c("SBS7b"))

summary(res.cut)

# 2. Plot cutpoint for SBS7b
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "SBS7b", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize

fit <- survfit(Surv(OS.time, OS) ~SBS7b, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)



# # calculate survival cutpoint for each signature by a loop
#
# signatures <- c("SBS7b" , "SBS7b",
#                 #    "SBS7d" , "SBS10b" ,"SBS6" ,  "SBS1" ,  "SBS30" , "SBS30",  "SBS11" ,  "SBS23",
#                   #"SBS26"  ,
#                   "SBS22"  )
#
# # 1. Determine the optimal cutpoint of variables
#
#   for (s in signatures){
#
#   results = surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "OS.time", event = "OS"
#                             ,variables = s )
# print(results)
#
# res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "OS.time", event = "OS",
#                          variables = s)
#
# summary(res.cut)
# print(res.cut)
# res.cat <- surv_categorize(res.cut)
# print(head(res.cat))
#
# # print curve
# pdf(paste0(s,"OS_survival_curve.pdf"), width = 10, height = 7)
#
# fit <- survfit(Surv(OS.time, OS) ~ s, data = res.cat)
#
#
# s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE,
#                     ggtheme = theme_classic2(base_size=20))
#
#
# print(s_plot)
# dev.off()
#
#
#
#
# }


##Determine overall survival cutpoints for SBS7d proportion

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "OS.time", event = "OS",
                         variables = c("SBS22"))

summary(res.cut)

# 2. Plot cutpoint for SBS22
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "SBS22", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize
pdf("survival_curve_SBS22.pdf", width = 10, height = 7)

fit <- survfit(Surv(OS.time, OS) ~SBS22, data = res.cat)


s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE,
                    ggtheme = theme_classic2(base_size=20))


print(s_plot)
dev.off()

##Determine survival cutpoints for Each signature ####

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "OS.time", event = "OS",
                         variables = c("SBS22"))

summary(res.cut)

# 2. Plot cutpoint for SBS22d
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "SBS22", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize

fit <- survfit(Surv(OS.time, OS) ~SBS22, data = res.cat)
OS_pvalue <- surv_pvalue(fit, res.cat)


# pdf("survival_curve_SBS22.pdf", width = 10, height = 7)
#
# s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE,
#                     ggtheme = theme_classic2(base_size=20))
#
#
# print(s_plot)
# dev.off()

##Determine DS survival cutpoints for SBS22 proportion

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "DSS.time", event = "DSS",
                         variables = c("SBS22"))

summary(res.cut)

# 2. Plot cutpoint for SBS22
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "SBS22", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize


fit <- survfit(Surv(DSS.time, DSS) ~SBS22, data = res.cat)
DSS_pvalue <- surv_pvalue(fit, res.cat)



# pdf("survival_curve_SBS22.pdf", width = 10, height = 7)
# s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE,
#                     ggtheme = theme_classic2(base_size=20))
#
#
# print(s_plot)
# dev.off()


##Determine PFI survival cutpoints for SBS22 proportion

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "PFI.time", event = "PFI",
                         variables = c("SBS22"))

summary(res.cut)

# 2. Plot cutpoint for SBS22
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "SBS22", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# # 4. Fit survival curves and visualize
# pdf("survival_curve_SBS22.pdf", width = 10, height = 7)

fit <- survfit(Surv(PFI.time, PFI) ~SBS22, data = res.cat)
PFI_pvalue <- surv_pvalue(fit, res.cat)




# s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE,
#                     ggtheme = theme_classic2(base_size=20))

#
# print(s_plot)
dev.off()

print(paste('Signature SBS22, : OS ',OS_pvalue ,', DSS ',DSS_pvalue,', PFI ', PFI_pvalue ))

signatures <- c("SBS7a" , "SBS7b",    "SBS7d" , "SBS10b" ,"SBS6" ,  "SBS1" ,
                "SBS30" , "SBS31",  "SBS11" ,  "SBS23",  "SBS26"  ,        "SBS22"  )


t.test(SBS22~ Cluster, data = SKCM_exp_info_clinical_primary_sig_sub)

#count number of non 0 sig proportions in each sample

SKCM_exp_info_clinical_primary_sig$CountSigs <- rowSums(SKCM_exp_info_clinical_primary_sig[478:527] != 0)

#t test difference in number of signatures present in each sample

t.test(CountSigs~ Cluster, data = SKCM_exp_info_clinical_primary_sig)

#count number of >0.12 sig proportions in each sample

SKCM_exp_info_clinical_primary_sig$CountSigs12 <- rowSums(SKCM_exp_info_clinical_primary_sig[478:527] > 0.12 )

#t test difference in number of signatures present in each sample

t.test(CountSigs12~ Cluster, data = SKCM_exp_info_clinical_primary_sig)

#Anova of 4 clusters
res.aov <- aov(SBS7a ~ cluster4, data = SKCM_exp_info_clinical_primary_sig)
# Summary of the analysis
summary(res.aov)

#Save SKCM_exp_info_clinical as Rdata

save(SKCM_exp_info_clinical_primary_sig, file = 'SKCM_exp_info_clinical_primary_signatures.RData')

#write csv
write.csv (SKCM_exp_info_clinical_primary_sig, file = "SKCM_exp_info_clinical_primary_signatures.csv")

##Determine survival cutpoints for Each signature ####

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "OS.time", event = "OS",
                         variables = c("CountSigs02"))

summary(res.cut)

# 2. Plot cutpoint for CountSigs02d
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "CountSigs02", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize

fit <- survfit(Surv(OS.time, OS) ~CountSigs02, data = res.cat)
OS_pvalue <- surv_pvalue(fit, res.cat)


pdf("survival_curve_CountSigs02.pdf", width = 10, height = 7)

s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE,
                    ggtheme = theme_classic2(base_size=20))


print(s_plot)
dev.off()

##Determine DS survival cutpoints for CountSigs02 proportion

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "DSS.time", event = "DSS",
                         variables = c("CountSigs02"))

summary(res.cut)

# 2. Plot cutpoint for CountSigs02
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "CountSigs02", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# 4. Fit survival curves and visualize


fit <- survfit(Surv(DSS.time, DSS) ~CountSigs02, data = res.cat)
DSS_pvalue <- surv_pvalue(fit, res.cat)



# pdf("survival_curve_CountSigs02.pdf", width = 10, height = 7)
# s_plot <-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = FALSE,
#                     ggtheme = theme_classic2(base_size=20))
#
#
# print(s_plot)
# dev.off()


##Determine PFI survival cutpoints for CountSigs02 proportion

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(SKCM_exp_info_clinical_primary_sig, time = "PFI.time", event = "PFI",
                         variables = c("CountSigs02"))

summary(res.cut)

# 2. Plot cutpoint for CountSigs02
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "CountSigs02", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# # 4. Fit survival curves and visualize
 pdf("survival_curve_CountSigs02.pdf", width = 10, height = 7)

fit <- survfit(Surv(PFI.time, PFI) ~CountSigs02, data = res.cat)
PFI_pvalue <- surv_pvalue(fit, res.cat)

print(s_plot)
dev.off()

print(paste('CountSigs02, : OS ',OS_pvalue ,', DSS ',DSS_pvalue,', PFI ', PFI_pvalue ))

####Metastatic ####


########Create balloon plot showing contribution of signatures to each patient ###

#Sort by cluster

signatures_and_sample_groups_met <-  signatures_and_sample_groups_met[order(signatures_and_sample_groups_met$Final_Sample_Group_Paper ),]

##### Prepare table for balloon plot
library(dplyr)
#Create copy
SKCM_exp_info_clinical_Metastatic_sig_sub <-signatures_and_sample_groups_met

#Remove signature columns with no values above 0
SKCM_exp_info_clinical_Metastatic_sig_sub <- SKCM_exp_info_clinical_Metastatic_sig_sub[, colSums(SKCM_exp_info_clinical_Metastatic_sig_sub != 0) > 0]

#calculate column means

mns <- (colMeans(SKCM_exp_info_clinical_Metastatic_sig_sub [,3:14]))
#convert to dataframe
mns<- as.data.frame(mns)
#convert rownames to column
mns$Signature <- rownames(mns)

#sort by mean
mns  <-  mns [order(-mns$mns),]

mnst <- t(mns)
mnst <- as.data.frame(mnst)
#bring columns SBS 7b, 7c and 7d after SBS7a
mnst <-
  mnst[moveme(names(mnst), "SBS7a first ")]

mnst <-
  mnst[moveme(names(mnst), "SBS7b after SBS7a ")]
mnst <-
  mnst[moveme(names(mnst), "SBS7c after SBS7b ")]
mnst <-
  mnst[moveme(names(mnst), "SBS7d after SBS7c ")]

### rename Row name column Sample
names(SKCM_exp_info_clinical_Metastatic_sig_sub)[1] <- "Sample"
names(SKCM_exp_info_clinical_Metastatic_sig_sub)[2] <- "Sample_Group"

### specify column order


col_order <- c('Sample', 'Sample_Group','SBS7a',
               'SBS7b',
               'SBS7c',
               'SBS7d',
               'SBS38',
               'SBS5',
               'SBS40',
               'SBS1',
               'SBS43',
               'SBS2',
               'SBS13',
               'SBS58')
SKCM_exp_info_clinical_Metastatic_sig_sub <- SKCM_exp_info_clinical_Metastatic_sig_sub[, col_order]  


### add in total UV exposure coliumns

SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals<-    SKCM_exp_info_clinical_Metastatic_sig_sub

SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals$SBS7_total <- rowSums(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals[ ,c('SBS7a','SBS7b','SBS7c', 'SBS7d')])

###read in suplementary clinical data

sup_clinical_df <- read.csv('TCGA_sup_clinical_data_with_sample_groups_and_main_clinical_data_no_dups.csv') 

###merge clinical and signature data

SKCM_exp_info_supplementary_clinical_Metastatic_sig_sub_with_totals_complete <- merge (SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals,sup_clinical_df, by = "Sample" )

###rename sample group column
names(SKCM_exp_info_supplementary_clinical_Metastatic_sig_sub_with_totals_complete)[2] <- "Sample_Group" 


###read in clinical data

clinical_df <- read.csv('SKCM_exp_info_clinical_November_22.csv') 

###merge clinical and signature data

SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete <- merge (SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals,clinical_df, by = "Sample" )

###rename sample group column
names(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete)[2] <- "Sample_Group" 

###remove duplicates by sample
SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete <- (SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete[!duplicated(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete$Sample), ])

####box plots of proportions  of signatures ##€##
sigs <- colnames (SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete[3:15])


j <- ggplot(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete , aes(x= Sample_Group , y= SBS7_total, fill= Sample_Group) )+ 
  geom_boxplot() +
  
  scale_fill_manual( name = "Sample Group", values=c(  "#7bd5d5" , "#7f4e97" )) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  #scale_fill_manual(values=c( '', '')) +
  # scale_x_discrete(labels=c('1'=' Cluster 1' ,  '2' ='Cluster 2'  )) +
  stat_compare_means(aes(group = Sample_Group), label = "p.format",
                     size = 8,
                     label.y =0.95
  )  +
  # geom_jitter() +
  
  labs(
    title= paste0 ("Proportion of ", 'SBS7_total', " Signature in Metastatic Samples"),
    x ="Sample Group", y = "Proportion") +
  theme(
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20) ,
    # panel.background = element_rect(fill = "white"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16,face="bold"),
    panel.background = element_rect(fill = 'white', colour = 'black')
    
  ) 
pdf(paste0("Proportion_of_SBS7_total_in_Metastatic_samples.pdf"), width = 8, height = 5)

print(j)
dev.off() 

# ### scatter plots ###
# 
# # Scatter plot by group 7 vs sbs38
# 
# t <-ggplot(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete, 
#            aes(x = SBS7_total, y = SBS38,  color = Sample_Group))+
#   geom_point()+
#   scale_color_manual(values=c("#296a8c","#fde725"))+
#   theme(
#     
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 20) ,
#     # panel.background = element_rect(fill = "white"),
#     axis.text=element_text(size=14),
#     axis.title=element_text(size=16,face="bold"),
#     panel.background = element_rect(fill = 'white', colour = 'black')
#     
#   ) 
# 
# 
# pdf("Scatter SBS7 vs SBS38.pdf", width = 8, height = 5)
# 
# print(t)
# dev.off() 
# 
# 
# # Scatter plot by group 40 vs sbs38
# 
# t <-ggplot(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete, 
#            aes(x = SBS40, y = SBS38,  color = Sample_Group))+
#   geom_point()+
#   scale_color_manual(values=c("#296a8c","#fde725"))+
#   theme(
#     
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 20) ,
#     # panel.background = element_rect(fill = "white"),
#     axis.text=element_text(size=14),
#     axis.title=element_text(size=16,face="bold"),
#     panel.background = element_rect(fill = 'white', colour = 'black')
#     
#   ) 
# 
# 
# pdf("Scatter SBS40 vs SBS38.pdf", width = 8, height = 5)
# 
# print(t)
# dev.off() 
# 
# # Scatter plot by group sbs38 vs age
# 
# t <-ggplot(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete, 
#            aes(x = age_at_initial_pathologic_diagnosis, y = SBS38,  color = Sample_Group))+
#   geom_point()+
#   scale_color_manual(values=c("#296a8c","#fde725"))+
#   theme(
#     
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 20) ,
#     # panel.background = element_rect(fill = "white"),
#     axis.text=element_text(size=14),
#     axis.title=element_text(size=16,face="bold"),
#     panel.background = element_rect(fill = 'white', colour = 'black')
#     
#   ) 
# 
# 
# pdf("Scatter age vs SBS38.pdf", width = 20, height = 5)
# 
# print(t)
# dev.off() 

# ###histogram of ratio of SBS7/SBS38
# SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete$SBS7_SBS38_ratio<-
#   (100*SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete$SBS7_total) /
#   (100*SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete$SBS38)
# 
# 
# 
# pdf("ratio of SBS7 to SBS38.pdf", width = 20, height = 5)
# 
# hist(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete$SBS7_SBS38_ratio)
# 
# dev.off() 
# 
# write.csv(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete, file = 'SKCM_mut_sigs_clinical.csv ')
# 
###balloon plot#### 
#melt into 2 column dataframe

### melt dataframe so ID comprises Patient Barcode and SampleTypeCateg for all


SKCM_exp_info_clinical_Metastatic_sig_sub_melt <- melt(SKCM_exp_info_clinical_Metastatic_sig_sub, id=c(1,2))

### Rename columns in melted dataframe

colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt)[colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt)=="variable"] <- "Signature"
colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt)[colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt)=="value"] <- "Proportion"

# This script forces plot to maintain the row order of the source dataframe
SKCM_exp_info_clinical_Metastatic_sig_sub_melt$Signature<-
  factor(SKCM_exp_info_clinical_Metastatic_sig_sub_melt$Signature,
         levels= c('SBS7a',
                   'SBS7b',
                   'SBS7c',
                   'SBS7d',
                   'SBS38',
                   'SBS5',
                   'SBS40',
                   'SBS1',
                   'SBS43',
                   'SBS2',
                   'SBS13',
                   'SBS58'))
SKCM_exp_info_clinical_Metastatic_sig_sub_melt$Sample <- factor(SKCM_exp_info_clinical_Metastatic_sig_sub_melt$Sample,
                                                             levels= SKCM_exp_info_clinical_Metastatic_sig_sub$Sample)

# ##remove rownames column
# SKCM_exp_info_clinical_Metastatic_sig_sub_melt<- SKCM_exp_info_clinical_Metastatic_sig_sub_melt[,-1]

# # Replace values == 0.000000000 with 0
#
# SKCM_exp_info_clinical_Metastatic_sig_sub_melt$Proportion[SKCM_exp_info_clinical_Metastatic_sig_sub_melt$Proportion == 0.000000000] <- NA
#

#create melted table with no patient Barcode


SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat <- melt(SKCM_exp_info_clinical_Metastatic_sig_sub[,2:14], id=c(1))
### Rename columns in melted dataframe
#
colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat)[colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat)=="variable"] <- "Signature"
colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat)[colnames(SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat)=="value"] <- "Proportion"


# This script forces plot to maintain the row order of the source dataframe
SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat$Signature<-
  factor(SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat$Signature,
         levels= c(levels= c('SBS7a',
                             'SBS7b',
                             'SBS7c',
                             'SBS7d',
                             'SBS38',
                             'SBS5',
                             'SBS40',
                             'SBS1',
                             'SBS43',
                             'SBS2',
                             'SBS13',
                             'SBS58') ))




#Balloon Plot showing proportion of each signature with samples in cluster overlaid


hplot <- ggballoonplot(
  SKCM_exp_info_clinical_Metastatic_sig_sub_melt_no_pat,
  x = 'Signature',
  y = 'Sample_Group',
  size = "Proportion",
  facet.by = NULL,
  size.range = c(-1, 20),
  shape = 21,
  color = "black",
  fill ='Sample_Group',
  palette = c("#7bd5d5" , "#7f4e97"),
  main = "Proportion of each signature in each sample, coloured by sample group",
  show.label = FALSE,
  font.label = list(size = 20, color = "black"),
  rotate.x.text = TRUE,
  text.size = 20,
  zlab = "Fold Change",
  font("title", size = 40) +
    font("x.text", size = 22) +
    font("y.text", size = 30)
  
  
)
pdf("Signature_Representation_Balloon_Plot_Metastatic_no_pat.pdf", width = 25, height = 5)

print(hplot)

dev.off()


#Balloon Plot showing proportion of each signature with samples in cluster separate


hplot <- ggballoonplot(
  SKCM_exp_info_clinical_Metastatic_sig_sub_melt,
  x = 'Sample',
  y = 'Signature',
  size = "Proportion",
  facet.by = NULL,
  size.range = c(-0.5, 8),
  shape = 21,
  color = 'Sample_Group',
  fill = 'Sample_Group',
  palette = c("#7bd5d5" , "#7f4e97"),
  main = "Proportion of each signature in each sample, coloured by sample group",
  show.label = FALSE,
  font.label = list(size = 12, color = "black"),
  rotate.x.text = TRUE,
  
  text.size = 20,
  font("title", size = 40) +
    font("x.text", size = 22) +
    font("y.text", size = 30)
  
) +
  (theme (panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(colour = "gray97"),
          legend.key.size  = unit(3,"cm"),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 35, face = "bold")
  ))+
  guides(fill = guide_legend(override.aes = list(size=10)))+
  labs(fill = "Sample_Group", size ='Proportion')


pdf("Signature_Representation_Balloon_Plot_sep_all_Samples_Met_Jan_23.pdf", width = 80, height = 20)

print(hplot)

dev.off()


# Scatter plot by group 7 vs sbs38

t <-ggplot(SKCM_exp_info_clinical_Metastatic_sig_sub_with_totals_complete, 
           aes(x = SBS7_total, y = SBS38,  color = Sample_Group))+
  geom_point()+
  scale_color_manual(values=c("#7bd5d5" , "#7f4e97"))+
  theme(
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20) ,
    # panel.background = element_rect(fill = "white"),
    axis.text=element_text(size=14),
    axis.title=element_text(size=16,face="bold"),
    panel.background = element_rect(fill = 'white', colour = 'black')
    
  ) 


pdf("Scatter SBS7 vs SBS38_met.pdf", width = 8, height = 5)

print(t)
dev.off() 

# ##Create histograms of proportion of each signature in samples coloured by cluster
# # Map smoke to fill, make the bars NOT stacked, and make them semitransparent
# ggplot(SKCM_exp_info_clinical_Metastatic_sig_sub, aes(x = SBS7a
#                                                    #, fill = cluster4
#                                                    )) +
#   geom_histogram(position = "identity", alpha = 0.4, binwidth  = 0.02)
# 

