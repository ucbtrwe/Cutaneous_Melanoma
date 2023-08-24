
###This code creates a barchart showing the Regenrich scores for the top 35 
###regulators identified by RegEnrich for 3 categories of genes and makes a scatterplot comparing 
## scores in primary and metastatic cohort.  


library('xlsx')

library(ggplot2)
library(dplyr)
library("data.table")    
library(plotly)
library(matrixStats)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}


#BiocManager::install("BioinformaticsFMRP/TCGAbiolinks") 

setwd("~/iCloud_Drive_Archive/Desktop/R_Proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()

####Primary Charts #####
#### load in amalgamated TF tables with minimum rank included #####
prim_TF_with_min <- read.xlsx('Combined_Regenrich_TF_Table_prim_cor_edited.xlsx', sheetName = "top_35")

### melt into 2 columns2
library(reshape2)
prim_TF_with_min_melt <- melt(prim_TF_with_min, id= 1)

### Rename columns in melted dataframe

colnames(prim_TF_with_min_melt)[colnames(prim_TF_with_min_melt)=="variable"] <- "Gene_Group"
colnames(prim_TF_with_min_melt)[colnames(prim_TF_with_min_melt)=="value"] <- "Score"

prim_TF_with_min_melt$Regulator <-factor(prim_TF_with_min_melt$Regulator,
                                  levels = prim_TF_with_min$Regulator )
#### Primary barchart 

# Create grouped bar chart
q<- ggplot(prim_TF_with_min_melt, aes(fill=Gene_Group, y= Score, x= Regulator)) + 
  scale_fill_manual(labels=c("PN Genes" ,"ATP-Independent HSP",
                             "ATP-Dependent HSP"),
                    name = '', values=c ("PN.Genes" = "#46044d","ATP.Independent.HSP" ="#9308a2",
                                         "ATP.Dependent.HSP" ="#e10df7")) +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'bold', size=11),
        axis.text.y = element_text(face = 'bold', size=15),
        axis.title.x = element_text(size=17, face="bold", vjust = 1.5),
        axis.title.y = element_text(size=17, face="bold"),
        panel.border = element_rect(color = "black", fill=NA ),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white" ),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        #legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=15, face = 'bold'),
        plot.margin = margin(1, 1, 1, 1, "cm"))
pdf ("Primary_TF_Regenrich_Score_barchart_August.pdf", width=15.5, height=6)
print(q)
dev.off()

#### Create Bar chart of individual gene groups ###

###Read in full table of Regenrich Scores
 reg_scores <- read.xlsx('Combined_Regenrich_TF_Table_edited_2.xlsx', sheetIndex = 1)



### PN Genes ##€#

### select PN column
prim_TF_PN <-  reg_scores[,c(1,2)]

#### sort by score ####
prim_TF_PN <- prim_TF_PN[order(desc(prim_TF_PN[,2])),]


###rename columns
colnames(prim_TF_PN)[1] <- "TF"
colnames(prim_TF_PN)[2] <- "Score"

#### set TFs as Factor by score ###

prim_TF_PN$TF <-factor(prim_TF_PN$TF,
                                  levels = prim_TF_PN$TF )


###select top 20 ###
prim_TF_PN_top <- prim_TF_PN[1:20,]

q<- ggplot(prim_TF_PN_top, aes(y=Score, x=TF)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
pdf ("Primary_PN_TF_Regenrich_Score_barchart.pdf", width=6, height=4)
print(q)
dev.off()

### ATP dependent HSP Genes ##€#

### select ATP_HSP column
prim_TF_ATP_HSP <-  reg_scores[,c(1,6)]

#### sort by score ####
prim_TF_ATP_HSP <- prim_TF_ATP_HSP[order(desc(prim_TF_ATP_HSP[,2])),]


###rename columns
colnames(prim_TF_ATP_HSP)[1] <- "TF"
colnames(prim_TF_ATP_HSP)[2] <- "Score"

#### set TFs as Factor by score ###

prim_TF_ATP_HSP$TF <-factor(prim_TF_ATP_HSP$TF,
                       levels = prim_TF_ATP_HSP$TF )


###select top 20 ###
prim_TF_ATP_HSP_top <- prim_TF_ATP_HSP[1:20,]

r<- ggplot(prim_TF_ATP_HSP_top, aes(y=Score, x=TF)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
pdf ("Primary_ATP_HSP_TF_Regenrich_Score_barchart.pdf", width=6, height=4)
print(r)
dev.off()

### sHSPs Genes ##€#

### select sHSP column
prim_TF_sHSP <-  reg_scores[,c(1,4)]

#### sort by score ####
prim_TF_sHSP <- prim_TF_sHSP[order(desc(prim_TF_sHSP[,2])),]


###rename columns
colnames(prim_TF_sHSP)[1] <- "TF"
colnames(prim_TF_sHSP)[2] <- "Score"

#### set TFs as Factor by score ###

prim_TF_sHSP$TF <-factor(prim_TF_sHSP$TF,
                            levels = prim_TF_sHSP$TF )


###select top 20 ###
prim_TF_sHSP_top <- prim_TF_sHSP[1:20,]

r<- ggplot(prim_TF_sHSP_top, aes(y=Score, x=TF)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
pdf ("Primary_sHSP_TF_Regenrich_Score_barchart.pdf", width=6, height=4)
print(r)
dev.off()

### proteasomes Genes ##€#

### select proteasome column
prim_TF_proteasome <-  reg_scores[,c(1,8)]

#### sort by score ####
prim_TF_proteasome <- prim_TF_proteasome[order(desc(prim_TF_proteasome[,2])),]


###rename columns
colnames(prim_TF_proteasome)[1] <- "TF"
colnames(prim_TF_proteasome)[2] <- "Score"

#### set TFs as Factor by score ###

prim_TF_proteasome$TF <-factor(prim_TF_proteasome$TF,
                         levels = prim_TF_proteasome$TF )


###select top 20 ###
prim_TF_proteasome_top <- prim_TF_proteasome[1:20,]

r<- ggplot(prim_TF_proteasome_top, aes(y=Score, x=TF)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
pdf ("Primary_proteasome_TF_Regenrich_Score_barchart.pdf", width=6, height=4)
print(r)
dev.off()

### UVs Genes ##€#

### select UV column
prim_TF_UV <-  reg_scores[,c(1,10)]

#### sort by score ####
prim_TF_UV <- prim_TF_UV[order(desc(prim_TF_UV[,2])),]


###rename columns
colnames(prim_TF_UV)[1] <- "TF"
colnames(prim_TF_UV)[2] <- "Score"

#### set TFs as Factor by score ###

prim_TF_UV$TF <-factor(prim_TF_UV$TF,
                               levels = prim_TF_UV$TF )


###select top 20 ###
prim_TF_UV_top <- prim_TF_UV[1:20,]

r<- ggplot(prim_TF_UV_top, aes(y=Score, x=TF)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
pdf ("Primary_UV_TF_Regenrich_Score_barchart.pdf", width=6, height=4)
print(r)
dev.off()

### Add column to say if TF is identified by Enrichr

#### read in Enrichr gene List ###
Enrichr_genes <- read.xlsx('Enrichr_all_PN.xlsx', sheetIndex = 1)

reg_scores$Enrichr <- ifelse(reg_scores$gene_name %in% Enrichr_genes$TF, "Yes", "No")

#### scatterplot PN vs shSP score

plot <- ggplot(reg_scores, aes(x=prim_PN_score, y=prim_sHSP_score, color = Enrichr)) +
  geom_point(aes(color = Enrichr))+  
 # scale_colour_manual(values=c ("Higher" = "#cf3b3b","Lower" ="#3b3bcf")) +
  ggtitle(''  ) +
  theme(plot.title = element_text(size = 10) ,
        legend.title = element_text(size = 10) ,
        panel.border = element_rect(color = "black", fill=NA ),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white" )
        
        #panel.grid.major = element_line(color = "gray") 
  )+
  labs(color = "TF identified by Enrichr ", 
       x ="PN", y = "sHSp") 
# geom_abline(intercept = 0,
#             slope = 1) + 
# geom_smooth(method='lm')
pdf ("PN_sHSP_score_primary.pdf",width=8, height=6)
print (plot)
dev.off()



### top ranked scatter plots ####

## read in table of top ranked TFs
reg_scores_ranked<- read.xlsx('Combined_Regenrich_TF_Table_edited_2.xlsx', sheetName = "top_table" )

reg_scores_ranked$Enrichr <- ifelse(reg_scores_ranked$gene_name %in% Enrichr_genes$TF, "Yes", "No")



reg_scores_top <-subset(reg_scores_ranked, reg_scores_ranked$min.rank <= 10)

#### scatterplot PN vs shSP score

plot <- ggplot(reg_scores_top, aes(x=prim_PN_score, y=prim_sHSP_score, color = Enrichr)) +
  geom_point(aes(color = Enrichr))+  
  # scale_colour_manual(values=c ("Higher" = "#cf3b3b","Lower" ="#3b3bcf")) +
  ggtitle(''  ) +
  theme(plot.title = element_text(size = 10) ,
        legend.title = element_text(size = 10) ,
        panel.border = element_rect(color = "black", fill=NA ),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white" )
        
        #panel.grid.major = element_line(color = "gray") 
  )+
  labs(color = "TF identified by Enrichr ", 
       x ="PN", y = "sHSp") 
# geom_abline(intercept = 0,
#             slope = 1) + 
# geom_smooth(method='lm')
pdf ("PN_sHSP_score_primary_top.pdf",width=8, height=6)
print (plot)
dev.off()

#### scatterplot PN vs shSP rank

plot <- ggplot(reg_scores_top, aes(x=prim_PN_rank, y=prim_sHSP_rank, color = Enrichr)) +
  geom_point(aes(color = Enrichr))+  
  # scale_colour_manual(values=c ("Higher" = "#cf3b3b","Lower" ="#3b3bcf")) +
  ggtitle(''  ) +
  theme(plot.title = element_text(size = 10) ,
        legend.title = element_text(size = 10) ,
        panel.border = element_rect(color = "black", fill=NA ),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white" )
        
        #panel.grid.major = element_line(color = "gray") 
  )+
  labs(color = "TF identified by Enrichr ", 
       x ="PN", y = "sHSp") 
# geom_abline(intercept = 0,
#             slope = 1) + 
# geom_smooth(method='lm')
pdf ("PN_sHSP_rank_primary_top.pdf",width=8, height=6)
print (plot)
dev.off()


####Create Matrix of Prim DEG targets of TFs identified by Regenrich ####
#### Create matrix showing which genes are regulated by each TF

##Read in list of significant PN Genes
Sig_PN_Gene_Table <- read.xlsx('TCGA_primary_PN_DSeq_Expression_Difference_final.xlsx', sheetName = 'DEG'  )

Sig_PN_Genes <-Sig_PN_Gene_Table$Gene
#Read in Enrichr TF Targets
TF_Targets  <- read.xlsx('Prim_DEG_Enrichr_TFs.xlsx', sheetIndex = 1)

### 

####amalgamate list of targets of TFs
TF_Targets$Targets <-'a'
### Concatenate rows with same marker

TF_Targets <- TF_Targets %>% 
  group_by(TF_Targets$Term) %>% 
  mutate(Targets = paste0(Genes, collapse = ";")) 

###remove extra columns ###

TF_Targets <-as.data.frame(TF_Targets[,c(1:10,82)])

###save methylation lists
write.xlsx(TF_Targets ,file = "Primary_DEG_Enrichr_TFs.xlsx" )


###remove duplicates ###

TF_Targets_no_dups <- TF_Targets[!duplicated(TF_Targets$Term), ]

TF_Targets_no_dups <- as.data.frame(TF_Targets_no_dups)
###Add marker as rowname

rownames(TF_Targets_no_dups) <- TF_Targets_no_dups$Term

###save Target lists
write.xlsx(TF_Targets_no_dups ,file = "Primary_DEG_Enrichr_TFs_no_dups.xlsx" )


##Read in list of Enrichr/Regenrich TFs

Regenrich_Enrichr_TFs_df <- 
read.xlsx('Combined_Regenrich_TF_Table_prim_cor_edited.xlsx', sheetName = 'Enrichr_Regenrich_TFs'  )

Regenrich_Enrichr_TFs<- Regenrich_Enrichr_TFs_df$TF
# ##Read in TF list
# 
# Enc_TI_df <- read.xlsx('Encode_Low_PN_TF_Table_eidted_2.xlsx', sheetIndex = 2)

x <- Regenrich_Enrichr_TFs
y <- Sig_PN_Genes

TF_target_df<- as.data.frame( matrix(ncol = length(x), nrow=length(y), dimnames = list(y,x)) )
TF_target_df$PN_Genes <- y
TF_target_df <- TF_target_df[moveme(names(TF_target_df), "PN_Genes first")]

#### create loop to write 'X'  if PN Gene is in the list of targets for the TF colname

## Read in csv of targets of C1C2TFs


t=0
for(t in 0:((length(x))-1))
  
{
  t <-t+1
  print(t)
  TF <- x[t]
  print(TF)
  TF_target_genes <- TF_Targets_no_dups [TF,12]
  TF_target_df[,TF]<- sapply(TF_target_df$PN_Genes, function(p)
  ifelse ((grepl(p,  TF_target_genes, fixed=TRUE) == TRUE) , 1,0))
}

### Delete PN Genes column

TF_target_df <- subset(TF_target_df[,-1])

TF_target_df_trans <-t(TF_target_df)

# Save rdata and write csv
save(TF_target_df_trans, file = 'Enrichr_TFs_matrix_primary.rdata')
write.xlsx  (TF_target_df_trans, file = 'Enrichr_TFs_matrix_primary.xlsx')


####make plotly interactive plot to identify genes
plot_ly(regulator_df_top, x=~met_sHSP_rank, y=~met_ATP_rank,
        text = ~gene_name
)

## plot of scores ###
plot <- ggplot(regulator_df_top, aes(x=met_sHSP_score, y=met_ATP_score, color = Enrichr)) +
  geom_point(aes(color = Enrichr))+  
  # scale_colour_manual(values=c ("Higher" = "#cf3b3b","Lower" ="#3b3bcf")) +
  ggtitle(''  ) +
  theme(plot.title = element_text(size = 10) ,
        legend.title = element_text(size = 10) ,
        panel.border = element_rect(color = "black", fill=NA ),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white" )
        
        #panel.grid.major = element_line(color = "gray") 
  )+
  labs(color = "TF identified by Enrichr ", 
       x ="sHSP", y = " ATP dependent HSP") +
  # scale_x_continuous(breaks = round(seq(0, max(regulator_df$met_sHSP_score), by = .2),1)) +
  #  scale_y_continuous(breaks = round(seq(0, max(regulator_df$met_ATP_score), by = .2),1))+
  xlim(0.89, 1.5) +
  ylim(0.89, 1.5 )+
  # 
  ggtitle("metastatic Regenrich Score of Regulators of\nATP dependent HSP genes vs sHSP Genes")+
  # geom_abline(intercept = 0, slope =1
  #, linetype, color, size
  #)
  pdf ("Prim_sHSP_ATPHSP_score.pdf",width=6, height=4)
print (plot)
dev.off()

####make plotly interactive plot to identify genes
plot_ly(regulator_df_top, x=~met_sHSP_score, y=~met_ATP_score,
        text = ~gene_name
)