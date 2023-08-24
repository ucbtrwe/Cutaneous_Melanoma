
#### This code creates a matrix showing which categories of Go Function each regulator is allocated to.

library(ggplot2)
 library(dplyr)
library(data.table)
library('xlsx')
# create moveme function to move columns
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

#load Rdata ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##### SET WORKING DIRECTORY , use getwd to designate new variable'current_dir' to use in sub directory names
setwd("~/iCloud_Drive_Archive/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()

#####Read in table of regulators and functions ####
Regulator_function_df <- read.csv("Regenrich_Regulators_functional_annotations.csv")

Regulator_function_df <-as.data.frame(Regulator_function_df)
####Create table of number of times each function appears

function_frequency_table <-  table(Regulator_function_df$term.description )

###Save table 

write.xlsx(function_frequency_table ,file = "Frequency_of_function_table.xlsx" )

###Concatenate rows of functions table ####


Regulator_function_df$Functions <-'a'

### Concatenate rows with same marker

Regulator_function_df_concatenated <-Regulator_function_df %>% 
  group_by(X.node) %>% 
  mutate(Functions = paste0(term.description, collapse = ";")) 



###remove duplicates ###

Regulator_function_df_no_dups <- as.data.frame(Regulator_function_df_concatenated[!duplicated(Regulator_function_df_concatenated$X.node), ])

demeth_Regulator_function_df_no_dups<- subset(Regulator_function_df_no_dups,Regulator_function_df_no_dups$Functions %like% "demeth")
write.csv(demeth_Regulator_function_df_no_dups,file ="demeth_regulators.csv")

###Add marker as rowname

rownames(Regulator_function_df_no_dups) <- Regulator_function_df_no_dups$X.node


####Read in table of GO functions
Functions_df <- read.xlsx ("Regulator Functions Table.xlsx",sheetIndex = 1 )

Functions_Categories <- names(Functions_df)
####Identify which Regulators functions include items from list

###https://stackoverflow.com/questions/44332906/finding-matches-between-each-element-of-a-list-of-strings-and-each-string-of-a-v
# df <- data.frame(V = c('Anameone','Bnametwo','Cnamethree'),
#                  L = I(list(c('name','asd'),c('dfg'),c('hey','C','hi'))))
# 
# sapply(1:nrow(df), function(x) any(sapply(df$L[[x]], function(y) grepl(y, df$V[x]))))
# 
 ####Matrix of Functionss #####
 #### Create matrix showing which genes are associated with each methylation marker
 ####Create Matrix of Marker genes #####
x<- rownames(Regulator_function_df_no_dups) 
  y <- Functions_Categories 

Regulator_GO_functions_df<- as.data.frame( matrix(ncol = length(y), nrow=length(x), dimnames = list(x,y)) )
Regulator_GO_functions_df$Regulators <- y
Regulator_GO_functions_df <- Regulator_GO_functions_df[moveme(names(Regulator_GO_functions_df), "Regulators first")]

 #### create loop to write '1'  if Met A Gene is in the list of targets for the Marker colname

 t=0
 for(t in 0:((length(x))-1))

 {
   t <-t+1
   print(t)
   reg <- x[t]
   print(reg)
   
   Regulator_GO_function_list <-Regulator_function_df_no_dups[reg,"Functions"]
     v=0   
     for(v in 0:((length(y))-1)) 
      { v <-v+1
     print(v)
   Functions_Category <- y[v]
   print(Functions_Category)
  
   members_of_func_cat <- Functions_df[,Functions_Category]
   members_of_func_cat<-members_of_func_cat[!is.na(members_of_func_cat)]
 
    members_of_func_cat <- paste(members_of_func_cat, collapse = "|")
    print(members_of_func_cat)
     Regulator_GO_functions_df[reg,Functions_Category]<-ifelse (( (grepl(members_of_func_cat,  Regulator_GO_function_list) == TRUE)), 1,0)
     }
     }
 
 ### order by regenrich minimum rank ###
Regenrich_rank_df <- read.xlsx("Regenrich_tables_met_and_prim_with_min.xlsx", sheetIndex = 1)
ranked_regulators <- Regenrich_rank_df$gene_name

Sorted_Regulator_GO_functions_df <- Regulator_GO_functions_df[match(ranked_regulators,row.names(Regulator_GO_functions_df)),]   


 # Save rdata and write csv
 save(Sorted_Regulator_GO_functions_df, file = 'GO_Functions_of_Regulators_inc_proteins.rdata')
 write.xlsx  (Sorted_Regulator_GO_functions_df, file = 'GO_Functions_of_Regulators_inc_proteins.xlsx')

 
 ### identify which regulators are differentially expressed in metastatic groups ### #

 # read csv of SKCM FKPM to sub folder
 FKPM_data<- loadRData( file =  'SKCM_FKPM.Rdata')
 
 # reduce table to include only genes from Regulators
 Reg_Gene_Expression <-  subset(FKPM_data,hgnc_symbol %in% ranked_regulators )
 
 #### Metastatic ####
 ###Select metastatic samples
 met_cols <- colnames(FKPM_data)[grep("06",colnames(FKPM_data))]
 
 Reg_Gene_Expression_met<- Reg_Gene_Expression[,c('hgnc_symbol',met_cols)]
 
 # remove duplicate rows by hgnc symbol
 Reg_Gene_Expression_met = Reg_Gene_Expression_met[!duplicated(Reg_Gene_Expression_met$hgnc_symbol),]
 
 
 # convert hgnc_symbol column to row name
 
 Reg_Gene_Expression_met <- data.frame(Reg_Gene_Expression_met, row.names = 1)
 
 ###save expression table ####
 save(Reg_Gene_Expression_met, file = "Reg_Gene_Expression_in_SKCM_met.Rdata")
 write.xlsx(Reg_Gene_Expression_met, file = "Reg_Gene_Expression_in_SKCM_met.xlsx")
 
 ###Identify TFs that are differentially expressed between sample groups ####
 ### read in list of sample group allocations ####
 sample_groups <- read.xlsx('Final_Expression_primary_and_metastatic_samples_with_groups.xlsx',sheetIndex = 4  )
 

 ### transpose expression table ###
 Reg_Gene_Expression_met_t <-  as.data.frame(t(Reg_Gene_Expression_met))
 
 ###alter rownames to match format of sample group allocation
 Reg_Gene_Expression_met_t$Sample<- row.names( Reg_Gene_Expression_met_t)
 Reg_Gene_Expression_met_t$Sample<-  substr(Reg_Gene_Expression_met_t$Sample,start=1,stop=15)
 ###replace dot by dash 
 Reg_Gene_Expression_met_t$Sample<-gsub(".", "-", Reg_Gene_Expression_met_t$Sample, fixed = TRUE)
 
 
 #### merge with  sample group allocation ####
 Reg_Gene_Expression_met_t<- merge(sample_groups,  Reg_Gene_Expression_met_t, by.x = 'Sample', by.y = 'Sample')
 
 ### calculate expression differences ####
 Expression_difference_df <- data.frame (Gene= character(0),
                                         Foldchange_A_over_B=numeric(0),P.value=numeric(0), mean_A =numeric(0), mean_B =numeric(0))
 
 
 
 Genes<- names(Reg_Gene_Expression_met_t[,3:591])
 n =0
 
 while(n<length(Genes)) try({
   n =n+1
   print(n)
   Gene <- Genes[n]
   print(Gene)
   A <-subset ( Reg_Gene_Expression_met_t[,Gene],(Reg_Gene_Expression_met_t$Sample_Group == "Metastatic (A)"))
   B<-subset ( Reg_Gene_Expression_met_t[,Gene],(Reg_Gene_Expression_met_t$Sample_Group == "Metastatic (B)"))
   # =Foldchange_A_over_B <- (mean(as.numeric(A))/mean(as.numeric(B)))
   mean_A<-mean(as.numeric(A))
   mean_A<-format(mean_A, round(mean_A, 3),scientific = FALSE)
   mean_A<- as.numeric(mean_A)
   mean_B <-  mean(as.numeric(B))
   mean_B <-  format(mean_B, round(mean_A, 3),scientific = FALSE)
   mean_B<- as.numeric(mean_B)
   Foldchange_A_over_B <- (mean_A/mean_B)
   results = t.test((as.numeric(A)),(as.numeric(B)))
   P.value<-results$p.value
   formatC(P.value, format = "e", digits = 3)
   P.value<-as.numeric(P.value)
   temp_df <- data.frame (Gene, P.value,Foldchange_A_over_B,mean_A,mean_B)
   # temp_df$TF <- TF
   # temp_df$P.value <- results$p.value
   # temp_df$Foldchange_A_over_B <- Foldchange
   Expression_difference_df<- rbind(Expression_difference_df, temp_df)
 })
 
 
 ##Add columnn listing higher meaen and adjusted p value
 Expression_difference_df$max_mean <- pmax(Expression_difference_df$mean_A,
                                           Expression_difference_df$mean_B)
 Expression_difference_df$P.value <- formatC(Expression_difference_df$P.value, format = "e", digits = 3)
 
 Expression_difference_df$P.value <- as.numeric(Expression_difference_df$P.value)
 
 Expression_difference_df$p_adjust <- p.adjust(Expression_difference_df$P.value, "hochberg")
 
 
 Expression_difference_df$p_adjust <- formatC(Expression_difference_df$p_adjust, format = "e", digits = 3)
 
 Expression_difference_df$p_adjust <- as.numeric(Expression_difference_df$p_adjust)
 
 ###Add expression column 
 Expression_difference_df$T_test_expression <- 
   ifelse(Expression_difference_df$p_adjust > 0.1 ,"No Significant Difference",
          ifelse(Expression_difference_df$Foldchange_A_over_B <= 0.8 , "Lower",
                 ifelse(Expression_difference_df$Foldchange_A_over_B >= 1.25 , "Higher", "No Significant Difference"))     )
 
 ####Save met TF expression table
 save(Expression_difference_df, file = "Met_Samples_Regulators_T_Test_Expression.Rdata")
 write.xlsx(Expression_difference_df, file = "Met_Samples_Regulators_T_Test_Expression.xlsx")
 
 ###merge T Test expression table with DESEq2 Expression list
 ###Read in met DESeq expression table
 DSEq_Exp <- read.csv('TCGA_metastatic_SKCM_raw_Gene_RNA_seq-reads_Dan_dds_DGE_results_with_gene_name.csv')
 
 ##
 T_Test_DESeq_exp <- merge ( Expression_difference_df,DSEq_Exp, by.x = 'Gene', by.y = 'hgnc_symbol')
 
 
 ###add final expression column
 T_Test_DESeq_exp$Final_Expression <- 
   ifelse( T_Test_DESeq_exp$T_test_expression == T_Test_DESeq_exp$DeSeq_Met_Expression,
           T_Test_DESeq_exp$T_test_expression, "No Significant Difference" )
 
#### sort by order of Regenrich Rank
 
 Sorted_T_Test_DESeq_exp <- T_Test_DESeq_exp[match(ranked_regulators,T_Test_DESeq_exp$Gene),]   
 
 ###add column with group with higher expression
 Sorted_T_Test_DESeq_exp$Higher_Expression_Metastatic <- 
   ifelse( Sorted_T_Test_DESeq_exp$Final_Expression == "Lower","Metastatic (B)",
           ifelse( Sorted_T_Test_DESeq_exp$Final_Expression == "Higher","Metastatic (A)", "No Significant Difference" ))
 
 
 
 
 save(Sorted_T_Test_DESeq_exp , file = "Regenrich_Regulators_expresison_in_metastatic_samples.Rdata")
 write.xlsx(Sorted_T_Test_DESeq_exp , file = "Regenrich_Regulators_expresison_in_metastatic_samples.xlsx")
 
 #### Primary ####
 ###Select Primary samples
 Prim_cols <- colnames(FKPM_data)[grep("01",colnames(FKPM_data))]
 
 Reg_Gene_Expression_Prim<- Reg_Gene_Expression[,c('hgnc_symbol',Prim_cols)]
 
 # remove duplicate rows by hgnc symbol
 Reg_Gene_Expression_Prim = Reg_Gene_Expression_Prim[!duplicated(Reg_Gene_Expression_Prim$hgnc_symbol),]
 
 
 # convert hgnc_symbol column to row name
 
 Reg_Gene_Expression_Prim <- data.frame(Reg_Gene_Expression_Prim, row.names = 1)
 
 ###save expression table ####
 save(Reg_Gene_Expression_Prim, file = "Reg_Gene_Expression_in_SKCM_Prim.Rdata")
 write.xlsx(Reg_Gene_Expression_Prim, file = "Reg_Gene_Expression_in_SKCM_Prim.xlsx")
 
 ###Identify TFs that are differentially expressed between sample groups ####

 
 ### transpose expression table ###
 Reg_Gene_Expression_Prim_t <-  as.data.frame(t(Reg_Gene_Expression_Prim))
 
 ###alter rownames to match format of sample group allocation
 Reg_Gene_Expression_Prim_t$Sample<- row.names( Reg_Gene_Expression_Prim_t)
 Reg_Gene_Expression_Prim_t$Sample<-  substr(Reg_Gene_Expression_Prim_t$Sample,start=1,stop=15)
 ###replace dot by dash 
 Reg_Gene_Expression_Prim_t$Sample<-gsub(".", "-", Reg_Gene_Expression_Prim_t$Sample, fixed = TRUE)
 
 ### read in list of sample group allocations ####
 sample_groups_prim <- read.xlsx('Final_Expression_primary_and_metastatic_samples_with_groups.xlsx',sheetIndex = 3  )
 
 #### merge with  sample group allocation ####
 Reg_Gene_Expression_Prim_t<- merge(sample_groups,  Reg_Gene_Expression_Prim_t, by.x = 'Sample', by.y = 'Sample')
 
 ### calculate expression differences ####
 Expression_difference_df_prim <- data.frame (Gene= character(0),
                                         Foldchange_A_over_B=numeric(0),P.value=numeric(0), mean_A =numeric(0), mean_B =numeric(0))
 
 
 
 Genes<- names(Reg_Gene_Expression_Prim_t[,3:591])
 n =0
 
 while(n<length(Genes)) try({
   n =n+1
   print(n)
   Gene <- Genes[n]
   print(Gene)
   A <-subset ( Reg_Gene_Expression_Prim_t[,Gene],(Reg_Gene_Expression_Prim_t$Sample_Group == "Primary (A)"))
   B<-subset ( Reg_Gene_Expression_Prim_t[,Gene],(Reg_Gene_Expression_Prim_t$Sample_Group == "Primary (B)"))
   # =Foldchange_A_over_B <- (mean(as.numeric(A))/mean(as.numeric(B)))
   mean_A<-mean(as.numeric(A))
   mean_A<-format(mean_A, round(mean_A, 3),scientific = FALSE)
   mean_A<- as.numeric(mean_A)
   mean_B <-  mean(as.numeric(B))
   mean_B <-  format(mean_B, round(mean_A, 3),scientific = FALSE)
   mean_B<- as.numeric(mean_B)
   Foldchange_A_over_B <- (mean_A/mean_B)
   results = t.test((as.numeric(A)),(as.numeric(B)))
   P.value<-results$p.value
   formatC(P.value, format = "e", digits = 3)
   P.value<-as.numeric(P.value)
   temp_df <- data.frame (Gene, P.value,Foldchange_A_over_B,mean_A,mean_B)
   # temp_df$TF <- TF
   # temp_df$P.value <- results$p.value
   # temp_df$Foldchange_A_over_B <- Foldchange
   Expression_difference_df_prim<- rbind(Expression_difference_df_prim, temp_df)
 })
 
 
 ##Add columnn listing higher meaen and adjusted p value
 Expression_difference_df_prim$max_mean <- pmax(Expression_difference_df_prim$mean_A,
                                           Expression_difference_df_prim$mean_B)
 Expression_difference_df_prim$P.value <- formatC(Expression_difference_df_prim$P.value, format = "e", digits = 3)
 
 Expression_difference_df_prim$P.value <- as.numeric(Expression_difference_df_prim$P.value)
 
 Expression_difference_df_prim$p_adjust <- p.adjust(Expression_difference_df_prim$P.value, "hochberg")
 
 
 Expression_difference_df_prim$p_adjust <- formatC(Expression_difference_df_prim$p_adjust, format = "e", digits = 3)
 
 Expression_difference_df_prim$p_adjust <- as.numeric(Expression_difference_df_prim$p_adjust)
 
 ###Add expression column 
 Expression_difference_df_prim$T_test_expression <- 
   ifelse(Expression_difference_df_prim$p_adjust > 0.1 ,"No Significant Difference",
          ifelse(Expression_difference_df_prim$Foldchange_A_over_B <= 0.8 , "Lower",
                 ifelse(Expression_difference_df_prim$Foldchange_A_over_B >= 1.25 , "Higher", "No Significant Difference"))     )
 
 ####Save Prim TF expression table
 save(Expression_difference_df_prim, file = "Prim_Samples_Regulators_T_Test_Expression.Rdata")
 write.xlsx(Expression_difference_df_prim, file = "Prim_Samples_Regulators_T_Test_Expression.xlsx")
 
 ###merge T Test expression table with DESEq2 Expression list
 ###Read in Prim DESeq expression table
 DSEq_Exp_prim <- read.csv('TCGA_primary_SKCM_raw_Whole_Genome_RNA_seq-reads_Dan_dds_DGE_results_with_gene_name_sig.csv')
 
 ##
 T_Test_DESeq_exp_prim <- merge ( Expression_difference_df_prim, DSEq_Exp_prim, by.x = 'Gene', by.y = 'hgnc_symbol')
 
 
 ###add final expression column
 T_Test_DESeq_exp_prim$Final_Expression <- 
   ifelse( T_Test_DESeq_exp_prim$T_test_expression == T_Test_DESeq_exp_prim$DSeq_Expression,
           T_Test_DESeq_exp_prim$T_test_expression, "No Significant Difference" )
 
 #### sort by order of Regenrich Rank
 
 Sorted_T_Test_DESeq_exp_prim <- T_Test_DESeq_exp_prim[match(ranked_regulators,T_Test_DESeq_exp_prim$Gene),]   
 
 ###add column with group with higher expression
 Sorted_T_Test_DESeq_exp_prim$Higher_Expression_Primary <- 
   ifelse( Sorted_T_Test_DESeq_exp_prim$Final_Expression == "Lower","Primary (B)",
           ifelse( Sorted_T_Test_DESeq_exp_prim$Final_Expression == "Higher","Primary (A)", "No Significant Difference" ))
 
 
 
 
 save(Sorted_T_Test_DESeq_exp_prim , file = "Regenrich_Regulators_expresison_in_Primary_samples.Rdata")
 write.xlsx(Sorted_T_Test_DESeq_exp_prim , file = "Regenrich_Regulators_expresison_in_Primary_samples.xlsx")
 
 
