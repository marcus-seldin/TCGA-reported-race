setwd('G:/My Drive/lab files/selma masri/race analysis TCGA')

#datasets from here: https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

#Next load the packages.  This you will run each time R is restarted
library(Rmisc)
library(survival)
library(ggplot2)
library(corrplot)
library(dplyr)
library(WGCNA)
library(MetBrewer)
library(pheatmap)
library(reshape2)
library(survminer)
library(RColorBrewer)
library(limma)

race_data = data.table::fread('./git analyses/panTCGA analysis datasets/Survival_SupplementalTable_S1_20171025_xena_sp')
cancer_annots = data.table::fread('./git analyses/panTCGA analysis datasets/TCGA_phenotype_denseDataOnlyDownload.tsv.gz')
colnames(race_data)
table(race_data$race)
#[Not Evaluated] 
#1160                                       160 
#[Unknown]          AMERICAN INDIAN OR ALASKA NATIVE 
#134                                        30 
#ASIAN                 BLACK OR AFRICAN AMERICAN 
#714                                      1020 
#NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER                                     WHITE 
#13                                      9360 

omit_cats = c('', '[Not Evaluated]', '[Unknown]')
sample_table = race_data[!race_data$race %in% omit_cats,]
sample_table$cancer_type = cancer_annots$`_primary_disease`[match(sample_table$sample, cancer_annots$sample)]

tt2 = sample_table

mut_data[1:20,1:10]

gene_expre = data.table::fread('./git analyses/panTCGA analysis datasets/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz')
gene_level = gene_expre[!is.na(gene_expre$sample),]
row.names(gene_level) = gene_level$sample
gene_level$sample=NULL

cnts_mat = as.data.frame(t(gene_level))
colnames(cnts_mat) = row.names(gene_level)
row.names(cnts_mat) = gsub('.', '-', row.names(cnts_mat), fixed = T)

cnts_mat = cnts_mat[row.names(cnts_mat) %in% tt2$sample,]
tt2 = tt2[tt2$sample %in% row.names(cnts_mat),]



### set up regression model for sex differences
table(tt2$treatment_outcome_first_course)
#FEMALE   MALE 
#5140   4556


cnts_mat[1:10,1:10]

driver_genes = c('KRAS', 'BRAF1', 'TP53', 'PIK3CA', 'PTEN')
##################################### #now divy by race and cancer type at the level of compare_set - first pan-cancer
perform_cor_call = function(race_call){
  cc1 = cnts_mat[,colnames(cnts_mat) %in% driver_genes]
  compare_set = tt2[tt2$race %in% race_call,]
  trt_set = compare_set[compare_set$OS==1,]
  cc1 = cc1[row.names(cc1) %in% trt_set$sample,]
cc1 = cc1[order(row.names(cc1)),]
trt_set = trt_set[order(trt_set$sample),]  

cors_1 = bicorAndPvalue(cc1, trt_set$OS.time, use = 'p')
  melted_cors = reshape2::melt(cors_1$bicor)
  head(melted_cors)
  melted_cors$Var2 = NULL
  colnames(melted_cors) = c('gene', 'bicor')
  melted_cors$pvalue = reshape2::melt(cors_1$p)$value
  melted_cors$race = paste0(race_call)
  return(melted_cors)
}
joined_data = as.data.frame(rbind(perform_cor_call('ASIAN'), perform_cor_call('WHITE'), perform_cor_call('BLACK OR AFRICAN AMERICAN'), perform_cor_call('AMERICAN INDIAN OR ALASKA NATIVE'), perform_cor_call('NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER')))
joined_data$bicor
cors_plot = dcast(joined_data, gene ~ race, value.var = 'bicor')
row.names(cors_plot) = cors_plot$gene
cors_plot$gene=NULL
cors_pval = dcast(joined_data, gene ~ race, value.var = 'pvalue')
row.names(cors_pval) = cors_pval$gene
cors_pval$gene=NULL
corrplot(as.matrix(cors_plot), p.mat = as.matrix(cors_pval), hc.order = T, 
         insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001), tl.col = "black")
