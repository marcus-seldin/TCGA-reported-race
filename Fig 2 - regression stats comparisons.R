setwd('G:/My Drive/lab files/selma masri/race analysis TCGA')

#datasets from here: https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

#Next load the packages.  This you will run each time R is restarted
library(Rmisc)
library(survival)
library(ggplot2)
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

gene_expre = data.table::fread('./git analyses/panTCGA analysis datasets/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')
gene_level = gene_expre[!duplicated(gene_expre$sample),]
row.names(gene_level) = gene_level$sample
gene_level$sample=NULL

cnts_mat = as.data.frame(t(gene_level))
colnames(cnts_mat) = row.names(gene_level)
row.names(cnts_mat) = gsub('.', '-', row.names(cnts_mat), fixed = T)

cnts_mat = cnts_mat[row.names(cnts_mat) %in% tt2$sample,]
tt2 = tt2[tt2$sample %in% row.names(cnts_mat),]


### set up regression model for sex differences
table(tt2$gender)
#FEMALE   MALE 
#5140   4556

female_sex = tt2$sample[tt2$gender=='FEMALE']

#############toggle here for category
cc1 = cnts_mat
cc1$dm = ifelse(row.names(cc1) %in% female_sex, 'Female', 'Male')
table(cc1$dm)
#Female   Male 
#5140   4556
design = model.matrix(~dm, data=cc1)
head(design)
table(cc1$dm)
dim(design)
new_cnts1 = as.data.frame(t(cc1[, !colnames(cc1)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)
row.names(fit)[1:10]

res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)
test1 = as.data.frame(cnts_mat[,colnames(cnts_mat) == 'KDM5D'])
row.names(test1) = row.names(cnts_mat)
test1$cond = ifelse(row.names(test1) %in% female_sex, 'Female', 'Male') 
res_table$race = 'combined'

full_de_table = res_table


##################################### #now divy by race and cancer type at the level of compare_set - first pan-cancer
perform_de_call = function(race_call){
cc1 = cnts_mat
compare_set = tt2[tt2$race %in% race_call,]
cc1 = cc1[row.names(cc1) %in% compare_set$sample,]
cc1$dm = ifelse(row.names(cc1) %in% female_sex, 'Female', 'Male')
table(cc1$dm)
#Female   Male 
#5140   4556
design = model.matrix(~dm, data=cc1)
head(design)
table(cc1$dm)
dim(design)
new_cnts1 = as.data.frame(t(cc1[, !colnames(cc1)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)
row.names(fit)[1:10]

res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)
test1 = as.data.frame(cnts_mat[,colnames(cnts_mat) == 'KDM5D'])
row.names(test1) = row.names(cnts_mat)
test1$cond = ifelse(row.names(test1) %in% female_sex, 'Female', 'Male') 
res_table$race = paste0(race_call)
res1 = res_table
  return(data.frame(res1))
}

asian_de = perform_de_call('ASIAN')
white_res = perform_de_call('WHITE')
black_de =  perform_de_call('BLACK OR AFRICAN AMERICAN')
nativeameri_de = perform_de_call('AMERICAN INDIAN OR ALASKA NATIVE')
pacific_isl_de = perform_de_call('NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER')

full_de_table = as.data.frame(rbind(full_de_table, asian_de, white_res, black_de, nativeameri_de, pacific_isl_de))
table(full_de_table$race)

ff1 = full_de_table[full_de_table$adj.P.Val<0.01,]
ff1$adp_Co = ifelse(ff1$adj.P.Val<0.001, 'adj_P < 0.001', 'adj_P < 0.01')
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-5, 'adj_P < 1e-5', paste0(ff1$adp_Co))
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-8, 'adj_P < 1e-8', paste0(ff1$adp_Co))
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-20, 'adj_P < 1e-20', paste0(ff1$adp_Co))
ff1 = na.omit(ff1)
ff1$adp_Co = factor(ff1$adp_Co, levels = c('adj_P < 0.01', 'adj_P < 0.001',  'adj_P < 1e-5', 'adj_P < 1e-8', 'adj_P < 1e-20' ))
pdf(file = 'sex differences TCGA by race - pan_cancer.pdf')
ggplot(ff1, aes(x=adp_Co, fill=race)) + geom_bar() + theme_minimal() + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(ff1$race))))) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('number if significant genes DE Male vs Female') + ylab('Number of DEGs M vs F - Limma') + xlab('')
dev.off()

############################################################################# do for cancer types
cc1 = cnts_mat
compare_set = tt2[tt2$cancer_type =='diffuse large B-cell lymphoma',]
cc1 = cc1[row.names(cc1) %in% compare_set$sample,]
cc1$dm = ifelse(row.names(cc1) %in% female_sex, 'Female', 'Male')
table(cc1$dm)
#Female   Male 
#5140   4556
design = model.matrix(~dm, data=cc1)
head(design)
table(cc1$dm)
dim(design)
new_cnts1 = as.data.frame(t(cc1[, !colnames(cc1)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)
row.names(fit)[1:10]

res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)
 
res_table$race = 'combined'

full_de_table = res_table

perform_de_call = function(race_call, cancer_type1){
  cc1 = cnts_mat
  compare_set = tt2[tt2$race %in% race_call,]
  compare_set = compare_set[compare_set$cancer_type %in% cancer_type1,]
  cc1 = cc1[row.names(cc1) %in% compare_set$sample,]
  cc1$dm = ifelse(row.names(cc1) %in% female_sex, 'Female', 'Male')
  table(cc1$dm)
  #Female   Male 
  #5140   4556
  design = model.matrix(~dm, data=cc1)
  head(design)
  table(cc1$dm)
  dim(design)
  new_cnts1 = as.data.frame(t(cc1[, !colnames(cc1)=='dm']))
  fit = lmFit(new_cnts1, design)
  fit = eBayes(fit)
  row.names(fit)[1:10]
  
  res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                       sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
  head(res_table)
  test1 = as.data.frame(cnts_mat[,colnames(cnts_mat) == 'KDM5D'])
  row.names(test1) = row.names(cnts_mat)
  test1$cond = ifelse(row.names(test1) %in% female_sex, 'Female', 'Male') 
  res_table$race = paste0(race_call)
  res1 = res_table
  return(data.frame(res1))
}

asian_de = perform_de_call('ASIAN', 'diffuse large B-cell lymphoma')
white_res = perform_de_call('WHITE', 'diffuse large B-cell lymphoma')

##black not possible
black_de =  perform_de_call('BLACK OR AFRICAN AMERICAN', 'diffuse large B-cell lymphoma')

##native american not possible
nativeameri_de = perform_de_call('AMERICAN INDIAN OR ALASKA NATIVE', 'diffuse large B-cell lymphoma')

##pacific_islander not possible
pacific_isl_de = perform_de_call('NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'diffuse large B-cell lymphoma')

full_de_table = as.data.frame(rbind(full_de_table, asian_de, white_res 
                                    #black_de, nativeameri_de, pacific_isl_de
                                    ))


ff1 = full_de_table[full_de_table$adj.P.Val<0.01,]
ff1$adp_Co = ifelse(ff1$adj.P.Val<0.001, 'adj_P < 0.001', 'adj_P < 0.01')
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-5, 'adj_P < 1e-5', paste0(ff1$adp_Co))
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-8, 'adj_P < 1e-8', paste0(ff1$adp_Co))
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-20, 'adj_P < 1e-20', paste0(ff1$adp_Co))
ff1 = na.omit(ff1)
ff1$adp_Co = factor(ff1$adp_Co, levels = c('adj_P < 0.01', 'adj_P < 0.001',  'adj_P < 1e-5', 'adj_P < 1e-8', 'adj_P < 1e-20' ))
pdf(file = 'sex differences TCGA by race - diffuse large B-cell lymphoma.pdf')
ggplot(ff1, aes(x=adp_Co, fill=race)) + geom_bar() + theme_minimal() + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(ff1$race))))) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('number if significant genes DE Male vs Female - diffuse large B-cell lymphoma') + ylab('Number of DEGs M vs F - Limma') + xlab('')
dev.off()


################################################
table(tt2$cancer_type)
cc1 = cnts_mat
compare_set = tt2[tt2$cancer_type =='colon adenocarcinoma',]
cc1 = cc1[row.names(cc1) %in% compare_set$sample,]
cc1$dm = ifelse(row.names(cc1) %in% female_sex, 'Female', 'Male')
table(cc1$dm)
#Female   Male 
#5140   4556
design = model.matrix(~dm, data=cc1)
head(design)
table(cc1$dm)
dim(design)
new_cnts1 = as.data.frame(t(cc1[, !colnames(cc1)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)
row.names(fit)[1:10]

res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)

res_table$race = 'combined'

full_de_table = res_table


asian_de = perform_de_call('ASIAN', 'colon adenocarcinoma')
white_res = perform_de_call('WHITE', 'colon adenocarcinoma')
black_de =  perform_de_call('BLACK OR AFRICAN AMERICAN', 'colon adenocarcinoma')

##native american not possible
nativeameri_de = perform_de_call('AMERICAN INDIAN OR ALASKA NATIVE', 'colon adenocarcinoma')

##pacific_islander not possible
pacific_isl_de = perform_de_call('NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'colon adenocarcinoma')

full_de_table = as.data.frame(rbind(full_de_table, asian_de, white_res, black_de
                                    #, nativeameri_de, pacific_isl_de
))


ff1 = full_de_table[full_de_table$adj.P.Val<0.01,]
ff1$adp_Co = ifelse(ff1$adj.P.Val<0.001, 'adj_P < 0.001', 'adj_P < 0.01')
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-5, 'adj_P < 1e-5', paste0(ff1$adp_Co))
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-8, 'adj_P < 1e-8', paste0(ff1$adp_Co))
ff1$adp_Co = ifelse(ff1$adj.P.Val<1e-20, 'adj_P < 1e-20', paste0(ff1$adp_Co))
ff1 = na.omit(ff1)
ff1$adp_Co = factor(ff1$adp_Co, levels = c('adj_P < 0.01', 'adj_P < 0.001',  'adj_P < 1e-5', 'adj_P < 1e-8', 'adj_P < 1e-20' ))
pdf(file = 'sex differences TCGA by race - colon adenocarcinoma.pdf')
ggplot(ff1, aes(x=adp_Co, fill=race)) + geom_bar() + theme_minimal() + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(ff1$race))))) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('number if significant genes DE Male vs Female - colon adenocarcinoma') + ylab('Number of DEGs M vs F - Limma') + xlab('')
dev.off()
