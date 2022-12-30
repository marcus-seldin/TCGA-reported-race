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

race_data = data.table::fread('./panTCGA phenotypes/Survival_SupplementalTable_S1_20171025_xena_sp')
cancer_annots = data.table::fread('./panTCGA phenotypes/TCGA_phenotype_denseDataOnlyDownload.tsv.gz')
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

binned_races= tt2 %>%
  dplyr::group_by(cancer_type, race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
pdf(file = 'Race by Cancer Type PanTCGA.pdf')
ggplot(binned_races, aes(x=cancer_type, y=log10(n), fill=race))+
  geom_col(position = 'dodge') +
  geom_text(aes(label = n), vjust = 0, hjust=0.5) + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(tt2$race))))) +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('individuals by race - panTCGA') + theme(legend.position = "none")
dev.off()

pdf(file = 'Race by Cancer Type PanTCGA - WITH LEGEND.pdf')
ggplot(binned_races, aes(x=cancer_type, y=log10(n), fill=race))+
  geom_col(position = 'dodge') +
  geom_text(aes(label = n), vjust = 0, hjust=0.5) + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(tt2$race))))) +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('individuals by race - panTCGA') 
dev.off()


#Pie chart of pan-TCGA race
binned_races= tt2 %>%
  dplyr::group_by(race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=25, face="bold")
  )
pdf(file = paste0('Pan-TCGA Race Distribution.pdf'))
ggplot(binned_races, aes(x = "", y = freq, fill =race)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=15, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=rev(met.brewer('Moreau', length(unique(tt2$race))))) + ggtitle(paste0('Pan-TCGA Race Distribution, n = ', length(row.names(tt2)))) + theme(plot.title = element_text(size=15)) +
  coord_polar(theta = "y") 
dev.off()





#Tumor vs surrounding tissue binning
tt2$tumor_vs_normal_number = gsub("(.*-\\s*(.*$))", "\\2", tt2$sample)
tt2$tumor_vs_normal_number = as.numeric(as.character(tt2$tumor_vs_normal_number))
table(tt2$tumor_vs_normal_number)
#01   02   03   05   06   07   11 
#9210   54  198   10  382    1 1282
tt2$tumor_vs_normal = ifelse(tt2$tumor_vs_normal_number<10, 'tumor', 'normal_tissue')
table(tt2$tumor_vs_normal)

tt3 = tt2[tt2$tumor_vs_normal =='tumor',]

binned_races= tt3 %>%
  dplyr::group_by(race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))


pdf(file = paste0('Tumor Sample Race Distribution.pdf'))
ggplot(binned_races, aes(x = "", y = freq, fill =race)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=15, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=rev(met.brewer('Moreau', length(unique(tt3$race))))) + ggtitle(paste0('Tumor Sample Race Distribution, n = ', length(row.names(tt3)))) + theme(plot.title = element_text(size=15)) +
  coord_polar(theta = "y") 
dev.off()

binned_races= tt3 %>%
  dplyr::group_by(cancer_type, race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
pdf(file = 'Tumor Sample by Cancer Type PanTCGA.pdf')
ggplot(binned_races, aes(x=cancer_type, y=log10(n), fill=race))+
  geom_col(position = 'dodge') +
  geom_text(aes(label = n), vjust = 0, hjust=0.5) + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(tt2$race))))) +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('individuals by race with tumor samples - panTCGA') + theme(legend.position = "none")
dev.off()


tt3 = tt2[tt2$tumor_vs_normal =='normal_tissue',]

binned_races= tt3 %>%
  dplyr::group_by(race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))


pdf(file = paste0('Normal Tissue Sample Race Distribution.pdf'))
ggplot(binned_races, aes(x = "", y = freq, fill =race)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=15, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=rev(met.brewer('Moreau', length(unique(tt3$race))))) + ggtitle(paste0('Normal Tissue Race Distribution, n = ', length(row.names(tt3)))) + theme(plot.title = element_text(size=15)) +
  coord_polar(theta = "y") 
dev.off()

binned_races= tt3 %>%
  dplyr::group_by(cancer_type, race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
pdf(file = 'Normal Tissue Sample by Cancer Type PanTCGA.pdf')
ggplot(binned_races, aes(x=cancer_type, y=log10(n), fill=race))+
  geom_col(position = 'dodge') +
  geom_text(aes(label = n), vjust = 0, hjust=0.5) + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(tt2$race))))) +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('individuals by race with normal tissue samples - panTCGA') + theme(legend.position = "none")
dev.off()



##############################################
#similar analysis but by staging
#clinical stage
table(tt2$clinical_stage)
unkown_key = c('', '[Discrepancy]')
tt2$staged_tumor = ifelse(tt2$clinical_stage %in% unkown_key, 'Not_Available', 'Available')
table(tt2$staged_tumor)
#Available Not_Available 
#2665          8472)

tt3 = tt2[tt2$staged_tumor =='Available',]

binned_races= tt3 %>%
  dplyr::group_by(race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))


pdf(file = paste0('Clinical Stage Availability Race Distribution.pdf'))
ggplot(binned_races, aes(x = "", y = freq, fill =race)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=15, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=rev(met.brewer('Moreau', length(unique(tt3$race))))) + ggtitle(paste0('Clinical Stage Availability Race Distribution, n = ', length(row.names(tt3)))) + theme(plot.title = element_text(size=15)) +
  coord_polar(theta = "y") 
dev.off()

binned_races= tt3 %>%
  dplyr::group_by(cancer_type, race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
pdf(file = 'Clinical Stage Availability by Cancer Type PanTCGA.pdf')
ggplot(binned_races, aes(x=cancer_type, y=log10(n), fill=race))+
  geom_col(position = 'dodge') +
  geom_text(aes(label = n), vjust = 0, hjust=0.5) + scale_fill_manual(values = rev(met.brewer('Moreau', length(unique(tt2$race))))) +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('Clinical Stage Availability by Cancer Type - panTCGA') + theme(legend.position = "none")
dev.off()


