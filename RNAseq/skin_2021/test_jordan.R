################################
# Power Analysis for Blanche 
################################

# package loading
library(RnaSeqSampleSize)
library(dplyr)
library(purrr)
library(ggplot2)

# Set seed 
set.seed(538734517)

# Power analysis for Blanche 

# This code contains the power analysis for Blanche's RNAseq grant. 
# Experiment design across all analyses is simple one factor 
# with heart - left ventricle vs. skin (three skin types); heart - left ventricle vs lung; and comparison of three skin types. Number of expected samples per group is 5-6.

# Set working directory to file locations
#setwd("/Users/jordan/Desktop/Blanche Power Analysis/Data")
#getwd()

# file parsing and cleaning 
counts2 <- readRDS("filtered_data.rds") # this is filtered data, done by August 

# filter rows according to proteome gene list short (will later need to also account for Reactome gene list .csv file)
library(xlsx)
gene_list_short <- read.xlsx("BI-ECM-proteome-gene-list-short.xlsx",1)
names(gene_list_short)
head(gene_list_short)
gene_list_short<-gene_list_short[,c(2, 1, 4)]
names(gene_list_short); head(gene_list_short)
gene_list_short<-gene_list_short[!is.na(gene_list_short$Ensembl),]
genelist<-gene_list_short$Ensembl
#counts2<-counts2[counts2$Name %in% gene_list_short$Ensembl,]
counts2<-data[data$Name %in% gene_list_short$Ensembl,]

##############################################################################################################################
# Counts2 is now your new filtered count matrix; now we select out columns based on which comparison we want to perform
##############################################################################################################################

# First analysis is heart - left ventricle vs skin - not sun exposed 

# Call in file for filtering out columns and filter columns needed for analysis 
columns<-read.delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",  sep="\t", header=TRUE)
heart<-columns[which(columns$SMTSD=="Heart - Left Ventricle"),]$SAMPID # control group 
skin_not_sun<-columns[which(columns$SMTSD=="Skin - Not Sun Exposed (Suprapubic)"),]$SAMPID # treatment group 
col_names<-c(as.character(heart), as.character(skin_not_sun)); col_names
sh_counts<-counts2 %>% select_if(colnames(counts2)=="Name" | colnames(counts2)=="Description" | colnames(counts2) %in% col_names)
sh_datamat <- data.matrix(sh_counts[,3:length(colnames(sh_counts))])
rownames(sh_datamat) <- sh_counts$Name
#head(sh_datamat)

full_counts <- data %>% select_if(colnames(data)=="Name" | colnames(data)=="Description" | colnames(data) %in% col_names)
full_datamat <- data.matrix(full_counts[,3:length(colnames(full_counts))])
rownames(full_datamat) <- full_counts$Name

set.seed(538734517)
subSampleNum <- 20
subSample<-sample(1:ncol(full_datamat),subSampleNum)
countsSub<-full_datamat[,subSample]
groupSub<-colnames(full_datamat) %in% skin_not_sun %>% as.numeric %>% .[subSample]
countsSub<-countsSub[which(rowMeans(full_datamat)>=1),]
y <- DGEList(countsSub, group=groupSub)
y<-calcNormFactors(y)
y<-estimateCommonDisp(y, verbose=TRUE)
y$common.dispersion

# estimate gene read count and dispersion distribution
distribution <- est_count_dispersion(sh_datamat,group=colnames(sh_datamat) %in% skin_not_sun %>% as.numeric)

full_dist <- est_count_dispersion(full_datamat,group=colnames(full_datamat) %in% skin_not_sun %>% as.numeric)

set.seed(123)
subSampleNum <- 20
subSample<-sample(1:ncol(full_datamat),subSampleNum)
countsSub<-full_datamat[,subSample]
groupSub<-colnames(full_datamat) %in% skin_not_sun %>% as.numeric %>% .[subSample]
countsSub<-countsSub[which(rowMeans(full_datamat)>=1),]
y <- DGEList(countsSub, group=groupSub)
y<-calcNormFactors(y)
y<-estimateCommonDisp(y, verbose=TRUE)
y$common.dispersion


packageVersion("edgeR")
y<-readRDS("dgelist.rds")

y <- calcNormFactors(y)
z <- estimateCommonDisp(y, verbose=TRUE)

z_colind <- which(colnames(full_datamat) %in% rownames(z$samples))
z_countsSub<-full_datamat[,z_colind]
z_groupSub<-colnames(full_datamat) %in% skin_not_sun %>% as.numeric %>% .[z_colind]
z_countsSub<-z_countsSub[which(rowMeans(full_datamat)>=1),]

# power estimation
# 6 selected Genes were not found in distributionObject, discarded
power_n5 <- est_power_distribution(n=5,f=0.01,rho=2,distributionObject=distribution,selectedGenes=genelist,storeProcess = TRUE)
mean(power_n5$power) # 0.0022
power_n10 <- est_power_distribution(n=30,f=0.05,rho=2,distributionObject=distribution,selectedGenes=genelist,storeProcess = TRUE)
mean(power_n10$power) #0.60

# power vs sample size plot for different coverage / FDR
coverage10_fdr1 <- est_power_curve(n=40,f=0.01,rho=2,lambda0=10,phi0=distribution$common.dispersion, m=56200, m1=158)
coverage20_fdr1 <- est_power_curve(n=40,f=0.01,rho=2,lambda0=20,phi0=distribution$common.dispersion, m=56200, m1=158)
coverage40_fdr1 <- est_power_curve(n=40,f=0.01,rho=2,lambda0=40,phi0=distribution$common.dispersion, m=56200, m1=158)
coverage60_fdr1 <- est_power_curve(n=40,f=0.01,rho=2,lambda0=60,phi0=distribution$common.dispersion, m=56200, m1=158)
coverage80_fdr1 <- est_power_curve(n=40,f=0.01,rho=2,lambda0=80,phi0=distribution$common.dispersion, m=56200, m1=158)
plot_power_curve(list(coverage10_fdr1,coverage20_fdr1,coverage40_fdr1,coverage60_fdr1,coverage80_fdr1))


coverage10_fdr5 <- est_power_curve(n=40,f=0.05,rho=2,lambda0=10,phi=distribution$common.dispersion, m=56200, m1=158)
coverage20_fdr5 <- est_power_curve(n=40,f=0.05,rho=2,lambda0=20,phi=distribution$common.dispersion, m=56200, m1=158)
coverage40_fdr5 <- est_power_curve(n=40,f=0.05,rho=2,lambda0=40,phi=distribution$common.dispersion, m=56200, m1=158)
coverage60_fdr5 <- est_power_curve(n=40,f=0.05,rho=2,lambda0=60,phi=distribution$common.dispersion, m=56200, m1=158)
coverage80_fdr5 <- est_power_curve(n=40,f=0.05,rho=2,lambda0=80,phi=distribution$common.dispersion, m=56200, m1=158)
plot_power_curve(list(coverage10_fdr5,coverage20_fdr5,coverage40_fdr5,coverage60_fdr5,coverage80_fdr5))

# optimization plots
result <- optimize_parameter(fun=est_power, opt1="n", opt2="lambda0", opt1Value=c(10,15,20,25,30,35,40),opt2Value = c(10,20,40,60,80), f=0.01, phi0=distribution$common.dispersion)

# optimization plots
result <- optimize_parameter(fun=est_power, opt1="n", opt2="lambda0", opt1Value=c(10,15,20,25,30,35,40),opt2Value = c(10,20,40,60,80), f=0.05, phi0=distribution$common.dispersion)


est_power(n=8, lambda0=20, phi0=0.07154, f=0.05, m=52000,m1=71)
power35 <- est_power_distribution(n=35,f=0.05,rho=2,distributionObject=distribution,selectedGenes=genelist,storeProcess = TRUE)
mean(power35$power)

gene_readcounts <- distribution$pseudo.counts.mean[which(names(distribution$pseudo.counts.mean) %in% genelist)]
gene_dispersions <- distribution$tagwise.dispersion[which(names(distribution$pseudo.counts.mean) %in% genelist)]
mean(gene_readcounts)
min(gene_readcounts)
mean(gene_dispersions)
ggplot(data=data.frame(counts=as.numeric(gene_readcounts),name=genenames[order(match(names(gene_readcounts), genelist))])) + geom_bar(aes(x=reorder(name,-counts),y=counts),stat="identity") + scale_y_log10() + coord_flip() + xlab("Gene name") + ylab("Mean count per gene across samples") + theme(axis.text.y = element_text(size=7)) + 
    ggtitle("Mean count per gene across samples from experimental data") + labs(caption="Median of ~82M reads/sample, see GTEx portal for more information.")