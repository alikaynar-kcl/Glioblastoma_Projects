# # rm(list = ls()) #clear all
# control+shift+c #multiple comment
############# oooo ---- PACKAGES LIBRARY ---- oooo #############
library(TCGAbiolinks)
library(DESeq2) # rna-seq
library(biomaRt)
library(rlist)
library(survival)
library(clusterProfiler)
library(limma)
library(heatmaps)
# install.packages(c("biomaRt", "rlist","survival", "limma", "clusterProfiler", "heatmaps"))
library(dplyr) # data wrangling
library(ggplot2) # plotting
library(readr) # Fast readr of files.
library(tidyverse)
library(magrittr)
library(tximport)
library(tximportData)
library(GenomicFeatures)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplotify)
library(plyr)
library(SummarizedExperiment)
library(RegParallel)

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
##########    TCGA GBM RNASeq dataset ###############
load("{PATH}/Data/TCGA_Meta.Rdata")
load("{PATH}/Data/infile_TCGA.Rdata") # expression matrix
load("{PATH}/Data/TCGA_TPM.Rdata") # Count matrix
TCGA_Meta <- as.data.frame(TCGA_Meta)
#################################################################################################################################################
##########    CGGA GBM RNASeq datasets ###############
load("{PATH}/Data/CGGA_693_GraterRSEM1.Rdata")
load("{PATH}/Data/CGGA_693_Count_Grater_1_RSEM.Rdata")
load("{PATH}/Data/CGGA_693_Primary_clinic.Rdata")
load("{PATH}/Data/CGGA_325_Count_Grater_1_RSEM.Rdata")
load("{PATH}/Data/CGGA_325_GraterRSEM1.Rdata")
load("{PATH}/Data/CGGA_325_Primary_clinic.Rdata")

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
############# SURVIVAL analysis for TCGA, CGGA_325 and CGGA_693  ###############

library(survival)
library(tidyverse)
library(imputeTS)

project='GBM_TCGA'
# set folder sot results
raw_path <- "{PATH}/Data/Cox/"
path_out <- "{PATH}/Data/Cox/"
setwd(path_out)
# 
# infile_CGGA_325_N<-t(CGGA_325_GraterRSEM1)
# #rownames(CGGA_693_Primary_clinic)<-CGGA_693_Primary_clinic$CGGA_ID
# infile_CGGA_325_N<-cbind(CGGA_325_Primary_clinic[,c(1, 8, 7)],infile_CGGA_325_N)
# save(file="{PATH}/Data/infile_CGGA_325_N.Rdata",infile_CGGA_325_N)

infile<- infile_TCGA_N

exp_index <- 4
cli_data_all <- infile[, 1:exp_index-1]
cli_data <- cli_data_all[, c('sample_id', 'status', 'survival_time')]
exp_data <- infile[, exp_index:ncol(infile)]
rownames(exp_data) <- cli_data$sample_id
# --- delete donors with NA tpm
ind <- apply(exp_data, 1, function(x) all(is.na(x)))
exp_data <- exp_data[!ind, ]
print(paste(project, 'has ', length(rownames(exp_data)[ind]), 'sample(s) with all NA TPM'))
print(rownames(exp_data)[ind])
sample_keep <- rownames(exp_data)
# --- replace missing expression value with 0
exp_data[is.na(exp_data)] <- 0
# --- keep FPKM >= 1
exp_data<- apply(exp_data, 2, as.numeric)
rownames(exp_data) <- sample_keep
exp_keep <- as.data.frame(exp_data[, which(colMeans(exp_data)>=1)]) ## keep expressed gene

infile <- merge(cli_data, exp_keep, by.x = 'sample_id', by.y = 0)
columns <- colnames(infile)
first_gene_col <- 4
gene_list <- colnames(infile)[first_gene_col:ncol(infile)]

infile_TCGA_N<-infile
survival_time <- as.numeric(infile_TCGA_N$survival_time)
status <- as.numeric(infile_TCGA_N$status)
sample_id <- as.factor(infile_TCGA_N$sample_id)

survival_time <- as.numeric(infile$survival_time)
status <- as.numeric(infile$status)
sample_id <- as.factor(infile$sample_id)

coef_each_gene<-NULL
p_each_gene<-NULL

for (j in first_gene_col:length(infile)){
  print(j)
  int_exp <- as.numeric(infile[,j])
  cox_result=summary(coxph(Surv(survival_time, status) ~ as.numeric(int_exp))) ##confirm label
  coef=as.matrix(cox_result$coefficients)[1,"coef"]
  p=as.matrix(cox_result$coefficients)[1,"Pr(>|z|)"]
  coef_each_gene<-cbind(coef_each_gene,coef)
  p_each_gene<-cbind(p_each_gene,p)
}
colnames(coef_each_gene)<-gene_list
colnames(p_each_gene)<-gene_list

out_file <- rbind(coef_each_gene, p_each_gene)
out_file_T <- t(out_file)
colnames(out_file_T) <- c('coef', 'p')
out_file_T <- as.data.frame(out_file_T)
out_file_T$ensembl <- rownames(out_file_T)
final <- out_file_T[,c('ensembl', 'coef', 'p')]
final$FDR_BH <- p.adjust(final$p, method = "BH")

## get average FPKM
exp_data<- apply(exp_keep, 2, as.numeric)
rownames(exp_data) <- infile$sample_id
mean_EXP <- as.data.frame(colMeans(exp_data))
colnames(mean_EXP) <- 'mean_exp' 

final_result2 <- merge(final, mean_EXP, by.x = 'ensembl', by.y = 0)
final_result2 <- final_result2 %>% relocate(mean_exp, .after=ensembl)
colnames(final_result2)[1] <- "genesymbol"
#file_out = paste0(path_out,  project, '_cox_p_GNSBL.txt')
file_out = paste0(project, '_cox_p_GNSBL.txt')
write.table(final_result2, file = file_out, row.names=F, col.names=T, sep = '\t', quote = F)
print(dim(final_result2))
print(paste0('Finsh cox for ', project))
print(paste0(project, ' including ', length(sample_keep), ' donors, ', 
             'with ', length(gene_list), ' expressed genes'))

Order_Final_result2<-final_result2[final_result2$p<0.05, ]
# save(file="{PATH}/Data/Order_Final_result2_FPKM.Rdata",Order_Final_result2)
# FPKM_cox<-Order_Final_result2
# TPM_cox<-Order_Final_result2
#################################################################################################################################################
### Favorable and Unfavorable gene selection for TCGA, CGGA_325 and CGGA_693

TCGA_Cox <- read.delim("GBM_TCGA_cox_p_GNSBL.txt")
# CGGA_693_Cox <- read.delim("GBM_CGGA_693_cox_p_GNSBL.txt")
# CGGA_325_Cox <- read.delim("GBM_CGGA_325_cox_p_GNSBL.txt")

Favorable_TCGA_Cox<-NULL
UnFavorable_TCGA_Cox<-NULL  
for (i in 1:nrow(TCGA_Cox)){
  if (TCGA_Cox$p[i] <0.05 && TCGA_Cox$coef[i] <0){Favorable_TCGA_Cox<-rbind(Favorable_TCGA_Cox,TCGA_Cox[i,])}
  else if (TCGA_Cox$p[i] <0.05 && TCGA_Cox$coef[i] >0){UnFavorable_TCGA_Cox<-rbind(UnFavorable_TCGA_Cox,TCGA_Cox[i,])}
}

write.table(Favorable_TCGA_Cox, file = "Favorable_TCGA_Cox.txt", row.names=F, col.names=T, sep = '\t', quote = F)
write.table(UnFavorable_TCGA_Cox, file = "UnFavorable_TCGA_Cox.txt", row.names=F, col.names=T, sep = '\t', quote = F)

#################################################################################################################################################
#################################################################################################################################################
##########  DEG Analysis for TCGA,  CGGA_325 and CGGA_693  ###############


TCGA_GBM_P<-infile_TCGA
TCGA_GBM_P<-TCGA_GBM_P[TCGA_GBM_P$status==1,] # Only samples from death patient used for DEG analysis
mode(TCGA_GBM_P$survival_time)="numeric"
ctoff_TCGA<-quantile(TCGA_GBM_P$survival_time, probs = c(0.33, 0.66)) # quantile setting

condition_TCGA<-NULL
for (i in 1:length(TCGA_GBM_P$survival_time)){
  if (TCGA_GBM_P$survival_time[i]>ctoff_TCGA[2])
  {condition_TCGA[i]<-"High"} # High survived patient 
  else if (TCGA_GBM_P$survival_time[i]<ctoff_TCGA[1])
  {condition_TCGA[i]<-"Low"}  # Low survived patient
  else{condition_TCGA[i]<-"Mediann"}
  condition_TCGA<-data.frame(condition_TCGA)
}
condition_TCGA<-t(condition_TCGA)
condition_TCGA<-data.frame(condition_TCGA)

TCGA_Count_1TPM<-TCGA_TPM# Count data of TCGA
rownames(TCGA_Count_1TPM)<-TCGA_Meta$sample
TCGA_Count_D<-TCGA_Count_1TPM[TCGA_GBM_P$sample_id,]

TCGA_Count_D_Clinic<-TCGA_GBM_P[, c(1:3)]
TCGA_Count_D_Clinic<-cbind(TCGA_Count_D_Clinic,condition_TCGA)

row_to_remove = which(TCGA_Count_D_Clinic$condition_TCGA=="Mediann") # median survived patient removed from data set
TCGA_Count_D_Clinic<-TCGA_Count_D_Clinic[-row_to_remove,]

rownames(TCGA_Count_D_Clinic)<-TCGA_Count_D_Clinic$sample_id
TCGA_Count_D<-TCGA_Count_D[TCGA_Count_D_Clinic$sample_id,]
TCGA_Count_D_Clinic$condition_TCGA<-as.factor(TCGA_Count_D_Clinic$condition_TCGA)
TCGA_Count_D<-t(TCGA_Count_D)
# DESeq data structure from matrix 
dds_TCGA_Count_D <- DESeqDataSetFromMatrix(countData=TCGA_Count_D,
                                           colData=TCGA_Count_D_Clinic, 
                                           design=~condition_TCGA)
dds_TCGA_Count_D$condition_TCGA <- relevel(dds_TCGA_Count_D$condition_TCGA, ref = "High") # High survived samples set as reference
dds_TCGA_Count_D_R = DESeq(dds_TCGA_Count_D)
result_TCGA_Count_D = results(dds_TCGA_Count_D_R )
summary(result_TCGA_Count_D)
View(as.data.frame(result_TCGA_Count_D))

nrow(result_TCGA_Count_D)
sum( is.na(result_TCGA_Count_D$pvalue) )
res_TCGA_Count_D = result_TCGA_Count_D[ ! is.na(result_TCGA_Count_D$pvalue), ]
res_TCGA_Count_D = res_TCGA_Count_D[ order(res_TCGA_Count_D$pvalue), ]

DESeq_Table_TCGA_Count_D<-as.data.frame(res_TCGA_Count_D)

DEG_TCGA_UP<-DESeq_Table_TCGA_Count_D[DESeq_Table_TCGA_Count_D$log2FoldChange>0,]
DEG_TCGA_UP<-DEG_TCGA_UP[DEG_TCGA_UP$pvalue<0.01,] # up-regulated genes

DEG_TCGA_DOWN<-DESeq_Table_TCGA_Count_D[DESeq_Table_TCGA_Count_D$log2FoldChange<0,]
DEG_TCGA_DOWN<-DEG_TCGA_DOWN[DEG_TCGA_DOWN$pvalue<0.01,] # down regulated genes

####################################################################################
# Same procedure for CGGA data sets 

CGGA_693_clinic<-CGGA_693_Primary_clinic[CGGA_693_Primary_clinic$Censor..alive.0..dead.1.==1,]# Only samples from death patient used for DEG analysis
rownames(CGGA_693_clinic)<-CGGA_693_clinic$CGGA_ID
CGGA_693_COUNT<-CGGA_693_Count_Grater_1_RSEM[,CGGA_693_clinic$CGGA_ID]

cutoffs_CGGA_693 = sort(unique(CGGA_693_clinic$OS))
Percentile33_CGGA_693 = quantile(CGGA_693_clinic$OS,0.33, na.rm=T)
Percentile66_CGGA_693 = quantile(CGGA_693_clinic$OS,0.66, na.rm=T)

condition_CGGA_693<-NULL
for (i in 1:length(CGGA_693_clinic$OS)){
  if (CGGA_693_clinic$OS[i]>Percentile66_CGGA_693)
  {condition_CGGA_693[i]<-"High"} # High survived patient 
  else if (CGGA_693_clinic$OS[i]<Percentile33_CGGA_693)
  {condition_CGGA_693[i]<-"Low"} # Low survived patient
  else{condition_CGGA_693[i]<-"Mediann"}
  condition_CGGA_693<-data.frame(condition_CGGA_693)
}
condition_CGGA_693<-t(condition_CGGA_693)
condition_CGGA_693<-data.frame(condition_CGGA_693)
CGGA_693_clinic<-cbind(condition_CGGA_693,CGGA_693_clinic)

row_to_remove = which(CGGA_693_clinic$condition_CGGA_693=="Mediann") # median survived patient removed from data set
CGGA_693_clinic<-CGGA_693_clinic[-row_to_remove,]
CGGA_693_COUNT<-CGGA_693_COUNT[,CGGA_693_clinic$CGGA_ID]
rownames(CGGA_693_clinic)<-CGGA_693_clinic$CGGA_ID
CGGA_693_clinic$condition_CGGA_693<-as.factor(CGGA_693_clinic$condition_CGGA_693)

dds_CGGA_693 <- DESeqDataSetFromMatrix(countData=CGGA_693_COUNT,
                                       colData=CGGA_693_clinic, 
                                       design=~condition_CGGA_693)

dds_CGGA_693$condition_CGGA_693 <- relevel(dds_CGGA_693$condition_CGGA_693, ref = "High") # High survived samples set as reference
dds_CGGA_693_R = DESeq(dds_CGGA_693)
result_CGGA_693_Count_D = results(dds_CGGA_693_R)
summary(result_CGGA_693_Count_D)
View(as.data.frame(result_CGGA_693_Count_D))

nrow(result_CGGA_693_Count_D)
sum( is.na(result_CGGA_693_Count_D$pvalue) )
res_CGGA_693_Count_D = result_CGGA_693_Count_D[ ! is.na(result_CGGA_693_Count_D$pvalue), ]
res_CGGA_693_Count_D = res_CGGA_693_Count_D[ order(res_CGGA_693_Count_D$pvalue), ]

DESeq_Table_CGGA_693_Count_D<-as.data.frame(res_CGGA_693_Count_D)

DEG_CGGA_693_UP<-DESeq_Table_CGGA_693_Count_D[DESeq_Table_CGGA_693_Count_D$log2FoldChange>0,]
DEG_CGGA_693_UP<-DEG_CGGA_693_UP[DEG_CGGA_693_UP$pvalue<0.01,] # up-regulated genes

DEG_CGGA_693_DOWN<-DESeq_Table_CGGA_693_Count_D[DESeq_Table_CGGA_693_Count_D$log2FoldChange<0,]
DEG_CGGA_693_DOWN<-DEG_CGGA_693_DOWN[DEG_CGGA_693_DOWN$pvalue<0.01,] # down regulated genes


#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#############   COX Hyper geometric test   ##############
deg_up_1 <- as.matrix(Unfavorable_TCGA_Cox$genesymbol)
#deg_up_1<-data.frame(deg_up_1)
deg_down_1<- as.matrix(Favorable_TCGA_Cox$genesymbol)
#deg_down_1<-data.frame(deg_down_1)
deg_up_2 <- as.matrix(UnFavorable_CGGA_693_Cox$genesymbol)
#deg_up_2<-data.frame(deg_up_2)
deg_down_2<- as.matrix(Favorable_CGGA_693_Cox$genesymbol)
#deg_down_2<-data.frame(deg_down_2)
bg_gene <- as.matrix(intersect(TCGA_Cox$genesymbol, CGGA_693_Cox$genesymbol))

TCGA_vs_CGGA_693_Cox<- DEGs_list_overlap(deg_up_1, deg_down_1, deg_up_2, deg_down_2, bg_gene)  
TCGA_vs_CGGA_693_Cox[["over_up"]] # Unfavorable overlaped  gene list
TCGA_vs_CGGA_693_Cox[["over_down"]] # Favorable overlaped  gene list

###### TWO gene list, Hyper geometric test, DEG
deg_up_1<- as.matrix(rownames(DEG_TCGA_UP))
deg_up_2<- as.matrix(rownames(DEG_CGGA_693_UP))
bg_num_d<-as.matrix(intersect(rownames(DESeq_Table_TCGA_Count_D), rownames(DESeq_Table_CGGA_693_Count_D)))
TCGA_vs_CGGA_693_CoxUnfav<-two_gene_list_overlap(deg_up_1,deg_up_2,bg_num_d)

#################################################################################################################################################
##################           WGCNA           #####################
library(WGCNA)
library(DESeq2)

path_raw<-"{PATH}/WGCNA_GBM/WGCNA_G/"
setwd(path_raw)
source(paste0(path_raw, 'WGCNA_topology_analysis.R'))

network_type<-"signed"#unsigned, signed

options(stringsAsFactors = FALSE)
# enableWGCNAThreads()
allowWGCNAThreads() # 

CGGA_325<-CGGA_325_Count_Grater_1_RSEM
CGGA_693<-CGGA_693_Count_Grater_1_RSEM
TCGA<- TCGA_TPM

GBM_Coh <- c('CGGA_325','CGGA_693','TCGA')  # They contain infile data
GBM_Coh_list = list(CGGA_325,CGGA_693,TCGA)
i=2
# for (i in 1:3){
exp_merge <- as.matrix.data.frame(GBM_Coh_list[[i]])
#  abc<-data.frame(is.na(exp_merge)), there is no missing value
# GOODSAMPL<-goodSamplesGenes(t(exp_merge), weights = NULL,minRelativeWeight = 0.1,verbose = 1, indent = 0)
# index<-which(GOODSAMPL[["goodGenes"]] %in% "FALSE")
# exp_merge<-exp_merge[-index,]
#count normalization
exp_vst<-vst(exp_merge)#vst row is gene, col is sample
powers<-c(1:25)
exp_vst<-t(exp_vst)

stand_dev = apply(exp_vst, 2, sd)
hist(stand_dev, breaks = 1000)

sum(stand_dev < 0.25)
abline(v = 0.25, col = "red")
exp_vst = exp_vst[, stand_dev >= 0.25]  # It selected depending on graphic (by eye)

sft = pickSoftThreshold(exp_vst, powerVector = powers,verbose = 5,blockSize = 14541,networkType = network_type)#pickSoftThreshold: row is sample, col is gene
FileName=paste("sft_selection_",GBM_Coh[1],".Rdata", sep="")
save(file=paste("sft_selection_",GBM_Coh[i],".Rdata", sep=""),sft)

#plot for the selection of soft threshold

pdf(file=paste("soft_threshold_selection_",GBM_Coh[i],".pdf",sep=""), width=12,height = 10)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.8,col="red");
abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.8,col="red")
dev.off()
###############################################################
#WGCNA: step-by-step network constructin
#softPower=sft$powerEstimate
adjacency = adjacency(exp_vst, power = 12,type = network_type)#-----------take long time-----
## # compute a matrix of gene-gene correlation values

TOM = TOMsimilarity(adjacency,TOMType = network_type);#-----------take long long time-----
# transform the adjacency into Topological Overlap Matrix (TOM)

dissTOM = 1-TOM  # calculate the corresponding dissimilarity (dissTOM)

# save(file=paste("adjacency_",GBM_Coh[i],".Rdata"),adjacency)
# Call the hierarchical clustering function (for gene-gene clustering)
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

minModuleSize = 20;#this parameter can be changed, e.g. minModuleSize=30
# Module identification using dynamic tree cut: method is hybrid?
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 1, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
#table(dynamicColors)
# Plot the dendrogram and colors underneath
#module primary
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
###################################################################################
## Calculate eigengenes and merge the modules whose expression profiles are very similar
# There is some NaN issues, it maybe related to some genes that have very small principal component contribution.
MEList = moduleEigengenes(exp_vst, colors = dynamicColors) #Principal component
### *** NA value, grey color >>> MEs<-na.omit(MEs)
### *** abc<-data.frame(is.na(exp_vst))
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);

METree = hclust(as.dist(MEDiss), method = "average");
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25 # it is about module size
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(exp_vst, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#sizeGrWindow(12, 9)
#module before after merging


png(paste("Cluster_Dendrogram_before_after_merge_",GBM_Coh[i],".png",sep=""), width = 600, height = 300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# Rename to moduleColors
moduleColors = mergedColors# final result: this is the module lable 
png(paste("Cluster_Dendrogram_Final_merge_",GBM_Coh[i],".png",sep=""), width = 600, height = 300)
plotDendroAndColors(geneTree, moduleColors,
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# Construct numerical labels corresponding to the colors
#colorOrder = c("grey", standardColors(50));
#moduleLabels = match(moduleColors, colorOrder)-1;
#MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts

save(mergedMEs, moduleColors, geneTree, file = paste("network-stepByStep_",GBM_Coh[i],".RData",sep=""))


File_NameD = degree_WGCNA(adjacency,moduleColors) 
File_NameD<-cbind(File_NameD,data.frame(moduleColors))
All_in_one<-NULL
all_M_color<-unique(moduleColors)
for (j in 1:length(all_M_color)){
  ind<-which(File_NameD$moduleColors %in% all_M_color[j])
  
  File_NameDegre = File_NameD[rownames(File_NameD)[ind],] 
  # assign(paste("Degree_",GBM_Coh[i], sep = ""), File_NameD )
  
  File_NameC =clustering_coefficient_WGCNA(adjacency,rownames(File_NameD)[ind])
  # assign(paste("Degree_",GBM_Coh[i], sep = ""), File_NameC )
  
  File_NameM=module_membership_WGCNA(exp_vst,MEs,all_M_color[j],rownames(File_NameD)[ind])
  
  # File_NameCor=corr_module_eigengene_vs_trait(MEs,trait,uni_cluster)
  Modul_top_Val<-cbind(File_NameDegre,File_NameC,File_NameM)
  All_in_one<-rbind(All_in_one,Modul_top_Val)
}

save(All_in_one,file=paste("All_in_one_",GBM_Coh[i],".RData",sep=""))

################################################################################

#Topology analysis for the modules identified by WGCNA
#Degree (intramodule connectivity) of genes in a specific module
#Clustering coefficient (cc) for genes in a specific module
#Module memership of genes in a specific module
library(WGCNA)
###############################################################################
###############################################################################
#Degrees of genes in a specific module
degree_WGCNA<-function(adjacency,moduleColors,moduleGenes){
  #input: 
  #adjacency: row and col names are gene symbols; the adjacency of whole network
  #moduleColors: its order should be mapped with the rowname or colnames (genes) of adjacency
  
  #output: 
  #degree: 4 colunms, total connectivity (kTotal), intramodular connectivity (kWithin), extra-modular connectivity (kOut), and the difference of the intra- and extra-modular connectivities (kDiff) for all genes; kTotal represents the degree of genes in the specific module.
  
  degree<-intramodularConnectivity(adjacency, as.matrix(moduleColors), scaleByMax = FALSE)
  degree<-degree[moduleGenes,]#order is same as moduleGenes
  return(degree)#mode is list
}
###############################################################################
###############################################################################
#Clustering coefficients of genes in a specific module
clustering_coefficient_WGCNA<-function(adjacency,moduleGenes){
  #adjacency can be all measured genes' or specific modules' adjacency
  int_adjacency<-adjacency[moduleGenes,moduleGenes]
  cc<-as.matrix(clusterCoef(int_adjacency))
  cc<-as.matrix(cc[moduleGenes,1])
  colnames(cc)<-"cc"
  cc<-as.data.frame(cc)
  return(cc)#mode is list
}
###############################################################################
###############################################################################
#correlation between module eigengenes and trait for meta_data
corr_module_eigengene_vs_trait<-function(MEs,trait,uni_cluster){
  #MEs: eigengenes matrix; row is sample, col is module first component
  #MEs = moduleEigengenes(exp_vst, moduleColors)$eigengenes
  #MEs = orderMEs(MEs)
  
  #trait: from meta_data; e.g. age; its order should be matched with row order of MEs
  #uni_cluster: cluster information
  nSamples<-nrow(MEs)
  moduleTraitCor = cor(MEs, trait, method="spearman",use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  loc<-match(uni_cluster$moduleColors,substring(rownames(moduleTraitCor),3))
  uni_cluster$moduleTraitCor=moduleTraitCor[loc]
  uni_cluster$moduleTraitPvalue=moduleTraitPvalue[loc]
  return(uni_cluster)#mode is list
}
###############################################################################
###############################################################################
#Module memberships of genes in a specific module
module_membership_WGCNA<-function(exp_vst,MEs,int_module,moduleGenes){
  #exp_vst: the vst normalization of count exp by DEseq2; row is sample, col is gene
  #MEs: eigengenes matrix; row is sample, col is module
  #MEs = moduleEigengenes(exp_vst, moduleColors)$eigengenes
  #MEs = orderMEs(MEs)
  #int_module: interested module's color
  #int_module and moduleGenes should be mapped (same module)
  
  nSamples=nrow(exp_vst)
  geneModuleMembership = as.data.frame(cor(exp_vst, MEs, method="spearman", use = "p"))#mode is list
  #geneModuleMembership: row is gene (13700), col is module (19)
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))#mode is list
  
  int_module=paste0("ME",int_module)
  MM<-as.matrix(geneModuleMembership[moduleGenes,int_module])
  MMP<-as.matrix(MMPvalue[moduleGenes,int_module])
  
  MM_matrix<-cbind(MM,MMP)
  colnames(MM_matrix)<-c("MM","MMPValue")
  rownames(MM_matrix)<-moduleGenes
  MM_matrix<-as.data.frame(MM_matrix)
  return(MM_matrix)#mode is list
}

#######################################################################################################################################
#######################################################################################################################################
# This part finds the prognostic modules 
load("{APTH}/All_in_one_WGCNA.RData") # All WGCNA results for three cohort 
TCGA_vs_CGGA_693_Cox[["over_up"]] 
list3 <- list.append(TCGA_vs_CGGA_693_Cox[["over_up"]] ,TCGA_vs_CGGA_325_Cox[["over_up"]], CGGA_325_vs_CGGA_693_Cox[["over_up"]])
gene_list_1<-as.matrix(unique(list3)) # Cox unfavorable gene pool
Backround_Genes<-length(intersect(TCGA_Cox$genesymbol,rownames(All_in_one_CGGA_325))) # Background gene list

result_All_in_one_CGGA_325<-NULL
Module_match_genes_CGGA_325<-NULL
Module_color<-unique(All_in_one_CGGA_325$moduleColors)# Modules , module numbers
for (i in 1:length(Module_color)){
  index_color<-which(All_in_one_CGGA_325$moduleColors %in% Module_color[i])# module index, which gene belongs to which module
  Module_genes<-as.matrix(rownames(All_in_one_CGGA_325[index_color,]))# Gathering module genes depending on module color
  result_each<-two_gene_list_overlap(gene_list_1, Module_genes,Backround_Genes) # Cox gene pool overlap with module gene list
  result_each<-cbind(result_each,"CGGA_325",Module_color[i]) # 
  result_All_in_one_CGGA_325<-rbind(result_All_in_one_CGGA_325,result_each)# result table
  
  match_Genes = matrix(0,length(Module_genes),1)# appointing module genes, that are also Cox genes, over the all gene list used for WGCNA
  index_ss<-which(Module_genes %in% gene_list_1)
  match_Genes[index_ss,]=1 # match gene marked as 1 in the table, so we can sort it over the all modules
  
  S1<-rep("CGGA_325", length(Module_genes))
  S2<-rep(Module_color[i], length(Module_genes))
  
  M_M_info<-cbind(Module_genes,match_Genes, S1, S2,All_in_one_CGGA_325[index_color,])
  Module_match_genes_CGGA_325<-rbind(Module_match_genes_CGGA_325,M_M_info)# final table formation with macth genes 
}
######
#############################################################################################################################
# Finding overlapped modules over the three data sets
# Example
####  MODULES Gene Comparison---All_in_one_CGGA_693   -vs-  All_in_one_TCGA
Backround_Genes<-length(intersect(rownames(All_in_one_TCGA),rownames(All_in_one_CGGA_693)))

Modul_comp_693_TCGA<-NULL
modul_color_693<-unique(All_in_one_CGGA_693$moduleColors)
modul_color_TCGA<-unique(All_in_one_TCGA$moduleColors)
for (i in 1:length(modul_color_693)){
  index_color1<-which(All_in_one_CGGA_693$moduleColors %in% modul_color_693[i])
  Module_genes1<-as.matrix(rownames(All_in_one_CGGA_693[index_color1,]))
  for (j in 1:length(modul_color_TCGA)){
    index_color2<-which(All_in_one_TCGA$moduleColors %in% modul_color_TCGA[j])
    Module_genes2<-as.matrix(rownames(All_in_one_TCGA[index_color2,]))
    result_each<-two_gene_list_overlap(Module_genes1, Module_genes2,Backround_Genes)
    result_each<-cbind(result_each,modul_color_693[i],modul_color_TCGA[j])
    Modul_comp_693_TCGA<-rbind(Modul_comp_693_TCGA,result_each)
    
  }
}
# write_xlsx(as.data.frame(Modul_comp_693_TCGA),"~/Desktop/UpgrdPostDoc/WGCNA_GBM/WGCNA_G/Modul_comp_693_TCGA.xlsx")
# save(Modul_comp_693_TCGA,Modul_comp_325_TCGA, Modul_comp_325_693, file="Module_Comp_3.RData")

############################################################################################################################
## Go Term enrichment for marker genes (overlapped prognostic genes)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)

source('deg_GoTerm_clusterProfiler.R')
mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl")

List_Modules_Genes<-list(rownames(All_in_one_CGGA_325),rownames(All_in_one_CGGA_693),rownames(All_in_one_TCGA))
background<-as.matrix(Reduce(intersect, List_Modules_Genes))

# Prognostic modules' GO results. Same procedures applied to other interested gene list

blue_325 <-rownames(All_in_one_CGGA_325[All_in_one_CGGA_325$moduleColors=="blue",])

blue_TCGA <-rownames(All_in_one_TCGA[All_in_one_TCGA$moduleColors=="blue",]) 
darkred_TCGA <- rownames(All_in_one_TCGA[All_in_one_TCGA$moduleColors=="darkred",])
salmon_TCGA <- rownames(All_in_one_TCGA[All_in_one_TCGA$moduleColors=="salmon",])

greenyellow_693 <- rownames(All_in_one_CGGA_693[All_in_one_CGGA_693$moduleColors=="greenyellow",])
lightcyan_693 <- rownames(All_in_one_CGGA_693[All_in_one_CGGA_693$moduleColors=="lightcyan",])

purple_693<-rownames(All_in_one_CGGA_693[All_in_one_CGGA_693$moduleColors=="purple",])


int_gene<-as.matrix(blue_325)
# int_gene<-as.matrix(salmon_TCGA)
# int_gene<-as.matrix(darkred_TCGA)
# int_gene<-as.matrix(blue_TCGA)
# int_gene<-as.matrix(lightcyan_693)
# int_gene<-as.matrix(greenyellow_693)
################################################################################

#int_gene<-as.matrix(read.csv("KIRC_markers.txt",header=T,sep="\t"))
heder_T<-paste("GO terms for Module blue_325 WGCNA", sep=" ")
ont1<-"CC"

pathway=enrichGO(gene=int_gene,OrgDb='org.Hs.eg.db',keyType = "SYMBOL",ont=ont1,universe = background ,pAdjustMethod = "BH",qvalueCutoff=0.05)
path_table_all<-deg_GoTerm_clusterProfiler(pathway)
Genelist_blue_salmon_lightcyan_goterm<-path_table_all


Headers_GO=paste(heder_T,ont1,sep=" ")
names_GF<-paste(Headers_GO, "png", sep=".")
names_GO<-paste(Headers_GO, "txt", sep=".")

png(paste("Bar plot", names_GF, sep=" "), width = 12, height = 12, units = 'in', res = 1000)
par(cex.main=3,cex.axis=1.5, cex.names=1.5)#cex.axis=1, cex.lab=1, , cex.sub=1
barplot(pathway, title=paste("Bar plot", Headers_GO, sep=" "),x = "Count", color = "p.adjust",font.size = 20,showCategory=20)
dev.off()
write.table(path_table_all,file=names_GO,sep="\t",row.names=F,col.names=T,quote=F)
# genesymbol to entrez gene ID 
int_gene_kegg = getBM(attributes = c("hgnc_symbol","entrezgene_id"),values =int_gene ,filters = "hgnc_symbol",mart = mart)
int_gene_kegg<-int_gene_kegg %>% drop_na()
background2 = getBM(attributes = c("hgnc_symbol","entrezgene_id"),values =background,filters = "hgnc_symbol",mart = mart)
background2<-background2 %>% drop_na()
#KEGG pathway
ekegg = enrichKEGG(gene = as.character(int_gene_kegg$entrezgene_id),
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   universe=as.character(background2$entrezgene_id),qvalueCutoff = 0.05,
                   use_internal_data = FALSE)

png(paste("Bar KEGG plot", names_GF, sep=" "), width = 12, height = 12, units = 'in', res = 1000)
barplot(ekegg, title=paste("Bar KEGG plot", names_GO, sep="_"), x = "Count", color = "p.adjust",font.size = 20, showCategory=20)
dev.off()

ekegg_print<-ekegg@result
write.table(ekegg_print ,file=paste("Bar KEGG plot",heder_T,"txt", sep="."),sep="\t",row.names=F,col.names=T,quote=F)

################################################################################
#### JAccard Index
############# Jaccard Function
Jaccard <- function(module_list1, module_list2) {
  i_len <- length(intersect(module_list1, module_list2))
  u_len <- length(union(module_list1, module_list2))
  return(i_len/u_len)
}
# changing module name, C1: CGGA_325, C2: CGGA_693, T: TCGA
CGGA_325_M_Color<-unique(All_in_one_CGGA_325$moduleColors)
CGGA_693_M_Color<-unique(All_in_one_CGGA_693$moduleColors)
TCGA_M_Color<-unique(All_in_one_TCGA$moduleColors)

for (i in 1:length(CGGA_325_M_Color)){
  M_name<-paste("C1.M",i,sep = "_")
  All_in_one_CGGA_325$moduleColors[All_in_one_CGGA_325$moduleColors == CGGA_325_M_Color[i]] <- M_name
}

for (i in 1:length(CGGA_693_M_Color)){
  M_name<-paste("C2.M",i,sep = "_")
  All_in_one_CGGA_693$moduleColors[All_in_one_CGGA_693$moduleColors == CGGA_693_M_Color[i]] <- M_name
}

for (i in 1:length(TCGA_M_Color)){
  M_name<-paste("T.M",i,sep = "_")
  All_in_one_TCGA$moduleColors[All_in_one_TCGA$moduleColors == TCGA_M_Color[i]] <- M_name
}
# Making Module Genes list
ModuleName_CGGA_693<-unique(All_in_one_CGGA_693$moduleColors)
cluster_genes_693<-rownames(All_in_one_CGGA_693[which(All_in_one_CGGA_693$moduleColors==ModuleName_CGGA_693[1]),])  

ModuleName_TCGA<-unique(All_in_one_TCGA$moduleColors)
cluster_genes_TCGA<-rownames(All_in_one_TCGA[which(All_in_one_TCGA$moduleColors==ModuleName_TCGA[1]),])  

nt<-length(ModuleName_CGGA_693)
nc<-length(ModuleName_TCGA)
out_Jacard_693_TCGA <- matrix(0, nrow = nt, ncol = nc, dimnames = list(ModuleName_CGGA_693, ModuleName_TCGA))
for (i in 1:nt) {
  for (j in 1:nc) {
    modulnum_t<-as.numeric(ModuleName_CGGA_693[i])
    GeneList_t<-rownames(All_in_one_CGGA_693[which(All_in_one_CGGA_693$moduleColors==ModuleName_CGGA_693[i]),])  
    modulnum_c<-as.numeric(ModuleName_TCGA[j])
    GeneList_c<-rownames(All_in_one_TCGA[which(All_in_one_TCGA$moduleColors==ModuleName_TCGA[j]),]) 
    
    out_Jacard_693_TCGA[i, j] <- Jaccard(GeneList_t, GeneList_c)
  }
}

png("Jaccard index of 693 and TCGA Modules.png", width = 16, height = 12, units = 'in', res = 400)
pheatmap::pheatmap(out_Jacard_693_TCGA, display_numbers = matrix(ifelse(out_Jacard_693_TCGA>=0.028, "*", ""), nrow(out_Jacard_693_TCGA)), main = "Jaccard index of 693 and TCGA Modules",
                   cellwidth = 30, cellheight = 20, fontsize = 16, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                   show_rownames = T, show_colnames = T)
# grid.text("693 Modules", y=-0.07, gp=gpar(fontsize=16))
# grid.text("325 Modules", x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()

