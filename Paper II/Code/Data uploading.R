## ---------------------------
## Script name: Data_Processing
## Purpose of script:Preparing dataset for statistical analysis
## Author: PhD Ali Kaynar
## Date Created: 01-03-2023
## Copyright (c) Ali Kaynar, 2023
## Email: 
## ---------------------------
## Notes:
## ---------------------------
## set working directory for Mac
setwd("~/Desktop/UpgrdPostDoc/Data_Coll_Pro")    # Ali's working directory 
## ---------------------------
## load up the packages we will need:
library(TCGAbiolinks)
library(sesameData)
library(sesame)
library(UpSetR)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(NMF)
library(data.table)
library(dplyr) 
library(org.Hs.eg.db)
library(biomaRt)
# source("functions/packages.R")       # loads up all the packages we need
## ---------------------------
## load up our functions into memory
#source("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/....R") 
## ---------------------------
projectG="TCGA-GBM"
query <- GDCquery(project = projectG,
                  access = "open",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

#samplesDown <- getResults(query,cols=c("sample.submitter_id"))
samplesDown <- getResults(query,cols=c("cases"))
# dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                   typesample = "TP")# TP: Primary Tumor, NT:Tissue Normal, TR: Recurrent Tumor
queryDown <- GDCquery(project = projectG, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      experimental.strategy = "RNA-Seq",
                      workflow.type = "STAR - Counts", 
                      barcode = samplesDown)
GDCdownload(query = queryDown)
dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename = paste(projectG,".rda", sep = ""))
GeneID<-data.frame(cbind(dataPrep1@rowRanges@elementMetadata@listData[["gene_id"]],dataPrep1@rowRanges@elementMetadata@listData[["gene_name"]], dataPrep1@rowRanges@elementMetadata@listData[["gene_type"]], dataPrep1@rowRanges@elementMetadata@listData[["hgnc_id"]],dataPrep1@rowRanges@elementMetadata@listData[["havana_gene"]]))
colnames(GeneID)<-c("gene_id","gene_name","gene_type","hgnc_id","havana_gene")
indexinG<-which(GeneID$gene_type %in% "protein_coding", arr.ind = TRUE) #
GeneID<-GeneID[indexinG,]
TCGA_GBM<-dataPrep1[indexinG,]

####################################################################################################################
# %%%%%%%%%%%%%%%----oooo----CGGA Datasets----oooo----%%%%%%%%%%%%%%%
####################################################################################################################
# Load Datasets
CGGA_Control_fpkm<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_RNAseq_Control_20_fpkm.txt")
CGGA_Control_count<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_RNAseq_Control_20_Counts.txt")

CGGA_693_clinical<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_mRNAseq_693_clinical.txt")
CGGA_693_fpkm<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_mRNAseq_693_fpkm.txt")
CGGA_693_count<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_mRNAseq_693_Counts.txt")

CGGA_325_clinical<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_mRNAseq_325_clinical.txt")
CGGA_325_fpkm<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_mRNAseq_325_fpkm.txt")
CGGA_325_count<-read.delim("~/Desktop/UpgrdPostDoc/Data_Coll_Pro/CGGA_mRNAseq_325_Counts.txt")

##########################################################################################################################
data_modified <- sapply(strsplit(CGGA_Control_fpkm$gene_name,"\\."), function(x) x[1])
annot <- select(org.Hs.eg.db,
                keys = data_modified,
                columns = c("GENETYPE",'SYMBOL','ENSEMBL', 'ENTREZID'),#,"GENETYPE", "ENZYME", 'ENTREZID', "GENENAME",
                keytype = 'SYMBOL')
# listDatasets(mart)
# Some gene have ensembl_gene_id but notexternal_gene_name
ensembl109 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=109)
HGNC_2_ENS = getBM(filters = "hgnc_symbol",attributes = c("hgnc_id","hgnc_symbol","ensembl_gene_id", "gene_biotype", "entrezgene_id", "description"), values = CGGA_Control_count$gene_name,mart = ensembl109)#, "entrezgene_id", "description"
HGNC_2_ENS_ProCod<-HGNC_2_ENS[ which(HGNC_2_ENS$gene_biotype=="protein_coding"),]
HGNC_2_ENS_ProCod<-HGNC_2_ENS_ProCod[!duplicated(HGNC_2_ENS_ProCod[,"hgnc_symbol"]),]
# ZHX1-C8ORF76
HGNC_2_ENS_ProCod[which(HGNC_2_ENS_ProCod$hgnc_symbol=="ZHX1-C8orf76"),"hgnc_symbol"]="ZHX1-C8ORF76"  
rownames(HGNC_2_ENS_ProCod)<-HGNC_2_ENS_ProCod$hgnc_symbol
HGNC_2_ENS_ProCod["CCL3L3","ensembl_gene_id"]="ENSG00000276085" # up-to-date
HGNC_2_ENS_ProCod["CCL3L3","entrezgene_id"]="414062"# up-to-date
HGNC_2_ENS_ProCod["ADAMTS13","ensembl_gene_id"]="ENSG00000160323"# up-to-date
save(file="~/Desktop/UpgrdPostDoc/Data_Coll_Pro/HGNC_2_ENS_ProCod.Rdata",HGNC_2_ENS_ProCod)
save(file="~/Desktop/UpgrdPostDoc/Data_Coll_Pro/HGNC_2_ENS.Rdata",HGNC_2_ENS)
##### HGNC_2_ENS_ProCod is the protein coding genes for CGGA database RNASeq datasets, 

##### Control
rownames(CGGA_Control_count)<-CGGA_Control_count$gene_name
CGGA_Control_count_F<-CGGA_Control_count[CGGA_Control_count$gene_name %in% HGNC_2_ENS_ProCod$hgnc_symbol,]
CGGA_Control_count_F<-cbind(HGNC_2_ENS_ProCod, CGGA_Control_count_F)

rownames(CGGA_Control_fpkm)<-CGGA_Control_fpkm$gene_name
CGGA_Control_fpkm_F<-CGGA_Control_fpkm[CGGA_Control_fpkm$gene_name %in% HGNC_2_ENS_ProCod$hgnc_symbol,]
CGGA_Control_fpkm_F<-cbind(HGNC_2_ENS_ProCod, CGGA_Control_fpkm_F)

##### 325
rownames(CGGA_325_count)<-CGGA_325_count$gene_name
CGGA_325_count_F<-CGGA_325_count[CGGA_325_count$gene_name %in% HGNC_2_ENS_ProCod$hgnc_symbol,]
CGGA_325_count_F<-cbind(HGNC_2_ENS_ProCod, CGGA_325_count_F)

rownames(CGGA_325_fpkm)<-CGGA_325_fpkm$Gene_Name
CGGA_325_fpkm_F<-CGGA_325_fpkm[CGGA_325_fpkm$Gene_Name %in% HGNC_2_ENS_ProCod$hgnc_symbol,]
CGGA_325_fpkm_F<-cbind(HGNC_2_ENS_ProCod[HGNC_2_ENS_ProCod$hgnc_symbol %in% CGGA_325_fpkm_F$Gene_Name,1:ncol(HGNC_2_ENS_ProCod)], CGGA_325_fpkm_F)

##### 693
rownames(CGGA_693_count)<-CGGA_693_count$gene_name
CGGA_693_count_F<-CGGA_693_count[CGGA_693_count$gene_name %in% HGNC_2_ENS_ProCod$hgnc_symbol,]
CGGA_693_count_F<-cbind(HGNC_2_ENS_ProCod, CGGA_693_count_F)

rownames(CGGA_693_fpkm)<-CGGA_693_fpkm$Gene_Name
CGGA_693_fpkm_F<-CGGA_693_fpkm[CGGA_693_fpkm$Gene_Name %in% HGNC_2_ENS_ProCod$hgnc_symbol,]
CGGA_693_fpkm_F<-cbind(HGNC_2_ENS_ProCod[HGNC_2_ENS_ProCod$hgnc_symbol %in% CGGA_693_fpkm_F$Gene_Name,1:ncol(HGNC_2_ENS_ProCod)], CGGA_693_fpkm_F)
