## ---------------------------
## Script name: Data_Processing
## Purpose of script:Preparing dataset for statistical analysis
## Author: PhD Ali Kaynar
## Date Created: 01-03-2023
## Copyright (c) Ali Kaynar, 2023
## ---------------------------
## Notes:
## ---------------------------
## set working directory for Mac
setwd("PATH")    # Ali's working directory 
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
#source("PATH/....R") 
## ---------------------------
# Data download process from GDC repository, step by step
projectG="TCGA-GBM" # define project name for glioblastoma as TCGA-GBM
query <- GDCquery(project = projectG,
                  access = "open",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

samplesDown <- getResults(query,cols=c("cases")) # case id is barcode ID. it includes all samples including duplicate
# dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                   typesample = "TP")# TP: Primary Tumor, NT:Tissue Normal, TR: Recurrent Tumor
queryDown <- GDCquery(project = projectG, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      experimental.strategy = "RNA-Seq",
                      workflow.type = "STAR - Counts", 
                      barcode = samplesDown)
GDCdownload(query = queryDown) # download files and store

dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename = paste(projectG,".rda", sep = "")) # prepare and integrate data structure
GeneID<-data.frame(cbind(dataPrep1@rowRanges@elementMetadata@listData[["gene_id"]],dataPrep1@rowRanges@elementMetadata@listData[["gene_name"]], dataPrep1@rowRanges@elementMetadata@listData[["gene_type"]], dataPrep1@rowRanges@elementMetadata@listData[["hgnc_id"]],dataPrep1@rowRanges@elementMetadata@listData[["havana_gene"]]))
colnames(GeneID)<-c("gene_id","gene_name","gene_type","hgnc_id","havana_gene")
indexinG<-which(GeneID$gene_type %in% "protein_coding", arr.ind = TRUE) # extract only protein coding genes
GeneID<-GeneID[indexinG,]
TCGA_GBM<-dataPrep1[indexinG,]# 