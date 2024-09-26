## ---------------------------
## Script name: Survival Analysis, KM and Cox
## Purpose of script:Perform survival analysis for Cancer 
## Author: PhD Ali Kaynar
## Date Created: 12-03-2023
## Copyright (c) Ali Kaynar, 2023
## Email: ali.kaynar@kcl.ac.uk
## ---------------------------
## Notes: 
# rm(list=ls())
## ---------------------------
## set working directory for Mac
setwd("{PATH}")    # Ali's working directory
## load up the packages we will need:
library(dplyr)
library(survival)
library(TCGAbiolinks)
GenerateSurvInputfromTCGA <- function(gene,cancerType,data) {
  #if (!hasArg("dataDir")) {
  #  dataDir = paste0("C:/Pan_cancer/TCGA_transcript_result_store/",cancerType,"/PKLR/")
  #} 
  LivingDays = as.numeric(data$survival_time)
  SurvInput= as.data.frame(LivingDays)
  SurvInput$EXP = as.numeric(data[,gene])
  SurvInput$DeadInd = data$status %in% c('1')
  SurvInput$SurvObj <- with(SurvInput, Surv(LivingDays, DeadInd))
  return(SurvInput)
}

raw_path <- "{PATH}/"
outPath <- "{PATH}/"
setwd(raw_path)
source(paste0(outPath, 'cheng_toolbox_beta_meng_no_plot.R'))

################################################################################
############################## DATA LOADINGG ###################################
load("{PATH}/TCGA_GBM.Rdata")
load("{PATH}/TCGA_GBM_Meta.RData")
################################################################################
##################### KAPLAN-MEIER (KM) ANALYSIS ###############################

# Create infile for survival, TCGA
cancer_name <- 'TCGA_GBM_tpm_Pur'
ens_gene_id<-TCGA_GBM@rowRanges@elementMetadata@listData[["gene_id"]]
ens_gene_id<- sapply(strsplit(ens_gene_id,"\\."), function(x) x[1])
gene_name<-TCGA_GBM@rowRanges@elementMetadata@listData[["gene_name"]]

infile<-TCGA_GBM@colData@rownames
infile<-data.frame(infile)
colnames(infile)<-c("sample_id")
infile$definition<-TCGA_GBM_Meta$definition
infile$purity<-TCGA_GBM_Meta$paper_ABSOLUTE.purity
infile$age<- TCGA_GBM_Meta$age_at_index
infile$status<- TCGA_GBM_Meta$vital_status
infile$survival_time<- TCGA_GBM_Meta$days_to_death
index_survival_NA<-which(is.na(infile$survival_time))
for (i in 1:length(index_survival_NA)){
  infile$survival_time[index_survival_NA[i]]<-TCGA_GBM_Meta$days_to_last_follow_up[index_survival_NA[i]]}



TCGA_GBM_tpm <- TCGAanalyze_Preprocessing(object = TCGA_GBM,
                                           cor.cut = 0.80,
                                           datatype = "tpm_unstrand")# "fpkm_unstrand"
# both tpm and fpkm version was used 
rownames(TCGA_GBM_tpm)<- ens_gene_id
TCGA_GBM_tpm<-t(TCGA_GBM_tpm)
infile<-cbind(infile,TCGA_GBM_tpm)
# that part for TCGA
index_survival_NA<-which(is.na(infile$survival_time))
infile<-infile[-index_survival_NA, ]

index_Pur_NA<-which(is.na(infile$purity))
infile<-infile[-index_Pur_NA, ]
#
index_status_NR<-which(infile$status=="Not Reported")
infile<-infile[-index_status_NR, ]

infile$status <- gsub('alive', 0, infile$status, ignore.case = T)
infile$status <- gsub('dead', 1, infile$status, ignore.case = T)
#
infile<-infile[infile$definition=="Primary solid Tumor",]
table(infile$status) 
################################################################################
# data arrangement for KM
exp_index <- min(grep('ENSG', colnames(infile)))
#cli_data_all <- infile[, 1:exp_index-1]
cli_data_all <- infile[, c(1, 5, 6)]
cli_data <- cli_data_all[, c('sample_id', 'status', 'survival_time')]

exp_data <- infile[, exp_index:ncol(infile)]
rownames(exp_data) <- infile$sample_id
### delete donors with NA tpm
ind <- apply(exp_data, 1, function(x) all(is.na(x)))
exp_data <- exp_data[!ind, ]
print(paste(cancer_name, 'has ', length(rownames(exp_data)[ind]), 'sample(s) with all NA TPM'))
print(rownames(exp_data)[ind])
sample_keep <- rownames(exp_data)
### replace missing expression value with 0
#exp_data[is.na(exp_data)] <- 0    ### ------- replace missing expression value with 0
### keep TPM >= 1
rownames(exp_data) <- sample_keep
exp_data<- apply(exp_data, 2, as.numeric)
exp_keep <- as.data.frame(exp_data[, which(colMeans(exp_data)>=1)]) ## keep expressed gene
rownames(exp_keep)<-sample_keep
infile <- merge(cli_data, exp_keep, by.x = 'sample_id', by.y = 0)
columns <- colnames(infile)
first_gene_col <- 4
gene_list <- colnames(infile)[first_gene_col:ncol(infile)]
## -------------------------end of infile format--------------------------------

path_out_pdf <- paste0(outPath, cancer_name, '/output_pdf/')
dir.create(path_out_pdf, recursive = TRUE)
path_out_stat <- paste0(outPath, cancer_name, '/output_stat/')
dir.create(path_out_stat, recursive = TRUE)

p_list<-NULL
cutoff_list<-NULL
coef_list<-NULL
for (j in 1:length(gene_list)){
  print(j)
  output<-GenerateSurvInputfromTCGA(gene_list[j],cancer_name,infile)
  setwd(path_out_pdf)
  result<-Cheng_generateKMplot(output,outFile=paste0(cancer_name,"_",gene_list[j]))
  log_rank_p<-as.matrix(result$logRankP)
  cut_off<-as.matrix(result$EXPcut)
  coef<-as.matrix(result$coef)
  
  p_list<-rbind(p_list,log_rank_p)
  cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,cut_off,coef)
}

final_result<-cbind(as.matrix(gene_list),coef_list,cutoff_list,p_list)
colnames(final_result)<-c("ensemble","coef","cutoff","p")
final_result <- as.data.frame(final_result)
final_result$FDR_BH <- p.adjust(final_result$p, method = "BH")
final_result2<-final_result
write.csv(final_result2, "final_result2_KM_TCGA.csv")
exp_data <- infile[, first_gene_col:ncol(infile)]
exp_data<- apply(exp_data, 2, as.numeric)
rownames(exp_data) <- infile$sample_id
median_EXP <- as.data.frame(colMeans(exp_data))
colnames(median_EXP) <- 'mean_exp'

final_result <- merge(final_result, median_EXP, by.x = 'ensemble', by.y = 0)
final_result <- final_result %>% relocate(mean_exp, .after=ensemble)
colnames(final_result)[1] <- "ensemble_gene_id"
setwd(path_out_stat)
write.table(final_result,file=paste0(cancer_name,"_KM_stat_result.txt"),sep="\t",row.names=F,col.names=T,quote=F)
print(dim(final_result))
print(paste0('Finsing Kaplan-Meier for ', cancer_name, '_TCGA'))
print(paste0(cancer_name, ' including ', length(sample_keep), ' donors, ', 
             'with ', length(gene_list), ' expressed genes'))

Order_Final_result2<-final_result[final_result$p<0.05, ]