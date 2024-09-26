############# SURVIVAL analysis for TCGA, CGGA_325 and CGGA_693  ###############

library(survival)
library(tidyverse)
library(imputeTS)

project='GBM_TCGA'
# set folder sot results
raw_path <- "{PATH}/Cox/"
path_out <- "{PATH}/Cox/"
setwd(path_out)
# 
# infile_CGGA_325_N<-t(CGGA_325_GraterRSEM1)
# #rownames(CGGA_693_Primary_clinic)<-CGGA_693_Primary_clinic$CGGA_ID
# infile_CGGA_325_N<-cbind(CGGA_325_Primary_clinic[,c(1, 8, 7)],infile_CGGA_325_N)
# save(file="{PATH}/infile_CGGA_325_N.Rdata",infile_CGGA_325_N)

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
# save(file="{PATH}/Order_Final_result2_FPKM.Rdata",Order_Final_result2)
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
