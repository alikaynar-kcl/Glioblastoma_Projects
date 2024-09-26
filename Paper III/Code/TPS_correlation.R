# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stats)
#library(p.adjust)
library(ggplot2)
library(ggpubr)
library(openxlsx)

TPS_GeneExp<-rbind(TCGA_TP@colData@listData[["paper_ABSOLUTE.purity"]],TCGA_TP@assays@data@listData[["tpm_unstrand"]])
colnames(TPS_GeneExp)<-TCGA_TP@colData@rownames
rownames(TPS_GeneExp)<-c("TPS_Score",TCGA_TP@rowRanges@elementMetadata@listData[["gene_name"]])

# define dataset
data <- as.data.frame(TPS_GeneExp)

# Separate TPS_Score and gene expression data
tps_score <- as.numeric(data[1,]) # make sure TPS_Scores be a numeric
gene_exp <- data[-1, ] # Exclude the TPS_Scores

exp_keep <- as.data.frame(gene_exp[which(rowMeans(gene_exp)>=1), ]) ## keep expressed gene
# Transpose gene_exp to make genes as columns for easier application
gene_exp_t <- t(as.matrix(exp_keep)) 
# Calculate Spearman's correlation
spearman_results <- apply(gene_exp_t, 2, function(gene) cor.test(tps_score, as.numeric(gene), method = "spearman"))
spearman_correlations <- sapply(spearman_results, function(x) x$estimate)
spearman_p_values <- sapply(spearman_results, function(x) x$p.value)

# Calculate  also Pearson's correlation (But spearman results used for further analysis)
pearson_results <- apply(gene_exp_t, 2, function(gene) cor.test(tps_score, as.numeric(gene), method = "pearson"))
pearson_correlations <- sapply(pearson_results, function(x) x$estimate)
pearson_p_values <- sapply(pearson_results, function(x) x$p.value)

# Create a dataframe for both correlations
cor_df <- data.frame(Gene = colnames(gene_exp_t), 
                     Spearman_Correlation = spearman_correlations, 
                     Spearman_P_Value = spearman_p_values, 
                     Pearson_Correlation = pearson_correlations, 
                     Pearson_P_Value = pearson_p_values)

# Adjust p-values for multiple testing using FDR for both Spearman and Pearson
cor_df$Adjusted_Spearman_P_Value <- p.adjust(cor_df$Spearman_P_Value, method = "fdr")
cor_df$Adjusted_Pearson_P_Value <- p.adjust(cor_df$Pearson_P_Value, method = "fdr")

# Filter for significant correlations (you can adjust this based on which p-values you're interested in)
significant_correlations_spearman <- cor_df %>% filter(Adjusted_Spearman_P_Value < 0.01)
significant_correlations_pearson <- cor_df %>% filter(Adjusted_Pearson_P_Value < 0.01)

# HPA_Info contains gene information from Human Protein Atlas
TPS_Comp_List<- cor_df; rownames(TPS_Comp_List)<-TPS_Comp_List$Gene
TPS_Comp_List <- merge(TPS_Comp_List,HPA_Info, by = 0, all.x = TRUE)

write_xlsx(TPS_Comp_List, "~/{PATH}/TPS_Comp_List.xlsx")

# TCGA_HPA_INFO <- merge(Model_input_T,HPA_Info, by = 0, all.x = TRUE)
# TCGA_HPA_INFO contains mean value of each group and also HPA information