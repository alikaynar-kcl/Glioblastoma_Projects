# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(readxl)

# Data extraction from the TCGA_GBM object
TCGA_TP<-cbind(high_tps, low_tps, mid_tps)
TCGA_TP_text <- TCGA_TP@assays@data@listData[["tpm_unstrand"]]

file_path <- "{PATH}/Supplementary File B.xlsx"  
df_figure_8 <- read_excel(file_path, sheet = "Figure_8", col_names = TRUE)

colnames(TCGA_TP_text) <- TCGA_TP@colData@rownames
rownames(TCGA_TP_text) <- TCGA_TP@rowRanges@elementMetadata@listData[["gene_name"]]
TCGA_TP_text <- TCGA_TP_text[rownames(TCGA_TP_text) %in% df_figure_8$Gene, ] # We focused on our significant gene list

# Sample information data extraction
listData <- TCGA_TP@colData@listData
TCGA_TP_text_Clin <- data.frame(
  barcode = as.character(listData[["barcode"]]),
  gender = as.character(listData[["gender"]]),
  paper_IDH_status = as.character(listData[["paper_IDH.status"]]),
  MGMT_promoter_status = as.character(listData[["paper_MGMT.promoter.status"]]),
  Chr7_gain_Chr10_loss = as.character(listData[["paper_Chr.7.gain.Chr.10.loss"]]),
  Chr19_20_co_gain = as.character(listData[["paper_Chr.19.20.co.gain"]]),
  TERT_promoter_status = as.character(listData[["paper_TERT.promoter.status"]]),
  ATRX_status = as.character(listData[["paper_ATRX.status"]]),
  Telomere_Maintenance = as.character(listData[["paper_Telomere.Maintenance"]])
)

# Transposing gene expression data for easier merging with TCGA_TP_text_Clin
expression_data <- t(TCGA_TP_text)
colnames(expression_data) <- rownames(TCGA_TP_text)
rownames(expression_data) <- NULL

# make sure expression data to be numeric
expression_data <- apply(expression_data, 2, function(x) as.numeric(as.character(x)))
rownames(expression_data) <- colnames(TCGA_TP_text)

# Merge with TCGA_TP_text_Clin info
data_merged <- merge(TCGA_TP_text_Clin, as.data.frame(expression_data), by.x = "barcode", by.y = "row.names", all.x = TRUE)

# Function to perform ANOVA for each gene across all measured features in TCGA_TP_text_Clin
perform_anova <- function(data, gene_name) {
  p_values <- c()
  for (clin_fea in names(data)[2:9]) {  # first column is 'barcode'
    # temporary dataset without NA values for the Sample information: TCGA_TP_text_Clin
    temp_data <- data %>% filter(!is.na(.[[clin_fea]]))
    clean_data <- temp_data[, c("barcode", clin_fea, gene_name), drop = FALSE]
    clean_data <- na.omit(clean_data)
    
    if (length(unique(clean_data[[clin_fea]])) > 1 && any(table(clean_data[[clin_fea]]) > 1)) {
      aov_result <- aov(as.formula(paste(gene_name, "~", clin_fea)), data = clean_data)
      p_value <- summary(aov_result)[[1]]["Pr(>F)"][1]
      p_values <- c(p_values, p_value)
    } else {
      p_values <- c(p_values, NA)
    }
  }
  names(p_values) <- names(data)[2:9]
  return(p_values)
}

# Apply ANOVA function across all genes and compile results into a data frame
anova_results <- t(sapply(colnames(data_merged)[10:ncol(data_merged)], function(gene) perform_anova(data_merged, gene)))
anova_results_df <- as.data.frame(anova_results)# Convert results to data frame

results_df <- matrix(ncol = dim(anova_results_df)[2], nrow = dim(anova_results_df)[1])

# for loops to fix anova_results_df by extracting only the first element from each list
for (j in 1:dim(anova_results_df)[1]) {
  for (i in 1:dim(anova_results_df)[2]){
    results_df[j, i] <- anova_results_df[i][[1]][[j]][1]
  }
}
colnames(results_df)<-colnames(anova_results_df)
rownames(results_df)<-rownames(anova_results_df)

Anova_Marker_Genes <- merge(data.frame(results_df),HPA_Info, by = 0, all.x = TRUE)
write.xlsx(as.data.frame(Anova_Marker_Genes), "{PATH}/Anova_Marker_Genes.xlsx")

# Function to create boxplots for a given gene
create_gene_plots <- function(gene_name, data_merged, results_df) {
  plots <- list()
  clin_feas <- colnames(results_df)
  
  for (i in 1:length(clin_feas)) {
    clin_fea <- clin_feas[i]
    p_value <- results_df[gene_name, clin_fea]
    
    # temporary dataset without NA values
    temp_data <- data_merged %>% filter(!is.na(.[[clin_fea]]))
    clean_data <- temp_data[, c("barcode", clin_fea, gene_name), drop = FALSE]
    clean_data <- na.omit(clean_data)
    
    # Generate the boxplot with colorful elements
    p <- ggplot(clean_data, aes_string(x = clin_fea, y = gene_name, fill = clin_fea)) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Set3") +  #a color palette
      labs(subtitle = paste("p-value:", format(p_value, digits = 4)),
           y = "TPM") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 10, hjust = 1))
    
    plots[[i]] <- p
  }
  return(plots)
}

# Example genes to plot
gene_to_plot <- c("SOX8", "NOVA1", "DSCAM", "LDLRAD3", "DPYSL4", "PCDHB9", "FGFBP3", "SALL3")

# Define the width ratios of the plotes for better visualisation
width_ratios <- c(1, 1, 1, 1.5, 1.5, 1, 1, 1)

# Create a list to store plots for each gene
all_gene_plots <- list()

# Loop through each gene and create plots
for (gene in gene_to_plot) {
  gene_plots <- create_gene_plots(gene, data_merged, results_df)
  combined_plots <- wrap_plots(gene_plots, ncol = 8, widths = width_ratios) + 
    plot_annotation(title = gene) +
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +  # Adjust plot margins
    theme(panel.spacing = unit(0.1, "cm"))  # Tiny space between plots
  all_gene_plots[[gene]] <- combined_plots
}

# Combine all gene plots into one figure
final_figure <- wrap_plots(all_gene_plots, ncol = 1) + 
  plot_layout(guides = 'collect')

# Save the final combined figure with specified dpi
ggsave(paste0("Marker_Anova_Genes_8",".tiff"), final_figure, width = 18, height = 24, dpi = 300)

print(final_figure)