library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(forcats) 

# Selected_genes is a vector containing the genes we're interested in
selected_genes <- c("KYNU","IDO1", "IDO2", "AFMID", "KMO", "HAAO", "QPRT", "NMNAT1", "NMNAT2", "NMNAT3", "NADSYN1", "NRK", "NAMPT", "TDO2")
Title_of_Interested_Genes<- "Essential Markers ATP salvage Hypoxanthine"
ncol_S=2;nrow_S=8

# Preparing the data frame with selected genes
selected_gene_correlations <- cor_df %>% # cor_df is a data frame containing TPS correlation results
  filter(Gene %in% selected_genes)

plots <- list() # The list to store plots

# Generate plots for each selected genes. This section include two separate plot then merge them into one.
for (i in 1:nrow(selected_gene_correlations)) {
  gene_name <- selected_gene_correlations$Gene[i]
  plot_data <- data.frame(TPS_Score = tps_score, Gene_Expression = as.numeric(gene_exp[gene_name, ]))
  
# Classify TPS score into High, Low, and remained, then assign colors
  plot_data$Color <- ifelse(plot_data$TPS_Score <= tps_quantile_25, 'blue', 
                            ifelse(plot_data$TPS_Score >= tps_quantile_75, 'red', 'black'))
  
# Rearrange positions of correlation info as label
  stat_cor_label_pos_x <- mean(range(plot_data$TPS_Score))
  stat_cor_label_pos_y <- max(plot_data$Gene_Expression, na.rm = TRUE) - 0.05 * diff(range(plot_data$Gene_Expression, na.rm = TRUE))
  
# Generate plot
  gene_name
  p1 <- ggplot(plot_data, aes(x = TPS_Score, y = Gene_Expression)) +
    geom_point(aes(color = Color)) +  
    scale_color_identity() + # Use the actual colors specified in the Color column
    geom_smooth(method = "lm", se = FALSE, color = "gray") + # Linear regression line
    stat_cor(method = "spearman", label.x = stat_cor_label_pos_x, label.y = stat_cor_label_pos_y) + # Spearman correlation label
    labs(title = gene_name, x = "TPS Score", y = "Gene Expression") +
    theme_minimal(base_size = 12) + 
    theme(
      plot.background = element_rect(fill = "white", colour = NA), # white background for the plot
      panel.background = element_rect(fill = "white", colour = NA), # white background for the panel
      legend.position = "none" # Hide legend for cleaner plots
    )
  
  # Extracting data for gene: C_H_L_TPM contain gene expression info for all genes (NT+high_tps+low_tps)
  gene_data2 <- data.frame(barcode = colnames(C_H_L_TPM), expression = C_H_L_TPM[gene_name, ])
  # Merging gene data with condition information
  merged_data2 <- merge(gene_data2, NT_high_low, by = "barcode")
  # Transforming data for plotting
  plot_data2 <- pivot_longer(merged_data2, cols = expression, names_to = "Gene", values_to = "Expression")
  
  plot_data2$condition <- fct_relevel(plot_data2$condition, "NT", "Low", "High")
  gene_name
  
  # Extract p-values for the gene from both DEG result data frames
  padj_low <- DESeq_Table_dds_NT_low_2$padj[DESeq_Table_dds_NT_low_2$gene_name == gene_name]
  padj_high <- DESeq_Table_dds_NT_high_2$padj[DESeq_Table_dds_NT_high_2$gene_name == gene_name]
  
  # Format the p-values to one decimal place
  formatted_padj_low <- sprintf("%.1e", padj_low)
  formatted_padj_high <- sprintf("%.1e", padj_high)
  
  # Create the two-layer title text by combining the descriptions and p-values
  title_text <- paste("DEG p-adj value", "\n", formatted_padj_low, "-", formatted_padj_high)
  
  # Plotting second piece with specific fill colors
  p2 <- ggplot(plot_data2, aes(x = condition, y = Expression, fill = condition)) +
    geom_boxplot() +
    scale_fill_manual(values = c("NT" = "green", "Low" = "blue", "High" = "red")) +
    scale_x_discrete(name = "Condition") +
    scale_y_continuous(name = "TPM Value") +
    #  ggtitle(gene_name) +
    
    ggtitle(title_text) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, family = "Times", color = "black",hjust = 0.5))
  
  plots[[i]]<- p1 + p2 + plot_layout(widths = c(4, 1)) # p1 and p2 portion 4 : 1
}

# Combine all plots into one figure using ggarrange
figure <- ggarrange(plotlist = plots,
                    ncol = ncol_S, nrow = nrow_S)
figure

# save the figure as tiff format
ggsave(paste0(Title_of_Interested_Genes,".tiff"), figure, width = 15, height = 24, dpi = 300)
#ggsave(paste0(Title_of_Interested_Genes,".tiff"), figure, width = 15, height = 8, dpi = 300)
