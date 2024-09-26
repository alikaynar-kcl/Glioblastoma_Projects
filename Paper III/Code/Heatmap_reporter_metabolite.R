library(readxl)
library(dplyr)
library(pheatmap)
library(grid)
library(gridExtra)
library(ggplot2)
library(viridis)  # For color-blind friendly palettes

# Define the metabolite IDs to be excluded, without their last character
Removed_Mets <- c('MAM02040', 'MAM01596', 'MAM02631', 'MAM02039', 'MAM02046', 'MAM02519', 
                  'MAM01597', 'MAM02751', 'MAM02759', 'MAM01334', 'MAM01285', 'MAM01371', 
                  'MAM02681', 'MAM02682', 'MAM01802', 'MAM01803', 'MAM02552', 'MAM02553', 
                  'MAM02554', 'MAM02555')

# Set the file path
file_path <- "{PATH}/Reporter_Met_Res.xlsx" # read the reporter metabolite tables

# Function to preprocess each sheet and remove unwanted metabolites
preprocess_sheet <- function(sheet_name, status_filter, pvalue_threshold, removed_ids) {
  df <- read_excel(file_path, sheet = sheet_name) %>%
    mutate(metID = gsub("[a-z]$", "", mets)) %>%  # Remove the last character
    filter(Status == status_filter, metPValues < pvalue_threshold, !metID %in% removed_ids) %>%
    select(metNames, metZScores, metPValues)
  
  df
}

# Status variable
Statu <- "UP" # Up regulated gene related reporter metabolites

# Reading and preprocessing each sheet according to specific filters
nt_vs_tp <- preprocess_sheet("NT_vs_TP", Statu, 0.05, Removed_Mets)
nt_vs_high <- preprocess_sheet("NT_vs_High", Statu, 0.05, Removed_Mets)
nt_vs_low <- preprocess_sheet("NT_vs_Low", Statu, 0.05, Removed_Mets)

# Combining the unique metNames from all datasets
combined_data <- data.frame(metNames = unique(c(nt_vs_tp$metNames, nt_vs_high$metNames, nt_vs_low$metNames)),
                            Z_Score_TP = 0, Z_Score_High = 0, Z_Score_Low = 0,
                            P_Value_TP = 1, P_Value_High = 1, P_Value_Low = 1)

nt_vs_tp_Ori <- read_excel(file_path, sheet = "NT_vs_TP") %>% filter(Status == Statu)
nt_vs_high_Ori <- read_excel(file_path, sheet = "NT_vs_High") %>% filter(Status == Statu)                           
nt_vs_low_Ori <- read_excel(file_path, sheet = "NT_vs_Low") %>% filter(Status == Statu)    

# Function to update Z-scores and p-values in combined_data based on original datasets (preserve original file)
update_Z_scores_and_p_values <- function(combined_df, original_df, score_column, pvalue_column) {
  combined_df[score_column] <- mapply(function(met) {
    original_score <- original_df %>% 
      filter(metNames == met) %>%
      pull(metZScores)
    
    if (length(original_score) == 0) {  # No match found in original data
      return(0)
    } else {
      return(original_score)
    }
  }, met = combined_df$metNames)
  
  combined_df[pvalue_column] <- mapply(function(met) {
    original_pvalue <- original_df %>% 
      filter(metNames == met) %>%
      pull(metPValues)
    
    if (length(original_pvalue) == 0) {  # No match found in original data
      return(1)
    } else {
      return(original_pvalue)
    }
  }, met = combined_df$metNames)
  
  combined_df
}

# Update Z-scores and p-values for each condition
combined_data <- update_Z_scores_and_p_values(combined_data, nt_vs_tp_Ori, "Z_Score_TP", "P_Value_TP")
combined_data <- update_Z_scores_and_p_values(combined_data, nt_vs_high_Ori, "Z_Score_High", "P_Value_High")
combined_data <- update_Z_scores_and_p_values(combined_data, nt_vs_low_Ori, "Z_Score_Low", "P_Value_Low")

# Prepare the matrix and generate the heatmap
heatmap_matrix <- combined_data %>%
  select(metNames, Z_Score_TP, Z_Score_High, Z_Score_Low) %>%
  column_to_rownames(var = "metNames") %>%
  as.matrix()

# Determine the range of positive and negative Z-scores for rescaling the color
max_z <- max(heatmap_matrix, na.rm = TRUE)
min_z <- min(heatmap_matrix, na.rm = TRUE)

# Generate breaks for positive and negative Z-scores separately
breaks <- c(seq(min_z, 0, length.out = 50), seq(0, max_z, length.out = 51)[-1])

# Generate color palette: 'magma' for negative and 'viridis' for positive
palette <- c(viridis::magma(50, direction = -1), viridis::viridis(50, direction = 1))

# Function to create significance labels based on p-values
create_significance_labels <- function(pvalues) {
  labels <- matrix("", nrow = nrow(pvalues), ncol = ncol(pvalues))
  labels[pvalues > 0.01 & pvalues <= 0.05] <- "*"
  labels[pvalues > 0.001 & pvalues <= 0.01] <- "**"
  labels[pvalues <= 0.001] <- "***"
  return(labels)
}

# Create a matrix of p-values
pvalue_matrix <- combined_data %>%
  select(metNames, P_Value_TP, P_Value_High, P_Value_Low) %>%
  column_to_rownames(var = "metNames") %>%
  as.matrix()

# Generate the significance labels
significance_labels <- create_significance_labels(pvalue_matrix)

# Adjusting the heatmap plot with specific cellwidth and cellheight
heatmap_plot <- pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = palette,
  breaks = breaks,  # Adjust the color scale based on actual range of Z-scores
  fontsize_row = 5,  # Adjust font size for row names
  fontsize_col = 5,  # Adjust font size for column names
  fontsize = 10,     # Adjust overall font size
  cellwidth = 15,    # Adjust cell width to make columns narrower
  cellheight = 5,    # Adjust cell height to make rows narrower if needed
  display_numbers = significance_labels, # Add significance labels
  number_color = "black" # Color of the significance labels
)

# Save the heatmap with adjusted size
ggsave("heatmap_Reporter_Map_Up.png", plot = heatmap_plot$gtable, width = 12, height = 44, dpi = 600)
