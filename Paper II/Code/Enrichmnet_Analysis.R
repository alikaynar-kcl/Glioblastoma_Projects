## Go Term enrichment for significantly altered genes
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(biomaRt)
source('{PATH}/deg_GoTerm_clusterProfiler.R')
# Interested gene list
## mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
# mart <- useMart("ensembl", host = "https://www.ensembl.org", dataset = "hsapiens_gene_ensembl")

# Set comparison and other initial parameters
Comparison <- "NT vs TR" # Define it accordingly
ont1 <- "BP"
heder_T <- paste0("Highlighting Significant GO Terms ",ont1, "; ", Comparison)
# Perform GO enrichment
background <- as.matrix(deg_results_df$gene_name) # deg_results_df is DEG  result table
int_gene <- as.matrix(deg_results_df1$gene_name) # Check if conversion to matrix is necessary; typically, a vector is used

pathway <- enrichGO(gene = int_gene, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = ont1, universe = background, pAdjustMethod = "BH", qvalueCutoff = 0.05)

# GeneExp<-NULL
# for (i in 1:length(pathway@result[["BgRatio"]])){
#   gene_ratio <- pathway@result[["GeneRatio"]][i]
#   s_gene_ratio <- as.numeric(strsplit(gene_ratio, "/")[[1]])
#   Bg_ratio <- pathway@result[["BgRatio"]][i]
#   s_Bg_ratio <- as.numeric(strsplit(Bg_ratio, "/")[[1]])
#   GeneExp[i]<-(s_Bg_ratio[1]*s_gene_ratio[2])/s_Bg_ratio[2]
# }
# pathway@result[["GeneExp"]]<-GeneExp

path_table_all<-deg_GoTerm_clusterProfiler(pathway)

# Generate and save GO enrichment plot and table
png(paste(heder_T, ".png", sep=" "), width = 12, height = 12, units = 'in', res = 300)
par(cex.main=3,cex.axis=1.5, cex.names=1.5)#cex.axis=1, cex.lab=1, , cex.sub=1
barplot(pathway, title=heder_T,x = "Count", color = "p.adjust",font.size = 12,showCategory=20)
dev.off()
# 
# Prepare plot file
png(paste(heder_T, ".png", sep=""), width = 12, height = 12, units = 'in', res = 300)

# Adjust margins (specifically increase the left margin for y-axis labels)
# par(mar=c(bottom, left, top, right))
# Default is usually c(5, 4, 4, 2) + 0.1
# increase the left margin for more space for y-axis variable (long names)
par(mar=c(5, 7, 4, 2) + 0.1)

barplot(pathway, main=heder_T, xlab = "Count", color = "p.adjust", cex.names = 1.5, las = 2, cex.axis = 1.5, cex.main = 3, showCategory=20)

# Close the plotting device
dev.off()

png(filename = paste(heder_T, ".png"), width = 12, height = 12, units = 'in', res = 300)
barplot(pathway, showCategory = 20, title = heder_T, col = rainbow(10), font.size = 12)
dev.off()

write.table(path_table_all, file = paste(heder_T,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Generate a dot plot for GO enrichment
pdf(file = paste("Dot_plot_GO_", Comparison, ".pdf"), width = 10, height = 8)
dotplot(pathway, showCategory = 20, split = ".significance", font.size = 12)
dev.off()

# For KEGG Enrichment, convert gene symbols to entrez ids
selected_rows <- deg_results_df[deg_results_df$padj < 1e-10 & abs(deg_results_df$log2FoldChange) > 2, ]

# Use the biomaRt package to map gene symbols to Entrez IDs
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
int_gene_kegg <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = rownames(selected_rows), mart = mart)
background2 <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = background, mart = mart)

# Perform KEGG enrichment
ekegg <- enrichKEGG(gene = int_gene_kegg$entrezgene_id, organism = "hsa", keyType = "entrez", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background2$entrezgene_id, qvalueCutoff = 0.05)

# Generate and save KEGG enrichment plot
png(filename = paste("Bar_plot_KEGG_", Comparison, ".png"), width = 12, height = 12, units = 'in', res = 300)
barplot(ekegg, showCategory = 20, title = paste("KEGG Enrichment", Comparison), col = rainbow(10), font.size = 12)
dev.off()

write.table(as.data.frame(ekegg), file = paste("KEGG_Enrichment_", Comparison, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# generate a dot plot for KEGG enrichment
pdf(file = paste("Dot_plot_KEGG_", Comparison, ".pdf"), width = 10, height = 8)
dotplot(ekegg, showCategory = 20, font.size = 12)
dev.off()