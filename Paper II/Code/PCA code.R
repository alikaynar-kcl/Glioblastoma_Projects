cancer_name<- "CGGA_693"
ddsSE <- CGGA_693_DESeq # DESeq object
###############################################################################
# Prepare data object for DESeq
ddsSE$definition<- relevel(ddsSE$definition, ref = "Solid Tissue Normal")
keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]
#Cook_Outliers<-assays(ddsSE)[["cooks"]]
#boxplot(log10(assays(ddsSE)[["cooks"]]), range=0, las=2)
#######################################
vst_d<-vst(ddsSE)
################################################################################
# PCA
path_out_2 <- paste0(outPath, cancer_name, '/')
setwd(path_out_2)
tiff(paste(cancer_name,"PCA.tiff", sep="_"), width = 8, height = 8, units = 'in', res = 300)
plotPCA(vst_d, intgroup=c("definition"))
dev.off()
