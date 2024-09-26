
############# Jaccard Function
Jaccard <- function(module_list1, module_list2) {
  i_len <- length(intersect(module_list1, module_list2))
  u_len <- length(union(module_list1, module_list2))
  return(i_len/u_len)
}
##########################################
# install.packages("vegan")
# library(vegan)
library("RColorBrewer")
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
library(VennDiagram)

########################## TCGA_vs_CGGA_693

nt<-result_end_TCGA_vs_CGGA_693$Module_No_TCGA
modules_TCGA<-unique(nt)
nt<-length(modules_TCGA)

nc<-result_end_TCGA_vs_CGGA_693$Module_No_693
modules_CGGA_693<-unique(nc)
nc<-length(modules_CGGA_693)

out_Jacard_TCGA_CGGA_693 <- matrix(0, nrow = nt, ncol = nc, dimnames = list(modules_TCGA, modules_CGGA_693))
for (i in 1:nt) {
  for (j in 1:nc) {
    modulnum_t<-as.numeric(modules_TCGA[i])
    GeneList_t<-moduleList_TCGA[modulnum_t][[1]]
    modulnum_c<-as.numeric(modules_CGGA_693[j])
    GeneList_c<-moduleList_CGGA_693[modulnum_c][[1]]
    
    out_Jacard_TCGA_CGGA_693[i, j] <- Jaccard(GeneList_t, GeneList_c)
  }
}

png("HighCC Modules of TCGA and CGGA_693.png", width = 6, height = 6, units = 'in', res = 600)
pheatmap::pheatmap(out_Jacard_TCGA_CGGA_693, display_numbers = matrix(ifelse(out_Jacard_TCGA_CGGA_693>0.03, "*", ""), nrow(out_Jacard_TCGA_CGGA_693)), main = "HighCC Modules of TCGA and CGGA_693",
                   cellwidth = 30, cellheight = 20, fontsize = 16, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                   show_rownames = T, show_colnames = T)
dev.off()

########################## TCGA_vs_CGGA_325

nt<-result_end_TCGA_vs_CGGA_325$Module_No_TCGA
modules_TCGA<-unique(nt)
nt<-length(modules_TCGA)

nc<-result_end_TCGA_vs_CGGA_325$Module_No_325
modules_CGGA_325<-unique(nc)
nc<-length(modules_CGGA_325)

out_Jacard_TCGA_CGGA_325 <- matrix(0, nrow = nt, ncol = nc, dimnames = list(modules_TCGA, modules_CGGA_325))
for (i in 1:nt) {
  for (j in 1:nc) {
    modulnum_t<-as.numeric(modules_TCGA[i])
    GeneList_t<-moduleList_TCGA[modulnum_t][[1]]
    modulnum_c<-as.numeric(modules_CGGA_325[j])
    GeneList_c<-moduleList_CGGA_325[modulnum_c][[1]]
    
    out_Jacard_TCGA_CGGA_325[i, j] <- Jaccard(GeneList_t, GeneList_c)
  }
}

png("HighCC Modules of TCGA and CGGA_325.png", width = 6, height = 6, units = 'in', res = 600)
pheatmap::pheatmap(out_Jacard_TCGA_CGGA_325, display_numbers = matrix(ifelse(out_Jacard_TCGA_CGGA_325>0.03, "*", ""), nrow(out_Jacard_TCGA_CGGA_325)), main = "HighCC Modules of TCGA and CGGA_325",
                   cellwidth = 30, cellheight = 20, fontsize = 16, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                   show_rownames = T, show_colnames = T)
dev.off()

########################## CGGA_693_vs_CGGA_325

nt<-result_end_CGGA_693_vs_CGGA_325$Module_No_CGGA_693
modules_CGGA_693<-unique(nt)
nt<-length(modules_CGGA_693)

nc<-result_end_CGGA_693_vs_CGGA_325$Module_No_325
modules_CGGA_325<-unique(nc)
nc<-length(modules_CGGA_325)

out_Jacard_CGGA_693_CGGA_325 <- matrix(0, nrow = nt, ncol = nc, dimnames = list(modules_CGGA_693, modules_CGGA_325))
for (i in 1:nt) {
  for (j in 1:nc) {
    modulnum_t<-as.numeric(modules_CGGA_693[i])
    GeneList_t<-moduleList_CGGA_693[modulnum_t][[1]]
    modulnum_c<-as.numeric(modules_CGGA_325[j])
    GeneList_c<-moduleList_CGGA_325[modulnum_c][[1]]
    
    out_Jacard_CGGA_693_CGGA_325[i, j] <- Jaccard(GeneList_t, GeneList_c)
  }
}

png("HighCC Modules of CGGA_693 and CGGA_325.png", width = 6, height = 6, units = 'in', res = 600)
pheatmap::pheatmap(out_Jacard_CGGA_693_CGGA_325, display_numbers = matrix(ifelse(out_Jacard_CGGA_693_CGGA_325>0.03, "*", ""), nrow(out_Jacard_CGGA_693_CGGA_325)), main = "HighCC Modules of CGGA_693 and CGGA_325",
                   cellwidth = 30, cellheight = 20, fontsize = 16, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                   show_rownames = T, show_colnames = T)
dev.off()

