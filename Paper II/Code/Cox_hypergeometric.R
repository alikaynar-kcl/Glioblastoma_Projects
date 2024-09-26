#rm(list=ls())
# deg_up_1,
# deg_down_1,
# deg_up_2,
# deg_down_2,
# bg_gene

source(paste0(outPath, 'DEGs_list_overlap.R'))
################################################################################
###### COX
# CGGA_Cox res table

Favorable_CGGA_Cox<-NULL
UnFavorable_CGGA_Cox<-NULL  
for (i in 1:nrow(CGGA_Cox)){
  if (CGGA_Cox$p[i] <0.05 && CGGA_Cox$coef[i] <0){Favorable_CGGA_Cox<-rbind(Favorable_CGGA_Cox,CGGA_Cox[i,])}
  else if (CGGA_Cox$p[i] <0.05 && CGGA_Cox$coef[i] >0){UnFavorable_CGGA_Cox<-rbind(UnFavorable_CGGA_Cox,CGGA_Cox[i,])}
}

##### COX
TCGA_Cox

Favorable_TCGA_Cox<-NULL
UnFavorable_TCGA_Cox<-NULL  
for (i in 1:nrow(TCGA_Cox)){
  if (TCGA_Cox$p[i] <0.05 && TCGA_Cox$coef[i] <0){Favorable_TCGA_Cox<-rbind(Favorable_TCGA_Cox,TCGA_Cox[i,])}
  else if (TCGA_Cox$p[i] <0.05 && TCGA_Cox$coef[i] >0){UnFavorable_TCGA_Cox<-rbind(UnFavorable_TCGA_Cox,TCGA_Cox[i,])}
}
i
################################################################################


##### COX
# CGGA_325_Cox<-na.omit(CGGA_325_Cox)
Favorable_CGGA_325_Cox<-NULL
UnFavorable_CGGA_325_Cox<-NULL  
for (i in 1:nrow(CGGA_325_Cox)){
  if (CGGA_325_Cox$p[i] <0.05 && CGGA_325_Cox$coef[i] <0){Favorable_CGGA_325_Cox<-rbind(Favorable_CGGA_325_Cox,CGGA_325_Cox[i,])}
  else if (CGGA_325_Cox$p[i] <0.05 && CGGA_325_Cox$coef[i] >0){UnFavorable_CGGA_325_Cox<-rbind(UnFavorable_CGGA_325_Cox,CGGA_325_Cox[i,])}
}

################################################################################
###### COX
deg_up_1 <- as.matrix(UnFavorable_TCGA_Cox$ensembl)
#deg_up_1<-data.frame(deg_up_1)
deg_down_1<- as.matrix(Favorable_TCGA_Cox$ensembl)
#deg_down_1<-data.frame(deg_down_1)
deg_up_2 <- as.matrix(UnFavorable_CGGA_Cox$ensembl_gene_id)
#deg_up_2<-data.frame(deg_up_2)
deg_down_2<- as.matrix(Favorable_CGGA_Cox$ensembl_gene_id)
#deg_down_2<-data.frame(deg_down_2)
bg_gene <- as.matrix(intersect(TCGA_Cox$ensembl, CGGA_325_Cox$ensembl_gene_id))


TCGA_vs_CGGA_693_Cox<- DEGs_list_overlap(deg_up_1, deg_down_1, deg_up_2, deg_down_2, bg_gene)  # 

################################################################################
###### COX
deg_up_1 <- as.matrix(UnFavorable_TCGA_Cox$ensembl)
#deg_up_1<-data.frame(deg_up_1)
deg_down_1<- as.matrix(Favorable_TCGA_Cox$ensembl)
#deg_down_1<-data.frame(deg_down_1)
deg_up_2 <- as.matrix(UnFavorable_CGGA_325_Cox$ensembl_gene_id)
#deg_up_2<-data.frame(deg_up_2)
deg_down_2<- as.matrix(Favorable_CGGA_325_Cox$ensembl_gene_id)
#deg_down_2<-data.frame(deg_down_2)
bg_gene <- as.matrix(intersect(CGGA_325_Cox$ensembl_gene_id, TCGA_Cox$ensembl))


TCGA_vs_CGGA_325_Cox<- DEGs_list_overlap(deg_up_1, deg_down_1, deg_up_2, deg_down_2, bg_gene)  # 

###############################################################################################################
## Comparison of CGGA_693 vs CGGA_325
###### COX
deg_up_1 <- as.matrix(UnFavorable_CGGA_Cox$ensembl_gene_id)
#deg_up_1<-data.frame(deg_up_1)
deg_down_1<- as.matrix(Favorable_CGGA_Cox$ensembl_gene_id)
#deg_down_1<-data.frame(deg_down_1)
deg_up_2 <- as.matrix(UnFavorable_CGGA_325_Cox$ensembl_gene_id)
#deg_up_2<-data.frame(deg_up_2)
deg_down_2<- as.matrix(Favorable_CGGA_325_Cox$ensembl_gene_id)
#deg_down_2<-data.frame(deg_down_2)
bg_gene <- as.matrix(intersect(CGGA_325_Cox$ensembl_gene_id, CGGA_Cox$ensembl_gene_id))


CGGA_vs_CGGA_325_Cox<- DEGs_list_overlap(deg_up_1, deg_down_1, deg_up_2, deg_down_2, bg_gene)  # 

List_UN_3_Cox<-list(UnFavorable_CGGA_325_Cox$ensembl_gene_id,UnFavorable_TCGA_Cox$ensembl_gene_id,UnFavorable_CGGA_Cox$ensembl_gene_id)
List_UN_3_Cox_Common<-list(CGGA_vs_CGGA_325_Cox[["over_up"]],TCGA_vs_CGGA_325_Cox[["over_up"]],TCGA_vs_CGGA_693_Cox[["over_up"]])
List_UN_3_Cox_Common_Down<-list(CGGA_vs_CGGA_325_Cox[["over_down"]],TCGA_vs_CGGA_325_Cox[["over_down"]],TCGA_vs_CGGA_693_Cox[["over_down"]])
List_DGE_UP_Common<-list(rownames(gene_DEG_CGGA_325_UP),rownames(gene_DEG_TCGA_Count_D_UP),rownames(gene_DEG_CGGA_693_UP))
List_F_3_Cox<-list(Favorable_CGGA_325_Cox$ensembl_gene_id,Favorable_TCGA_Cox$ensembl_gene_id,Favorable_CGGA_Cox$ensembl_gene_id)

####################

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
library(VennDiagram)
#############################

venn.diagram(List_UN_3_Cox,
             category.names = c("CGGA_325","TCGA","CGGA_693"),
             filename = "COX_3_UNFAV_COM.png",
             output=TRUE,
             imagetype="png" ,
             height = 900 ,
             width = 900 ,
             resolution = 400,
             compression = "lzw",
             lwd = 1,
             lty = 'blank',
             fill = myCol,
             cex = .8,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.8,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)

venn.diagram(List_F_3_Cox,
             category.names = c("CGGA_325","TCGA","CGGA_693"),
             filename = "COX_FAV_COM.png",
             output=TRUE,
             imagetype="png" ,
             height = 900 ,
             width = 900 ,
             resolution = 400,
             compression = "lzw",
             lwd = 1,
             lty = 'blank',
             fill = myCol,
             cex = .8,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.8,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)
