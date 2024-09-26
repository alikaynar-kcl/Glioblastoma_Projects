readRDS(high_tps, file = "{PATH}/Data/high_tps.rds")
readRDS(low_tps, file = "{PATH}/Data/low_tps.rds")
readRDS(mid_tps, file = "{PATH}/Data/mid_tps.rds")
readRDS(TCGA_NT, file = "{PATH}/Data/TCGA_NT.rds")

TCGA_GBM<-cbind(TCGA_NT, high_tps, low_tps, mid_tps)