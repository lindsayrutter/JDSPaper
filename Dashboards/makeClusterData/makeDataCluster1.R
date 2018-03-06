library(bigPint)

ID <- readRDS("../../Bioinformatics/Pictures/FilterNotSig/Clustering_data_FDR_05/Sig_4_1.Rds")
load("../../Bioinformatics/Pictures/FilterNotSig/soybean_ir_noFilt_metrics.rda")

allMetrics = soybean_ir_noFilt_metrics[["N_P"]]

metricsCluster <- allMetrics[which(allMetrics$ID %in% ID),]
metricsCluster$ID <- as.character(metricsCluster$ID)

metrics <- list()
metrics[["N_P"]] <- metricsCluster

save(metrics, file = "cluster1_metrics.rda")
