library(bigPint)

ID <- readRDS("../../Bioinformatics/Pictures/Clustering_data/FocusOn/Sig_4_3.Rds")
load("../../Bioinformatics/Pictures/FilterNotSig/soybean_ir_noFilt_metrics.rda")

allMetrics = soybean_ir_noFilt_metrics[["N_P"]]
allMetrics$PValue = 1

metricsCluster <- allMetrics[which(allMetrics$ID %in% ID),]
metricsCluster$ID <- as.character(metricsCluster$ID)

metrics <- list()
metrics[["N_P"]] <- metricsCluster

save(metrics, file = "cluster3_metrics.rda")
