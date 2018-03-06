library(bigPint)

ID <- readRDS("../../Bioinformatics/Pictures/liverKidney/Clustering_data_FDR_001_TMMvRaw_Keep/Sig_8_1.Rds")
load("../../Bioinformatics/Pictures/liverKidney/KL_TMM_Metrics.rda")

allMetrics = metricList[["K_L"]]

metricsCluster <- allMetrics[which(allMetrics$ID %in% ID),]
metricsCluster$ID <- as.character(metricsCluster$ID)

metrics <- list()
metrics[["K_L"]] <- metricsCluster

save(metrics, file = "keep1_metrics.rda")
