library(bigPint)
data(soybean_cn_metrics)
metrics <- soybean_cn_metrics

length(which(metrics[["S1_S2"]]$FDR<0.05)) #86
length(which(metrics[["S1_S2"]]$FDR<0.1)) #166
length(which(metrics[["S1_S3"]]$FDR<0.05)) #451
length(which(metrics[["S1_S3"]]$FDR<0.1)) #732
min(metrics[["S1_S2"]]$PValue) #4.388051e-08
min(metrics[["S1_S3"]]$PValue) #2.060093e-08

rm(list=ls())
load("../soybean_cn_switch12_metrics.rda")
metrics <- soybean_cn_switch12_metrics
length(which(metrics[["S1_S2"]]$FDR<0.05)) #0
length(which(metrics[["S1_S2"]]$FDR<0.1)) #0
length(which(metrics[["S1_S2"]]$FDR<0.9)) #0
min(metrics[["S1_S2"]]$PValue) #0.007008889
min(metrics[["S1_S3"]]$PValue) #1.559117e-07

rm(list=ls())
load("../soybean_cn_switch13_metrics.rda")
metrics <- soybean_cn_switch13_metrics
length(which(metrics[["S1_S3"]]$FDR<0.05)) #0
length(which(metrics[["S1_S3"]]$FDR<0.1)) #0
length(which(metrics[["S1_S2"]]$FDR<0.9)) #6
min(metrics[["S1_S2"]]$PValue) #2.280505e-05
min(metrics[["S1_S3"]]$PValue) #0.01209505