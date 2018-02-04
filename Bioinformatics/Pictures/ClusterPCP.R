library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

data(soybean_ir)
data(soybean_ir_metrics)
metrics <- soybean_ir_metrics
metrics <- metrics[["N_P"]]
metrics <- metrics[with(metrics, order(PValue)), ]
sbLog = soybean_ir
sbLog[,-1] <- log(sbLog[,-1]+1)
soybean_ir_metrics[["N_P"]]$ID <- as.character(soybean_ir_metrics[["N_P"]]$ID)

outDir = "Clustering_data2"

#13771 DEGs p-value<0.05
plotDEG(sbLog, soybean_ir_metrics, option="scatterPoints", threshVar = "PValue", fileName=paste(getwd(), "/", outDir, "/", "SM_DEG_PVal_0.05.jpg", sep=""))

# 2050 DEGs FDR<0.05
plotDEG(sbLog, soybean_ir_metrics, option="scatterPoints", threshVar = "PValue", threshVal = 0.05/nrow(sbLog), fileName=paste(getwd(), "/", outDir, "/", "SM_DEG_FDR_0.05.jpg", sep=""))

# 1681 DEGs FDR<0.01
plotDEG(sbLog, soybean_ir_metrics, option="scatterPoints", threshVar = "PValue", threshVal = 0.01/nrow(sbLog), fileName=paste(getwd(), "/", outDir, "/", "SM_DEG_FDR_0.01.jpg", sep=""))

#### Determine the five streak genes (all have PValues between 0.58 and 0.90)
streakIDs = which(metrics$ID %in% c("Glyma.18G057100.Wm82.a2.v1", "Glyma.18G092200.Wm82.a2.v1", "Glyma.09G164000.Wm82.a2.v1", "Glyma.08G269900.Wm82.a2.v1", "Glyma.07G079300.Wm82.a2.v1"))
metrics$PValue = 1e-10
streakMetrics = metrics[streakIDs,]
metrics = list()
metrics[["N_P"]] = streakMetrics

plotDEG(data=sbLog, dataMetrics=metrics, threshVar = "PValue", outDir=outDir, option="parallelCoord", fileName=paste(getwd(), "/", outDir, "/", "N_P_Streak_SM.jpg", sep=""), lineSize=0.5)

soyScale = soybean_ir
soyScale[,-1] <- t(apply(as.matrix(soybean_ir[,-1]), 1, scale))

plotDEG(data=soyScale, dataMetrics=metrics, threshVar = "PValue", outDir=outDir, option="parallelCoord", fileName=paste(getwd(), "/", outDir, "/", "N_P_Streak_SM_Scaled.jpg", sep=""), lineSize=0.5)
