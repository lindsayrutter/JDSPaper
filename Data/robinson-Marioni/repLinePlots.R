# The raw and TMM normalized versions both look good

#Y1_Y4 raw has 6158 DEGs (which(dataMetrics$K_L$PValue < (0.05/nrow(data))))
data = readRDS("dataRaw.rds")
data[,-1] <- log(data[,-1]+1)
dataMetrics = readRDS("metricListRaw.rds")
plotRepLine(data, dataMetrics)

#Y1_Y4 between has 4818 DEGs (which(dataMetrics$K_L$PValue < (0.05/nrow(data))))
data = readRDS("dataTMM.rds")
data[,-1] <- log(data[,-1]+1)
dataMetrics = readRDS("metricListTMM.rds")
plotRepLine(data, dataMetrics)
