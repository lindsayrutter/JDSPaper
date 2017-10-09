#Y1_Y4 raw has 173 DEGs
data = readRDS("dataRaw.rds")
data = data[,c(1,2,7,8)]
data = setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","Y1.1","Y1.2","Y4.1","Y4.2")
data = as.data.frame(data)
data[,-1] <- log(data[,-1]+1)
dataMetrics = readRDS("metricListRaw.rds")
plotRepLine(data, dataMetrics)

#Y1_Y4 between has 412 DEGs
data = readRDS("dataBetween.rds")
data = data[,c(1,2,7,8)]
data = setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","Y1.1","Y1.2","Y4.1","Y4.2")
data = as.data.frame(data)
data[,-1] <- log(data[,-1]+1)
dataMetrics = readRDS("metricListBetween.rds")
plotRepLine(data, dataMetrics)
