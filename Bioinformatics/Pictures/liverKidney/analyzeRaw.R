# This is the data from reference 6 of TMM Robinson paper

library(data.table)
load("LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]

baseOutDir = getwd()

# Obtain R1 values
outDir = paste0(baseOutDir, "/R1")
dataSel <- data[,c(1:4,7:9)]
plotMDS(dataSel[,-1])
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel, outDir = outDir)
boxSel = dataSel[,-1] %>% gather(Sample,Count)
bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
png(filename=paste0(outDir = outDir,"/boxplot.jpg"))
bPlot
dev.off()
plotScatterStatic(dataSel, threshOrth = 0.5, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 1, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 2, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 3, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 4, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
plotScatterInteractive(dataSel, xbins=20)

# Obtain R2 values
outDir = paste0(baseOutDir, "/R2")
dataSel <- data[,c(1,5:6,11)]
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel, outDir = outDir)
boxSel = dataSel[,-1] %>% gather(Sample,Count)
bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
png(filename=paste0(outDir = outDir,"/boxplot.jpg"))
bPlot
dev.off()
plotScatterStatic(dataSel, threshOrth = 0.5, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 1, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 2, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 3, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 4, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
plotScatterInteractive(dataSel, xbins=20)

# Obtain K values for R1 and R2
outDir = paste0(baseOutDir, "/K")
dataSel <- data[,c(1:6)]
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel, outDir = outDir)
boxSel = dataSel[,-1] %>% gather(Sample,Count)
bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
png(filename=paste0(outDir = outDir,"/boxplot.jpg"))
bPlot
dev.off()
plotScatterStatic(dataSel, threshOrth = 0.5, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 1, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
plotScatterInteractive(dataSel, xbins=15)

# Obtain L values for R1 and R2
outDir = paste0(baseOutDir, "/L")
dataSel <- data[,c(1,7:11)]
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)

plotScatterStatic(dataSel, outDir = outDir)
boxSel = dataSel[,-1] %>% gather(Sample,Count)
bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
png(filename=paste0(outDir = outDir,"/boxplot.jpg"))
bPlot
dev.off()
plotScatterStatic(dataSel, threshOrth = 0.5, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 1, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
plotScatterInteractive(dataSel, xbins=15)
