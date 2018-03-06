# This is the data from edgeR section 4.3

data = readRDS("dataLib.rds")
baseOutDir = "/Users/lindz/JDSPaper/Data/edgeR-4.3/libNorm"

# Obtain M_F values
outDir = paste0(baseOutDir, "/M_F")
dataSel <- data
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel, outDir = outDir)
boxSel = dataSel[,-1] %>% gather(Sample,Count)
bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
png(filename=paste0(outDir = outDir,"/boxplot.jpg"))
bPlot
dev.off()
plotScatterStatic(dataSel, threshOrth = 0.5, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 2, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, threshOrth = 4, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99, outDir = outDir, option="prediction")
plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
plotScatterInteractive(dataSel, xbins=20)
