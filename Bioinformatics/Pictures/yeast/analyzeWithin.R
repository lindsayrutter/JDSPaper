# This is the data from EDASeq vignette

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
data = readRDS("dataWithin.rds")
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","Y1.1","Y1.2","Y2.1","Y2.2","Y7.1","Y7.2","Y4.1","Y4.2","D.1","D.2","D.7","G.1","G.2","G.3")
data = as.data.frame(data)

outDir = "/Users/lindz/JDSPaper/Bioinformatics/Pictures/yeast/within"

dataSel <- data
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel, outDir = outDir, option="point")
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
