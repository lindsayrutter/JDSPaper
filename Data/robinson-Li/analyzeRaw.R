# This is the data from reference 12 of TMM Robinson paper

library(data.table)
library(readxl)
data1 <- read_excel("~/JDSPaper/Data/robinson-Cloonan/Grimmond_lengths.xls")

data = data1[,c(3,5,6)]
data = as.data.frame(data)
colnames(data) = c("ID","ES","EB")

baseOutDir = "/Users/lindz/JDSPaper/Data/robinson-Cloonan/raw"

# Obtain ES-EB values
outDir = paste0(baseOutDir, "/ESEB")
dataSel <- data
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
