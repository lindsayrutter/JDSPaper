library(edgeR)
library(data.table)
library(bigPint)
library("RColorBrewer")
library("gplots")
library(bigPint)

data(soybean_ir)
data = soybean_ir
rownames(data) = data[,1]
data[,-1] = log(data[,-1]+1)

ind <- which(apply(data[,-1], 1, function(X) any(X<2)))
data = data[-(ind),]

set.seed(1)
keep = sample(1:nrow(data), 40)
data2 <- data[keep,]
ggpairs(data2, columns = 2:7)



select <- order(rowMeans(data2[,-1]),decreasing=TRUE)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(300)
heatmap.2(as.matrix(data[,-1]), col=hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6))
