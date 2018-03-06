library(bigPint)
library("RColorBrewer")
library("gplots")

data("soybean_cn")
data=soybean_cn

select <- order(rowMeans(dat[,-1]),decreasing=TRUE)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(300)
heatmap.2(as.matrix(dat[,-1]), col=hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6))

