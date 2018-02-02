library(bigPint)
library(dplyr)

data(soybean_ir)
data(soybean_ir_metrics)
metrics <- soybean_ir_metrics
metrics <- metrics[["N_P"]]
metrics <- metrics[with(metrics, order(PValue)), ]

L120 = soybean_ir
RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
L120 = mutate(L120, mean = (N.1+N.2+N.3+P.1+P.2+P.3)/6, stdev = RowSD(cbind(N.1,N.2,N.3,P.1,P.2,P.3)))
L120q1 = L120[,2:7]
L120q1s = t(apply(as.matrix(L120q1), 1, scale))
data = L120q1s
data = as.data.frame(data)
data$ID = L120$ID
data = data[,c(7,1,2,3,4,5,6)]
colnames(data) = c("ID","N.1","N.2","N.3","P.1","P.2","P.3")

baseOutDir = "Clustering_N_P_Top100"
id100 = metrics$ID[1:1000]
id100df = which(data$ID %in% id100)
data = data[id100df,]

plotName = "N_P"
outDir= paste0(getwd(),"/",baseOutDir)

dendo = abs(data[,-1])
rownames(data[,-1]) = NULL
# Euclidean distance between rows of matrix
d = dist(as.matrix(dendo))
# Hierarchical clustering using ward.D linkage
hc = hclust(d, method="ward.D")
plotName = "N_P_Dendogram"
jpeg(file = paste(outDir, "/", plotName, ".jpg", sep=""))
plot(hc,main=plotName, xlab=NA, sub=NA)
invisible(dev.off())

getPCP <- function(nC, hc, data, outDir){
  
  boxDat <- data %>% gather(key, val, c(-ID))
  colnames(boxDat) <- c("ID", "Sample", "Count")
  
  colList = c("red","magenta", "darkorange","darkgreen","blue","purple")
  k = cutree(hc, k=nC)
  ###########################
  plot_clusters = lapply(1:nC, function(i){
    plotName = "N_P"
    x = as.data.frame(data[which(k==i),])
    nGenes = nrow(x)
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    write.table(xNames, file = paste(outDir, "/", plotName, "_", nC, "_", i, ".txt", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)  
    
    xPCP <- x[,c(1:7)]
    pcpDat2 <- xPCP %>% gather(key, val,c(-ID))
    colnames(pcpDat2) <- c("ID", "Sample", "Count")
    p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat2, aes_string(x = 'Sample', y = 'Count', group = 'ID'), color = colList[i]) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(legend.position = "none", axis.title=element_text(size=12), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=15), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))
    fileName = paste(outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
    jpeg(fileName)
    plot(p)
    invisible(dev.off())
    p    
  })
  ###########################
  jpeg(file = paste(outDir, "/N_P_", nC, ".jpg", sep=""), height = 200 * ceiling((nC+1)/3), width = min(200 * (nC+1), 800))
  # We allow up to 4 (now 3) plots in each column
  p = do.call("grid.arrange", c(plot_clusters, ncol=3)) #change from 4 to 3 ceiling(nC/3)
  invisible(dev.off())
}

for (i in c(2,3,4,5,6)){
  getPCP(nC=i, hc=hc, data=data, outDir=outDir)
}

