data(soybean_cn)
data(soybean_cn_metrics)
#plotDEG(soybean_cn, soybean_cn_metrics)

data=soybean_cn; dataMetrics=soybean_cn_metrics; outDir=getwd(); pointSize=0.5; bluePointSize=0.1; redPointSize=0.1; greyPointSize=0.1; lineSize=0.1; xbins=10; piLevel=0.95; threshFC=3; threshOrth=3; threshVar="FDR"; threshVal=0.05; logFC="logFC"; PValue="PValue"; option="scatterHexagon"

colNames <- colnames(data)
myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- myPairs[-which(myPairs=="ID")]
colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])

ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)

maxVal = max(data[,-1])
minVal = min(data[,-1])
maxRange = c(minVal, maxVal)
xbins=xbins
buffer = maxRange[2]/xbins

i=1;j=3
group1 = myPairs[i]
group2 = myPairs[j]
datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
rowDEG <- c(rowDEG1, rowDEG2)
degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
degID <- c(degID1, degID2)
degData <- datSel[which(datSel$ID %in% degID),]

my_fn <- function(data, mapping, ...){
  xChar = as.character(mapping$x)
  yChar = as.character(mapping$y)
  x = data[,c(xChar)]
  y = data[,c(yChar)]
  # removed shape=1
  h <- hexbin(x=x, y=y, xbins=xbins, IDs=TRUE, xbnds=maxRange, ybnds=maxRange) 

  # Added
  # counts <- hexTapply(h, bindata$factor, table)
  # counts <- t(simplify2array(counts))
  # counts <- melt(counts)
  # colnames(counts) <- c("ID", "factor", "counts")
  counts <- t(simplify2array(h@count))
  counts <- melt(counts)
  colnames(counts) <- c("ID", "factor", "counts")
  
  hexdf <- data.frame(hcell2xy(h), hexID=h@cell, counts=counts)
  hexdf <- merge(counts, hexdf)
  attr(hexdf, "cID") <- h@cID
  
  hexdf$countColor <- Hmisc::cut2(hexdf$counts, g=6)
  hexdf$countColor2 <- as.factor(unlist(lapply(as.character(hexdf$countColor), function(x) substring(strsplit(gsub(" ", "", x, fixed = TRUE), ",")[[1]][1], 2))))
  hexdf$countColor2 <- factor(hexdf$countColor2, levels = as.character(sort(as.numeric(levels(hexdf$countColor2)))))
  
  for (i in 1:(length(levels(hexdf$countColor2))-1)){
    levels(hexdf$countColor2)[i] <- paste0(levels(hexdf$countColor2)[i],"-",levels(hexdf$countColor2)[i+1])
  }
  
  levels(hexdf$countColor2)[length(levels(hexdf$countColor2))] <- paste0(levels(hexdf$countColor2)[length(levels(hexdf$countColor2))], "+")
  
  my_breaks = levels(hexdf$countColor2)
  clrs <- RColorBrewer::brewer.pal(length(my_breaks)+3, "Blues")
  clrs <- clrs[3:length(clrs)]
  
  hexdf$counts [hexdf$counts == 0] <- NA
  
  p <- ggplot(hexdf, aes(x=x, y=y, hexID=hexID, counts=counts, fill=countColor2)) + geom_hex(stat="identity") + scale_fill_manual(labels = as.character(my_breaks), values = rev(clrs), name = "Cases count") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer)) + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = "orange", size = pointSize)
  
  #coord_fixed(ratio=1)
  #+ coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
  
  #p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer)) + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = "orange", size = pointSize)
  p
}

# Problems with hexagon sizes if adjacent hexagon sizes not populated
p1 <- ggpairs(datSel[,-1], lower = list(continuous = my_fn))
ggplotly(p1) # Problem fixed in Plotly
