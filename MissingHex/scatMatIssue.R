set.seed(1)
dat <- data.frame(ID = paste0("ID",1:20), A.1 = abs(rnorm(20)), A.2 = abs(rnorm(20)), B.1 = abs(rnorm(20)), B.2 = abs(rnorm(20)))
dat$ID <- as.character(dat$ID)

maxVal = max(dat[,-1])
minVal = min(dat[,-1])
maxRange = c(minVal, maxVal)
binSize=3
buffer = maxRange[2]/(binSize)

my_fn <- function(data, mapping, ...){
  xChar = as.character(mapping$x)
  yChar = as.character(mapping$y)
  x = data[,c(xChar)]
  y = data[,c(yChar)]
  h <- hexbin(x=x, y=y, xbins=binSize, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
  hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
  attr(hexdf, "cID") <- h@cID
  
  hexdf$countColor <- cut2(hexdf$counts, g=6)
  hexdf$countColor2 <- as.factor(unlist(lapply(as.character(hexdf$countColor), function(x) substring(strsplit(gsub(" ", "", x, fixed = TRUE), ",")[[1]][1], 2))))
  hexdf$countColor2 <- factor(hexdf$countColor2, levels = as.character(sort(as.numeric(levels(hexdf$countColor2)))))
  for (i in 1:(length(levels(hexdf$countColor2))-1)){
    levels(hexdf$countColor2)[i] <- paste0(levels(hexdf$countColor2)[i],"-",levels(hexdf$countColor2)[i+1])
  }
  levels(hexdf$countColor2)[length(levels(hexdf$countColor2))] <- paste0(levels(hexdf$countColor2)[length(levels(hexdf$countColor2))], "+")
  my_breaks = levels(hexdf$countColor2)
  clrs <- brewer.pal(length(my_breaks)+3, "Blues")
  clrs <- clrs[3:length(clrs)]
  
  p <- ggplot(hexdf, aes(x=x, y=y, hexID=hexID, counts=counts, fill=countColor2)) + geom_hex(stat="identity") + scale_fill_manual(labels = as.character(my_breaks), values = rev(clrs), name = "Cases count") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_fixed(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
  p
}

# Hexbins from static plot do not look consistent
p <- ggpairs(dat[,2:ncol(dat)], lower = list(continuous = my_fn))

gP <- ggplotly(p)
for (i in 1:(length(gP$x$data)-1)){
  info <- gP$x$data[i][[1]]$text
  info2 <- strsplit(info,"[<br/>]")
  myIndex <- which(startsWith(info2[[1]], "counts:"))
  gP$x$data[i][[1]]$text <- info2[[1]][myIndex]
}
gP$x$data[length(gP$x$data)][[1]]$text <- NULL

# Now, interactive plot is missing certain hexbins
gP
