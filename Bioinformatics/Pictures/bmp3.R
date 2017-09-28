# Script demonstrating that boxplots and MDS plots do not show information about individual observations. In contrast, parallel coordinate plots do show individual observations and allow for us to better see whether replications are consistent for individual observations.

library(ggplot2)
library(cowplot) # Need for combining ggplot2 output into one plot aesthetically
library(data.table)

# This function creates a boxplot, MDS plot, and parallel coordinate plot for five replications
makePlots <- function(A.1, A.2, A.3, B.1, B.2, B.3, i){
  dat <- data.frame(ID = paste0("ID", 1:50), A.1, A.2, A.3, B.1, B.2, B.3)
  datM <- melt(dat, id.vars = "ID")
  colnames(datM) <- c("ID", "Sample", "Count")
  
  boxPlots[[i]] <<- ggplot(datM, aes(Sample, Count)) + geom_boxplot() + theme(text = element_text(size=12))
  
  # Convert DF from scatterplot to PCP
  datt <- data.frame(t(dat))
  names(datt) <- as.matrix(datt[1, ])
  datt <- datt[-1, ]
  datt[] <- lapply(datt, function(x) type.convert(as.character(x)))
  setDT(datt, keep.rownames = TRUE)[]
  dat_long <- melt(datt, id.vars ="rn" )
  colnames(dat_long) <- c("Sample", "ID", "Count")
  
  pcpPlots[[i]] <<- ggplot(dat_long) + geom_line(aes(x = Sample, y = Count, group = ID, color = ID)) + theme(legend.position="none", text = element_text(size=12))
  
  tDat <- t(dat[,2:7]) #orig 2:6
  datD <- as.matrix(dist(tDat))
  fit <- cmdscale(datD, eig = TRUE, k = 2)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  mdsPlots[[i]] <<- qplot(x,y) + geom_text(label = names(x), nudge_y = 0.35) + labs(x = "Dim 1", y = "Dim 2") + theme(text = element_text(size=12))
}

set.seed(30) #30, 50
boxPlots <- vector('list', 2)
pcpPlots <- vector('list', 2)
mdsPlots <- vector('list', 2)

# In the first case, we purposely create individual observations that will be inconsistent across their replications
A.1=c(rnorm(25,10),rnorm(25,6))
A.2=c(rnorm(25,10),rnorm(25,6))
A.3=c(rnorm(25,10),rnorm(25,6))
B.1=c(rnorm(25,6),rnorm(25,10))
B.2=c(rnorm(25,6),rnorm(25,10))
B.3=c(rnorm(25,6),rnorm(25,10))

makePlots(A.1, A.2, A.3, B.1, B.2, B.3, 1)

# In the second case, we purposely create individual observations that will be more consistent across their replications
A.1=c(rnorm(50,8))
A.2=c(rnorm(50,8))
A.3=c(rnorm(50,8))
B.1=c(rnorm(50,8))
B.2=c(rnorm(50,8))
B.3=c(rnorm(50,8))

makePlots(A.1, A.2, A.3, B.1, B.2, B.3, 2)

# View the first case
plot_grid(boxPlots[[1]], mdsPlots[[1]], pcpPlots[[1]], labels=c("A", "B", "C"), ncol = 1, nrow = 3)

# View the second case
plot_grid(boxPlots[[2]], mdsPlots[[2]], pcpPlots[[2]], labels=c("A", "B", "C"), ncol = 1, nrow = 3)
