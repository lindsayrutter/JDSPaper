library(devtools)
library(yeastRNASeqRisso2011)
library(EDASeq)
library(DESeq)
library(RCurl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(bigPint)
library(gridExtra)
library(cowplot)
install_github("drisso/yeastRNASeqRisso2011")

# Helper function
formatDF <- function(df){
  setDT(df, keep.rownames = TRUE)[]
  colnames(df) = c("ID","Y1.1","Y1.2","Y2.1","Y2.2","Y7.1","Y7.2","Y4.1","Y4.2","D.1","D.2","D.7","G.1","G.2","G.3")
  df = as.data.frame(df)
  df[,c(2:ncol(df))] = log(df[,c(2:ncol(df))]+1)
  df
}

# Read in three .rda files
githubURL <- "https://github.com/drisso/yeastRNASeqRisso2011/blob/master/data/"
load(url(paste0(githubURL, "geneLevelCounts.rda?raw=true")))
load(url(paste0(githubURL, "laneInfo.rda?raw=true")))
load(url(paste0(githubURL, "geneInfo.rda?raw=true")))

# Prepare data frame based on parameters from .rda files
data(yeastGC)
colnames(laneInfo)[2] <- "conditions"
means <- rowMeans(geneLevelCounts)
filter <- means >= 10
geneLevelCounts <- geneLevelCounts[filter,]
sub <- intersect(rownames(geneLevelCounts), names(yeastGC))
mat <- as.matrix(geneLevelCounts[sub, ])
data <- newSeqExpressionSet(mat, phenoData=laneInfo, featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))

# Run within and between normalization on data frame
dataWithin <- withinLaneNormalization(data, "gc", which="full", offset=FALSE)
dataNorm <- betweenLaneNormalization(dataWithin, which="median")
counts <- as(dataWithin,"CountDataSet")
dataWithin <- as.data.frame(counts@assayData$counts)
counts <- as(dataNorm,"CountDataSet")
dataBetween <- as.data.frame(counts@assayData$counts)
dataWithin <- formatDF(dataWithin)
dataWithin <- select(dataWithin, ID, Y1.1, Y1.2, Y4.1, Y4.2)
dataBetween <- formatDF(dataBetween)
dataBetween <- select(dataBetween, ID, Y1.1, Y1.2, Y4.1, Y4.2)

# Use bigPint function to create scatterplot matrices showing within and between normalization
retWithin <- plotScatterStatic(dataWithin, option="point", pointSize = 0.1, saveFile = FALSE)
retBetween <- plotScatterStatic(dataBetween, option="point", pointSize = 0.1, saveFile = FALSE)

# Arrange two scatterplot matrices into grid
plot1 <- plot_grid(ggmatrix_gtable(retWithin[["Y1_Y4"]]), labels=c("A"), ncol = 1, nrow = 1, label_size=9) + theme(plot.background = element_rect(size=0.1,linetype="solid",color="black"))
plot2 <- plot_grid(ggmatrix_gtable(retBetween[["Y1_Y4"]]), labels=c("B"), ncol = 1, nrow = 1, label_size=9) + theme(plot.background = element_rect(size=0.1,linetype="solid",color="black"))
grid.arrange(plot1, plot2, ncol=1)

# Save figure (Optional)
# jpeg(filename=paste0(getwd(), "/", "yeastWithinBetween.jpg"), height=600, width=300)
# grid.arrange(plot1, plot2, ncol=1)
# dev.off()
