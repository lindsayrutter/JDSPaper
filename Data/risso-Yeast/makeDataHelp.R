# This vignette for EDASeq has an error for the withinLaneNormalization, so here I am combining parts of it with the help function example for withinLaneNormalization (which by itself also has an error)
source("https://bioconductor.org/biocLite.R")
biocLite("yeastRNASeq")
library(yeastRNASeq)

data(geneLevelCounts)
data(yeastGC)

sub <- intersect(rownames(geneLevelCounts), names(yeastGC))
mat <- as.matrix(geneLevelCounts[sub, ])
data <- newSeqExpressionSet(mat, phenoData=laneInfo, featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))
norm <- withinLaneNormalization(data, "gc", which="full", offset=FALSE)
