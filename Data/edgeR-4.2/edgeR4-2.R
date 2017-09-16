# This is the data from Section 4.2 of edgeR vignette
# https://mayoclinic.pure.elsevier.com/en/publications/tumor-transcriptome-sequencing-reveals-allelic-expression-imbalan

library(NBPSeq)
library(edgeR)
library(data.table)
data(arab)
data <- as.data.frame(arab)

setDT(data, keep.rownames = TRUE)[]
colnames(data) <- c("ID","m.1","m.2","m.3","h.1","h.2","h.3")
data <- as.data.frame(data)
data[,c(2:ncol(data))] = log(data[,c(2:ncol(data))]+1)
plotScatterStatic(data)
