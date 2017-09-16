# This is the data from Section 4.1 of edgeR vignette (https://www.bioconductor.org/packages/3.3/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
# https://mayoclinic.pure.elsevier.com/en/publications/tumor-transcriptome-sequencing-reveals-allelic-expression-imbalan

library(readxl)
data <- read_excel("~/JDSPaper/Data/journal.pone.0009317.s009.xls")
data <- data[3:nrow(data),]
dataSel <- data[,c(1,4:9)]
colnames(dataSel) <- c("ID","N.8","T.8","N.33","T.33","N.51","T.51")
datNums = dataSel[,2:ncol(dataSel)]
datNums2 = lapply(datNums, function(x) as.numeric(as.character(x)))
dataSel[,2:ncol(dataSel)] = datNums2
dataSel = dataSel[,c(1,2,4,6,3,5,7)]
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel)
