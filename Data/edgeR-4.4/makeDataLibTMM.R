# This is the data from Section 4.3 of edgeR vignette
# https://mayoclinic.pure.elsevier.com/en/publications/tumor-transcriptome-sequencing-reveals-allelic-expression-imbalan

library(readxl)
data <- read_excel("~/JDSPaper/Data/edgeR-4.4/GSE60450_Lactation-GenewiseCounts.xls")
data = data[,-c(1,2)]
data = data[,c(1,3,5,7,9,11)]
colnames(data) = c("L2.1", "L2.2", "L2.3", "L1.1", "L1.2", "L1.3")
data=as.data.frame(data)
saveRDS(data, "dataRaw.rds")

# Make dataLib object
libSizes <- colSums(data)
dataLib = data.frame('a'=data[,1]/libSizes[1])
for (i in 2:ncol(data)){
  dataLib[,i] = data[,i]/libSizes[i]
}
colnames(dataLib) = colnames(data)
rownames(dataLib) = rownames(data)
setDT(dataLib, keep.rownames = TRUE)[]
colnames(dataLib)[1] = "ID"
dataLib = as.data.frame(dataLib)
saveRDS(dataLib, "dataLib.rds")

# Make dataTMM object
tmmSizes <- calcNormFactors(data, logratioTrim=.3)
dataTMM = data.frame('a'=data[,1]/tmmSizes[1])
for (i in 2:ncol(data)){
  dataTMM[,i] = data[,i]/tmmSizes[i]
}
colnames(dataTMM) = colnames(data)
rownames(dataTMM) = rownames(data)
setDT(dataTMM, keep.rownames = TRUE)[]
colnames(dataTMM)[1] = "ID"
dataTMM = as.data.frame(dataTMM)
saveRDS(dataTMM, "dataTMM.rds")

# Make dataTMM2 object (see if get same TMM normalization factors)
y <- DGEList(counts=data, group= c(rep(1,3),rep(2,3)))
y <- calcNormFactors(y, method="TMM")
identical(tmmSizes, y[[2]]$norm.factors)
identical(data, as.data.frame(y[[1]]))

