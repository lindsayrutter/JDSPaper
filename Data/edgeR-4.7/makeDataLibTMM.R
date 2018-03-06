# This is the data from Section 4.7 of edgeR vignette
# https://mayoclinic.pure.elsevier.com/en/publications/tumor-transcriptome-sequencing-reveals-allelic-expression-imbalan


library(readxl)
data <- read_excel("~/JDSPaper/Data/edgeR-4.7/GSE86297_RNA-Seq_Read_Counts.xls")
data <- data[,6:13]
data = as.data.frame(data)
colnames(data) = c("A.1","A.2","B.1","B.2","C.1","C.2","D.1","D.2")
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
y <- DGEList(counts=data, group= c(rep(1,2),rep(2,2),rep(3,2),rep(4,2)))
y <- calcNormFactors(y, method="TMM")
identical(tmmSizes, y[[2]]$norm.factors)
identical(data, as.data.frame(y[[1]]))

