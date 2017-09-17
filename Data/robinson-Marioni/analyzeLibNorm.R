# This is the data from reference 6 of TMM Robinson paper.
# It is the data that has been normalized by library size, created in makeDataLibTMM.R

data = readRDS("dataLib.rds")

baseOutDir = "/Users/lindz/JDSPaper/Data/robinson-Marioni/libNorm"

# Obtain R1 values
outDir = paste0(baseOutDir, "/R1")
dataSel <- data[,c(1:4,7:9)]
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel, outDir = outDir)

