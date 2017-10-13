# This is derived from explicitprogram.R
data <- read.table("Bottomly.Data/bottomly_count_table.txt",header=T)
count.data <- data[, -1]
rownames(count.data) <- data[, 1]
count.data <- count.data[which(rowSums(count.data) >= 10), ]
# Lindsay added next line
colnames(count.data) <- c("C.1","C.2","C.3","C.4","C.5","C.6","C.7","C.8","C.9","C.10","D.1","D.2","D.3","D.4","D.5","D.6","D.7","D.8","D.9","D.10","D.11")
data2 <- read.table("Bottomly.Data/bottomly_phenodata.txt", header=T)
condition <- factor(as.character(c(rep("C",10),rep("D",11)))) # Changed by Lindsay
cond <- condition # Added by Lindsay

# This is derived from edgeRsimul.R
mode.norm <- "TMM"
dispersion.type <- "auto"
cds <- DGEList(count.data, group=cond)
cds <- calcNormFactors(cds, method=mode.norm )
cds <- estimateCommonDisp(cds)
cds <- estimateTagwiseDisp(cds)
et  <- exactTest(cds, dispersion=dispersion.type)

#536
length(which(et[[1]]$PValue<10^-4))
p4i <- which(et[[1]]$PValue<10^-4)
#302
length(which(et[[1]]$PValue<10^-6))
p6i <- which(et[[1]]$PValue<10^-6)

cdsDF <- as.data.frame(cds[[1]])
edgeRSig4 <- cdsDF[p4i,]
edgeRSig6 <- cdsDF[p6i,]

saveRDS(edgeRSig4, "edgeRSig4.rds")
saveRDS(edgeRSig6, "edgeRSig6.rds")


