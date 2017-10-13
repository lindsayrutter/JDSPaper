# This is derived from explicitprogram.R
data <- read.table("Bottomly.Data/bottomly_count_table.txt",header=T)
count.data <- data[, -1]
rownames(count.data) <- data[, 1]
count.data <- count.data[which(rowSums(count.data) >= 10), ]
# Lindsay added next line
colnames(count.data) <- c("C.1","C.2","C.3","C.4","C.5","C.6","C.7","C.8","C.9","C.10","D.1","D.2","D.3","D.4","D.5","D.6","D.7","D.8","D.9","D.10","D.11")
data2 <- read.table("Bottomly.Data/bottomly_phenodata.txt", header=T)
condition <- factor(as.character(c(rep("C",10),rep("D",11)))) # changed by Lindsay
cond <- condition #added by Lindsay

# This is derived from limma.voomsimul.R
method.norm <- "TMM"
dge <- DGEList(counts=count.data, genes=rownames(count.data), group=factor(cond))
dge <- calcNormFactors(dge, method.norm)
design <- model.matrix(~factor(cond))
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
temp <- topTreat(fit, n=Inf)

#9771
length(which(temp$P.Value<10^-4))
p4i <- which(temp$P.Value<10^-4)
#8977
length(which(temp$P.Value<10^-6))
p6i <- which(temp$P.Value<10^-6)

genes <- rownames(fit[[1]])

cdsDF <- as.data.frame(cds[[1]])
edgeRSig4 <- cdsDF[p4i,]
edgeRSig6 <- cdsDF[p6i,]

saveRDS(edgeRSig4, "edgeRSig4.rds")
saveRDS(edgeRSig6, "edgeRSig6.rds")





