data(soybean_ir)
soybean_ir <- soybean_ir
soybean_ir[,-1] <- log(soybean_ir[,-1]+1)
plotScatterStatic(soybean_ir)

# Looks terrible because taking ceiling of non-integers
# Also tried rlog on original data of integers, and not good
# rownames(soybean_ir) <- soybean_ir[,1]
# soybean_ir <- soybean_ir[,-1]
# soybean_ir <- as.matrix(soybean_ir)
# coldata <- data.frame(row.names = colnames(soybean_ir), treatment = unlist(lapply(colnames(countdata), function (x) unlist(strsplit(x, "[.]"))[1])))
# dds <- DESeqDataSetFromMatrix(countData = as.data.frame(ceiling(soybean_ir)), colData = coldata, design = ~ treatment)
# dds <- DESeq(dds)
# rld <- rlog(dds)
# rData = as.data.frame(assay(rld))
# setDT(rData, keep.rownames = TRUE)[]
# colnames(rData)[1] <- "ID"
# rData = as.data.frame(rData)
# plotScatterStatic(rData)
