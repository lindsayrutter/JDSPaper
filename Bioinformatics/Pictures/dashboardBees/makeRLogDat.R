dat <- readRDS("/Users/lindz/bigPint/tblshoot/AllPairs/data_limma.Rds")
dat <- dat$counts
dat <- as.data.frame(dat)

# rlog() method requires dds object
coldata = data.frame(row.names = colnames(dat), virus = unlist(lapply(colnames(dat), function (x) substring(unlist(strsplit(x, "[.]"))[1],1,1))), diet = unlist(lapply(colnames(dat), function (x) substring(unlist(strsplit(x, "[.]"))[1],2,2))), treatment = unlist(lapply(colnames(dat), function (x) unlist(strsplit(x, "[.]"))[1])))
dds = DESeqDataSetFromMatrix(countData = dat, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)
logDat <- rlog(dds)
dat <- as.data.frame(assay(logDat))
setDT(dat, keep.rownames = TRUE)[]
colnames(dat)[1] <- "ID"
dat <- as.data.frame(dat)
saveRDS(dat, "RLogDat.Rds")
