dat <- read.delim(file="/Users/lindz/BeeVirusDiet/AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(dat) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
countdata <- dat[ , order(names(dat))]
countdata <- as.matrix(countdata)
coldata <- data.frame(row.names = colnames(countdata), virus = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],1,1))), diet = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],2,2))), treatment = unlist(lapply(colnames(countdata), function (x) unlist(strsplit(x, "[.]"))[1])))
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)
rld <- rlog(dds)

saveRDS(list(dds, rld), "beeDataDDSRLD.rds")
