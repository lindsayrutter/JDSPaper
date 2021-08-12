dds <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[1]]
rld <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[2]]

dat <- as.data.frame(assay(rld))
setDT(dat, keep.rownames = TRUE)[]
colnames(dat)[1] <- "ID"
dat$ID <- as.factor(dat$ID)
dat <- as.data.frame(dat)

myLevels <- unique(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- list()

# Runs exact test on all pairs of groups and saves in list
k=1
for (i in 1:(length(myLevels)-1)){
  for (j in (i+1):(length(myLevels))){
    res <- results(dds, contrast=c("treatment",myLevels[i],myLevels[j]))
    dat[[paste(i,j,"FC",sep="-")]] <- res@listData$log2FoldChange
    dat[[paste(i,j,"pval",sep="-")]] <- -1*log10(res@listData$pvalue)
    myPairs[[k]] <- paste(myLevels[i], " and ", myLevels[j])
    k=k+1
  }
}

saveRDS(dat, "beeVolcanoData.rds")
