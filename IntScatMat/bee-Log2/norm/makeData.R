dat <- read.delim(file="AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(dat) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
countdata <- dat[ , order(names(dat))]
countdata <- as.matrix(countdata)

y <- DGEList(counts=countdata)

filter <- apply(y, 1, function(x) length(x[x>5])>=8)
filtered <- y[filter,]
filtered <- as.matrix(filtered)

cpm.new <- cpm(filtered, TRUE, TRUE)
cpm.norm <- betweenLaneNormalization(cpm.new, which="full", round=FALSE)
L120 = cpm.L120.norm

setDT(filtered, keep.rownames = TRUE)[]
colnames(filtered)[1] <- "ID"
filtered$ID <- as.character(filtered$ID)
filtered <- as.data.frame(filtered)

beeData <- filtered

save(beeData,file="beeData.Rda")
