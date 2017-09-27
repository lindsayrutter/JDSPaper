library(edgeR)

load("All_leaves040615.rda")

ct <- countTable
ct2 <- assays(ct)[[1]]
rownames(ct2) <- ct@rowRanges@elementMetadata@listData$ID
colnames(ct2) <- unlist(strsplit(colnames(ct2), "\\."))[seq(1, 17*3, 3)]
leaves.all <- ct2
load("All_roots.rda")
ct.roots <- countTable
ct2 <- assays(ct.roots)[[1]]
rownames(ct2) <- ct.roots@rowRanges@elementMetadata@listData$ID
colnames(ct2) <- unlist(strsplit(colnames(ct2), "\\."))[seq(1, 18*3, 3)]
roots.all <- ct2
all <- cbind(leaves.all, roots.all[,c(1:15,17:18)])
# (54044, 34)
y <- DGEList(counts=all)

# Now try to threshold count number and normalize on the six samples
#L120 is (56044, 6)
L120 = y[,c("ML08R","ML14R","ML22R","ML11R","ML27R","ML33R")]
colnames(L120)=c("s8","s14","s22","s11","s27","s33")

#L120 is (39120, 6)
# Make sure each gene has at least one count in at least half of the six samples
L120 <- L120[rowSums(L120$counts>1)>=ncol(L120)/2,]
L120 <- as.data.frame(L120[[1]])
colnames(L120) <- c("N.1","N.2","N.3","P.1","P.2","P.3")
soydata <- L120

save(soydata,file="soyData.Rda")
