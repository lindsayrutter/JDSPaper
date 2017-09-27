load("../All_leaves040615.rda")

ct <- countTable
ct2 <- assays(ct)[[1]]
rownames(ct2) <- ct@rowRanges@elementMetadata@listData$ID
colnames(ct2) <- unlist(strsplit(colnames(ct2), "\\."))[seq(1, 17*3, 3)]
leaves.all <- ct2
load("../All_roots.rda")
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
#L120a <- calcNormFactors(L120, method="none") # does nothing
# Now positive and negative
cpm.L120.new <- cpm(L120, TRUE, TRUE)
cpm.L120.norm <- betweenLaneNormalization(cpm.L120.new, which="full", round=FALSE)
L120 = cpm.L120.norm

RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

L120t = L120
L120 = as.data.frame(L120t)
L120 = mutate(L120, mean = (s8+s14+s22+s11+s27+s33)/6, stdev = RowSD(cbind(s8,s14,s22,s11,s27,s33)))
rownames(L120)=rownames(L120t)
# L120 = (39120, 8)

# The first quartile threshold of mean counts across the 5 samples
q1T = as.numeric(summary(L120$mean)["1st Qu."])
L120q1 = subset(L120,mean>q1T)
# The first quartile threshold of standard deviation across the 5 samples
q1Ts = as.numeric(summary(L120q1$stdev)["1st Qu."])
# L120q1 (22004, 8)
L120q1 = subset(L120q1,stdev>q1Ts)
# filt (17116, 8)
filt = subset(L120,mean<=q1T|stdev<=q1Ts)
ind = seq(1, nrow(L120), by=10)
L120Plot=L120[ind, ]

model = loess(mean ~ stdev, data=L120q1)
# L120q1 (9809, 8)
L120q1 = L120q1[which(sign(model$residuals) == 1),]

L120q1 = L120q1[,1:6]
L120q1s = t(apply(as.matrix(L120q1), 1, scale))
colnames(L120q1s)=c("N.1","N.2","N.3","P.1","P.2","P.3")
colnames(L120q1)=c("N.1","N.2","N.3","P.1","P.2","P.3")
filt = filt[,1:6]
colnames(filt)=c("N.1","N.2","N.3","P.1","P.2","P.3")
# filt (17116, 8)
filt = rbind(filt,L120q1[which(sign(model$residuals) == -1),])
# filt (29311, 6)
filts = t(apply(as.matrix(filt), 1, scale))
colnames(filts)=c("N.1","N.2","N.3","P.1","P.2","P.3")
colnames(filt)=c("N.1","N.2","N.3","P.1","P.2","P.3")

setDT(L120q1, keep.rownames = TRUE)[]
colnames(L120q1)[1] <- "ID"
bindata <- L120q1[,-1]
bindata <- apply(bindata, 2, function(x) {ifelse(x < 0, 0, x)})
bindata <- cbind(L120q1[,1], bindata)
bindata <- as.data.frame(bindata)
bindata$ID <- as.character(bindata$ID)

save(bindata,file="bindataL120.Rda")
