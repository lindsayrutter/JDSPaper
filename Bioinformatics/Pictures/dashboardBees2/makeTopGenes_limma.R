library(Glimma)
library(GGally)
library(ggplot2)
library(edgeR)

beeCounts <- read.delim(file="AllLaneCount.txt", row.names=1, stringsAsFactors = FALSE)
colnames(beeCounts) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
beeCounts <- beeCounts[ , order(names(beeCounts))]
x <- DGEList(counts=beeCounts)

#samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
#colnames(x) <- samplenames
group <- as.factor(colnames(x))
group <- factor(c(rep("NC",6), rep("NP",6), rep("NR",6), rep("NS",6), rep("VC",6), rep("VP",6), rep("VR",6), rep("VS",6)))
x$samples$group <- group
lane <- as.factor(c("L12","L12","L12","L12","L34","L34","L12","L12","L12","L12","L34","L34","L12","L12","L34","L34","L34","L34","L12","L12","L34","L34","L34","L34","L12","L12","L12","L12","L34","L34","L12","L12","L12","L12","L34","L34","L12","L12","L34","L34","L34","L34","L12","L12","L34","L34","L34","L34"))
x$samples$lane <- lane
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
keep.exprs <- rowSums(cpm>1)>=8 # tried filtering up to 24 and not much difference
x <- x[keep.exprs,, keep.lib.sizes=FALSE] # 15,314 to 10,626
dim(x)

x <- calcNormFactors(x, method = "TMM")

library(RColorBrewer)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(
   NCvsNP = NC-NP,
   NCvsNR = NC-NR,
   NCvsNS = NC-NS,
   NCvsVC = NC-VC,
   NCvsVP = NC-VP,
   NCvsVR = NC-VR,
   NCvsVS = NC-VS,
   NPvsNR = NP-NR,
   NPvsNS = NP-NS,
   NPvsVC = NP-VC,
   NPvsVP = NP-VP,
   NPvsVR = NP-VR,
   NPvsVS = NP-VS,
   NRvsNS = NR-NS,
   NRvsVC = NR-VC,
   NRvsVP = NR-VP,
   NRvsVR = NR-VR,
   NRvsVS = NR-VS,
   NSvsVC = NS-VC,
   NSvsVP = NS-VP,
   NSvsVR = NS-VR,
   NSvsVS = NS-VS,
   VCvsVP = VC-VP,
   VCvsVR = VC-VR,
   VCvsVS = VC-VS,
   VPvsVR = VP-VR,
   VPvsVS = VP-VS,
   VRvsVS = VR-VS,
   levels = colnames(design))

v <- voom(x, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)

pairNames <- colnames(contr.matrix)
topGenes <- list()
genePval <- list()
for (i in 1:length(pairNames)) {
  temp <- topTreat(efit, coef=i, n=Inf)
  setDT(temp, keep.rownames = TRUE)[]
  colnames(temp)[1] <- "ID"
  temp = as.data.frame(temp)
  #sigRows <- which(temp$adj.P.Val<0.05)
  #topGenes[[ pairNames[i] ]] <- temp[sigRows,]
  topGenes[[ pairNames[i] ]] <- temp
  genePval[[ pairNames[i] ]] <- temp
}
saveRDS(topGenes, file="topGenes_limma.Rds")
#saveRDS(genePval, file="genePval_limma.Rds")
#saveRDS(x, file="data_limma.Rds")





