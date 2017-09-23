
rawDEG = readRDS("metricListRaw.rds")
withinDEG = readRDS("metricListWithin.rds")
betweenDEG = readRDS("metricListBetween.rds")

pairs = names(rawDEG)

raw = length(which(rawDEG[[names(rawDEG)[1]]]$FDR<0.05))
within = length(which(withinDEG[[names(withinDEG)[1]]]$FDR<0.05))
between = length(which(betweenDEG[[names(betweenDEG)[1]]]$FDR<0.05))
C = c(raw, within, between)

for (i in 2:length(pairs)){
  raw = length(which(rawDEG[[names(rawDEG)[i]]]$FDR<0.05))
  within = length(which(withinDEG[[names(withinDEG)[i]]]$FDR<0.05))
  between = length(which(betweenDEG[[names(betweenDEG)[i]]]$FDR<0.05))
  X = c(raw, within, between)
  C <- rbind(C, X)
}
C = as.data.frame(C)
colnames(C) = c("raw","within","between")
rownames(C) = pairs
degDF = C
saveRDS(degDF, "compareDEG.rds")
