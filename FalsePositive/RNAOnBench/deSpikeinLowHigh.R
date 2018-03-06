library(devtools)
install_git("https://github.com/plger/RNAontheBENCH", build_vignettes=TRUE)

library(RNAontheBENCH)

data(exampledata)

jpeg('high-AUC.jpg')
resHigh=deSpikein(exampleGeneLevel, method="voom", norm="TMM", quantification="Tophat-featureCounts")
dev.off()

jpeg('low-AUC.jpg')
resLow=deSpikein(exampleGeneLevel, method="t", norm="linear", quantification="Tophat-featureCounts")
dev.off()


