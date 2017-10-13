library(devtools)
install_git("https://github.com/plger/RNAontheBENCH", build_vignettes=TRUE)

library(RNAontheBENCH)
vignette("RNAontheBENCH")

data(exampledata)
dir.create("example")
write.table(exampleTranscriptLevel,"example/w12.transcripts.quant",sep="\t",quote=F)
write.table(exampleGeneLevel,"example/w12.genes.quant",sep="\t",quote=F)

# Normalization = TMM
deNanostring(exampleGeneLevel, method="edgeR", norm="TMM", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq", norm="TMM", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq2", norm="TMM", quantification="Tophat-featureCounts")
# Error
# deNanostring(exampleGeneLevel, method="EBSeq", norm="TMM", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="voom", norm="TMM", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="t", norm="TMM", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="logt", norm="TMM", quantification="Tophat-featureCounts")

# Normalization = Linear
deNanostring(exampleGeneLevel, method="edgeR", norm="linear", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq", norm="linear", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq2", norm="linear", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="voom", norm="linear", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="t", norm="linear", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="logt", norm="linear", quantification="Tophat-featureCounts")

# Normalization = RLE
deNanostring(exampleGeneLevel, method="edgeR", norm="RLE", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq", norm="RLE", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq2", norm="RLE", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="voom", norm="RLE", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="t", norm="RLE", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="logt", norm="RLE", quantification="Tophat-featureCounts")

# Normalization = Upperquartile
deNanostring(exampleGeneLevel, method="edgeR", norm="upperquartile", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq", norm="upperquartile", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq2", norm="upperquartile", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="voom", norm="upperquartile", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="t", norm="upperquartile", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="logt", norm="upperquartile", quantification="Tophat-featureCounts")

# Normalization = None
deNanostring(exampleGeneLevel, method="edgeR", norm="none", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq", norm="none", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="DESeq2", norm="none", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="voom", norm="none", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="t", norm="none", quantification="Tophat-featureCounts")
deNanostring(exampleGeneLevel, method="logt", norm="none", quantification="Tophat-featureCounts")

