library(RUVSeq)
library(zebrafishRNASeq)
library(GGally)

data(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]

genes <- rownames(filtered)[grep("ˆENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("ˆERCC", rownames(filtered))]

logDF <- log2(filtered)+1

scatmat(logDF)
ggplot(stack(logDF), aes(x = factor(ind, levels = names(logDF)), y = values)) + geom_boxplot()
