library(ggplot2)
library(hexbin)
library (reshape2)
set.seed(1) #Create data
bindata <- data.frame(x=rnorm(100), y=rnorm(100))
fac_probs <- dnorm(seq(-3, 3, length.out=26))
fac_probs <- fac_probs/sum(fac_probs)
bindata$factor <- sample(letters, 100, replace=TRUE, prob=fac_probs)
bindata$factor <- as.factor (bindata$factor)

h <- hexbin(bindata, xbins=5, IDs=TRUE, xbnds=range(bindata$x), ybnds=range(bindata$y))

# Calculate counts for each ID in each facet
counts <- hexTapply(h, bindata$factor, table)
counts <- t(simplify2array(counts))
counts <- melt(counts)
colnames (counts)  <- c ("ID", "factor", "counts")

hexdf <- data.frame(hcell2xy(h), ID=h@cell)
hexdf <- merge(counts, hexdf)




