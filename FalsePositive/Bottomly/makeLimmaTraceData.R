# This is derived from explicitprogram.R
data <- read.table("Bottomly.Data/bottomly_count_table.txt",header=T)
count.data <- data[, -1]
rownames(count.data) <- data[, 1]
count.data <- count.data[which(rowSums(count.data) >= 10), ]
# Lindsay added next line
colnames(count.data) <- c("C.1","C.2","C.3","C.4","C.5","C.6","C.7","C.8","C.9","C.10","D.1","D.2","D.3","D.4","D.5","D.6","D.7","D.8","D.9","D.10","D.11")
data2 <- read.table("Bottomly.Data/bottomly_phenodata.txt", header=T)
condition <- factor(as.character(c(rep("C",10),rep("D",11)))) # changed by Lindsay
cond <- condition #added by Lindsay

# Helper function
tran.est <- function(count.data,cond){
  tmp <- uniroot(tran.reg,interval=c(0.1,1),tmat=as.matrix(count.data),cond=cond)
  return(tmp$root)
}
# Helper function
tran.reg <- function(x,tmat,cond){
  #  start.time <- proc.time()
  require(matrixStats)
  levs <- unique(cond)
  tmat <- as.matrix(tran1(tmat,x))
  means1 <- rowMeans(tmat[,cond==levs[1]])
  means2 <- rowMeans(tmat[,cond==levs[2]])
  vars1 <- rowVars(tmat[,cond==levs[1]])
  vars2 <- rowVars(tmat[,cond==levs[2]])
  means <- c(means1,means2)
  vars <- c(vars1,vars2)
  #  plot(means,vars)
  #  print(x)
  #  print(proc.time()-start.time)
  slope <- (coef(lm(vars~means))[2])
  #  print(slope)
  return(slope)
}
# Helper function
tran1 <- function(y,c1){
  ty <- log(y+sqrt(y^2+y/c1^2)+1/(2*c1^2))
  return(ty)
}
# Helper function
mednorm <- function(mat){
  cmeds <- apply(mat,2,median)
  medmed <- median(cmeds)
  appnorm <- function(x){
    return(medmed*x/cmeds)
  }
  norm.mat <- t(apply(mat,1,appnorm))
  return(norm.mat)
}

# This is derived from limma.transSimul.R
c1 <- tran.est(count.data,cond)
count.data <- mednorm(tran1(count.data,c1))
#  dge <- DGEList(counts=count.data, genes=rownames(count.data), group=factor(cond))
#  dge <- calcNormFactors(dge, method.norm)
design <- model.matrix(~factor(cond))
#  v <- voom(dge, design, plot = FALSE)
fit <- lmFit(count.data, design)
fit <- eBayes(fit)
temp <- topTreat(fit, n=Inf)


#9771
length(which(temp$P.Value<10^-4))
p4i <- which(temp$P.Value<10^-4)
#8977
length(which(temp$P.Value<10^-6))
p6i <- which(temp$P.Value<10^-6)

genes <- rownames(fit[[1]])

cdsDF <- as.data.frame(cds[[1]])
edgeRSig4 <- cdsDF[p4i,]
edgeRSig6 <- cdsDF[p6i,]

saveRDS(edgeRSig4, "edgeRSig4.rds")
saveRDS(edgeRSig6, "edgeRSig6.rds")





