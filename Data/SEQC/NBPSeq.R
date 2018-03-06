library(seqc)
options(width=110, digits=2)
ls(2)

require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
norm.data.matrix <- normalize.nbpseq(data.matrix,sample.list)
#p <- stat.nbpseq(norm.data.matrix,sample.list,contrast)

#estimate.norm.factors (estimate normalization factors)
#prepare.nbp (adjust library sizes)
#estimate.disp (fit a dispersion model)
#exact.nb.test exact.nb.test(obj, grp1, grp2, print.level = 1)

norm.factors = estimate.norm.factors(data.matrix, lib.sizes = colSums(data.matrix), method = "AH2010")










## Load Arabidopsis data
data(arab);
## Specify treatment groups
## grp.ids = c(1, 1, 1, 2, 2, 2); # Numbers or strings are both OK
grp.ids = rep(c("mock", "hrcc"), each=3);
## Estimate normalization factors
norm.factors = estimate.norm.factors(arab);
print(norm.factors);
## Prepare an NBP object, adjust the library sizes by thinning the
## counts. For demonstration purpose, only use the first 100 rows of
## the arab data.
set.seed(999);
obj = prepare.nbp(arab[1:100,], grp.ids, lib.size=colSums(arab), norm.factors=norm.factors);
print(obj);
## Fit a dispersion model (NBQ by default)
obj = estimate.disp(obj);
plot(obj);
## Perform exact NB test
## grp1 = 1;
## grp2 = 2;
grp1 = "mock";
grp2 = "hrcc";
obj = exact.nb.test(obj, grp1, grp2);
