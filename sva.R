library(sva)
library(bladderbatch)
library(pamr)
library(limma)
data(bladderdata)
# contains batch numbers and more for the 57 samples
pheno = pData(bladderEset)
# count table
edata = exprs(bladderEset)

# full model (variable of interest because no adjustment variables)
mod = model.matrix(~as.factor(cancer), data=pheno)
# null model (no adjustment variables, so just have intercept)
mod0 = model.matrix(~1,data=pheno)

# Now that the model matrices have been created, we can apply the sva function to estimate batch and other artifacts.
# Identifies the number of latent factors
n.sv = num.sv(edata,mod,method="leek")
# Apply sva to estimate the surrogate variables
svobj = sva(edata,mod,mod0,n.sv=n.sv)

############################################################### 
# ADJUSTING FOR SURROGATE VARIABLES USING THE F.PVALUE FUNCTION

# Ftest compared models mod and mod0
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")

# Almost 70% of qValues had FDR<0.05 which seems artifically high. So, we now do the same analysis taking into account the surrogate variables. To do so, we include the surrogate variables as adjustment variables in both models. Now these are the adjusted P-values and Q-values accounting for surrogate variables.
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

############################################################### 
# APPLYING COMBAT TO ADJUST FOR KNOWN BATCHES

batch = pheno$batch
# Create a model matrix for the adjustment variables, including the variable of interest. Note that you do not include batch in creating this model matrix - it will be included later in the ComBat function. In this case there are no other adjustment variables so we simply fit an intercept term.
modcombat = model.matrix(~1, data=pheno)

# This returns an expression matrix, with the same dimensions as your original dataset. This new expression matrix has been adjusted for batch
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# These P-values and Q-values now account for the known batch e↵ects included in the batch variable.
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

