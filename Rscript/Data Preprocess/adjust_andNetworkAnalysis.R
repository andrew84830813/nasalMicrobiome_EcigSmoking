## written by: Andrew Hinton (andrew84@email.unc.edu)


## Install/Load Required Pacakages
library(compositions)
library(data.table)
library(reshape2)
library(vegan)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)



###########################
## Remove Sex effect fom by treatment data
###########################
i =  3# [1 : 4]
df = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv") ) [,-1]
df2  = data.frame( read_csv(file = "Output/allSubjectGroup_nasalMbiome_byGenus.csv") )[,-1]
dat = data.frame(df)
## Remove Sparse Features and impute Zeroes
keep = getFeats_bySparsity(dat[,-1],mxSparsePercent = .8)
fcol = colnames(dat)[1]
dat = subset(dat,select = c(fcol,keep))
dat = dat[rowSums(dat[,-1])>0,]
#percent sparse
sapply(1:nrow(dat), function(x) sum(dat[x,-1]==0) )/ncol(dat[,-1])
factor = 1
ph = compositions::clo(dat[,-1])
impFact = min(ph[ph>0]) / factor
dat = data.frame(Status = factor(dat[,1]), fastImputeZeroes( compositions::clo(dat[,-1]) , impFactor  = impFact ) )
# compute clr matrix
dat.clr = data.frame(clr(dat[,-1]))
rda_tree = rda(dat.clr ~ . , data=data.frame(Sex = factor(df[,1])))
adjPLR = residuals(rda_tree)
#dat =  data.frame(Status = df[,1],clrInv(adjPLR))
dat =  data.frame(Status = df2[,1],clrInv(adjPLR))
#comparison
c = combinat::combn2( unique(dat$Status) )
dat = dat %>%
  filter(Status%in%c[i,])
write_csv(file = dat,"Output/byTreatment_adjSex_byGenus_80sparse.csv")
#######################################


