library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(glmnet)
library(vegan)
# library(keras)
# library(PRROC)
library(energy)
library(stringr)
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(easyCODA)

library(selEnergyPermR)



# Detect Cores and Start Cluster (optional)
# ensure there is enough memory for parallel computing as multiple cores scales the amount of memory required
num_cores = parallel::detectCores()-2
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)



taxa_level = 1

switch (taxa_level,
        {
          df = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byLevel7.csv") ) [,-1] ### ****
          df2  = data.frame( read_csv(file = "Output/allSubjectGroup_nasalMbiome_byLevel7.csv") )[,-1]
          dat = data.frame(df)
          cnames = read_csv("Output/colNamesTable_byLevel7.csv")

        },

        {
          df = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv") ) [,-1]
          dat = data.frame(df)
          cnames = read_csv("Output/colNamesTable_byGenus.csv")


        }
)


## get colnames
## Remove Sparse Features and impute Zeroes
useALR = F
df.metadata = read.csv("Output/metadata_all.csv")
procData = selEnergyPermR::processCompData(df,minPrevalence = .85)
glm.train = procData$processedData
#dat = zCompositions::cmultRepl(glm.train[,-1])
dat =fastImputeZeroes( compositions::clo(procData$processedData[,-1]),impFactor = procData$impFactor)
df = data.frame(Status = df.metadata$SubjectGroup,dat)

if(useALR){
  ## find best ref
  ff = FINDALR(dat)
  ## with feature selection
  ref = ff$procrust.ref
  dat.plr = data.frame(alr(df[,-1],ivar = ref))
  # dat.plr = data.frame(alr(df[,-1]))

}else{
  dat.plr = selEnergyPermR::calcLogRatio(df)[,-1]
}
y_label = df[,1]




# compute plr matrix
adjustFor = T ## sex and age

if(adjustFor){
  dat.clr = data.frame(clr(dat))
  rda_tree = rda(dat.clr ~ . , data=data.frame(
    Age = df.metadata$Age,
    Sex = factor(df.metadata$Sex)
  )
  )
  adjPLR = residuals(rda_tree)
  df =  data.frame(Status = df$Status,clrInv(adjPLR))

  if(useALR){
    ## find best ref
    ff = FINDALR(df[,-1])
    ## with feature selection
    ref = ff$procrust.ref
    dat.plr = data.frame(alr(df[,-1],ivar = ref))
    ref_name = colnames(df[,-1])[ref]
    # dat.plr = data.frame(alr(df[,-1]))
  }else{
    dat.plr = selEnergyPermR::calcLogRatio(df)[,-1]
  }

  glm.train = df
  y_label = df[,1]
}
#write_csv(df,file = "Results/raData_sepImpute_genusLevel_adjAgeSex_85sparse.csv")




## pca
pc = prcomp(clr(glm.train[,-1]))
pc.df = data.frame(Status  = df$Status,pc$x)
ggplot(pc.df,aes(PC1,PC2,col = Status))+
  geom_point()



## compute logratios
nfolds = 1
f = c(5,10,15,20,25,50,100)
h = 0
## compute scores
dcv_scores = DiCoVarML::computeDCV(train_data = glm.train[,-1],y_train = glm.train$Status,num_repeats = 1,num_folds = nfolds)
scr = dcv_scores$dcv
dcv_str = DiCoVarML::computeDCVStrength(dcv_scores)
res_df = data.frame()


for(num_feats in f){

  ## select Features
  topN = dplyr::top_n(dcv_str,n = num_feats,wt = Str)
  tbl = glm.train %>%
    dplyr::select(Status,one_of(topN$Node))

  r = selEnergyPermR::selectionEnergy.scaled(inputData = tbl,
                                                    patience = 500,
                                                    dcv_nfold = nfolds,
                                                    optimizationMetric = "combinedF"
  )

  ff = r$finalSubset
  x = selEnergyPermR::normalizedEnergy(ff[,-1],ff$Status)$H
  i = data.frame(nfeats = num_feats, h = x,num_ratios = ncol(ff)-1 )
  res_df = rbind(res_df,i)
  View(res_df)

  if(x>h){
    sep.true = r
    h = x
  }
}






kr = colnames(sep.true$finalSubset[,-1])
sc = sep.true$rawDCV
scores = sc
sep.true$optimResult
res = selEnergyPermR::featureSlectionPerformance(tbl = sep.true$finalSubset)
res
targetRatios = res$performance$NumRatios

feats = sep.true$finalSubset
#write_csv(feats,file = "Results/keyLogRatiosfeatures_top15_combF_sepImpute_genusLevel_adjAgeSex.csv")

selEnergyPermR::normalizedEnergy(feats[,-1],feats$Status)$H



## nukk
## compute logratios
null.df = data.frame()
ntrials = 10

for(i in 1:ntrials){
  suppressMessages({suppressWarnings({
    set.seed(i)
    permy = sample(glm.train$Status,length(glm.train$Status));permy
    null_df = data.frame(Status = permy,glm.train[,-1])


    dcv_scores = DiCoVarML::computeDCV(train_data = glm.train[,-1],y_train = permy,num_repeats = 1,num_folds = nfolds)
    scr = dcv_scores$dcv
    dcv_str = DiCoVarML::computeDCVStrength(dcv_scores)
    topN = dplyr::top_n(dcv_str,n = num_feats,wt = Str)
    tbl = null_df %>%
      dplyr::select(Status,one_of(topN$Node))

    ## selper
    sep_null = selEnergyPermR::selectionEnergy.scaled(inputData = tbl,targetFeats = targetRatios,
                                                      patience = 500,
                                                       dcv_nfold = nfolds,
                                                      optimizationMetric = "combinedF"
    )

    null_res = selEnergyPermR::featureSlectionPerformance(tbl = sep_null$finalSubset)
    null_res$performance

    null.df = rbind(null.df,data.frame(Seed = i,null_res$performance))
  })})

  print(i)
  View(null.df)

}

p = (sum(null.df$combinedF>res$performance$combinedF)+1) / (ntrials+1)
p
null.df$obsCombinedF = res$performance$combinedF
null.df$obsEnergyF = res$performance$EnergyF
null.df$p = p
write_csv(null.df,file = "Results/updatedSEP_top15_n100_combF_sepImpute_genusLevel.csv")

hist(null.df$combinedF)




# PCA ---------------------------------------------------------------------
##
pc = prcomp(feats[,-1])
pc.df = data.frame(Status  = feats$Status,pc$x)
ggplot(pc.df,aes(PC1,PC2,col = Status))+
  geom_point()



# PLS ---------------------------------------------------------------------



