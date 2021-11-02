## written by: Andrew Hinton (andrew84@email.unc.edu)

rm(list=ls())
gc()

### Load Install/Required Packages  ####
library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(vegan)
library(keras)
library(PRROC)
library(energy)
library(Rfast)
library(dplyr)
library(foreach)
library(parallel)
library(readr)
library(tidyr)
library(stringr)


## Start Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)

## Load Functions
Rcpp::sourceCpp("Functions/cpp_ksTest3.cpp")
source(file = "Functions/functions1.R")


## Model Training Setup
train_control <- trainControl(method="cv",
                              repeats = 5,
                              number = 10,
                              #sampling = "rose",
                              #number=100,
                              seeds = NULL,
                              classProbs = TRUE,
                              savePredictions = T,
                              allowParallel = TRUE,
                              summaryFunction = caret::multiClassSummary
                              #summaryFunction = prSummary
)

tc1  <- trainControl(method="cv",
                     repeats = 5,
                     number = 5,search = "random",
                     #sampling = "rose",
                     #number=100,
                     seeds = NULL,
                     classProbs = TRUE,
                     savePredictions = T,
                     allowParallel = TRUE,
                     summaryFunction = caret::multiClassSummary
                     #summaryFunction = prSummary
)



######===========================================================-
### Read Data Scenario  ####
######===========================================================-
dataset_ = 2
seed_ = 10
set.seed(seed_)


switch(dataset_,

       ## By Sex
       {dat = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv") ) [,-1];
       fname = "bySex_NasalMicrobiomeGenus_80sparse"}, ## - 1
       ## By Treatment (File derived after adjustment)
       {dat = data.frame( read_csv(file = "Output/byTreatment_adjSex_byGenus_80sparse.csv") );
       fname = "byTreatGenusadjSex_NasalMicrobiome_80sparse"} ## - 2

)

dat = data.frame(dat)


###------------------------------------------------------------------------------------------------------------------------###
## Process Zeroes/Sparisity
###------------------------------------------------------------------------------------------------------------------------###
hist(rowSums(dat[,-1]))
mnV=min(dat[,-1])
sum(dat[,-1]==mnV)
clo(table(dat[,1]))
minorityClass = names(which.min(table(dat[,1])))
majorityClass = as.character(unique(dat[dat[,1]!=minorityClass,1]))
keep = getFeats_bySparsity(dat[,-1],mxSparsePercent = .8)
fcol = colnames(dat)[1]
dat = subset(dat,select = c(fcol,keep))
dat = dat[rowSums(dat[,-1])>0,]
factor = 1
ph = compositions::clo(dat[,-1])
impFact = min(ph[ph>0]) / factor
dat = data.frame(Status = factor(dat[,1]), dat[,-1]   )


###------------------------------------------------------------------------------------------------------------------------###
## MST Parms
###------------------------------------------------------------------------------------------------------------------------###
energyMax_Parms = data.frame(patienceInterations = 25,scale_Estat = "scale",
                             nreps.Energy = 1e5,
                             useIG = T,
                             maxFeatures = 500,
                             dcvMethod = 2
)
## Feature Set Sizes to consider
bootstrap = F
featureSet = c(5,10,20,40,60,80,100)
nr = 1
use_tarFeatSelection = T
mlModel = "pls"
###------------------------------------------------------------------------------------------------------------------------###
## Empirical Results ####
###------------------------------------------------------------------------------------------------------------------------###
repeats = 1
perf_ = data.frame()
for(r in 1:repeats){

  mst.empirical = mstTest(trainData1 = dat[,-1],num_Repeats = nr,seed_ = r,
                          ytrain1 = dat[,1],energyMaxParms = energyMax_Parms,
                          permLabels = F,testModels = mlModel,dcv_Folds = 5,dcv_Repeats = 1,
                          useTarSelection = use_tarFeatSelection,featureSetSize = featureSet )
  ## Visualize
  allFeature = unlist(mst.empirical$allFeature_Energy[[1]])
  names(allFeature) = paste0("AF_",names(allFeature))
  allFeature = data.frame(t(allFeature),row.names = r)
  ## Output
  np = mst.empirical$performance
  np = cbind(np,allFeature)
  keyMetrics = c("TrainAUC","estat","raw_estat", "norm_esat","numRatios","numParts",names(allFeature))
  perf1 = subset(np,select = c("seed","feats","method","labels",keyMetrics))
  perf1 = gather(perf1,"Metric","value",5:ncol(perf1))
  perf_ = rbind(perf_,perf1)
  View(perf_)
}
perf.ss1 = spread(perf_,"Metric","value")
perf.ss = cbind(perf.ss1[,!colnames(perf.ss1)%in%keyMetrics],perf.ss1[,keyMetrics]/perf.ss1$numRatios)
perf.s = gather(perf.ss,"Metric","value",5:ncol(perf.ss))
perf = aggregate(value ~ feats + method + Metric, data = perf.s,
                 FUN = function(x) mean(x) )
perf.sd = aggregate(value ~ feats + method + Metric, data = perf.s,
                    FUN = function(x) sd(x) )

perf = aggregate(value ~ feats + method + Metric, data = perf_,
                 FUN = function(x) mean(x) )
perf.sd = aggregate(value ~ feats + method + Metric, data = perf_,
                    FUN = function(x) sd(x) )


ggplot(perf,aes(feats,value,label = round(value,digits = 3) ))+
  geom_line()+
  geom_point(alpha = 1,size = 3,fill = "white",pch = 21,col  ="black")+
  geom_label(col = "black",size = 3)+3
  facet_grid(Metric~method,scales  = "free")

## Select Test Feature Size using estat
estats = unique(mst.empirical$performance$estat)
mxE =  max(perf.ss1$estat)
fts = perf.ss1$feats[which.max(perf.ss1$norm_esat)]#unique(mst.empirical$performance$feats[wE])
fid = which(featureSet==fts)
model1 = perf.ss1$method[which.max(perf.ss1$norm_esat)]
model = strsplit(model1, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)[[1]][1]
mdl = mst.empirical$resampledMetricList[[fid]][[model1]]
bestTune = data.frame( mdl$bestTune)
probs = left_join(bestTune,mdl$pred)
aucVec.df = mdl$resample
mm = str_split(aucVec.df$Resample,pattern = "\\.",simplify = T)[,2]
logLossVec = data.frame()
for(mc in unique(mm)){
  ph = mean(aucVec.df$logLoss[mm==mc])
  logLossVec = rbind(logLossVec,data.frame(Rep = mc,MeanLL = ph))
}
feats.df = mst.empirical$features[[fid]]

## Visualize
rt = prcomp(feats.df[,-1])
coords.df = data.frame(Group = feats.df[,1],rt$x)
library(ggsci)
ggplot(coords.df,aes(PC1,PC2,Fill = Group))+
  scale_alpha_continuous(range = c(0, 1))+
  scale_shape_manual(values = c(21,24,22))+
  scale_color_lancet()+
  scale_fill_lancet()+
  geom_point(aes(shape = Group,fill  = Group),alpha = .6,size = 3)+
  theme_bw()+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 20),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))

targetRatios = ncol(feats.df[,-1])

## Write Table
write_csv(mst.empirical$performance,path = paste0(fname,"_",fts,"_repeatedCV_MST_empirical.csv"))
write_csv(feats.df,paste0(fname,"_",fts,"_keyFeatures.csv"))



###------------------------------------------------------------------------------------------------------------------------###
## Null Results ####
###------------------------------------------------------------------------------------------------------------------------###
nullPerf = data.frame()
resampledAUC.null = list()

nreps = 150
for(r in 1:nreps){
  ## Perform MST Test
  mst.null = mstTest(trainData1 = dat[,-1],
                     ytrain1 = dat[,1],mxNumRatios = targetRatios,
                     permLabels = T,seed = r,num_Repeats = nr,
                     useTarSelection = use_tarFeatSelection,featureSetSize = fts,energyMaxParms = energyMax_Parms,
                     testModels = model )

  allFeature = unlist(mst.null$allFeature_Energy[[1]])
  names(allFeature) = paste0("AF_",names(allFeature))
  allFeature = data.frame(t(allFeature),row.names = r)
  ## Output
  np = mst.null$performance
  np = cbind(np,allFeature)
  nullPerf = rbind(nullPerf,np)
  View(nullPerf)
  ## Save Results
  write_csv(nullPerf,path = paste0(fname,"_",fts,"_repeatedCV_MST_permuted.csv"))
}

