rm(list=ls())
gc()

### Load Required Packages  ####
library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(tidyverse)
library(glmnet)
library(vegan)
library(keras)
library(PRROC)
library(energy)

## Start Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)


source(file = "Functions/functions1.R")


## Model Training Setup
train_control <- trainControl(method="cv",
                              repeats = 5,
                              number = 10,
                              seeds = NULL,
                              classProbs = TRUE,
                              savePredictions = T,
                              allowParallel = TRUE,
                              summaryFunction = caret::multiClassSummary
)

tc1  <- trainControl(method="repeatedcv",
                     repeats = 20,
                     number = 10,search = "random",
                     #sampling = "rose",
                     #number=100,
                     seeds = NULL,
                     classProbs = TRUE,
                     savePredictions = T,
                     allowParallel = TRUE,
                     summaryFunction = caret::multiClassSummary
)



##-----------------------------------------*
### Nasal By Sex and by Genus ####
##-----------------------------------------*
colNames = read_csv("Output/colNamesTable_byGenus.csv")
colNames = separate(colNames,col = 1,into = c('Family','Genus'),sep = "g__")
colNames$Family = str_replace(colNames$Family,"f__","")
colNames$Family = str_replace_all(colNames$Family,"\\.","")
Name = if_else(stringr::str_length(colNames$Genus )==0,colNames$Family,colNames$Genus)
colNames$names = Name
colNames = data.frame(colNames)
## Load Data
dat = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv") ) [,-1]
table(dat$Sex)
sampleID = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv")[,1] )
prot = data.frame(read_csv(file = "Output/NLF_bySex.csv"));colnames(prot) = str_replace_all(colnames(prot),"..ug.mL.",replacement = "")
procData = processCompData(dat,minPrevalence = .80)
dat = procData$processedData
impFact = procData$impFactor
minorityClass = procData$minClss
majorityClass = procData$majClass


## load empirical Results
results.emp = read_csv("Output/bySex_NasalMicrobiomeGenus_80sparse_10_repeatedCV_MST_empirical.csv")
results.perm = read_csv("Output/bySex_NasalMicrobiomeGenus_80sparse_10_repeatedCV_MST_permuted.csv")
## Load Features
feature.df = read_csv(file = "Output/bySex_NasalMicrobiomeGenus_80sparse_10_keyFeatures.csv")
feature.df = data.frame(feature.df)

## Integrate  Data
lrs.df = getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat[,-1],Class = dat[,1])
lrs.df = data.frame(ID  = sampleID$X1,lrs.df)
dat.intg = left_join(prot,lrs.df)
bool = !is.na(rowSums(dat.intg[,-2:-1]))
dat.intg = dat.intg[bool,]
cormat = cor(dat.intg[,-2:-1],method = "spearman")



###--------------------------------------*
## Visualizer Microbiome before after Sex adj. ####
##----------------------------------------*
## By Genus
df = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv") ) [,-1]
df2  = data.frame( read_csv(file = "Output/allSubjectGroup_nasalMbiome_byGenus.csv") )[,-1]
dat = data.frame(df)
lrs.unadj = getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat[,-1],Class = dat[,1])
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
dat.clr = data.frame(clr(dat[,-1]))
## rda
rda_tree = rda(dat.clr ~ . , data=data.frame(Sex = factor(df[,1])))
adjPLR = fitted(rda_tree)
adjPLR = residuals(rda_tree)
RsquareAdj(rda_tree)
coefficients(rda_tree)
anova(rda_tree,permutations = 5000)
#dat =  data.frame(Status = df[,1],clrInv(adjPLR))
dat =  data.frame(Status = df[,1],clrInv(adjPLR))
lrs.adj = getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat[,-1],Class = dat[,1])
## Visualize Sex after adj
pc = cmdscale(d = dist(lrs.adj[,-1]))
pc.df = data.frame(labels = dat[,1],pc)
ggplot(pc.df,aes(X1,X2,col = labels))+
  geom_point(size = 3,alpha = .75)+
  theme_bw()+
  ggtitle("After Adj.")+
  stat_ellipse()
#before adj
pc = cmdscale(d = dist(feature.df[,-1]))
pc.df = data.frame(labels = feature.df[,1],pc)
ggplot(pc.df,aes(X1,X2,col = labels))+
  geom_point(size = 3,alpha = .75)+
  theme_bw()+
  ggtitle("Before Adj.")+
  stat_ellipse()


##---------------------------*
## permanova mbiome feats ####
##-------------------------****
#dispersion test
d = dist(feature.df[,-1])
mod = vegan::betadisper(d,group = feature.df[,1])
anova(mod)
permutest(mod)
mod.HSD <- TukeyHSD(mod)
plot(mod,)
boxplot(mod)
plot(TukeyHSD(mod))
plot(mod, ellipse = F, hull = T,label = T) # 1 sd data ellipse
a.df = data.frame(Type = feature.df[,1])
pmv = vegan::adonis2(d~Type,data = a.df,permutations = 5000)
(pmv)
vegan::anosim(x = d,grouping = df[,1])


##-------------------------*
## PErmutation Test ####
##-------------------------*
## Estat
dt = results.perm$estat
hist(dt)
testStat =  results.emp$estat[which.max(results.emp$norm_esat)]
1-ecdf(dt)(testStat)
mm = mean(dt)
ss = sd(dt)
pval = 1-pnorm(testStat,mean = mm,sd = ss)
x = results.perm$estat
permRep = length(x)
x = rnorm(1000,mean = mean(x),sd(x))
t_stat = results.emp$estat[which.max(results.emp$norm_esat)]
tstat_p = pval
xden = density(x)$x
yden = density(x)$y
xx.df = data.frame(X = xden)
x.df = data.frame(X = xden,Y = yden)
jitter.df = data.frame(X = x,jitter = rnorm(length(x),quantile(yden,probs = .3),sd = 0.005))
library(ggthemes)
tiff(filename = "Figures/bySex_selPermEnergyPermTest.tiff",width = 4,height = 3,units = "in",res = 300)
ggplot(x.df,aes(X,Y))+
  geom_histogram(data = jitter.df,aes(y = ..density..),colour = "white") +
  geom_line(col = "darkgreen",size = 1)+
  geom_point(data = jitter.df,aes(X,jitter),col = "royalblue4",alpha = .1)+
  ggthemes::theme_tufte()+
  ggthemes::geom_rangeframe() +
  geom_vline(xintercept = t_stat,lty = "dashed",col = "red",size = 1)+
  xlab("E-statistic")+
  ylab("Density")+
  labs(caption = paste0("scaled-Estatistic = ",round(t_stat,digits = 4)," , pval = ",round(tstat_p,digits = 4), " (permutations = ",permRep,")"))+
  theme(axis.text = element_text(hjust = 1,size = 8),panel.grid = element_blank(),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 8),
        #legend.direction = "horizontal",
        legend.title = element_text(size = 8),legend.margin =  margin(t = 0, b = 0,unit='cm'),
        #legend.box.background = element_rect(colour = "black"),legend.key = element_rect(colour = "black"),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 8),plot.caption = element_text(size = 10,face = "bold"))
dev.off()

tiff(filename = "Figures/bySex_selPermEnergyKeyFeats.tiff",width = 4,height = 3,units = "in",res = 300)
ggplot(results.emp,aes(feats,norm_esat))+
  geom_point()+
  geom_line()+
  ggthemes::theme_tufte()+
  ggthemes::geom_rangeframe() +
  xlab("# Taxa")+
  ylab("Normalized Energy")+
  geom_vline(xintercept = results.emp$feats[which.max(results.emp$norm_esat)], col = "red", lty = "dashed")+
  scale_x_continuous(breaks = results.emp$feats)+
  theme(axis.text = element_text(hjust = 1,size = 8),panel.grid = element_blank(),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 8),
        #legend.direction = "horizontal",
        legend.title = element_text(size = 8),legend.margin =  margin(t = 0, b = 0,unit='cm'),
        #legend.box.background = element_rect(colour = "black"),legend.key = element_rect(colour = "black"),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 8),plot.caption = element_text(size = 8))
dev.off()



##-----------------------------------------*
## PCA ####
##-----------------------------------------*
library(ggsci)
rt = prcomp(feature.df[,-1])
coords.df = data.frame(Group = feature.df[,1],rt$x)
ggplot(coords.df,aes(PC1,PC2,Fill = Group))+
  stat_density_2d(geom = "polygon", aes(alpha = (..level..)^3, fill = Group),na.rm = T,show.legend = F)+
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

##-----------------------------------------*
## Visualize Signature ####
##-----------------------------------------*
colNames = data.frame(colNames)
feature.df = data.frame(feature.df)
sign = feature.df %>%
  gather("Ratio","Value",2:ncol(feature.df)) %>%
  group_by(Status,Ratio) %>%
  summarise(mnValue = mean(Value),CI = ((1.96*sd(Value))/sqrt(n())),lb = mnValue - CI,ub = mnValue + CI)
sign = data.frame(sign)
order = sign %>%
  filter(Status=='Female') %>%
  arrange(desc(mnValue))
sign$Ratio = factor(sign$Ratio,levels = order$Ratio)
sign = data.frame(sign)
lcolor = if_else(sign$mnValue>0,"blue","red")
## parse ratio
sign = separate(sign,col =2,into = c("numID","denomID"),sep = "___",remove = F )
cn.num = colNames;colnames(cn.num)[2:3] = c("num","numID")
cn.denom = colNames;colnames(cn.denom)[2:3] = c("denom","denomID")
sign = left_join(sign,cn.denom[,2:3],by = "denomID")
sign = left_join(sign,cn.num[,2:3],by = "numID")
sign$ratio_name = paste0("frac(",sign$num, ",", sign$denom,")")


tiff(filename = "Figures/bySex_lrSign_colChart.tiff",width = 3.5,height = 3.5,units = "in",res = 300)
ggplot(sign,aes(mnValue,Status,col = Status,fill = Status))+
  geom_col(alpha = .9,width = .75,col = "black")+
  scale_color_lancet()+
  scale_fill_lancet()+
  geom_errorbar(aes(xmin = lb, xmax = ub),width = .5,size = .5,alpha = 1,col = "black")+
  geom_point(aes(fill = Status),shape = 21,size = 2,col = "black")+
  #geom_vline(xintercept = 0,col = "black",size = 1,lty = "dashed")+
  facet_grid(ratio_name~.,labeller = label_parsed)+
  xlab("Mean logratio")+
  theme_bw()+
  theme(
    axis.line = element_line(),strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0,face = "bold",size = 8,hjust = .5),
    axis.ticks = element_line(colour = "black",size = 1),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8,face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.margin = margin(0,0,0,0,unit = "cm"),
    legend.box.spacing = unit(0,units = "in"),
    legend.box.margin = margin(0,0,0.1,0,unit = "cm"),
    legend.key.height = unit(.1,units = "in"),
    axis.text.x =  element_text(size = 8),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    plot.caption = element_text(size = 8,face = "italic"))
dev.off()



  ##-----------------------------------------*
## Visulaize CEntroids ####
##-----------------------------------------*
## Isomap
phd = parallelDist::parallelDist(as.matrix(feature.df[,-1]),method = "euclidean")
ismp = vegan::isomap(phd ,ndim = 2,k=5)
coords.df = data.frame(Treatment = feature.df[,1],ismp$points)
## Classical MDS
mds = cmdscale(phd,eig = T)
eg = mds$eig
eg[eg<0] = 0
eg = 100*clo(eg)
coords.df = data.frame(Treatment = feature.df[,1],mds$points)
colnames(coords.df)[-1] = paste0("Dim",1:ncol(coords.df[,-1]))
## Compute Centroids
centroid = coords.df %>%
  gather(key = "Dim",value = "Val",2:ncol(coords.df)) %>%
  group_by(Treatment,Dim) %>%
  summarise(mn = mean(Val)) %>%
  spread(key = "Dim",value = "mn")
disp = left_join(coords.df,centroid,by = 'Treatment')
## Plot
tiff(filename = "Figures/bySex_centroids_MDS.tiff",width = 3.5,height = 3.5,units = "in",res = 300)
ggplot(data = coords.df,aes(Dim1,Dim2,col = Treatment,fill = Treatment))+
  geom_hline(yintercept = 0,alpha = .5)+
  geom_vline(xintercept = 0,alpha = .5)+
  geom_segment(aes(x = Dim1.y , y = Dim2.y , xend = Dim1.x , yend = Dim2.x , colour = Treatment),
               data = disp,alpha = .5,size = .5)+
  geom_point(alpha = 1,size = 1)+
  geom_point(data = centroid,aes(Dim1,Dim2,fill = Treatment),size = 3,alpha = 1,pch = 21,col = "black")+
  scale_color_lancet()+
  scale_fill_lancet()+
  xlab(paste0("PCo1 (",round(eg[1],2),"% of Variation)"))+
  ylab(paste0("PCo2 (",round(eg[2],2),"% of Variation)"))+
  theme_bw()+
  labs(fill = "Group",color = "Group")+
  theme(axis.text = element_text(hjust = 1,size = 8),panel.grid = element_blank(),
        axis.title.x = element_text(size = 8,face = "bold"),
        axis.title.y = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 8),
        #legend.direction = "horizontal",
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        #legend.box.background = element_rect(colour = "black"),legend.key = element_rect(colour = "black"),
        legend.position = "top",legend.key.size = unit(.25,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 8),plot.caption = element_text(size = 8))
dev.off()





##-----------------------------------------*
## LR Network ####
##-----------------------------------------*
#### Kruskall Test
krus.test_df = data.frame()
# feats.df =  mst.empirical$features[[1]]
lrs = feature.df[,-1]
cnames = colnames(lrs)
for(i in 1:ncol(lrs)){
  ph = kruskal.test(x = lrs[,i],g  = factor(feature.df[,1]))
  ph = data.frame(Ratio =cnames[i],pval = ph$p.value,Statistic = ph$statistic )
  krus.test_df = rbind(krus.test_df,ph)
}
#False Discovery Rate estimation
krus.test_df$p.adjust = p.adjust(krus.test_df$pval,method = "BH")
pval_level = 0.05
fdrLevel =  max(krus.test_df$p.adjust[krus.test_df$pval <= pval_level])
krus.test_df$signf = if_else(krus.test_df$p.adjust>0.05,F,T)
signRatios = krus.test_df %>%
  filter(signf==T)
#Importance df
colNames = data.frame(ID = colNames$ID,Taxa = colNames$names)

imp.df = data.frame(Ratio = krus.test_df$Ratio,Imp = krus.test_df$Statistic)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)

num = left_join(data.frame(keyRats),data.frame(Num = colNames$ID,Tx = colNames$Taxa))$Tx
denom = left_join(data.frame(keyRats),data.frame(Denom = colNames$ID,Tx = colNames$Taxa))$Tx
keyRats = data.frame(Ratio = paste(num,denom,sep = " / "),Num = num,Denom = denom,Imp = keyRats$Imp)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = if_else((imp.df$Imp)<0,0,(imp.df$Imp))
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")
plot(g,layout = layout_with_fr)
toGephi(g,"bySexEcig")





### LDA Projection
y_train = feature.df[,1]
mdls = trainML_Models(trainLRs =  data.frame(feature.df[,-1]),
                      testLRs =  data.frame(feature.df[,-1]),
                      ytrain = y_train, y_test = y_train,
                      cvMethod = "repeatedcv",mtry_ = 1,numFolds = 10,numRepeats = 20,
                      testIDs = NULL,
                      bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                      models = "pls")
ph = mdls$performance

## Test Predictions ####
predsAll = mdls$predictionMatrix
mdlNames = unique(predsAll$model)
levels = as.character(unique(y_train))
mroc = pROC::multiclass.roc(as.character(y_train),predsAll[,levels])

## Get Predictions Cross Validation ####
mdlName =mdlNames[1]
mdl = mdls$models[[mdlName]]
bestTune = data.frame( mdl$bestTune)
probs = left_join(bestTune,mdl$pred)
probs = probs %>%
  group_by(pred,obs,rowIndex) %>%
  summarise_all(.funs = mean)
rocobj5 = pROC::roc(probs$obs,probs[,minorityClass])
x = pROC::ci.auc(rocobj5)
roc.df = data.frame(sens = rocobj5$sensitivities,spec = 1-rocobj5$specificities)

tiff(filename = "Figures/bySex_AUC.tiff",width = 3,height = 3,units = "in",res = 300)
ggplot(roc.df,aes(spec,sens))+
  geom_line(size = 1)+
  cowplot::theme_cowplot()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  annotate("text", x=.55, y=.45,size = 3,
           label = paste0("AUC: ",round(pROC::auc(rocobj5),digits = 3)," (95% CI: 0.842-0.883)")  ) +
  geom_abline(slope = 1,col = "grey")+
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
       )
dev.off()



# tiff(filename = "Figures/bySex_AUC.tiff",width = 3.5,height = 3.5,units = "in",res = 300)
# par(mar=c(0,0,0,0))
#  pROC::plot.roc(probs$obs,probs[,minorityClass],col= "black",
#                 cex.axis = .75,cex.lab = .75,
#                     print.auc.x =.65, print.auc.y = .3,print.auc.cex = 0.6,
#                     #main = "LOOCV GLM w/FS",
#                     #print.thres=TRUE,print.thres.best.method="youden",
#                     percent=F,
#                     ci = TRUE,                  # compute AUC (of AUC by default)
#                     print.auc = TRUE)
#  dev.off()


 mdl = mdls$models$pls1
loadings = data.frame(Ratio = rownames(glm.mdl1$finalModel$scaling),glm.mdl1$finalModel$scaling)
loadings =  loadings[order(loadings$LD1),]
loadings$Ratio

predictions <-predict(glm.mdl1$finalModel,feature.df[,-1])
lda.df = data.frame(Label = feature.df[,1],predictions$x)
ggplot(lda.df,aes(Label,LD1,col = Label))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = .1,alpha = 0.15,size = 3)+
  theme_bw()+
  scale_color_aaas()+
  theme(legend.position = "none")+
  ggtitle("Training Set LDA Scores")
ggplot(lda.df,aes(Label,LD2,col = Label))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = .1,alpha = 0.15,size = 3)+
  theme_bw()+
  scale_color_aaas()+
  theme(legend.position = "none")+
  ggtitle("Training Set LDA Scores")
ggplot(loadings,aes(reorder(Ratio,LD1),LD1))+
  geom_col(col = "black")+
  #geom_jitter(alpha = 0.5,width = .1)+
  scale_alpha_continuous(range = c(0, 1))+
  scale_shape_manual(values = c(21,24,22))+
  scale_color_lancet()+
  scale_fill_lancet()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  xlab("log-ratio")+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 20),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))







### PLS Analysis
glm.mdl1=train(x = feature.df[,-1],
               y = feature.df[,1],
               method = "pls",tuneGrid = expand.grid(ncomp = 1),
               metric = "AUC",
               trControl = tc1
)
glm.mdl1$finalModel$Xvar
## Get laodings
loading = as.matrix(glm.mdl1$finalModel$loadings)[,1]
loading = data.frame(Ratio = names(loading),Weight = as.vector(loading))
## Get Scores
scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("X",1:ncol(scores))
glm.mdl1$finalModel$loadings
coords.df = data.frame(Group = feature.df[,1],scores)
# loading$Ratio  =str_replace(loading$Ratio ,"___"," / ")


loading = separate(loading,col =1,into = c("numID","denomID"),sep = "___",remove = F )
cn.num = colNames;colnames(cn.num)[] = c("numID","num")
cn.denom = colNames;colnames(cn.denom)[3:4] = c("denomID","denom")
loading = left_join(loading,cn.denom[,3:4],by = "denomID")
loading = left_join(loading,cn.num[3:4],by = "numID")
loading$ratio_name = paste0("frac(",loading$num, ",", loading$denom,")")
loading = loading %>%
  arrange(desc(Weight))
loading$Ratio = factor(loading$Ratio,levels = loading$Ratio)



tiff(filename = "Figures/bySex_plsComps.tiff",width = 3,height = 3,units = "in",res = 300)
ggplot(coords.df,aes(Group,X1,Fill = Group,col = Group))+
  geom_boxplot()+
  geom_jitter(alpha = 0.5,width = .1)+
  scale_alpha_continuous(range = c(0, 1))+
  scale_shape_manual(values = c(21,24,22))+
  scale_color_lancet()+
  scale_fill_lancet()+
  theme_bw()+
  ylab('Comp.1')+
  theme(legend.position = "none",panel.background = element_blank(),
        axis.title.y = element_text(face = "bold",size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8,face = "bold"),
        axis.text.y = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),
        legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))
dev.off()







tiff(filename = "Figures/bySex_plsCompsWeights.tiff",width = 3,height = 3,units = "in",res = 300)
ggplot(loading,aes(Ratio,Weight))+
  geom_col(fill = if_else(loading$Weight>0,pal_lancet()(2)[2],pal_lancet()(2)[1]),col = "black",width = .6)+
  #geom_jitter(alpha = 0.5,width = .1)+
  scale_alpha_continuous(range = c(0, 1))+
  scale_shape_manual(values = c(21,24,22))+
  scale_color_lancet()+
  scale_fill_lancet()+
  theme_bw()+
  coord_flip()+
  scale_x_discrete(labels = parse(text = loading$ratio_name))+
  xlab('Comp.1 Weights')+
  geom_hline(yintercept = 0)+
  theme(legend.position = "top",
        axis.title.x = element_text(face = "bold",size = 8),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        #axis.text.x = element_text(size = 8,angle = 45,vjust = 1,hjust = 1,face = "bold"),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))
dev.off()


ggplot(coords.df,aes(X1,X2,Fill = Group))+
  #stat_density_2d(geom = "polygon", aes(alpha = (..level..)^3, fill = Group),na.rm = T,show.legend = F)+
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





##-----------------------------------------*
## Data integration analysis ####
##-----------------------------------------*
prot.lr = data.frame(ID = prot$ID,calcLogRatio(prot[,-1]))
dat.intg = left_join(prot.lr,lrs.df)
bool = !is.na(rowSums(dat.intg[,-2:-1]))
dat.intg = dat.intg[bool,]

corData = dat.intg[,-2:-1]
cc = combinat::combn2(1:ncol(corData))
cn =colnames(corData)
cor.df = data.frame()
for(i in 1:nrow(cc)){
 x = cc[i,1]
 y = cc[i,2]
 ph = cor.test(corData[,x],corData[,y])
 ph = data.frame(X = cn[x],Y = cn[y],Cor = ph$estimate,p = ph$p.value)
 cor.df = rbind(cor.df,ph)
}
cor.df$m1 = if_else(str_detect(cor.df$X,"__V"),1,0)
cor.df$m2 = if_else(str_detect(cor.df$Y,"__V"),1,0)
cor.df$keep = if_else(cor.df$m1==1 | cor.df$m2==1,1,0)
cor.df$sum = cor.df$m1+cor.df$m2
cor.df = cor.df[cor.df$keep==1,]
cor.df = cor.df[cor.df$sum!=2,]
cor.df$p.adj = p.adjust(cor.df$p,method = "BH")

# Select Data
intFeats = dat.intg[,-2:-1]
### Mbiome Only
y_train = dat.intg$Status
mbiomePositions = str_detect(colnames(intFeats),"_V")
intFeats1 = intFeats[,mbiomePositions]
mdls = trainML_Models(trainLRs =  data.frame(intFeats1),
                      testLRs =  data.frame(intFeats1),
                      ytrain = y_train, y_test = y_train,
                      cvMethod = "repeatedcv",mtry_ = 1,numFolds = 10,numRepeats = 50,
                      testIDs = NULL,
                      bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                      models = "pls")
mdls$performance
res.mbiom = mdls$models$pls1$resample
res.mbiom$Rep = str_split(res.mbiom$Resample,"\\.",simplify = T,n = 2)[,2]
res.mbiom =  aggregate(AUC ~ Rep, data = res.mbiom,
                       FUN = function(x) mean(x) )
res.mbiom$Model = "Mbiome Only"


### Protein and Mbiome Integration
## DCV Feature Selection
cc = unique(dat.intg[,-1]$Status)
cc = combinat::combn2(cc)
dcvScores = data.frame()
for(i in 1:nrow(cc)){
  pwComp = c(cc[i,1],cc[i,2])
  ph = dat.intg[dat.intg$Status%in%pwComp,-1]
  ph = finalDCV(logRatioMatrix = ph,includeInfoGain = T)$lrs
  ph$comb =paste0(cc[i,1],"_",cc[i,2])
  dcvScores = rbind(dcvScores,ph)
}
dcvScores = aggregate(rowmean ~ Ratio, data = dcvScores,
                      FUN = function(x) mean(x) )
dcvScores = dcvScores[order(dcvScores$rowmean,decreasing = T),]
dcvScores = dcvScores[dcvScores$rowmean>0,]
intFeats = subset(intFeats,select = dcvScores$Ratio)
mbiomePositions = str_detect(colnames(intFeats),"_V");intFeats2 = intFeats[,!mbiomePositions]
## Protein + Mbiome
mdls = trainML_Models(trainLRs =  data.frame(intFeats),
                      testLRs =  data.frame(intFeats),
                      ytrain = y_train, y_test = y_train,
                      cvMethod = "repeatedcv",mtry_ = 1,numFolds = 10,numRepeats = 50,
                      testIDs = NULL,
                      bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                      models = "rf")
mdls$performance
varImp(mdls$models[[1]])
res = mdls$models$pls1$resample
res$Rep = str_split(res$Resample,"\\.",simplify = T,n = 2)[,2]
res =  aggregate(AUC ~ Rep, data = res,
                 FUN = function(x) mean(x) )
res$Model = "Mbiome+Protein"
## Protein Only
mdls = trainML_Models(trainLRs =  data.frame(intFeats2),
                      testLRs =  data.frame(intFeats2),
                      ytrain = y_train, y_test = y_train,
                      cvMethod = "repeatedcv",mtry_ = 1,numFolds = 10,numRepeats = 50,
                      testIDs = NULL,
                      bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                      models = "pls")
mdls$performance
varImp(mdls$models$pls1)
res.prot = mdls$models$pls1$resample
res.prot$Rep = str_split(res.prot$Resample,"\\.",simplify = T,n = 2)[,2]
res.prot =  aggregate(AUC ~ Rep, data = res.prot,
                      FUN = function(x) mean(x) )
res.prot$Model = "Protein Only"
wilcox.test(res$AUC,res.mbiom$AUC)
auc.df = rbind(res.mbiom,res,res.prot)
ggplot(auc.df,aes(reorder(Model,AUC),AUC))+
  geom_boxplot()+
  geom_jitter(width = .1)+
  theme_bw()+
  theme(axis.title.x = element_blank())
mean(res$AUC)
mean(res.mbiom$AUC)






##-----------------------------------------*
## Data integration analysis ####
##-----------------------------------------*
prot.lr = data.frame(ID = prot$ID,calcLogRatio(prot[,-1]))
dat.intg = left_join(prot.lr,lrs.df)
bool = !is.na(rowSums(dat.intg[,-2:-1]))
dat.intg = dat.intg[bool,]
cormat = cor(dat.intg[,-2:-1],method = "spearman")
intFeats = dat.intg[,-2:-1]

### Mbiome Only
y_train = dat.intg$Status
mbiomePositions = str_detect(colnames(intFeats),"_V")
intFeats1 = intFeats[,mbiomePositions]
mdls = trainML_Models(trainLRs =  data.frame(intFeats1),
                      testLRs =  data.frame(intFeats1),
                      ytrain = y_train, y_test = y_train,
                      cvMethod = "repeatedcv",mtry_ = 1,numFolds = 10,numRepeats = 50,
                      testIDs = NULL,
                      bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                      models = "pls")
mdls$performance
res.mbiom = mdls$models$pls1$resample
res.mbiom$Rep = str_split(res.mbiom$Resample,"\\.",simplify = T,n = 2)[,2]
res.mbiom =  aggregate(AUC ~ Rep, data = res.mbiom,
                       FUN = function(x) mean(x) )
res.mbiom$Model = "Mbiome Only"


### Protein and Mbiome Integration
## DCV Feature Selection
cc = unique(dat.intg[,-1]$Status)
cc = combinat::combn2(cc)
dcvScores = data.frame()
for(i in 1:nrow(cc)){
  pwComp = c(cc[i,1],cc[i,2])
  ph = dat.intg[dat.intg$Status%in%pwComp,-1]
  ph = finalDCV(logRatioMatrix = ph,includeInfoGain = T)$lrs
  ph$comb =paste0(cc[i,1],"_",cc[i,2])
  dcvScores = rbind(dcvScores,ph)
}
dcvScores = aggregate(rowmean ~ Ratio, data = dcvScores,
                      FUN = function(x) mean(x) )
dcvScores = dcvScores[order(dcvScores$rowmean,decreasing = T),]
dcvScores = dcvScores[dcvScores$rowmean>0,]
intFeats = subset(intFeats,select = dcvScores$Ratio)
mbiomePositions = str_detect(colnames(intFeats),"_V");intFeats2 = intFeats[,!mbiomePositions]
## Protein + Mbiome
mdls = trainML_Models(trainLRs =  data.frame(intFeats),
                      testLRs =  data.frame(intFeats),
                      ytrain = y_train, y_test = y_train,
                      cvMethod = "repeatedcv",mtry_ = 1,numFolds = 10,numRepeats = 50,
                      testIDs = NULL,
                      bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                      models = "pls")
mdls$performance
varImp(mdls$models$pls1)
res = mdls$models$pls1$resample
res$Rep = str_split(res$Resample,"\\.",simplify = T,n = 2)[,2]
res =  aggregate(AUC ~ Rep, data = res,
                 FUN = function(x) mean(x) )
res$Model = "Mbiome+Protein"
## Protein Only
mdls = trainML_Models(trainLRs =  data.frame(intFeats2),
                      testLRs =  data.frame(intFeats2),
                      ytrain = y_train, y_test = y_train,
                      cvMethod = "repeatedcv",mtry_ = 1,numFolds = 10,numRepeats = 50,
                      testIDs = NULL,
                      bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                      models = "pls")
mdls$performance
varImp(mdls$models$pls1)
res.prot = mdls$models$pls1$resample
res.prot$Rep = str_split(res.prot$Resample,"\\.",simplify = T,n = 2)[,2]
res.prot =  aggregate(AUC ~ Rep, data = res.prot,
                      FUN = function(x) mean(x) )
res.prot$Model = "Protein Only"
wilcox.test(res$AUC,res.mbiom$AUC)
auc.df = rbind(res.mbiom,res,res.prot)
ggplot(auc.df,aes(reorder(Model,AUC),AUC))+
  geom_boxplot()+
  geom_jitter(width = .1)+
  theme_bw()+
  theme(axis.title.x = element_blank())
mean(res$AUC)
mean(res.mbiom$AUC)


### relationship to pls comp
y_train = dat.intg$Status
glm.mdl1=train(x = intFeats1,
               y = y_train,
               method = "pls",tuneGrid = expand.grid(ncomp = 1),
               metric = "AUC",
               trControl = tc1
)
glm.mdl1
## Get laodings
loading = as.matrix(glm.mdl1$finalModel$loadings)[,1]
loading = data.frame(Ratio = names(loading),Weight = as.vector(loading))
## Get Scores
scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)
coords.df = cbind(coords.df,prot.lr[bool,-2:-1])
yy = coords.df[,2]
xx = coords.df[,-2:-1]
tcr  <- trainControl(method="repeatedcv",
                     repeats = 10,
                     number = 10,
                     seeds = NULL,
                     classProbs = F,
                     savePredictions = T,
                     allowParallel = TRUE
                     #summaryFunction = caret::
                     #summaryFunction = prSummary
)
glm.mdl1=train(x = xx,
               y = yy,
               method = "nnet",
               #metric = "AUC",
               trControl = tcr
)
glm.mdl1

ph = data.frame(tar = yy,xx)
cm = cor(ph)
