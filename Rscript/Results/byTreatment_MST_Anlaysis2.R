rm(list=ls())
gc()

### Install/Load Required Packages  ####
library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(glmnet)
library(vegan)
library(keras)
library(PRROC)
library(energy)
library(dplyr)
library(foreach)
library(parallel)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggsci)

## Start Cluster
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)

source(file = "Functions/functions1.R")


## Model Traing Setup
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
## Load Data
colNames = read_csv("Output/colNamesTable_byGenus.csv")
colNames = separate(colNames,col = 1,into = c('Family','Genus'),sep = "g__")
colNames$Family = str_replace(colNames$Family,"f__","")
colNames$Family = str_replace_all(colNames$Family,"\\.","")
Name = if_else(stringr::str_length(colNames$Genus )==0,colNames$Family,colNames$Genus)
colNames$names = Name
## data
dat = data.frame( read_csv(file = "Output/byTreatment_adjSex_byGenus_80sparse.csv") )
table(dat[,1])
sampleID = data.frame( ID = data.frame(read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv")[,1] )[,1])
prot = data.frame(read_csv(file = "Output/NLF_bySubjectGroup.csv"));colnames(prot) = str_replace_all(colnames(prot),"..ug.mL.",replacement = "")
procData = processCompData(dat)
dat = procData$processedData
impFact = procData$impFactor
minorityClass = procData$minClss
majorityClass = procData$majClass
impClass = "Smoker"


## Protein COncentration Variablility
ph = data.frame(Concetration = rowSums(prot[,-2:-1]))
tiff(filename = "Figures/byTreatment_proteinConcentration.tiff",width = 3,height = 3,units = "in",res = 300)
ggplot(ph,aes(Concetration))+
  geom_histogram(col = "white") +
  ggthemes::theme_tufte()+
  ggthemes::geom_rangeframe() +
  xlab("Total Concentration (ug/mL)")+
  ylab("Count")+
  #labs(caption = paste0("scaled-Estatistic = ",round(t_stat,digits = 4)," , pval = ",round(tstat_p,digits = 4), " (permutations = ",permRep,")"))+
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



## load empirical Results
results.emp = read_csv("Output/byTreatGenusadjSex_NasalMicrobiome_80sparse_20_repeatedCV_MST_empirical.csv")
results.perm = read_csv("Output/byTreatGenusadjSex_NasalMicrobiome_80sparse_20_repeatedCV_MST_permuted.csv")
## Load Features
feature.df = read_csv(file = "Output/byTreatGenusadjSex_NasalMicrobiome_80sparse_20_keyFeatures.csv")
feature.df = data.frame(feature.df)

## INtegrate  Data
lrs.df = getLogratioFromList(Ratio = colnames(feature.df[,-1]),raMatrix = dat[,-1],Class = dat[,1])
#feature.df = data.frame(Status = dat[,1],lrs.df)
lrs.df = data.frame(ID  = sampleID$ID,lrs.df)
dat.intg = left_join(prot,lrs.df)
bool = !is.na(rowSums(dat.intg[,-2:-1]))
keepID = data.frame(ID = dat.intg$ID[bool])
metadata = read.csv("Data/2019_09_13 ELFmetadata.csv")
colnames(metadata)[1] = "ID"
metadata = left_join(keepID,metadata)
dat.intg = dat.intg[bool,]
cormat = cor(dat.intg[,-2:-1],method = "spearman")


##-------------------------*
## PErmutation Test
##-------------------------*
library("ggthemes")## Estat
hist(results.perm$estat)
testStat =  results.emp$estat[which.max(results.emp$norm_esat)]
dt = results.perm$estat
1-ecdf(dt)(testStat)
mm = mean(results.perm$estat)
ss = sd(results.perm$estat)
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

tiff(filename = "Figures/byTreatment_selPermEnergyPermTest.tiff",width = 4,height = 3,units = "in",res = 300)
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

tiff(filename = "Figures/byTreatment_selPermEnergyKeyFeats.tiff",width = 4,height = 3,units = "in",res = 300)
ggplot(results.emp,aes(feats,norm_esat))+
  geom_point()+
  geom_line()+
  theme_tufte()+
  geom_rangeframe() +
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
sign = feature.df %>%
  gather("Ratio","Value",2:ncol(feature.df)) %>%
  group_by(Status,Ratio) %>%
  summarise(mnValue = mean(Value),CI = ((1.96*sd(Value))/sqrt(n())),lb = mnValue - CI,ub = mnValue + CI)
order = sign %>%
  filter(Status==impClass) %>%
  arrange(desc(mnValue))
sign$Ratio = factor(sign$Ratio,levels = order$Ratio)
lcolor = if_else(sign$mnValue>0,"blue","red")
heatmap.dat = sign
## parse ratio
sign = separate(sign,col =2,into = c("numID","denomID"),sep = "___",remove = F )
cn.num = colNames;colnames(cn.num)[3:4] = c("numID","num")
cn.denom = colNames;colnames(cn.denom)[3:4] = c("denomID","denom")
sign = left_join(sign,cn.denom[,3:4],by = "denomID")
sign = left_join(sign,cn.num[3:4],by = "numID")
sign$ratio_name = paste0("frac(",sign$num, ",", sign$denom,")")

write_csv(sign,path = "Output/suppTables/byTreatment_indvAnlaysisofLogRatios.csv")
tiff(filename = "Figures/byTreatment_lrSign_colChart.tiff",width = 3.5,height = 7.1,units = "in",res = 300)
ggplot(sign,aes(mnValue,Status,col = Status,fill = Status))+
  geom_col(alpha = .9,width = .75,col = "black")+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  geom_errorbar(aes(xmin = lb, xmax = ub),width = .5,size = .5,alpha = 1,col = "black")+
  geom_point(aes(fill = Status),shape = 21,size = 2,col = "black")+
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
## Visualize CEntroids ####
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
tiff(filename = "Figures/byTreatment_centroids_MDS.tiff",width = 3,height = 3,units = "in",res = 300)
ggplot(data = coords.df,aes(Dim1,Dim2,col = Treatment,fill = Treatment))+
  geom_hline(yintercept = 0,alpha = .5)+
  geom_vline(xintercept = 0,alpha = .5)+
  geom_segment(aes(x = Dim1.y , y = Dim2.y , xend = Dim1.x , yend = Dim2.x , colour = Treatment),
               data = disp,alpha = .5,size = .5)+
  geom_point(alpha = 1,size = 1)+
  geom_point(data = centroid,aes(Dim1,Dim2,fill = Treatment),size = 3,alpha = 1,pch = 21,col = "black")+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
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
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 8),plot.caption = element_text(size = 8))
dev.off()

##-----------------------------------------*
### PLS Projection
##-----------------------------------------*

##---------------------------------*
### PLS ####
##---------------------------------*
num_comps = 2
glm.mdl1=train(x = feature.df[,-1],
               y = feature.df[,1],probMethod = "bayes",
               method = "pls",tuneGrid = expand.grid(ncomp = num_comps),
               metric = "AUC",
               trControl = tc1
)
vi = varImp(glm.mdl1,scale = F)
vi = data.frame(Ratio  =rownames(vi$importance),clo(vi$importance))
vi = gather(vi,"Group",Value,2:ncol(vi))
vi.smoker = vi[vi$Group==impClass,]
vi.smoker = vi.smoker[order(vi.smoker$Value,decreasing = T),]
vi$Ratio = factor(vi$Ratio,levels = vi.smoker$Ratio)
##parse ratio
vi = separate(vi,col =1,into = c("numID","denomID"),sep = "___",remove = F )
cn.num = colNames;colnames(cn.num)[3:4] = c("numID","num")
cn.denom = colNames;colnames(cn.denom)[3:4] = c("denomID","denom")
vi = left_join(vi,cn.denom[,3:4],by = "denomID")
vi = left_join(vi,cn.num[3:4],by = "numID")
vi$ratio_name = paste0("frac(",vi$num, ",", vi$denom,")")
nn = left_join(vi.smoker,vi)

tiff(filename = "Figures/byTreatment_varImp.tiff",width = 3.35,height = 4.1,units = "in",res = 300)
ggplot(vi,aes(Ratio,Value,fill  = Group))+
  geom_col(col = "black",width = .5)+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  coord_flip()+
  scale_x_discrete(labels = parse(text = nn$ratio_name))+
  ylab("Relative Importance")+
  labs(fill = "Group",color = "Group")+
  theme(axis.text.y  = element_text(hjust = .5,size = 7),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 8),
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        axis.text.x = element_text(size = 8,face = "bold",hjust = .5)
        ,plot.caption = element_text(size = 8))
dev.off()

## Get laodings
loading = as.matrix(glm.mdl1$finalModel$loading.weights)[,1:num_comps]
loading = data.frame(Ratio = rownames(loading),(loading))
loading = gather(loading,'Component',"Value",2:ncol(loading))
## Get Scores
scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("X",1:ncol(scores))
glm.mdl1$finalModel$loadings
coords.df = data.frame(Group = feature.df[,1],scores)
colnames(coords.df)[-1] =paste("Comp",1:ncol(coords.df[,-1]))
coords.df1 = gather(coords.df,'Component',"Value",2:ncol(coords.df))
coords.df1$Group = factor(coords.df1$Group)

ggplot(coords.df1,aes(Group,Value))+
  geom_hline(yintercept = 0,lty = "dashed",col = "grey40")+
  geom_boxplot(aes( fill = Group),outlier.shape = NA)+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  geom_jitter(alpha = 0.5,width = .075,size = 1)+
  # scale_fill_viridis_d(end = .75,option = "C")+
  # scale_color_viridis_d(end = .75,option = "C")+
  #scale_fill_d3()+
  theme_bw()+
  facet_grid(.~Component)+
  theme(legend.position = "none",
        axis.title =element_blank(),
        axis.text.x =  element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))


ggplot(loading,aes(reorder(Ratio,Value),Value))+
  geom_hline(yintercept = 0,lty = "dashed",col = "grey40")+
  geom_col()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  # scale_fill_viridis_d(end = .75,option = "C")+
  # scale_color_viridis_d(end = .75,option = "C")+
  #scale_fill_d3()+
  theme_bw()+
  facet_grid(.~Component,scales = "free")+
  theme(legend.position = "none",
        axis.title =element_blank(),
        axis.text.x =  element_text(size = 8,angle = 45,hjust = 1),
        axis.text.y = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))






###----------------------------------------*
### Heatmap of signature
###----------------------------------------*
feats = heatmap.dat[,1:3]
feats = spread(feats,"Status","mnValue")
feats = data.frame(feats)

# log ratio matrix
mat = as.matrix(feats[,-1])
rownames(mat) = feats[,1]
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(1000)
# heatmap microbiome composition
tiff(filename = "Figures/byTreatment_lrSign.tiff",width = 2.5,height = 3.25,units = "in",res = 300)
gplots::heatmap.2( mat,
                   col = viridis::viridis(n = 1000,option = "D"),
                   Rowv = TRUE,Colv = F,
                   margins = c(4,6), lwid = c(1,4),lhei = c(0.1,4),
                   sepwidth=c(0.01,0.01),keysize = 3,
                   sepcolor="white",colsep=1:ncol(mat),rowsep=1:nrow(mat),
                   hclustfun = function(x) hclust(x, method = "average"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=F, symkey=TRUE, srtCol = 45,
                   scale = "row",
                   density.info="density", trace="none",
                   #ColSideColors = sc,
                   #main = "Presensitization Gut Mircobiome",
                   #key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",reference$i,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = .75,cexCol = .75,
                   #adjCol = c(1,1)
)
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
toGephi(g,"byTreatmentEcig")






##-----------------------------------------*
## ROC Analysis ####
##-----------------------------------------*
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
mroc = pROC::multiclass.roc(as.character(probs$obs),probs[,levels])
pROC::auc(mroc)
x = pROC::ci.auc(mroc)
pROC::ci(mroc)
rs <-mroc[['rocs']]
nms = names(rs)

## Bootstrap Confdidcence Interval multiclass
mcAUC = c()
nreps = 1000
for(i in 1:nreps){
  ph = sample_n(probs,size = nrow(probs),replace = T)
  mroc = pROC::multiclass.roc(as.character(ph$obs),ph[,levels])
  mcAUC[i] = pROC::auc(mroc)
}
quantile(mcAUC,probs = .025)
quantile(mcAUC,probs = .975)
mean(mcAUC)
mcx = paste0("Multiclass AUC:",round(mean(mcAUC),digits = 3)," (95% CI:",round(quantile(mcAUC,probs = .025),digits = 3),"-",round(quantile(mcAUC,probs = .975),digits = 3),")")

roc.df = data.frame()
for(i in 1:length(nms)){
  rc = rs[[i]]
  x = pROC::ci.auc(rc[[1]])
  x = paste0(" AUC=",round(x[2],digits = 3)," (95% CI:",round(x[1],digits = 3),"-",round(x[3],digits = 3),")")
  ph = data.frame(sens = rc[[1]]$sensitivities,spec = 1 - rc[[1]]$specificities)
  ph$Comp = paste0(nms[i],x)
  roc.df = rbind(roc.df,ph)
}
## Figure
tiff(filename = "Figures/byTreatment_AUC.tiff",width = 3,height = 4.1,units = "in",res = 300)
ggplot(roc.df,aes(spec,sens,col = Comp,group = Comp))+
  geom_step(size = .75)+
  #geom_line(size = .75)+
  cowplot::theme_cowplot()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  annotate("text", x=.03, y=1.1,size = 2.5, hjust = 0, fontface = "bold",
           label = paste(mcx) ) +
  geom_abline(slope = 1,col = "grey")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.line = element_line(size = .5),
        legend.margin=margin(-1,-10,-10,-10),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.175, "cm"),
        legend.key.size = unit(0.25, "cm"),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )+
  guides(color=guide_legend(nrow=3,byrow=TRUE))
dev.off()







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
table(y_train)
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
mdls$models$pls1

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
  ph = finalDCV(logRatioMatrix = ph,includeInfoGain = T,nfolds = 10,numRepeats = 20)$lrs
  ph$comb =paste0(cc[i,1],"_",cc[i,2])
  dcvScores = rbind(dcvScores,ph)
}
dcvScores1 = aggregate(rowmean ~ Ratio, data = dcvScores,
          FUN = function(x) mean(x) )

dcvScores2 = dcvScores1[order(dcvScores1$rowmean,decreasing = T),]
dcvScores2 = dcvScores2[dcvScores2$rowmean>0,]
intFeats = subset(intFeats,select = dcvScores2$Ratio)

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
x = data.frame(intFeats2)
mdls = trainML_Models(trainLRs =  x,
                      testLRs =  x,
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

wilcox.test(res.prot$AUC,res.mbiom$AUC)
wilcox.test(res.prot$AUC,res$AUC)
wilcox.test(res$AUC,res.mbiom$AUC)
auc.df = rbind(res.mbiom,res,res.prot)
p.adjust(c(wilcox.test(res.prot$AUC,res.mbiom$AUC)$p.value,
           wilcox.test(res$AUC,res.prot$AUC)$p.value,
           wilcox.test(res$AUC,res.mbiom$AUC)$p.value),method = "BH")

##-------------------------*
### Mediator Response ####
##--------------------------*
## plot scores
ph = dcvScores1
ph = ph[str_detect(ph$Ratio,"V",negate = T),]
ph$Ratio = str_replace(ph$Ratio,"___"," / ")
tiff(filename = "Figures/byTreatment_mediatorDCV.tiff",width = 3,height = 4,units = "in",res = 300)
ggplot(ph,aes(reorder(Ratio,rowmean),rowmean))+
  geom_col()+
  coord_flip()+
  ylab("DCV")+
  ggthemes::theme_tufte()+
  ggthemes::geom_rangeframe() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8,face = "bold"))
dev.off()
##-----------------------------------------*
## Visualize Mediator Signature ####
##-----------------------------------------*
ph = data.frame(Status = y_train,intFeats2)
sign =  ph %>%
  gather("Ratio","Value",2:ncol(ph))
sign$Ratio = str_replace_all(sign$Ratio,pattern = "___",replacement = " / ")
my_comparisons = list(c("Ecig","Nonsmoker"),c("Ecig","Smoker"),c("Nonsmoker","Smoker"))
rts = unique(sign$Ratio)
for( i in 1:length(rts)){
  ph1 = sign[sign$Ratio==rts[i],]
  fname = paste0("Figures/byTreatment_mediatorBoxplot",i,".tiff")
  #kruskall walis test
  tiff(filename = fname,width = 3,height = 2,units = "in",res = 300,bg = "transparent")
  p = ggpubr::ggboxplot(ph1, x = "Status", y = "Value",
                        outlier.shape = NA,
                        color = "Status",title = rts[i],
                        palette = RColorBrewer::brewer.pal(3,name = "Set2"),
                        xlab = "Ratio",ylab = "Value",
                        add = "jitter",add.params = list(size = 3,shape = 21,linetype = "black",alpha = .5,width = .01),
                        ggtheme = theme(legend.position = "none",
                                        axis.line = element_line(),axis.title = element_text(size = 8),text = element_text(size = 8),
                                        axis.text = element_text(face = "bold",size = 8),
                                        axis.title.x = element_blank(),
                                        axis.title.y = element_text(size = 8,face = "bold"),
                                        panel.background = element_blank()
                        ))+
    ggpubr::stat_compare_means(size = 2.5,label.y = 4)  +
    ggpubr::stat_compare_means(comparisons = my_comparisons,size = 2.5,
                               label = "p",
                               #label.y = 3.3,
                               tip.length = 0)# Add pairwise comparisons p-value
  idx <- which(sapply(p$layers, function(l) "PositionJitter" %in% class(l$position)))
  p$layers[[idx]]$aes_params$size <- 1
  p
  dev.off()
}



##-----------------------------------------*
## Medaitor ROC Analysis ####
##-----------------------------------------*
y_train = dat.intg$Status
mdls = trainML_Models(trainLRs =  intFeats2,
                      testLRs =  intFeats2,
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
mroc = pROC::multiclass.roc(as.character(probs$obs),probs[,levels])
pROC::auc(mroc)
x = pROC::ci.auc(mroc)
pROC::ci(mroc)
rs <-mroc[['rocs']]
nms = names(rs)
## Bootstrap Confdidcence Interval multiclass
mcAUC = c()
nreps = 1000
for(i in 1:nreps){
  ph = sample_n(probs,size = nrow(probs),replace = T)
  mroc = pROC::multiclass.roc(as.character(ph$obs),ph[,levels])
  mcAUC[i] = pROC::auc(mroc)
}
quantile(mcAUC,probs = .025)
quantile(mcAUC,probs = .975)
mean(mcAUC)
mcx = paste0("Multiclass AUC:",round(mean(mcAUC),digits = 3)," (95% CI:",round(quantile(mcAUC,probs = .025),digits = 3),"-",round(quantile(mcAUC,probs = .975),digits = 3),")")

roc.df = data.frame()
for(i in 1:length(nms)){
  rc = rs[[i]]
  x = pROC::ci.auc(rc[[1]])
  x = paste0(" AUC=",round(x[2],digits = 3)," (95% CI:",round(x[1],digits = 3),"-",round(x[3],digits = 3),")")
  ph = data.frame(sens = rc[[1]]$sensitivities,spec = 1 - rc[[1]]$specificities)
  ph$Comp = paste0(nms[i],x)
  roc.df = rbind(roc.df,ph)
}
## Figure
tiff(filename = "Figures/byTreatment_mediatorAUC.tiff",width = 3,height = 2.5,units = "in",res = 300)
ggplot(roc.df,aes(spec,sens,col = Comp,group = Comp))+
  geom_step(size = .75)+
  #geom_line(size = .75)+
  cowplot::theme_cowplot()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  annotate("text", x=.03, y=.2,size = 2.25, hjust = 0, fontface = "bold",
           label = paste(mcx) ) +
  geom_abline(slope = 1,col = "grey")+
  theme(legend.position = "top",legend.title = element_blank(),
        axis.line = element_line(size = .5),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.35, "cm"),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
  )+
  guides(color=guide_legend(nrow=3,byrow=TRUE))
dev.off()



#### Kruskall Test
# ph2 = data.frame(Status = y_train,intFeats1)
# ph2 = gather(ph2,"Ratio",)
# #ph2 = data.frame(Status = y_train,Value = intFeats2[,1])
# pairwise.test = ph2 %>%
#   rstatix::wilcox_test(Value ~ Status) %>%
#   adjust_pvalue(method = 'BH')
# pairwise.test$comp = colnames(intFeats2)[1]
#
# krusTest = data.frame()
# for(i in 1:ncol(intFeats2))
# pairwise.test = ph2 %>%
#   rstatix::kruskal_test(Value ~ Status) %>%
#   adjust_pvalue(method = 'BH')
# pairwise.test$comp = colnames(intFeats2)[1]


krus.test_df = data.frame()
pairwise.test  =data.frame()
# feats.df =  mst.empirical$features[[1]]
lrs = intFeats1
cnames = colnames(lrs)
for(i in 1:ncol(lrs)){
  ph = kruskal.test(x = lrs[,i],g  = factor(y_train))
  ph = data.frame(Ratio =cnames[i],pval = ph$p.value,Statistic = ph$statistic )
  pwtest = pairwise.wilcox.test(lrs[,i],g = y_train,p.adjust.method = "BH")
  pd = data.frame(Ratio = cnames[i],g1 = rownames(pwtest$p.value),pwtest$p.value)
  pd = gather(pd,"g2","p",3:4)
  pd   = pd[!is.na(pd$p),]
  pairwise.test = rbind(pairwise.test,pd)
  krus.test_df = rbind(krus.test_df,ph)
}
write_csv(pairwise.test,path = "Output/byTreatment_pairwiseWilcoxonTest_Mbiome.csv")


#False Discovery Rate estimation
krus.test_df$p.adjust = p.adjust(krus.test_df$pval,method = "BH")
pval_level = 0.1
fdrLevel =  max(krus.test_df$p.adjust[krus.test_df$pval <= 0.05])
### suplemental table mediators #####
krus.test_df$signf = if_else(krus.test_df$p.adjust>pval_level,F,T)
signRatios = krus.test_df %>%
  filter(signf==T)
imp.df = data.frame(Ratio = krus.test_df$Ratio,Imp = krus.test_df$Statistic)
keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g)$weight = if_else((imp.df$Imp)<0,0,(imp.df$Imp))
#g = minimum.spanning.tree(graph = g,weights = -E(g)$weight)
w_degree  = strength(g,mode = "total")

tiff(filename = "Figures/byTreatment_mediatorNetwork.tiff",width = 3,height = 3,units = "in",res = 300)
par(mar=c(1,1,1,1)+.5)
set.seed(1)
plot(g,layout = layout_with_fr,vertex.size = 10,vertex.label.cex = .85,edge.curved = .2,
     edge.arrow.size = .5,edge.arrow.width = 1)
dev.off()


##---------------------------*
## permanova protein feats ####
##-------------------------****
#dispersion test
d = dist(intFeats2)
mod = vegan::betadisper(d,group = y_train)
anova(mod)
permutest(mod)
mod.HSD <- TukeyHSD(mod)
plot(mod,)
boxplot(mod)
plot(TukeyHSD(mod))
plot(mod, ellipse = F, hull = F,label = T,label.cex = .5) # 1 sd data ellipse
#permanova
a.df = data.frame(Type = y_train)
pmv = vegan::adonis2(d~Type,data = a.df,permutations = 1000)
(pmv)

##---------------------------*
## permanova mbiome feats ####
##-------------------------****
#dispersion test
d = dist(intFeats1)
mod = vegan::betadisper(d,group = y_train)
anova(mod)
permutest(mod)
mod.HSD <- TukeyHSD(mod)
plot(mod,)
boxplot(mod)
plot(TukeyHSD(mod))
tiff(filename = "Figures/byTreatment_permanovaProteins.tiff",width = 3,height = 2,units = "in",res = 300)
plot(mod, ellipse = F, hull = T,label = T) # 1 sd data ellipse
dev.off()
#permanova
a.df = data.frame(Type = y_train)
pmv = vegan::adonis2(d~Type,data = a.df,permutations = 5000)
(pmv)

#vegan::anosim(x = d,grouping = labels)


###----------------------------------------*
### Heatmap of NLF Protein signature ####
###----------------------------------------*
# log ratio matrix
mat = as.matrix(intFeats2)
#mat = mat[dat.intg$Status=="Ecig",]
col
col  =RColorBrewer::brewer.pal(3,"Set2")
cc = dat.intg$Status
cc[cc=="Smoker"] = col[1]
cc[cc=="Nonsmoker"] = col[2]
cc[cc=="Ecig"] = col[3]
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(1000)
# heatmap microbiome composition
tiff(filename = "Figures/byTreatment_lrSignalProteins.tiff",width = 7,height = 7,units = "in",res = 300)
gplots::heatmap.2( mat,
                   col = col,#viridis::viridis(n = 1000,option = "A"),
                   Rowv = TRUE,Colv = T,
                   margins = c(10,6), lwid = c(1,4),lhei = c(0.1,4),
                   sepwidth=c(0.01,0.01),keysize = 3,
                   sepcolor="white",colsep=1:ncol(mat),rowsep=1:nrow(mat),
                   hclustfun = function(x) hclust(x, method = "average"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=F, symkey=TRUE, srtCol = 45,
                   scale = "col",
                   density.info="density", trace="none",
                   #ColSideColors = cc,
                   #main = "Presensitization Gut Mircobiome",
                   #key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",reference$i,")",sep = "" ),
                   RowSideColors = cc,
                   cexRow = .85,cexCol = .85,
                   #adjCol = c(1,1)
)
dev.off()








##------------------------------------------*
### Model Performance Comparisions ####
##--------------------------------------------*
auc.df1 =  aggregate(AUC ~ Model, data = auc.df,
                    FUN = function(x) mean(x) )
auc.df.sd =  aggregate(AUC ~ Model, data = auc.df,
                     FUN = function(x) sd(x) )
auc.df1$CI = (1.96*auc.df.sd$AUC) / (sqrt(nrow(res.mbiom)))
auc.df1$lb = auc.df1$AUC-auc.df1$CI
auc.df1$ub = auc.df1$AUC+auc.df1$CI


my_comparisons = list(c("Mbiome Only","Mbiome+Protein"),c("Mbiome Only","Protein Only"),c("Protein Only","Mbiome+Protein"))
auc.df$Model = factor(auc.df$Model,levels = c("Protein Only","Mbiome Only","Mbiome+Protein"))
#kruskall walis test
tiff(filename = "Figures/byTreatment_aucIntegration.tiff",width = 2.5,height = 3.25,units = "in",res = 300,bg = "transparent")
p = ggpubr::ggboxplot(auc.df, x = "Model", y = "AUC",
                      outlier.shape = NA,
                      color = "Model",
                      palette = rep("black",3),
                      xlab = "AUROC",ylab = "AUC",
                      add = "jitter",add.params = list(size = 3,shape = 21,linetype = "black",alpha = .5,width = .01),
                      ggtheme = theme(legend.position = "none",
                                      axis.line = element_line(),
                                      axis.text = element_text(face = "bold",size = 8,angle = 45,hjust = 1),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_text(size = 8,face = "bold"),
                                      panel.background = element_blank()
                      ))+
  ggpubr::stat_compare_means(comparisons = my_comparisons,size = 2.5,
                             label = "p",
                             #label.y = 3.3,
                             tip.length = 0)# Add pairwise comparisons p-value
idx <- which(sapply(p$layers, function(l) "PositionJitter" %in% class(l$position)))
p$layers[[idx]]$aes_params$size <- 1
p
dev.off()




###-----------------------------*
### relationship to pls comp
###-----------------------------*
y_train = dat.intg$Status
glm.mdl1=train(x =  data.frame(intFeats1),
               y = y_train,
               method = "pls",tuneGrid = expand.grid(ncomp = 2),
               metric = "AUC",
               trControl = tc1
)
glm.mdl1
## Get laodings
loading = as.matrix(glm.mdl1$finalModel$loading.weights)[,1]
loading = data.frame(Ratio = names(loading),Weight = as.vector(loading))
loading$Ratio  =str_replace(loading$Ratio,"___"," / ")

tiff(filename = "Figures/byTreatment_loadingDatIntegration.tiff",width = 3.25,height = 2,units = "in",res = 300,bg = "transparent")
ggplot(loading,aes(reorder(Ratio,Weight),Weight))+
  geom_col()+
  theme_classic()+
  theme(axis.text.y = element_text(hjust = ,size = 8),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8,face = "bold",angle = 45,hjust = 1),
        axis.title.y = element_text(size = 8,face = "bold"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        plot.caption = element_text(size = 8))
dev.off()

## Yhat
glm.mdl1=train(x =  data.frame(intFeats1),
               y = y_train,
               method = "pls",tuneGrid = expand.grid(ncomp = 2),
               metric = "AUC",
               trControl = tc1
)
ldaScaling = glm.mdl1$finalModel$scaling
preProc  <- preProcess(intFeats1,method = "center")
A <- as.matrix(predict(preProc, intFeats1))
biplot.Loading = glm.mdl1$finalModel$loading.weights[,]
biplot.Loading = data.frame(Ratio = rownames(biplot.Loading),biplot.Loading )
loading_ =glm.mdl1$finalModel$projection
loading_ = data.frame(Ratio = rownames(loading_),Comp1 = loading_[,1],Comp2 = loading_[,2],row.names = 1:nrow(loading_))
x = as.matrix(loading_[,-1]);colnames(x) = paste0("Comp",1:2)
T.test = A %*% x
Q = as.matrix(glm.mdl1$finalModel$Yloadings[1,])
Yhat = (T.test %*% (Q))
yy = data.frame(Status = y_train,yhat = Yhat)
ggplot(yy,aes(Status,yhat))+
  geom_boxplot()

glm.mdl1$finalModel$Yscores

## Get Scores
scores = glm.mdl1$finalModel$Yscores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)
yy = coords.df[,2:3]

scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)
yy = coords.df[,2:3]




# coords.df = cbind(coords.df,prot.lr[bool,-2:-1])
# yy = coords.df[,2:3]
# xx = coords.df[,-3:-1]
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
corData =  data.frame( y = intFeats2$DEFB4A.2___Neutrophil.Elastase,yy)
t1=train(x = corData[,-1],
               y = (corData[,1]),
               method = "pls",
               #metric = "AUC",
               trControl = tcr
)
t1$fin
# ph = data.frame(tar = yy[,1],intFeats2)
# cm = cor(ph)


###------------------------------------------*
## data integrated correlation analysis ####
###-------------------------------------------*
## use microbial signature only
glm.mdl1=train(x =  data.frame(intFeats1),
               y = y_train,
               method = "pls",tuneGrid = expand.grid(ncomp = 2),
               metric = "AUC",
               trControl = tc1
)
plsVariation = glm.mdl1$finalModel$Xvar/glm.mdl1$finalModel$Xtotvar
## Yscores
scores = glm.mdl1$finalModel$Yscores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)
yy = coords.df[,2:3]
## Scores
scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)
yy = coords.df[,2:3]
myrda <- rda(yy,intFeats2)
plot(myrda,scaling=2)

corData =  data.frame( yy,intFeats2)
md = lm(Lactoferrin___IL.8~Comp.1,data = corData)
md = lm(Lactoferrin___IL.8~poly(Comp.1, 2, raw = TRUE),data = corData)
summary(md)
cc = combinat::combn2(1:ncol(corData))
cn =colnames(corData)
cor.df = data.frame()
for(i in 1:nrow(cc)){
  x = cc[i,1]
  y = cc[i,2]
  ph = cor.test(corData[,x],corData[,y])
  ci = t(as.vector(ph$conf.int))
  ph = data.frame(X = cn[x],Y = cn[y],Cor = ph$estimate,ci,p = ph$p.value)
  cor.df = rbind(cor.df,ph)
}
cor.df$m1 = dplyr::if_else(str_detect(cor.df$X,"Comp.1"),1,0)
cor.df$m2 = dplyr::if_else(str_detect(cor.df$Y,"Comp.1"),1,0)
cor.df$keep = dplyr::if_else(cor.df$m1==1 | cor.df$m2==1,1,0)
cor.df$sum = cor.df$m1+cor.df$m2
cor.df = cor.df[cor.df$keep==1,]
cor.df = cor.df[cor.df$sum!=2,]
cor.df = cor.df[-1,]
cor.df$p.adj = p.adjust(cor.df$p,method = "fdr")
write_csv(cor.df,path = "Output/byTreatmentnasalMbiome_MediatorCorrelationTable.csv")

corData1 = data.frame(y_train, yy,intFeats2)
colnames(corData1) = str_replace(colnames(corData1),"___"," / ")
## Figure <6B> ####
md = lm(Neutrophil.Elastase___IL.8~Comp.1,data = corData);summary(md)
mymat <- as.data.frame(predict(md, newdata = corData, interval = "c"))
mymat <- mutate(mymat, xvals = corData$Comp.1)
tiff(filename = "Figures/byTreatment_comp1Corr_wProteinRatio.tiff",width = 2.75,height = 2,units = "in",res = 300,bg = "transparent")
ggplot(corData1,aes(Comp.1,`Neutrophil.Elastase / IL.8`))+
  stat_smooth(method = "lm", se = F, size = 1)+
  geom_point(aes(col = y_train,shape = y_train),alpha = .95,size= 2)+
  theme_bw()+
  xlab(paste0("PLSDA-Comp.1 (",round(plsVariation[1]*100,2),"% of variation)"))+
  geom_line(data = mymat, mapping = aes(x = xvals, y = lwr), linetype = "dashed") +
  geom_line(data = mymat, mapping = aes(x = xvals, y = upr), linetype = "dashed") +
  #geom_text(aes(-5,0, label=(paste(expression(" R "^2*" = 0.3445,p =  ,FDR = 0.0261")))),parse = TRUE,size = 3)+
  geom_text(aes(-4,7, label=(paste("PCC = 0.3445,p = 0.0065,FDR = 0.0261"))),size = 2.5)+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.y = element_text(hjust = 0,size = 8),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8,face = "bold"),
        axis.title = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 8),
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        plot.caption = element_text(size = 8))
dev.off()

md = lm(Total.IgA___IL.8~Comp.1,data = corData);summary(md)
mymat <- as.data.frame(predict(md, newdata = corData, interval = "c"))
mymat <- mutate(mymat, xvals = corData$Comp.1)
tiff(filename = "Figures/byTreatment_comp1Corr_wProteinRatio2.tiff",width = 2.75,height = 2,units = "in",res = 300,bg = "transparent")
ggplot(corData1,aes(Comp.1,`Total.IgA / IL.8`))+
  stat_smooth(method = "lm", se = F, size = 1)+
  geom_point(aes(col = y_train,shape = y_train),alpha = .95,size= 2)+
  theme_bw()+
  xlab(paste0("PLSDA-Comp.1 (",round(plsVariation[1]*100,2),"% of variation)"))+
  geom_line(data = mymat, mapping = aes(x = xvals, y = lwr), linetype = "dashed") +
  geom_line(data = mymat, mapping = aes(x = xvals, y = upr), linetype = "dashed") +
  #geom_text(aes(-5,0, label=(paste(expression(" R "^2*" = 0.3445,p =  ,FDR = 0.0261")))),parse = TRUE,size = 3)+
  geom_text(aes(-4,1.5, label=(paste("PCC = 0.2904,p = 0.0231,FDR = 0.0463"))),size = 2.5)+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.y = element_text(hjust = 0,size = 8),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8,face = "bold"),
        axis.title = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 8),
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        plot.caption = element_text(size = 8))
dev.off()

tiff(filename = "Figures/byTreatment_comp1.tiff",width = 3.25,height = 2,units = "in",res = 300,bg = "transparent")
ggplot(coords.df,aes(Group,Comp.1,fill = Group))+
  geom_hline(yintercept = 0,lty = "dashed")+
  geom_violin(alpha = .5,draw_quantiles = .5,trim = F,size = .5) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.75,col = "black")+
  coord_flip()+
  theme_classic()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  theme(axis.text.y = element_text(hjust = 1,size = 8,face = "bold"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8,face = "bold"),
        axis.title.x = element_text(size = 8,face = "bold"),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        legend.position = "none",legend.key.size = unit(.15,units = "in"),
        plot.caption = element_text(size = 8))
dev.off()


ggplot(corData1,aes(Comp.1,Neutrophil.Elastase___IL.8))+
  stat_smooth(method = "lm", se = F, size = 1)+
  geom_point()+
  theme_classic()
ggplot(corData1,aes(Comp.1, Lactoferrin___IL.8))+
  #stat_ellipse()+
  geom_smooth(method = "lm",se  =F )+
  #stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1)+
  geom_point()+
  theme_classic()

ggplot(corData1,aes(Neutrophil.Elastase___IL.8,fill = y_train))+
  geom_density()+
  theme_classic()



### cor prot and mbiome
corData =  data.frame( intFeats)
cc = combinat::combn2(1:ncol(corData))
cn =colnames(corData)
cor.df = data.frame()
for(i in 1:nrow(cc)){
  x = cc[i,1]
  y = cc[i,2]
  ph = cor.test(corData[,x],corData[,y])
  ci = t(as.vector(ph$conf.int))
  ph = data.frame(X = cn[x],Y = cn[y],Cor = ph$estimate,ci,p = ph$p.value)
  cor.df = rbind(cor.df,ph)
}
metadata.all =   data.frame( metadata)
metaData.all = cbind(metadata.all,corData,scores)
plot(metaData.all$Comp.1,metaData.all$Age)
# Check Correlation with metadata
# discete
summary(lm(Comp.1~factor(Visit_Season),data = metaData.all))
summary(lm(Comp.1~factor(NoseSide),data = metaData.all))
summary(lm(Comp.1~factor(Sex),data = metaData.all))
summary(lm(Comp.1~factor(Race),data = metaData.all))
## Continous
summary(lm(Comp.1~SerumCotinineELISA,data = metaData.all))
summary(lm(Comp.1~BMI,data = metaData.all))
summary(lm(Comp.1~BMI,data = metaData.all))

## SerumConintine
bool = metaData.all$SerumCotinineELISA!=0
cor.test(metaData.all$Comp.1[bool],metaData.all$SerumCotinineELISA[bool])
plot(metaData.all$Comp.1[bool],metaData.all$SerumCotinineELISA[bool])
## SerumCon to Neutrophil
cor.test(metaData.all$Neutrophil.Elastase___IL.8[bool],metaData.all$SerumCotinineELISA[bool])
plot(metaData.all$Neutrophil.Elastase___IL.8[bool],metaData.all$SerumCotinineELISA[bool])




## By Group
corMethod = "pearson"
betweenTreatmentCorrelations = data.frame()
impComparistions = data.frame()

## Smoker
levels_ = unique(y_train)
corData =  data.frame( intFeats[y_train==levels_[1],])
colnames(corData) = str_replace_all(colnames(corData),"___"," / ")
metadata1 =   data.frame( metadata[y_train==levels_[1],])
metaData.df = cbind(metadata1,corData,scores[y_train==levels_[1],])
summary(lm(Comp.1~factor(Visit_Season),data = metaData.df))
plot(metaData.df$Comp.1,metaData.df$Age)
allCor = cor(corData,method = corMethod)
colnames(corData) = str_replace_all(colnames(corData),"___"," / ")
cc = combinat::combn2(1:ncol(corData))
cn =colnames(corData)
cor.df = data.frame()
for(i in 1:nrow(cc)){
  x = cc[i,1]
  y = cc[i,2]
  ph = cor.test(corData[,x],corData[,y],method = corMethod)
  ci = t(as.vector(ph$conf.int))
  ph = data.frame(X = cn[x],Y = cn[y],Cor = ph$estimate,ci,p = ph$p.value)
  cor.df = rbind(cor.df,ph)
}
cor.df$p.adj = p.adjust(cor.df$p,method = "BH")
betweenTreatmentCorrelations = rbind(betweenTreatmentCorrelations,data.frame(cor.df,Group = levels_[1]))
c1 = cor.df
fdr = 0.1
smokerCor = cor.df[cor.df$p<=fdr,]
impComparistions = rbind(impComparistions,data.frame(smokerCor,Group = levels_[1]))
mn = min(abs(smokerCor$Cor))
keep = unique(c(smokerCor$X,smokerCor$Y))
cn = colnames(allCor)
rn = rownames(allCor)
allCor1 = allCor[rn %in% keep,cn %in% keep]
allCor1[abs(allCor1)<mn]=0
cmat = abs(allCor1)
g1 = graph_from_adjacency_matrix(as.matrix(cmat),weighted = T,mode = "undirected")
g1 = igraph::simplify(g1, remove.loops = TRUE,
                      edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g1)$weight
plot(g1)
toGephi(g1,"smokerCorMat")

corData = subset(corData,select = sort(colnames(corData)))
p.mat <- cor.mtest(corData)
lt = lower.tri(p.mat,diag = F)
p.mat[lt] = p.adjust(p.mat[lt],method = "fdr")
lt = upper.tri(p.mat,diag = F)
p.mat[lt] = p.adjust(p.mat[lt],method = "fdr")
M = cor(corData)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
tiff(filename = "Figures/byTreatent_smokerCor.tiff",width = 6,height = 6,units = "in",res = 300,bg = "transparent")
corrplot::corrplot(M, method="circle", col=col(200),
         type="upper", mar = c(0,0,0,0),
         addrect = 10,outline = T,rect.col = "black",addshade = "all",
         #order="hclust",
         tl.cex = .65,
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.1, insig = "label_sig", pch.cex = 1.25,pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE
)
dev.off()

### nonSmoker
levels_ = unique(y_train)
corData =  data.frame( intFeats[y_train==levels_[2],])
colnames(corData) = str_replace_all(colnames(corData),"___"," / ")
allCor = cor(corData,method = corMethod)
cc = combinat::combn2(1:ncol(corData))
cn =colnames(corData)
cor.df = data.frame()
for(i in 1:nrow(cc)){
  x = cc[i,1]
  y = cc[i,2]
  ph = cor.test(corData[,x],corData[,y],method = corMethod)
  ci = t(as.vector(ph$conf.int))
  ph = data.frame(X = cn[x],Y = cn[y],Cor = ph$estimate,ci,p = ph$p.value)
  cor.df = rbind(cor.df,ph)
}
cor.df$p.adj = p.adjust(cor.df$p,method = "BH")
betweenTreatmentCorrelations = rbind(betweenTreatmentCorrelations,data.frame(cor.df,Group = levels_[2]))
fdr = 0.1
nonsmokerCor = cor.df[cor.df$p<=fdr,]
impComparistions = rbind(impComparistions,data.frame(nonsmokerCor,Group = levels_[2]))
mn = min(abs(nonsmokerCor$Cor))
keep = unique(c(nonsmokerCor$X,nonsmokerCor$Y))
cn = colnames(allCor)
rn = rownames(allCor)
allCor1 = allCor[rn %in% keep,cn %in% keep]
allCor1[abs(allCor1)<mn]=0
cmat = abs(allCor1)
g1 = graph_from_adjacency_matrix(as.matrix(cmat),weighted = T,mode = "undirected")
g1 = igraph::simplify(g1, remove.loops = TRUE,
                      edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g1)$weight
plot(g1)
toGephi(g1,"nonsmokerCorMat")



corData = subset(corData,select = sort(colnames(corData)))
p.mat <- cor.mtest(corData)
lt = lower.tri(p.mat,diag = F)
p.mat[lt] = p.adjust(p.mat[lt],method = "fdr")
lt = upper.tri(p.mat,diag = F)
p.mat[lt] = p.adjust(p.mat[lt],method = "fdr")
M = cor(corData)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
tiff(filename = "Figures/byTreatent_NonsmokerCor.tiff",width = 6,height = 6,units = "in",res = 300,bg = "transparent")
corrplot::corrplot(M, method="circle", col=col(200),
                   type="upper", mar = c(0,0,0,0),
                   addrect = 10,outline = T,rect.col = "black",addshade = "all",
                   #order="hclust",
                   tl.cex = .65,
                   #addCoef.col = "black", # Add coefficient of correlation
                   tl.col="black", tl.srt=45, #Text label color and rotation
                   # Combine with significance
                   p.mat = p.mat, sig.level = 0.1, insig = "label_sig", pch.cex = 1.25,pch.col = "black",
                   # hide correlation coefficient on the principal diagonal
                   diag=FALSE
)
dev.off()




### ecig
levels_ = unique(y_train)
corData =  data.frame( intFeats[y_train==levels_[3],])
colnames(corData) = str_replace_all(colnames(corData),"___"," / ")
allCor = cor(corData,method = corMethod)
cc = combinat::combn2(1:ncol(corData))
cn =colnames(corData)
cor.df = data.frame()
for(i in 1:nrow(cc)){
  x = cc[i,1]
  y = cc[i,2]
  ph = cor.test(corData[,x],corData[,y],method = corMethod)
  ci = t(as.vector(ph$conf.int))
  ph = data.frame(X = cn[x],Y = cn[y],Cor = ph$estimate,ci,p = ph$p.value)
  cor.df = rbind(cor.df,ph)
}
cor.df$p.adj = p.adjust(cor.df$p,method = "BH")
betweenTreatmentCorrelations = rbind(betweenTreatmentCorrelations,data.frame(cor.df,Group = levels_[3]))
fdr = 0.1
nonsmokerCor = cor.df[cor.df$p<=fdr,]
impComparistions = rbind(impComparistions,data.frame(nonsmokerCor,Group = levels_[3]))
mn = min(abs(nonsmokerCor$Cor))
keep = unique(c(nonsmokerCor$X,nonsmokerCor$Y))
cn = colnames(allCor)
rn = rownames(allCor)
allCor1 = allCor[rn %in% keep,cn %in% keep]
allCor1[abs(allCor1)<mn]=0
cmat = abs(allCor1)
g1 = graph_from_adjacency_matrix(as.matrix(cmat),weighted = T,mode = "undirected")
g1 = igraph::simplify(g1, remove.loops = TRUE,
                      edge.attr.comb = igraph_opt("edge.attr.comb"))
E(g1)$weight
toGephi(g1,"ecigCorMat")

corData = subset(corData,select = sort(colnames(corData)))
p.mat <- cor.mtest(corData)
lt = lower.tri(p.mat,diag = F)
p.mat[lt] = p.adjust(p.mat[lt],method = "fdr")
lt = upper.tri(p.mat,diag = F)
p.mat[lt] = p.adjust(p.mat[lt],method = "fdr")
M = cor(corData)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
tiff(filename = "Figures/byTreatent_EcigCor.tiff",width = 6,height = 6,units = "in",res = 300,bg = "transparent")
corrplot::corrplot(M, method="circle", col=col(200),
                   type="upper", mar = c(0,0,0,0),
                   addrect = 10,outline = T,rect.col = "black",addshade = "all",
                   #order="hclust",
                   tl.cex = .65,
                   #addCoef.col = "black", # Add coefficient of correlation
                   tl.col="black", tl.srt=45, #Text label color and rotation
                   # Combine with significance
                   p.mat = p.mat, sig.level = 0.1, insig = "label_sig", pch.cex = 1.25,pch.col = "black",
                   # hide correlation coefficient on the principal diagonal
                   diag=FALSE
)
dev.off()



### between Treatmen Correlations
impComparistions$Comp = paste0("(",impComparistions$X,") - (",impComparistions$Y,")")
compar = unique(impComparistions$Comp)
betweenTreatmentCorrelations1 = spread(betweenTreatmentCorrelations[,c(1:3,8)],"Group","Cor")
betweenTreatmentCorrelations1$Comp = paste0("(",betweenTreatmentCorrelations1$X,") - (",betweenTreatmentCorrelations1$Y,")")
keepComp = betweenTreatmentCorrelations1[betweenTreatmentCorrelations1$Comp%in%compar,]

ref =data.frame()
g = unique(betweenTreatmentCorrelations$Group)
for(gr in g){
  ph = betweenTreatmentCorrelations %>%
    filter(Group==gr)
  keepComp.ecig = left_join(keepComp,ph)
  phh = keepComp.ecig
  #phh$signf = if_else(phh$p.adj<0.1,"*","")
  phh$Ratio = str_replace(string = phh$X,pattern = " / ",replacement = "___")
  phh = left_join(phh,loading_)
  keepComp.ecig$X = paste0("(",phh$num,"/",phh$denom,")")
  phh = keepComp.ecig
  phh$Ratio = str_replace(string = phh$Y,pattern = " / ",replacement = "___")
  phh = left_join(phh,loading_)
  keepComp.ecig$Y = paste0("(",phh$num,"/",phh$denom,")")
  keepComp.ecig$p.adj = p.adjust(keepComp.ecig$p,method = "BH")
  keepComp.ecig$signf = if_else(keepComp.ecig$p.adj<0.1,"*","")
  keepComp.ecig = keepComp.ecig %>%
    select(X,Y,signf,Group)
  ref = rbind(ref,keepComp.ecig)
}
ref$comp = paste0(ref$X," - ",ref$Y)
ref = spread(ref,"Group","signf")

ref.cor =data.frame()
g = unique(betweenTreatmentCorrelations$Group)
for(gr in g){
  ph = betweenTreatmentCorrelations %>%
    filter(Group==gr)
  keepComp.ecig = left_join(keepComp,ph)
  phh = keepComp.ecig
  #phh$signf = if_else(phh$p.adj<0.1,"*","")
  phh$Ratio = str_replace(string = phh$X,pattern = " / ",replacement = "___")
  phh = left_join(phh,loading_)
  keepComp.ecig$X = paste0("(",phh$num," / ",phh$denom,")")
  phh = keepComp.ecig
  phh$Ratio = str_replace(string = phh$Y,pattern = " / ",replacement = "___")
  phh = left_join(phh,loading_)
  keepComp.ecig$Y = paste0("(",phh$num," / ",phh$denom,")")
  keepComp.ecig$signf = if_else(keepComp.ecig$p.adj<0.1,"*","")
  keepComp.ecig = keepComp.ecig %>%
    select(X,Y,Cor,Group)
  ref.cor = rbind(ref.cor,keepComp.ecig)
}
ref.cor$comp = paste0(ref.cor$X," - ",ref.cor$Y)
ref.cor = spread(ref.cor,"Group","Cor")



mat = as.matrix(ref.cor[,levels_])
rownames(mat) = ref.cor$comp
mat.signf = ref[,levels_]

tiff(filename = "Figures/byTreatment_corrHeatmap.tiff",width = 5,height = 7,units = "in",res = 300,bg = "transparent")
gplots::heatmap.2( mat,
                   col = gplots::redblue(n = 100), #viridis::viridis(n = 1000,direction = 1,option = "C"),#gplots::redblue(n = 1000) ,#viridis::viridis(n = 1000,option = "D"),#,
                   Rowv = TRUE,Colv = T,cellnote = mat.signf,notecex = 1.2,notecol = "black",
                   margins = c(2,20),
                   sepwidth=c(0.01,0.01),lwid = c(.1,2),lhei = c(.02,1),
                   #sepcolor="white",colsep=1:ncol(mat),rowsep=1:nrow(mat),
                   hclustfun = function(x) hclust(x, method = "average"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=F, symkey=TRUE, srtCol = 0,
                   scale = "none",adjCol = c(.5,-.5),adjRow = c(0.05,0),
                   density.info="density", trace="none",
                   #ColSideColors = sc,
                   #main = "Presensitization Gut Mircobiome",
                   #key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",reference$i,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = .65,cexCol = .75,
                   #adjCol = c(1,1)
)
dev.off()

tiff(filename = "Figures/byTreatment_corrHeatmapLegend.tiff",width = 5,height = 7,units = "in",res = 300,bg = "transparent")
gplots::heatmap.2( mat,
                   col = gplots::redblue(n = 100), #viridis::viridis(n = 1000,direction = 1,option = "C"),#gplots::redblue(n = 1000) ,#viridis::viridis(n = 1000,option = "D"),#,
                   Rowv = TRUE,Colv = T,cellnote = mat.signf,notecex = 1.2,notecol = "black",
                   # margins = c(2,20),
                   #sepwidth=c(0.01,0.01),lwid = c(.1,2),
                   lhei = c(.25,1),
                   # #sepcolor="white",colsep=1:ncol(mat),rowsep=1:nrow(mat),
                   hclustfun = function(x) hclust(x, method = "average"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=T, symkey=TRUE, srtCol = 0,keysize = 2,key.title = "",key.xlab = "PCC",key.par = list(cex=.75),
                   scale = "none",adjCol = c(.5,-.5),adjRow = c(0.05,0),
                   density.info="density", trace="none",
                   #ColSideColors = sc,
                   #main = "Presensitization Gut Mircobiome",
                   #key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",reference$i,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = .65,cexCol = .75,
                   #adjCol = c(1,1)
)
dev.off()




### ci
ph.all  = betweenTreatmentCorrelations
ph.all$comp = paste0("(",ph.all$X,") - (",ph.all$Y,")")
ph = ph.all[ph.all$p.adj<=0.1,]


keyComparisons = unique(ph$comp)
kc = ph.all[ph.all$comp %in% keyComparisons,]

kc1 = spread(kc[,c(1:3,8:9)],"Group","Cor")
mat = as.matrix(kc1[,levels_])
rownames(mat) = kc1$comp



#tiff(filename = "Figures/byTreatment_corrHeatmap.tiff",width = 5,height = 7,units = "in",res = 300,bg = "transparent")
gplots::heatmap.2( mat,
                   col = gplots::redblue(n = 100), #viridis::viridis(n = 1000,direction = 1,option = "C"),#gplots::redblue(n = 1000) ,#viridis::viridis(n = 1000,option = "D"),#,
                   Rowv = TRUE,Colv = T,
                   cellnote = round(mat,2),notecex = .80,notecol = "black",
                   margins = c(5,15),
                   sepwidth=c(0.01,0.01),lwid = c(.1,2),lhei = c(.1,1),
                   sepcolor="white",colsep=1:ncol(mat),rowsep=1:nrow(mat),
                   hclustfun = function(x) hclust(x, method = "average"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=F, symkey=TRUE, srtCol = 0,
                   scale = "none",
                   density.info="density", trace="none",
                   #ColSideColors = sc,
                   #main = "Presensitization Gut Mircobiome",
                   #key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",reference$i,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = .65,cexCol = .75,
                   #adjCol = c(1,1)
)
dev.off()

##-------------------------*
## Network by subjectGroup ####
##------------------------------*
h = 3
w = 3
ewScale = 2.5
nodeSize = 25

##========*
## ecig
##=======*
net = kc[kc$Group=="Ecig",]
colNodes = c()
g = 1
fdr = 0.1
for(i in 1:nrow(net)){

  if(net$p.adj[i]<=fdr){
   colNodes[g] = as.vector(net[i,1])
   g = g+1
   colNodes[g] = as.vector(net[i,2])
   g= g+1
  }

}
g1 = graph_from_edgelist(as.matrix(net[,c(1:2)]),directed = F)
g1 = igraph::simplify(g1, remove.loops = TRUE,
                      edge.attr.comb = igraph_opt("edge.attr.comb"))
el = data.frame(X = get.edgelist(g1)[,1],Y =  get.edgelist(g1)[,2])
el  =left_join(el,net)
edgeCol = if_else(el$p.adj<=fdr,if_else(el$Cor<0,"red","blue") ,"grey")
eweights = abs(el$Cor)
vnames = names(V(g1))
vcol = if_else(vnames %in% colNodes,"gold","grey")
vLabelCol = if_else(vnames %in% colNodes,"royalblue","black")

tiff(filename = "Figures/byTreatment_ecigNet.tiff",width = w,height = h,units = "in",res = 300)
par(mar=c(0,2.5,0,2))
set.seed(3)
plot.igraph(g1,layout = layout_with_kk,vertex.color = vcol,vertex.label.color = vLabelCol,vertex.size = nodeSize,
            edge.color = edgeCol,edge.width = ewScale * eweights,
            edge.arrow.size = .5, edge.arrow.width = 1.25,edge.curved = .2,vertex.label.cex = .55,
)
dev.off()

##========*
## smoker
##=======*
net = kc[kc$Group=="Smoker",]
colNodes = c()
g = 1
fdr = 0.1
for(i in 1:nrow(net)){

  if(net$p.adj[i]<=fdr){
    colNodes[g] = as.vector(net[i,1])
    g = g+1
    colNodes[g] = as.vector(net[i,2])
    g= g+1
  }

}
g1 = graph_from_edgelist(as.matrix(net[,c(1:2)]),directed = F)
g1 = igraph::simplify(g1, remove.loops = TRUE,
                      edge.attr.comb = igraph_opt("edge.attr.comb"))
el = data.frame(X = get.edgelist(g1)[,1],Y =  get.edgelist(g1)[,2])
el  =left_join(el,net)
edgeCol = if_else(el$p.adj<=fdr,if_else(el$Cor<0,"red","blue") ,"grey")
eweights = abs(el$Cor)
vnames = names(V(g1))
vcol = if_else(vnames %in% colNodes,"gold","grey")
vLabelCol = if_else(vnames %in% colNodes,"royalblue","black")

tiff(filename = "Figures/byTreatment_smokerNet.tiff",width = w,height = h,units = "in",res = 300)
par(mar=c(0,2.5,0,2))
set.seed(3)
plot.igraph(g1,layout = layout_with_kk,vertex.color = vcol,vertex.label.color = vLabelCol,vertex.size = nodeSize,
            edge.color = edgeCol,edge.width = ewScale * eweights,
            edge.arrow.size = .5, edge.arrow.width = 1.25,edge.curved = .2,vertex.label.cex = .55,
)
dev.off()


##========*
## Nonsmoker
##=======*
net = kc[kc$Group=="Nonsmoker",]
colNodes = c()
g = 1
fdr = 0.1
for(i in 1:nrow(net)){

  if(net$p.adj[i]<=fdr){
    colNodes[g] = as.vector(net[i,1])
    g = g+1
    colNodes[g] = as.vector(net[i,2])
    g= g+1
  }

}
g1 = graph_from_edgelist(as.matrix(net[,c(1:2)]),directed = F)
g1 = igraph::simplify(g1, remove.loops = TRUE,
                      edge.attr.comb = igraph_opt("edge.attr.comb"))
el = data.frame(X = get.edgelist(g1)[,1],Y =  get.edgelist(g1)[,2])
el  =left_join(el,net)
edgeCol = if_else(el$p.adj<=fdr,if_else(el$Cor<0,"red","blue") ,"grey")
eweights = abs(el$Cor)

vnames = names(V(g1))
vcol = if_else(vnames %in% colNodes,"gold","grey")
vLabelCol = if_else(vnames %in% colNodes,"royalblue","black")


tiff(filename = "Figures/byTreatment_NonsmokerNet.tiff",width = w,height = h,units = "in",res = 300)
par(mar=c(0,2.5,0,2))
set.seed(3)
plot.igraph(g1,layout = layout_with_kk,vertex.color = vcol,vertex.label.color = vLabelCol,vertex.size = nodeSize,
            edge.color = edgeCol,edge.width = ewScale * eweights,
            edge.arrow.size = .5, edge.arrow.width = 1.25,edge.curved = .2,vertex.label.cex = .55,
)
dev.off()






##---------------------------------------*
### Score plot ####
##----------------------------------------*
### Proteins
glm.mdl1=train(x =  data.frame(intFeats2),
               y = y_train,
               method = "pls",tuneGrid = expand.grid(ncomp = 2),
               metric = "AUC",
               trControl = tc1
)
## Get Scores
plsdaVariation = glm.mdl1$finalModel$Xvar  / glm.mdl1$finalModel$Xtotvar
loading_ =glm.mdl1$finalModel$projection
loading_ = data.frame(Ratio = rownames(loading_),Comp1 = loading_[,1],Comp2 = loading_[,2],row.names = 1:nrow(loading_))
scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)
### Scores Plot
scalingFactor = 3
tiff(filename = "Figures/byTreatment_scorePlot_Proteins.tiff",width = 5,height = 4,units = "in",res = 300)
ggplot(coords.df,aes(Comp.1,Comp.2))+
  geom_segment(data = loading_,aes(x = 0,y=0,xend = scalingFactor*Comp1, yend = scalingFactor*Comp2),size = 1,alpha = 1,
               lineend = "round", linejoin = "bevel",
               arrow = arrow(length = unit(0.1, "inches")))+
  stat_ellipse(geom = "polygon",aes(fill = Group),alpha = .2,type = "t",level = .9)+
  geom_point(aes(fill = Group),alpha = 1,size = 2,shape = 21)+
  geom_text(data = loading_,aes(label=Ratio,x = scalingFactor*Comp1,y = scalingFactor*Comp2), col = "red",
             size=2,
             hjust = .5,
             nudge_x = .1,nudge_y = .25,
             check_overlap = T)+
  xlab(paste0("PLSDA-Comp.1 (",round(plsdaVariation[1]*100,2),"% of variation)"))+
  ylab(paste0("PLSDA-Comp.2 (",round(plsdaVariation[2]*100,2),"% of variation)"))+
  theme_bw()+
  ggtitle("Nasal Proteins Only")+
   geom_vline(xintercept = 0,col = "grey40")+
  geom_hline(yintercept = 0,col = "grey40")+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "top",legend.title = element_text(face = "bold",size = 8),
        legend.box = "vertical",legend.margin = margin(0,0,0,-.1),
        axis.text = element_text(size = 8),axis.title = element_text(size = 8,face = "bold"),
        legend.spacing.y = unit(-.1,"cm"),
        legend.text = element_text(size = 8))
dev.off()
##--------*
### Mbiome
##----------*
glm.mdl1=train(x =  data.frame(intFeats1),
               y = y_train,
               method = "pls",tuneGrid = expand.grid(ncomp = 2),
               metric = "AUC",
               trControl = tc1
)

comp1Var = glm.mdl1$finalModel$Xvar[1]/glm.mdl1$finalModel$Xtotvar
comp2Var = glm.mdl1$finalModel$Xvar[2]/glm.mdl1$finalModel$Xtotvar

## Get Loading and Plot
loading_ = glm.mdl1$finalModel$projection
loading_ = data.frame(Ratio = rownames(loading_),Comp1 = loading_[,1],Comp2 = loading_[,2],row.names = 1:nrow(loading_))
loading_ = separate(loading_,col =1,into = c("numID","denomID"),sep = "___",remove = F )
cn.num = colNames;colnames(cn.num)[3:4] = c("numID","num")
cn.denom = colNames;colnames(cn.denom)[3:4] = c("denomID","denom")
loading_ = left_join(loading_,cn.denom[,3:4],by = "denomID")
loading_ = left_join(loading_,cn.num[3:4],by = "numID")
loading_$ratio_name = paste0("frac(",loading_$num, ",", loading_$denom,")")
loading_ = loading_ %>%
  arrange(desc(Comp1))
loading_$Ratio  = factor(x = loading_$Ratio,levels = loading_$Ratio)


tiff(filename = "Figures/byTreatment_loadingsC1_Mbiome.tiff",width = 4,height = 5,units = "in",res = 300)
ggplot(loading_,aes(Ratio,Comp1))+
  geom_col()+
  coord_flip()+
  ylab("Comp.1 Loadings")+
  theme_bw()+
  scale_x_discrete(labels = parse(text = loading_$ratio_name))+
  theme(legend.position = "top",legend.title = element_text(face = "bold",size = 8),
        legend.box = "vertical",legend.margin = margin(0,0,0,-.1),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8,face = "bold"),
        axis.title.y = element_blank(),
        legend.spacing.y = unit(-.1,"cm"),
        legend.text = element_text(size = 8))
dev.off()


ggplot(loading_,aes(reorder(Ratio,Comp2),Comp2))+
  geom_col()+
  xlab("logratio")+
  theme(legend.position = "top",legend.title = element_text(face = "bold",size = 8),
        legend.box = "vertical",legend.margin = margin(0,0,0,-.1),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(size = 8,face = "bold"),
        legend.spacing.y = unit(-.1,"cm"),
        legend.text = element_text(size = 8))



## Get Scores
plsdaVariation = glm.mdl1$finalModel$Xvar  / glm.mdl1$finalModel$Xtotvar
scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)
scalingFactor = 13
loading_ = glm.mdl1$finalModel$projection
loading_ = data.frame(Ratio = rownames(loading_),Comp1 = loading_[,1],Comp2 = loading_[,2],row.names = 1:nrow(loading_))

##parse ratio
loading_ = separate(loading_,col =1,into = c("numID","denomID"),sep = "___",remove = F )
cn.num = colNames;colnames(cn.num)[3:4] = c("numID","num")
cn.denom = colNames;colnames(cn.denom)[3:4] = c("denomID","denom")
loading_ = left_join(loading_,cn.denom[,3:4],by = "denomID")
loading_ = left_join(loading_,cn.num[3:4],by = "numID")
loading_$ratio_name = paste0("frac(",loading_$num, ",", loading_$denom,")")
nn = left_join(vi.smoker,vi)




### Scores Plot
tiff(filename = "Figures/byTreatment_scorePlot_Mbiome.tiff",width = 5,height = 4,units = "in",res = 300)
ggplot(coords.df,aes(Comp.1,Comp.2))+
  stat_ellipse(geom = "polygon",aes(fill = Group),alpha = .2,type = "t",level = .9)+
  geom_point(aes(fill = Group),alpha = 1,size = 2,shape = 21)+
  xlab(paste0("PLSDA-Comp.1 (",round(plsdaVariation[1]*100,2),"% of variation)"))+
  ylab(paste0("PLSDA-Comp.2 (",round(plsdaVariation[2]*100,2),"% of variation)"))+
  theme_bw()+
  #ggtitle("Nasal Microbiome Only")+
  geom_vline(xintercept = 0,col = "grey40")+
  geom_hline(yintercept = 0,col = "grey40")+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "top",legend.title = element_text(face = "bold",size = 8),
        legend.box = "vertical",legend.margin = margin(.5,0,-5,-.1),
        axis.text = element_text(size = 8),axis.title = element_text(size = 8,face = "bold"),
        legend.spacing.y = unit(-.1,"cm"),
        legend.text = element_text(size = 8))
dev.off()
##--------*
### Mbiome+Proteins
##----------*
glm.mdl1=train(x =  data.frame(intFeats),
               y = y_train,
               method = "pls",tuneGrid = expand.grid(ncomp = 2),
               metric = "AUC",
               trControl = tc1
)
plsdaVariation = glm.mdl1$finalModel$Xvar  / glm.mdl1$finalModel$Xtotvar
loading_ =glm.mdl1$finalModel$projection
loading_ = data.frame(Ratio = rownames(loading_),Comp1 = loading_[,1],Comp2 = loading_[,2],row.names = 1:nrow(loading_))
## Get Scores
scores = glm.mdl1$finalModel$scores
scores = foreach(ii = 1:ncol(scores),.combine = cbind)%do%{
  data.frame(scores[,ii])
}
colnames(scores) = paste0("Comp.",1:ncol(scores))
coords.df = data.frame(Group = y_train,scores)

## Parse
## parse ratio
loading_ = separate(loading_,col =1,into = c("numID","denomID"),sep = "___",remove = F )
loading_$parse = if_else(str_detect(string = loading_$Ratio,pattern = "___V"),1,0)
## mbiome
ph = loading_ %>%
  filter(parse==1)
cn.num = colNames;colnames(cn.num)[3:4] = c("numID","num")
cn.denom = colNames;colnames(cn.denom)[3:4] = c("denomID","denom")
ph = left_join(ph,cn.denom[,3:4],by = "denomID")
ph = left_join(ph,cn.num[3:4],by = "numID")
ph$ratio_name = paste0("frac(",ph$num, ",", ph$denom,")")
## proteins
ph2 = loading_ %>%
  filter(parse!=1)
ph2$num = ph2$numID
ph2$denom = ph2$denomID
ph2$ratio_name = paste0("frac(",ph2$num, ",", ph2$denom,")")
##combine
loading_ = rbind(ph,ph2)
loading_$mag =sqrt((loading_$Comp1^2+loading_$Comp2^2))
loading_$rank = rank(-loading_$mag )
loading_$label = if_else(loading_$rank%in%1:8,loading_$ratio_name,"")
loading_$col = if_else(loading_$rank%in%1:8,"darkblue","grey")
loading1 = loading_ %>%
  filter(rank %in%1:8)


### Scores Plot
scalingFactor = 13
tiff(filename = "Figures/byTreatment_scorePlot_MbiomeProteins.tiff",width = 5,height = 3.9,units = "in",res = 300)
ggplot(coords.df,aes(Comp.1,Comp.2))+
  stat_ellipse(geom = "polygon",aes(fill = Group),alpha = .2,type = "t",level = .9)+
  geom_point(aes(col = Group,shape = Group),alpha = 1,size = 2)+
  geom_vline(xintercept = 0,col = "grey40")+
  geom_hline(yintercept = 0,col = "grey40")+
  geom_segment(data = loading_,
               aes(x = 0,y=0,xend = scalingFactor*Comp1, yend = scalingFactor*Comp2),
               size = .5,alpha = 1,col = loading_$col,
               lineend = "round", linejoin = "bevel",
               arrow = arrow(length = unit(0.1, "inches")))+
  ggrepel::geom_label_repel(data = loading1,aes(label=label ,
                                               x = scalingFactor*Comp1,
                                               y = scalingFactor*Comp2),size = 2,
                            box.padding = 1,parse = T,force = 2,max.overlaps = 20,col  ="black")+
  xlab(paste0("PLSDA-Comp.1 (",round(plsdaVariation[1]*100,2),"% of variation)"))+
  ylab(paste0("PLSDA-Comp.2 (",round(plsdaVariation[2]*100,2),"% of variation)"))+
  theme_bw(base_rect_size = 1)+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.y = element_text(hjust = 0,size = 8),
        #axis.line = element_line(size = 1),
        #line = element_line(size = 1),
        #axis.line.x.top = element_line(size = 5),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8,face = "bold"),
        axis.title = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 8),
        legend.margin = margin(0,0,0,0,unit = "cm"),legend.box.spacing = unit(0,units = "in"),
        legend.title = element_blank(),
        legend.position = "top",legend.key.size = unit(.15,units = "in"),
        plot.caption = element_text(size = 8))
dev.off()



