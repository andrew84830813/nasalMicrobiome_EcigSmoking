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



taxa_level = 2

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
useALR = T
df.metadata = read.csv("Output/metadata_all.csv")
procData = selEnergyPermR::processCompData(df,minPrevalence = .85)
glm.train = procData$processedData
dat = zCompositions::cmultRepl(glm.train[,-1])
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
  df =  data.frame(Status = df2[,1],clrInv(adjPLR))

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


# Regression Parm
type_family = dplyr::if_else(length(classes)>2,"multinomial","binomial")
scale_ = F
classes = unique(y_label)
lam = 0.001
alp = 1


if(useALR){
  colnames(dat.plr)=  c(paste0(colnames(dat.plr),"___",ref_name))
}




# Select lamda and alpha --------------------------------------------------

# alpha_ = seq(0,1,by = .1)
# lambda_ = seq(0.001,.2,length.out = 40)

# res  = data.frame()
#
# nreps = 100
#
# for(a in alpha_){
#   compTime2 = system.time({
#
#     dev = foreach::foreach(i=1:nreps,.combine = rbind)%dopar%{
#       cv.clrlasso <- glmnet::cv.glmnet(as.matrix(dat.plr),y_label,
#                                        standardize=scale_,
#                                        alpha=a,lambda = lambda_,parallel = F,
#                                        family=type_family,
#                                        #nfolds = 5,
#                                        #type.measure = "auc"
#       )
#       min(cv.clrlasso$cvm)
#     }
#
#   })
#
#   ph = data.frame(Alpha = a,lambda_,dev = mean(dev))
#   res = rbind(res,ph)
#   message(a)
# }
#
#
# # get solution
# i = which.min(res$dev)
# lam = res$lambda_[i]
# alp = res$Alpha[i]



# Compute Test Statistic -------------------------------------------------


mdl = glmnet::glmnet(as.matrix(dat.plr),y_label,
                     standardize=scale_, alpha=alp,
                     family = type_family,
                     lambda = lam,
                     #penalty.factor = c(0,rep(1,ncol(dat.plr)-1)),
                     #type.multinomial = "grouped"
)
## make predictions
response = stats::predict(mdl, newx = as.matrix(dat.plr),type = "response")
link = stats::predict(mdl, newx = as.matrix(dat.plr),type = "link")
scores = pROC::auc(pROC::multiclass.roc(y_label,response[,,1]))
loss = MLmetrics::MultiLogLoss(y_true = factor(y_label),y_pred = response[,,1])



link_classes = colnames(link)
link = data.frame(link);colnames(link)=link_classes
md.df = data.frame()
for(c in link_classes){
  bool = (y_label==c)
  m = abs(mean(link[bool,c]) - mean(link[!bool,c]))
  ph =data.frame(class = c,md = m)
  md.df = rbind(md.df,ph)
}
totDiff = sum(md.df$md)



# Compute Null DIstribution -----------------------------------------------

nreps = 1000
res.null = foreach::foreach(i=1:nreps,.combine = rbind)%dopar%{
  set.seed(i)
  yt = sample(y_label)
  mdl.null = glmnet::glmnet(as.matrix(dat.plr),yt,
                            standardize=scale_, alpha=alp,family = type_family,
                            lambda = lam,
                            #penalty.factor = c(0,rep(1,ncol(dat.plr)-1))
  )

  link.null = stats::predict(mdl.null, newx = as.matrix(dat.plr),type = "link")
  lc = colnames(link.null)
  link.null = data.frame(link.null);colnames(link.null)=lc

  md.null =data.frame()
  for(c in link_classes){
    bool = (yt==c)
    m = abs(mean(link.null[bool,c]) - mean(link.null[!bool,c]))
    ph =data.frame(class = c,md = m,rep = i)
    md.null = rbind(md.null,ph)
  }

  md.null

}




# Compute Pval and visulaize result ---------------------------------------

null = spread(res.null,key = "class","md")
null$totDiff = rowSums(null[,-1])
hist(null$totDiff )
pval = (sum(totDiff <=null$totDiff )+1)/(nreps+1)
dat = data.frame(mean_diff = null$totDiff)



#pdf(file = "Figures/revision_diProPerm_NasalMicrobiome_bySex.pdf",width = 4,height = 4)
ggplot(dat,aes(mean_diff))+
  geom_histogram(col = "black",fill = "grey")+
  geom_vline(xintercept = totDiff,col = "red",lty = "dashed")+
  labs(caption = paste0("t-stat = ",round(totDiff,5),"; Empirical pval = ",round(pval,5),"; permutations = ",nreps))+
  theme_classic()+
  ggtitle("Nasal Microbiome")+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 12,hjust = .5,face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+

  )
dev.off()








# Select feature and Visualize --------------------------------------------

features = (stats::coef(mdl))
feat.df = data.frame()
coef.df = data.frame()
for(o in 1:length(features)){
  ph = as.matrix(features[[o]])
  feat = ph[-1,]
  keep = feat[abs(feat)>0]
  coef.df = rbind(coef.df,data.frame(Group = link_classes[o],Ratio = names(keep),coef = as.numeric(keep)))
  feat.df = rbind(feat.df,data.frame(Ratio = names(keep),coef = as.numeric(keep)))
}
feat.df =feat.df %>%
  group_by(Ratio) %>%
  summarise(coef = sum(coef)) %>%
  filter(coef!=0)
features = feat.df$coef[-1]



## Ecig
expGroup = "Ecig"
cc =  coef.df %>%
  filter(Group==expGroup)
imp.df = data.frame(Ratio = cc$Ratio,Imp = abs(cc$coef),coef = cc$coef)
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
## stack ratio for consisitnet interpretation
keyRats2 = keyRats
## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b)
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$coef)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$coef))
weights.df = weights.df %>%
  group_by(Part) %>%
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,expGroup,"Other")
weights.df$col = factor(weights.df$col,levels = c(expGroup,"Other"))
logcontrast = weights.df
logcontrast$expGroup = expGroup

## Smoker
expGroup = "Smoker"
cc =  coef.df %>%
  filter(Group==expGroup)
imp.df = data.frame(Ratio = cc$Ratio,Imp = abs(cc$coef),coef = cc$coef)
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
## stack ratio for consisitnet interpretation
keyRats2 = keyRats
## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b)
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$coef)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$coef))
weights.df = weights.df %>%
  group_by(Part) %>%
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,expGroup,"Other")
weights.df$col = factor(weights.df$col,levels = c(expGroup,"Other"))
weights.df$expGroup = expGroup
logcontrast = rbind(logcontrast,weights.df)


## Smoker
expGroup = "Nonsmoker"
cc =  coef.df %>%
  filter(Group==expGroup)
imp.df = data.frame(Ratio = cc$Ratio,Imp = abs(cc$coef),coef = cc$coef)
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
## stack ratio for consisitnet interpretation
keyRats2 = keyRats
## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b)
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$coef)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$coef))
weights.df = weights.df %>%
  group_by(Part) %>%
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,expGroup,"Other")
weights.df$col = factor(weights.df$col,levels = c(expGroup,"Other"))
weights.df$expGroup = expGroup
logcontrast = rbind(logcontrast,weights.df)

## Rename taxa
colnames(cnames)[2] = "Part"
logcontrast = left_join(logcontrast,cnames)
logcontrast$cn = str_replace(logcontrast$cn,pattern = "\\.\\.",replacement = ".")
logcontrast1 = separate(logcontrast,
                       col = 5,
                       into = c("Family","Genus"),
                       sep = "\\.")
highest_taxa = logcontrast1$Genus
logcontrast1$taxa = if_else(logcontrast1$Genus=="g__",logcontrast1$Family,logcontrast1$Genus)
logcontrast1$taxa = str_replace(logcontrast1$taxa,pattern = "g__",replacement = "")
logcontrast1$taxa = str_replace(logcontrast1$taxa,pattern = "f__",replacement = "")
write_csv(logcontrast1,file = "Output/dpp_byExpGroup_importantTaxa.csv")


## visualize
lc = logcontrast1 %>%
  filter(expGroup=="Nonsmoker")
ggplot(lc,aes(reorder(taxa,Coef),Coef,fill = col))+
  geom_col(width = .7)+
  scale_fill_manual(values = c("grey",ggsci::pal_jama(alpha = 1)(5)[3]))+
  facet_grid(.~expGroup)+
  coord_flip()+
  #scale_fill_manual(values = ggsci::pal_aaas()(2)[2:1])+
  theme_bw()+
  # facet_grid(.~expGroup)+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+

  )




lc = logcontrast1 %>%
  filter(expGroup=="Smoker") %>%
  arrange(desc(Coef)) %>%
  group_by(taxa,col,expGroup) %>%
  summarise(Coef = sum(Coef))
ggplot(lc,aes(reorder(taxa,Coef),Coef,fill = col))+
  geom_col(width = .7)+
  scale_fill_manual(values = c("grey",ggsci::pal_jama(alpha = 1)(5)[2]))+
  facet_grid(.~expGroup)+
  coord_flip()+
  #scale_fill_manual(values = ggsci::pal_aaas()(2)[2:1])+
  theme_bw()+
  # facet_grid(.~expGroup)+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+

  )


lc = logcontrast1 %>%
  filter(expGroup=="Ecig")
ggplot(lc,aes(reorder(taxa,Coef),Coef,fill = col))+
  geom_col(width = .7)+
  scale_fill_manual(values = c(ggsci::pal_jama(alpha = 1)(7)[6],"grey"))+
  facet_grid(.~expGroup)+
  coord_flip()+
  #scale_fill_manual(values = ggsci::pal_aaas()(2)[2:1])+
  theme_bw()+
  # facet_grid(.~expGroup)+
  theme(legend.position = "none",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+

  )
