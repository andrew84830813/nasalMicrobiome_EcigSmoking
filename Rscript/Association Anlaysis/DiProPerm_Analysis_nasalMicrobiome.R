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
        },

        {
          df = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv") ) [,-1]
          df2  = data.frame( read_csv(file = "Output/allSubjectGroup_nasalMbiome_byGenus.csv") )[,-1]
          dat = data.frame(df)

        }
)


## get colnames
cnames = read_csv("Output/colNamesTable_byLevel7.csv")





## Remove Sparse Features and impute Zeroes
useALR = F
df.metadata = read.csv("Output/metadata_all.csv")
procData = selEnergyPermR::processCompData(df,minPrevalence = .85)
glm.train = procData$processedData
dat = zCompositions::cmultRepl(glm.train[,-1])
df = data.frame(Status = df2[,1],dat)

if(useALR){
  dat.plr = alr(df[,-1])
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
    dat.plr = alr(df[,-1])
  }else{
    dat.plr = selEnergyPermR::calcLogRatio(df)[,-1]
  }



  glm.train = df
  y_label = df[,1]
}






# Select lamda and alpha --------------------------------------------------

alpha_ = seq(0,1,by = .1)
lambda_ = seq(0.001,.2,length.out = 40)
scale_ = F
classes = unique(y_label)
res  = data.frame()

for(a in alpha_){
  type_family = dplyr::if_else(length(classes)>2,"multinomial","binomial")
  compTime2 = system.time({
    cv.clrlasso <- glmnet::cv.glmnet(as.matrix(dat.plr),y_label,
                                     standardize=scale_,
                                     alpha=a,lambda = lambda_,parallel = T,
                                     family=type_family,
                                     #nfolds = 5,
                                     #type.measure = "auc"
    )
  })

  ph = data.frame(Alpha = a,lambda_,dev = (cv.clrlasso$cvm))
  res = rbind(res,ph)
  message(a)
}


# get solution
i = which.min(res$dev)
lam = res$lambda_[i]
alp = res$Alpha[i]



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

nreps = 500
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


pdf(file = "Figures/revision_diProPerm_NasalMicrobiome.pdf",width = 4,height = 4)
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
names(features)=  feat.df$Ratio[-1]






# Compute Log Contrast ----------------------------------------------------

expGroup = "Ecig"
cc =  coef.df %>%
  filter(Group==expGroup)
imp.df = data.frame(Ratio = cc$Ratio,Imp = abs(cc$coef),raw = cc$coef)
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)

## stack ratio for consisitnet interpretation
keyRats2 = keyRats
keyRats2$Num = keyRats$Denom
keyRats2$Denom= keyRats$Num
keyRats2$Ratio= paste0(keyRats$Denom,"___",keyRats$Num)
keyRats2$raw = -keyRats$raw
keyRats2 = rbind(keyRats,keyRats2)
### keep negative egdes (more abundace more likely Female)
keyRats = keyRats2 %>%
  filter(raw<0)

## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b)
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$raw)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$raw))
weights.df = weights.df %>%
  group_by(Part) %>%
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,expGroup,"Other")
weights.df$col = factor(weights.df$col,levels = c(expGroup,"Other"))
logcontrast = weights.df
logcontrast$expGroup = expGroup


library(igraph)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)


vertices = data.frame(Part = V(g)$name,Label =  V(g)$name)
vertices = left_join(vertices,weights.df)
v_color = if_else(vertices$Coef>0,ggsci::pal_lancet(alpha = .6)(2)[2],ggsci::pal_lancet(alpha = .6)(2)[1])
vertices$abs_coef = abs(vertices$Coef)

el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)
col = if_else(el_act$raw>0,ggsci::pal_lancet(alpha = 1)(2)[1],ggsci::pal_lancet(alpha = 1)(2)[2])
el_act$col = col




expGroup = "Smoker"
cc =  coef.df %>%
  filter(Group==expGroup)
imp.df = data.frame(Ratio = cc$Ratio,Imp = abs(cc$coef),raw = cc$coef)
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)

## stack ratio for consisitnet interpretation
keyRats2 = keyRats
keyRats2$Num = keyRats$Denom
keyRats2$Denom= keyRats$Num
keyRats2$Ratio= paste0(keyRats$Denom,"___",keyRats$Num)
keyRats2$raw = -keyRats$raw
keyRats2 = rbind(keyRats,keyRats2)
### keep negative egdes (more abundace more likely Female)
keyRats = keyRats2 %>%
  filter(raw>0)

## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b)
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$raw)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$raw))
weights.df = weights.df %>%
  group_by(Part) %>%
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,expGroup,"Other")
weights.df$col = factor(weights.df$col,levels = c(expGroup,"Other"))
weights.df$expGroup = expGroup
logcontrast = rbind(logcontrast,weights.df)



library(igraph)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)


vertices = data.frame(Part = V(g)$name,Label =  V(g)$name)
vertices = left_join(vertices,weights.df)
v_color = if_else(vertices$Coef>0,ggsci::pal_lancet(alpha = .6)(2)[2],ggsci::pal_lancet(alpha = .6)(2)[1])
vertices$abs_coef = abs(vertices$Coef)

el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)
col = if_else(el_act$raw>0,ggsci::pal_lancet(alpha = 1)(2)[1],ggsci::pal_lancet(alpha = 1)(2)[2])
el_act$col = col

plot(g,
     vertex.color = v_color,#rep(ggsci::pal_lancet(alpha = .75)(9)[8],vcount(g)),
     layout = igraph::layout.kamada.kawai,
     vertex.frame.color = "white",vertex.frame.width = .5,
     vertex.label.cex = .8,
     vertex.label.color = "black",
     #edge.color  =col,
     edge.width = el_act$Imp ,
     vertex.size = abs(vertices$Coef) * 15+2,
     edge.curved = .2,
     edge.arrow.size = .51)





expGroup = "Nonsmoker"
cc =  coef.df %>%
  filter(Group==expGroup)
imp.df = data.frame(Ratio = cc$Ratio,Imp = abs(cc$coef),raw = cc$coef)
keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)

## stack ratio for consisitnet interpretation
keyRats2 = keyRats
keyRats2$Num = keyRats$Denom
keyRats2$Denom= keyRats$Num
keyRats2$Ratio= paste0(keyRats$Denom,"___",keyRats$Num)
keyRats2$raw = -keyRats$raw
keyRats2 = rbind(keyRats,keyRats2)
### keep negative egdes (more abundace more likely Female)
keyRats = keyRats2 %>%
  filter(raw<0)

## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b)
weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$raw)
weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$raw))
weights.df = weights.df %>%
  group_by(Part) %>%
  summarise_all(.funs = sum)
weights.df$col = if_else(weights.df$Coef>0,expGroup,"Other")
weights.df$col = factor(weights.df$col,levels = c(expGroup,"Other"))
weights.df$expGroup = expGroup
logcontrast = rbind(logcontrast,weights.df)

library(igraph)
el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
g = igraph::simplify(g, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))
el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)


vertices = data.frame(Part = V(g)$name,Label =  V(g)$name)
vertices = left_join(vertices,weights.df)
v_color = if_else(vertices$Coef>0,ggsci::pal_lancet(alpha = .6)(2)[2],ggsci::pal_lancet(alpha = .6)(2)[1])
vertices$abs_coef = abs(vertices$Coef)

el_act = data.frame(get.edgelist(g))
el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
el_act = left_join(el_act,imp.df)
col = if_else(el_act$raw>0,ggsci::pal_lancet(alpha = 1)(2)[1],ggsci::pal_lancet(alpha = 1)(2)[2])
el_act$col = col




# Visualize log contrast by exp group ----------------------------------------


# logcontrast_ = logcontrast %>%
#   filter(col != "Other")

colnames(cnames)[2] = "Part"
logcontrast_ = left_join(logcontrast_,cnames)
logcontrast = left_join(logcontrast,cnames)

logcontrast = separate(logcontrast,col = 5,into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep = "\\.\\.")
logcontrast = data.frame(logcontrast)

highest_taxa = logcontrast$Species
highest_taxa = if_else(is.na(logcontrast$Species),"s__",logcontrast$Species)

ph = logcontrast[,11]
locs = stringr::str_locate_all(ph,"_")
ids = sapply(locs, max)
lt = stringr::str_sub(ph,start = ids+1,end = ids+1)
lt = if_else(is.na(lt),"",lt)
bool = lt==""
ph.taxa = highest_taxa
ph.taxa[bool] = ""


taxa_lab = 10:6

for(c in taxa_lab){
  ph = logcontrast[,c]
  locs = stringr::str_locate_all(ph,"_")
  ids = sapply(locs, max)
  lt = stringr::str_sub(ph,start = ids+1,end = ids+1)
  lt = if_else(is.na(lt),"",lt)
  bool = lt!="" & ph.taxa==""
  ph.taxa[bool] = ph[bool]
}


ph.taxa = str_replace_all(ph.taxa,pattern = "\\.",replacement = "")
logcontrast.df = data.frame(Taxa = ph.taxa,logcontrast[,1:6])
logcontrast.df = left_join(logcontrast.df,cnames)
logcontrast.df = logcontrast.df %>%
  arrange(desc(cn))
logcontrast.df$Taxa = factor(logcontrast.df$Taxa,levels = unique(logcontrast.df$Taxa))



lc = logcontrast.df %>%
  filter(expGroup=="Nonsmoker") %>%
  arrange(desc(Coef)) %>%
  group_by(Taxa,col,expGroup) %>%
  summarise(Coef = sum(Coef))

ggplot(lc,aes(reorder(Taxa,Coef),Coef,fill = col))+
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






lc = logcontrast.df %>%
  filter(expGroup=="Smoker") %>%
  arrange(desc(Coef)) %>%
  group_by(Taxa,col,expGroup) %>%
  summarise(Coef = sum(Coef))

ggplot(lc,aes(reorder(Taxa,Coef),Coef,fill = col))+
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


lc = logcontrast.df %>%
  filter(expGroup=="Ecig") %>%
  arrange(desc(Coef)) %>%
  group_by(Taxa,col,expGroup) %>%
  summarise(Coef = sum(Coef))

ggplot(lc,aes(reorder(Taxa,Coef),Coef,fill = col))+
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



# Integrate Immun Mediators -----------------------------------------------

## Load Data
colNames = read_csv("Output/colNamesTable_byGenus.csv")
colNames = separate(colNames,col = 1,into = c('Family','Genus'),sep = "g__")
colNames$Family = str_replace(colNames$Family,"f__","")
colNames$Family = str_replace_all(colNames$Family,"\\.","")
Name = if_else(stringr::str_length(colNames$Genus )==0,colNames$Family,colNames$Genus)
colNames$names = Name
## data
#dat = data.frame( read_csv(file = "Output/byTreatment_adjSex_byGenus_80sparse.csv") )
#table(dat[,1])
sampleID = data.frame( SampleID = data.frame(read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv")[,1] )[,1])
prot = data.frame(read_csv(file = "Output/NLF_bySubjectGroup.csv"));colnames(prot) = str_replace_all(colnames(prot),"..ug.mL.",replacement = "")
colnames(prot)[1] = "SampleID"
md = df.metadata[,1:26]


## Join Data
mediators_links = inner_join(link.df,prot)
mediators_links = inner_join(mediators_links,md)
pwr = selEnergyPermR::calcLogRatio(df = data.frame(Status = "d",mediators_links[,8:14]))[,-1]




#  Visualize log contrast  -----------------------------------------------

library(ggsci)

link.df =data.frame(SampleID = mediators_links$SampleID,Group = mediators_links$SubjectGroup,mediators_links[,link_classes])
link.df$totalScore = rowSums(link.df[,-2:-1])


pdf(file = "Figures/revision_nasalMicrobiome_logContrastScatterplot.pdf",width = 3,height = 4)
ggplot(link.df,aes(Ecig,Smoker,col = Group))+
  geom_point(size = 3)+
  theme_bw()+
  scale_color_manual(values =  c(ggsci::pal_jama(alpha = 1)(7)[c(2,3,6)]))+
  xlab("Log Contrast 1 (Ecig Signature)")+
  ylab("Log Contrast 2 (Smoker Signature)")+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 8),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))
dev.off()


links.gat = gather(link.df,"group_link","value",3:5)
links.gat$multinomial = if_else(links.gat$group_link==links.gat$Group,links.gat$Group,"Other")

pdf(file = "Figures/revision_nasalMicrobiome_logContrastBoxplot.pdf",width = 4,height = 4)
ggplot(links.gat,aes(multinomial,value,fill = multinomial))+
  scale_fill_manual(values = c(ggsci::pal_jama(alpha = 1)(7)[6],ggsci::pal_jama(alpha = 1)(7)[3],"grey",ggsci::pal_jama(alpha = 1)(7)[2]))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = .1)+
  facet_grid(.~group_link,scales = "free")+
  theme_bw()+
  #scale_fill_lancet()+
  ylab("Log Contrast Score")+
  xlab("Pairwise Comparsion (Multinomial Logistics Regression)")+
  theme(legend.position = "none",strip.background = element_blank(),strip.text = element_text(face = "bold",size = 8),
        axis.title = element_text(face = "bold",size = 10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))
dev.off()




## Pairwsie linear modeling of links
res.ovr = data.frame()
medRatios = colnames(pwr)
for(c in link_classes){
  res = data.frame()
  for(r in link_classes){

    ## remove outlier
    if(r!="Ecig"){
      if(r=="Smoker"){
        outlier = which.min(link.proteins[,r])
      }else{
        outlier = which.max(link.proteins[,r])
      }
      ph = data.frame(Y = mediators_links[-outlier,c],X = link.proteins[-outlier,r]
                      #exp = factor(mediators_links$Status,levels = c("Nonsmoker","Ecig","Smoker"))
      )
    }else{
      ph = data.frame(Y = mediators_links[,c],X = link.proteins[,r]
                      #exp = factor(mediators_links$Status,levels = c("Nonsmoker","Ecig","Smoker"))
      )
    }

    mdl = lm(Y~.,data = ph)
    an = anova(mdl)
    xx = summary(mdl)
    xx = xx$coefficients
    rr = data.frame(Y = c,X =paste(r,"_proteins") ,effect = mdl$coefficients[2],std_error = xx[2,2],p_effect = xx[2,4],F = an$`F value`[1],p = an$`Pr(>F)`[1])
    res = rbind(res,rr)
  }

  res$p.adj = p.adjust(res$p,method = "BH")
  res.ovr = rbind(res.ovr,res)

}





## Pairwsie spearman correl modeling of links
res.ovr = data.frame()
medRatios = colnames(pwr)
for(c in link_classes){
  res = data.frame()
  for(r in link_classes){

    ## remove outlier
    if(r!="Ecig"){
      if(r=="Smoker"){
        outlier = which.min(link.proteins[,r])
      }else{
        outlier = which.max(link.proteins[,r])
      }
      ph = data.frame(Y = mediators_links[-outlier,c],X = link.proteins[-outlier,r]
                      #exp = factor(mediators_links$Status,levels = c("Nonsmoker","Ecig","Smoker"))
      )
    }else{
      ph = data.frame(Y = mediators_links[,c],X = link.proteins[,r]
                      #exp = factor(mediators_links$Status,levels = c("Nonsmoker","Ecig","Smoker"))
      )
    }

    mdl = cor.test(ph$Y,ph$X,method = "spearman")
    rr = data.frame(Y = c,X =paste(r,"_proteins") ,cor = mdl$estimate,
                    p = mdl$p.value)
    res = rbind(res,rr)
  }

  res$p.adj = p.adjust(res$p,method = "BH")
  res.ovr = rbind(res.ovr,res)

}



tbl = data.frame(mediators_links,
                 Smoker.proteins = link.proteins$Smoker,
                 Nonsmoker.proteins = link.proteins$Nonsmoker,
                 Ecig.proteins = link.proteins$Ecig)

outlier = which.min(link.proteins$Smoker)
ggplot(tbl[-outlier,],aes(Smoker,Smoker.proteins))+
  geom_point(aes(col = Group))+
  geom_smooth(method = "lm")+
  theme_bw()+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 12),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))
ggplot(tbl[-outlier,],aes(Nonsmoker,Smoker.proteins))+
  geom_point(aes(col = Group))+
  geom_smooth(method = "lm")+
  theme_bw()+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 12),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))


outlier = which.max(link.proteins$Nonsmoker)
ggplot(tbl[-outlier,],aes(Smoker,Nonsmoker.proteins))+
  geom_point(aes(col = Group))+
  geom_smooth(method = "lm")+
  theme_bw()+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 12),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))
ggplot(tbl[-outlier,],aes(Nonsmoker,Nonsmoker.proteins))+
  geom_point(aes(col = Group))+
  geom_smooth(method = "lm")+
  theme_bw()+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold",size = 12),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.05,units = "in"),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        legend.background = element_rect(colour = "black"))

