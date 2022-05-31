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
        },

        {
          df = data.frame( read_csv(file = "Output/bySex_nasalMbiome_byGenus.csv") ) [,-1]
          dat = data.frame(df)

        }
)


## get colnames
cnames = read_csv("Output/colNamesTable_byLevel7.csv")



## Remove Sparse Features and impute Zeroes
useALR = T
df.metadata = read.csv("Output/metadata_all.csv")
procData = selEnergyPermR::processCompData(df,minPrevalence = .85)
glm.train = procData$processedData
dat = zCompositions::cmultRepl(glm.train[,-1])
df = data.frame(Status = df[,1],dat)

if(useALR){
  ## find best ref
  ff = FINDALR(dat)
  ## with feature selection
  ref = ff$procrust.ref
  dat.plr = data.frame(alr(df[,-1],ivar = ref))
  dat.plr = data.frame(alr(df[,-1]))

}else{
  dat.plr = selEnergyPermR::calcLogRatio(df)[,-1]
}
y_label = df[,1]







# Select lamda and alpha --------------------------------------------------

alpha_ = seq(0,1,by = .1)
lambda_ = seq(0.001,.2,length.out = 40)
scale_ = F
classes = unique(y_label)
res  = data.frame()
type_family = dplyr::if_else(length(classes)>2,"multinomial","binomial")
nreps = 100

for(a in alpha_){
  compTime2 = system.time({

    dev = foreach::foreach(i=1:nreps,.combine = rbind)%dopar%{
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(dat.plr),y_label,
                                       standardize=scale_,
                                       alpha=a,lambda = lambda_,parallel = F,
                                       family=type_family,
                                       #nfolds = 5,
                                       #type.measure = "auc"
      )
      min(cv.clrlasso$cvm)
    }

  })

  ph = data.frame(Alpha = a,lambda_,dev = mean(dev))
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
scores = pROC::auc(pROC::roc(y_label,(response)))
loss = MLmetrics::MultiLogLoss(y_true = factor(y_label),y_pred = response[,1])



link_classes = unique(y_label)
md.df = data.frame()
bool = (y_label==link_classes[1])
m = abs(mean(link[bool,1]) - mean(link[!bool,1]))
totDiff = m



## projection
pr = as.numeric(link[,1])
m1 = mean(pr[df[,1]!=link_classes[1]])
m2 = mean(pr[df[,1]==link_classes[1]])
md = m1-m2



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

  md.null =data.frame()
  bool = (yt==link_classes[1])
  m = abs(mean(link.null[bool,1]) - mean(link.null[!bool,1]))
  ph =data.frame(class = c,md = m,rep = i)
  md.null = rbind(md.null,ph)
  md.null

}



# Compute Pval and visulaize result ---------------------------------------


null = res.null
hist(null$md )
pval = (sum(totDiff <=null$md )+1)/(nreps+1)
dat = data.frame(mean_diff = null$md)


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

features = as.matrix(stats::coef(mdl))
features =  features[abs(features)>0,]
features =features[-1]

