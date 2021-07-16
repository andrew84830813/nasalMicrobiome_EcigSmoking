library(compositions)
library(ggplot2)


source(file = "Functions/functions1.R")


## Read Data
df = data.frame(read_csv("Output/bySex_nasalMbiome_byGenus.csv"))
md = data.frame(read_csv(file = "Data/2019_09_13 ELFmetadata.csv"))
dat = data.frame(Status = df[,2],df[,-2:-1],row.names = df[,1])
md = left_join(data.frame(SampleID = df$X1),md)


## Process Zeroes/Sparisity
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

## by Nose Side
pc = prcomp(clr(dat[,-1]))
pc.df = data.frame(md,pc$x)
ggplot(pc.df,aes(PC1,PC2,col = NoseSide))+
  geom_point(size = 2)+
  stat_ellipse()

##  by Seq
pc.df = data.frame(md,pc$x)
ggplot(pc.df,aes(PC1,PC2,col = factor(SeqBatch)))+
  geom_point(size = 2)+
  stat_ellipse()

##  by Season
pc.df = data.frame(md,pc$x)
ggplot(pc.df,aes(PC1,PC2,col = Visit_Season))+
  geom_point(size = 2)+
  stat_ellipse()

##  by Sex
pc.df = data.frame(md,pc$x)
ggplot(pc.df,aes(PC1,PC2,col = Sex))+
  geom_point(size = 2)+
  stat_ellipse()

##  by Treatment
pc.df = data.frame(md,pc$x)
ggplot(pc.df,aes(PC1,PC2,col = SubjectGroup))+
  geom_point(size = 2)+
  stat_ellipse()

