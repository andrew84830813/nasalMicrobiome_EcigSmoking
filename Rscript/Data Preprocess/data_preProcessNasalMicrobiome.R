## written by: Andrew Hinton (andrew84@email.unc.edu)


rm(list=ls())
gc()


### Load Required Packages  ####
library(tidyverse)



## initalize parallel compute cluster ####
clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)

## Load Helper Functions ####
functions_Path = "Functions/"
setwd(functions_Path)
file_names <- dir()
for(i in 1:length(file_names)){
  source(file = file_names[i])
}
setwd("..")



############################################################-
## Read Data ####
############################################################-
metadata = read_csv("Data/2019_09_13 ELFmetadata.csv")
test = read.delim("Data/OTU_Table.txt")
taxo = read_tsv("Data/taxonomy.tsv")
colnames(taxo)[1] = "OTUID"
test = right_join(taxo,test)
test = test %>%
  filter(!is.na(Taxon))
############################################################-






#By overall level 7
test = gather(data = test,key = "SampleID",value = "Counts",3:ncol(test))
test = test[,c(2:4)]
test = test %>%
  group_by(SampleID,Taxon) %>%
  summarise(Counts = sum(Counts))
test = test %>%
  spread(key = "Taxon",value = "Counts",fill = 0)

col_Counts = data.frame(Taxon = names(colSums(test[,-1])),colSums(test[,-1]))


############################################################-
## Pre Process OTU Table by aggregated on Genus ####
############################################################-
# test = separate(test,col = 2,into =c("Kingdom","Phylum","Class","Order","Family","Genus","Species") ,sep = ";")
# test$Family = str_trim(test$Family,side = "both")
# test$Genus = str_trim(test$Genus,side = "both")
# test = test %>%
#   mutate(Genus = if_else(is.na(Genus),paste(Family,"|g__",sep = ""),paste(Family,"|",Genus,sep = ""))) %>%
#   filter(!is.na(Family)) %>%
#   filter(Family!="f__")
# test = gather(data = test,key = "SampleID",value = "Counts",9:ncol(test))
# test = test[,c(7,9,10)]
# test = test %>%
#   group_by(SampleID,Genus) %>%
#   summarise(Counts = sum(Counts))
# test = test %>%
#   spread(key = "Genus",value = "Counts",fill = 0)
############################################################-



############################################################-
## Sample Metadata ####
############################################################-
allGroups = right_join(metadata,test)
rownames(allGroups) = allGroups$SampleID
## Output  -  << Metadata >>
write_csv(df.metadata[,1:25],"Output/sampleMetaData.csv")
############################################################-



############################################################-
### AlL Subject Groups ####
############################################################-
featurePosition = which(str_detect(colnames(allGroups),pattern = "k__")==T)
df = data.frame(Status=allGroups$SubjectGroup,allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
write_csv(colnamesTable,"Output/colNamesTable.csv")
colnames(df) = c("SubjectGroup",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
df.metadata = allGroups
## Output  - << Nasal Microbiome By Subject Group >>
write.csv(df,"Output/allSubjectGroup_nasalMbiome.csv",row.names = T)
############################################################-



############################################################-
### Nasal Micrbiome by Sex ####
############################################################-
df = data.frame(Sex=allGroups$Sex,allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("Sex",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
## Output  - << Nasal Microbiome By Subject Group >>
write.csv(df,"Output/bySex_nasalMbiome.csv",row.names = T)
############################################################-


############################################################-
### Nasal Micrbiome by Sex*treatment ####
############################################################-
df = data.frame(Sex_Treatment=paste0(allGroups$Sex,"_",allGroups$SubjectGroup),allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("Sex",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
## Output  - << Nasal Microbiome By Subject Group >>
write.csv(df,"Output/bysexTreat_nasalMbiome.csv",row.names = T)
############################################################-




############################################################-
### Nasal Micrbiome Ecig vs. NonSmoker ####
############################################################-
ecig_ns = right_join(metadata,test) %>%
  filter(SubjectGroup!="Smoker") %>%
  dplyr::rename(Status = SubjectGroup)
rownames(ecig_ns) = ecig_ns$SampleID
df = data.frame(Status=ecig_ns$Status,ecig_ns[,26:(ncol(ecig_ns)-1)])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("Status",as.character(colnamesTable$ID))
rownames(df) = ecig_ns$SampleID
## Output  - << Nasal Microbiome By Subject Group >>
write.csv(df,"Output/ecig-vs-nonSmoker_nasalMbiome.csv",row.names = T)
############################################################-




############################################################-
### Nasal Micrbiome Ecig vs. Smoker ####
############################################################-
ecig_ns = right_join(metadata,test) %>%
  filter(SubjectGroup!="Nonsmoker") %>%
  dplyr::rename(Status = SubjectGroup)
rownames(ecig_ns) = ecig_ns$SampleID
df = data.frame(Status=ecig_ns$Status,ecig_ns[,26:(ncol(ecig_ns)-1)])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("Status",as.character(colnamesTable$ID))
rownames(df) = ecig_ns$SampleID
## Output  - << Nasal Microbiome By Subject Group >>
write.csv(df,"Output/ecig-vs-Smoker_nasalMbiome.csv",row.names = T)
############################################################-




############################################################-
### Nasal Micrbiome NonSmoker vs. Smoker ####
############################################################-
ecig_ns = right_join(metadata,test) %>%
  filter(SubjectGroup!="Ecig") %>%
  dplyr::rename(Status = SubjectGroup)
rownames(ecig_ns) = ecig_ns$SampleID
df = data.frame(Status=ecig_ns$Status,ecig_ns[,26:(ncol(ecig_ns)-1)])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("Status",as.character(colnamesTable$ID))
rownames(df) = ecig_ns$SampleID
## Output  - << Nasal Microbiome By Subject Group >>
write.csv(df,"Output/nonSmoker-vs-Smoker_nasalMbiome.csv",row.names = T)
############################################################-



