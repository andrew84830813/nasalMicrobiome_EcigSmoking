## written by: Andrew Hinton (andrew84@email.unc.edu)


rm(list=ls())
gc()


### Install/Load Required Packages  ####
library(dplyr)
library(readr)
library(tidyr)
library(stringr)



# Read Data  --------------------------------------------------------------

metadata = read_csv("Data/2019_09_13 ELFmetadata.csv")
test = read.delim("Data/OTU_Table_contaminantsFilteredOut.txt")
taxo = read_tsv("Data/taxonomy.tsv")
colnames(taxo)[1] = "OTUID"
test = right_join(taxo,test)
test = test %>%
  filter(!is.na(Taxon))



# Process by Taxa Level ---------------------------------------------------

processBy = 2
switch (processBy,

  ## Pre Process OTU Table by Level 7:  (1)
  {
    test = gather(data = test,key = "SampleID",value = "Counts",3:ncol(test));
    test = test[,c(2:4)];
    test = test %>%
      group_by(SampleID,Taxon) %>%
      summarise(Counts = sum(Counts));
    test = test %>%
      spread(key = "Taxon",value = "Counts",fill = 0)
    col_Counts = data.frame(Taxon = names(colSums(test[,-1])),colSums(test[,-1]));
    ptrn = "k__";
  },

  ## Pre Process OTU Table by aggregated on Family/Genus: (2)
  {
    test = separate(test,col = 2,into =c("Kingdom","Phylum","Class","Order","Family","Genus","Species") ,sep = ";");
    test$Family = str_trim(test$Family,side = "both");
    test$Genus = str_trim(test$Genus,side = "both");
    test = test %>%
      mutate(Genus = if_else(is.na(Genus),paste(Family,"|g__",sep = ""),paste(Family,"|",Genus,sep = ""))) %>%
      filter(!is.na(Family)) %>%
      filter(Family!="f__");
    test = gather(data = test,key = "SampleID",value = "Counts",9:ncol(test));
    test = test[,c(7,9,10)];
    test = test %>%
      group_by(SampleID,Genus) %>%
      summarise(Counts = sum(Counts));
    test = test %>%
      spread(key = "Genus",value = "Counts",fill = 0);
    ptrn = "f__";
  }
)




## Sample Metadata ####
allGroups = right_join(metadata,test)
rownames(allGroups) = allGroups$SampleID



### AlL Subject Groups ####
featurePosition = which(str_detect(colnames(allGroups),pattern = ptrn)==T)
df = data.frame(Status=allGroups$SubjectGroup,allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
write_csv(colnamesTable,"Output/colNamesTable_byGenus.csv")
colnames(df) = c("SubjectGroup",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
df.metadata = allGroups
## Output  - << Nasal Microbiome By Subject Group >>
write.csv(df,"Output/allSubjectGroup_nasalMbiome_byGenus.csv",row.names = T)



### Nasal Microbiome by Sex ####
df = data.frame(Sex=allGroups$Sex,allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("Sex",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
## Output  - << Nasal Microbiome By Sex>>
write.csv(df,"Output/bySex_nasalMbiome_byGenus.csv",row.names = T)



### Nasal Microbiome by Season ####
df = data.frame(Season=allGroups$Visit_Season,allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("Season",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
## Output  - << Nasal Microbiome By Sex>>
write.csv(df,"Output/bySeason_nasalMbiome_byGenus.csv",row.names = T)


### Nasal Microbiome by Nostril ####
df = data.frame(Sex=allGroups$NoseSide,allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("NoseSide",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
## Output  - << Nasal Microbiome By Sex>>
write.csv(df,"Output/byNoseSide_nasalMbiome_byGenus.csv",row.names = T)



### Nasal Microbiome by Nostril ####
df = data.frame(Sex=allGroups$SeqBatch,allGroups[,featurePosition])
colnamesTable = data.frame(cn = colnames(df[,-1]),ID = sapply(1:ncol(df[,-1]), function(x) paste("V",x,sep = "")))
colnames(df) = c("SeqBatch",as.character(colnamesTable$ID))
rownames(df) = allGroups$SampleID
## Output  - << Nasal Microbiome By Sex>>
write.csv(df,"Output/bySeqBatch_nasalMbiome_byGenus.csv",row.names = T)



