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
NLF_protein = read_csv("Data/2020_06_11 NLF Data.csv")
nlf_df = left_join(metadata,NLF_protein)
############################################################-


############################################################-
## NLF by Subject Group ####
############################################################-
df = nlf_df %>% 
  dplyr::select(SubjectGroup,one_of(colnames(NLF_protein)[-1])) #Sex
df = data.frame(df)
rownames(df) = nlf_df$SampleID
#remove samples with missing data
remove_row = which(is.na(df$Total.IgA..pg.mL.))
df = df[-remove_row,]
#Convert units so that everything is ug/ML
df$DEFB4A.2..pg.mL. = df$DEFB4A.2..pg.mL.*1e-6
df$DEFB1..pg.mL. = df$DEFB1..pg.mL.*1e-6
df$Neutrophil.Elastase..ng.mL. = df$Neutrophil.Elastase..ng.mL.*0.001
df$Total.IgA..pg.mL. = df$Total.IgA..pg.mL.*1e-6
df$IL.8..pg.mL. = df$IL.8..pg.mL.*1e-6 
colnames(df) = str_replace(string = colnames(df),pattern = "pg",replacement = "ug")
colnames(df) = str_replace(string = colnames(df),pattern = "ng",replacement = "ug")
df = df %>% 
  dplyr::rename(Status = SubjectGroup)
## Output  -  << NLF by Subject Group >>
write.csv(df,"Output/NLF_bySubjectGroup.csv",row.names = T)
############################################################-


############################################################-
## NLF by Sex ####
############################################################-
df = nlf_df %>% 
  dplyr::select(Sex,one_of(colnames(NLF_protein)[-1])) #Sex
df = data.frame(df)
rownames(df) = nlf_df$SampleID
#remove samples with missing data
remove_row = which(is.na(df$Total.IgA..pg.mL.))
df = df[-remove_row,]
#Convert units so that everything is ug/ML
df$DEFB4A.2..pg.mL. = df$DEFB4A.2..pg.mL.*1e-6
df$DEFB1..pg.mL. = df$DEFB1..pg.mL.*1e-6
df$Neutrophil.Elastase..ng.mL. = df$Neutrophil.Elastase..ng.mL.*0.001
df$Total.IgA..pg.mL. = df$Total.IgA..pg.mL.*1e-6
df$IL.8..pg.mL. = df$IL.8..pg.mL.*1e-6 
colnames(df) = str_replace(string = colnames(df),pattern = "pg",replacement = "ug")
colnames(df) = str_replace(string = colnames(df),pattern = "ng",replacement = "ug")
df = df %>% 
  dplyr::rename(Status = Sex)
## Output  -  << NLF by Subject Group >>
write.csv(df,"Output/NLF_bySex.csv",row.names = T)
## NLF with Nasal Microbiome by Sex ####
df = data.frame(X1 = rownames(df),df)
nasalMicrobiome = read_csv(file = "Output/bySex_nasalMbiome.csv")
############################################################-




############################################################-
## NLF with Nasal Microbiome by Sex ####
############################################################-
df = nlf_df %>% 
  dplyr::select(Sex,one_of(colnames(NLF_protein)[-1])) #Sex
df = data.frame(df)
rownames(df) = nlf_df$SampleID
#remove samples with missing data
remove_row = which(is.na(df$Total.IgA..pg.mL.))
df = df[-remove_row,]
#Convert units so that everything is ug/ML
df$DEFB4A.2..pg.mL. = df$DEFB4A.2..pg.mL.*1e-6
df$DEFB1..pg.mL. = df$DEFB1..pg.mL.*1e-6
df$Neutrophil.Elastase..ng.mL. = df$Neutrophil.Elastase..ng.mL.*0.001
df$Total.IgA..pg.mL. = df$Total.IgA..pg.mL.*1e-6
df$IL.8..pg.mL. = df$IL.8..pg.mL.*1e-6 
colnames(df) = str_replace(string = colnames(df),pattern = "pg",replacement = "ug")
colnames(df) = str_replace(string = colnames(df),pattern = "ng",replacement = "ug")
df = df %>% 
  dplyr::rename(Status = Sex)
## Output  -  << NLF by Subject Group >>
write.csv(df,"Output/NLF_bySex.csv",row.names = T)
############################################################-



