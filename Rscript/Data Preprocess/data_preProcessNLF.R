## written by: Andrew Hinton (andrew84@email.unc.edu)


rm(list=ls())
gc()


### Install/Load Required Packages  ####
library(dplyr)
library(readr)
library(tidyr)
library(stringr)




############################################################-
## Read Data ####
############################################################-
metadata = read_csv("Data/2019_09_13 ELFmetadata.csv")
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
############################################################-


