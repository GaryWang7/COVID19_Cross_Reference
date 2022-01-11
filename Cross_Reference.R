## This is the R script to perform cross-reference
#### Peptides ####
# MART1: ELAGIGILTV
# M1:   GILGFVFTL
# LLL:  LLLDRLNQL
library(dplyr)
library(stringr)
library(tidyr)
library(patchwork)
library(immunarch)
library(circlize)
library(networkD3)
library(htmlwidgets)
library(hrbrthemes)
library(ggalluvial)
library(tidyverse)
library(ggpubr)
library(viridis)
#### Functions ####
rearrange_ImmuneCODE<-function(df,select = F){
  df<-separate(df,col = "TCR.BioIdentity", 
               sep = "\\+",
               into = c("amino_acid","v_resolved","j_resolved"),
               remove = T)
  names(df)[which(names(df)=="Amino.Acids")]<-"Epitope"
  if(select == T) df<-select(df,"amino_acid","v_resolved","j_resolved","Epitope")
  df<-subset(df, amino_acid!= "unproductive" )
  unique(df)
}
rearrange_VDJdb<-function(df,peptide = NA){
  df<-df[grep(peptide,df$Epitope,fixed = T),]
  df<-subset(df,df$Gene=="TRB")
  names(df)[which(names(df)=="CDR3")]<-"amino_acid"
  names(df)[which(names(df)=="V")]<-"v_resolved"
  names(df)[which(names(df)=="J")]<-"j_resolved"
  df<-select(df,2:13)
}
reformat_iSMART<-function(df,output_dir = NA,split = F,sample = "All"){ 
  # Split will split the df according to sample_name. Here no sample_name is needed.
  # Input file format (in the order of columns):
  # CDR3 amino acid sequence (Starting from C, ending with the first F/L in motif [FL]G.)
  # Variable gene name in Imgt format: TRBVXX-XX*XX
  # Joining gene name (optional)
  # Frequency (optional)
  # Other information (optional)
  if(split==F) {
    df_new <-
      select(df,
             amino_acid,
             v_resolved,
             j_resolved,
             productive_frequency,
             sample_name)
    if (!is.na(output_dir)) {
      write.table(
        df_new,
        file = paste0(output_dir, "/", sample, "_iSMART.tsv"),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        col.names=TRUE
      )
    }
  }else if(split==T){
    sample_list<-unique(df$sample_name)
    for(i in 1:length(sample_list)){
      df_temp<-subset(df,sample_name==sample_list[i])
      reformat_iSMART(df,output_dir = output_dir,split = F,sample = sample_list[i])
    }
  }
}
#### Read and process raw data #### 
setwd("Y:/Files/JHU/Labs/Git_Repos/COVID19_Cross_Reference")
# Read and rearrange M1 and MART1 data
M1<-read.csv("Data/VDJdb_M1_TCR.tsv",sep = "\t")
MART1<-read.csv("Data/VDJdb_MART1_TCR.tsv",sep = "\t")
M1<-unique(rearrange_VDJdb(M1,peptide = "GILGFVFTL"))
MART1<-unique(rearrange_VDJdb(MART1,peptide = "ELAGIGILTV"))
# Read ImmuneCODE data
ImmuneCODE_raw<-rbind(read.csv("Data/peptide-detail-ci.csv"),
                      read.csv("Data/peptide-detail-cii.csv"))
ImmuneCODE<-rearrange_ImmuneCODE(ImmuneCODE_raw)
LLL<-unique(ImmuneCODE[grep("LLLDRLNQL",ImmuneCODE$Epitope,fixed = T),])

#### Cross reference by exact match ####
MART1_M1 <-subset(MART1,MART1$amino_acid%in% M1$amino_acid)
MART1_LLL<-subset(MART1,MART1$amino_acid%in% LLL$amino_acid)
SARS_M1 <-subset(ImmuneCODE, amino_acid %in% M1$amino_acid)
## Output peptide list ##
Peptide<-as.data.frame(table(SARS_M1$Epitope))
#### Cross reference by Similarity Match ####
LLL$Epitope<-"LLL"
MART1$Epitope<-"MART1"
M1$Epitope<-"M1"
All<-rbind(LLL,
           select(M1,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(MART1,c("amino_acid","v_resolved","j_resolved","Epitope")))
write.table(
  All,
  file = paste0("Data/Output/iSMART", "/", "all", "_iSMART.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names=TRUE
)
