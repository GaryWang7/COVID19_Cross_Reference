## This is the R script to perform cross-reference
#### Peptides ####
# MART1: ELAGIGILTV
# M1:   GILGFVFTL
# COVID TSCAN sequences:  LLLDRLNQL, KLWAQCVQL, YLQPRTFLL, ALWEIQQVV, YLFDESGEFKL, LLYDANYFL
library(Biostrings)
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
find_cross_reactive <-function(query_df,target_epitope,cluster_result){
  # Function to find target-specific TCRs that are clustered with TCRs from query_df
  shared_cluster_IDs <- subset(cluster_result, amino_acid %in% query_df$amino_acid)$clusterID
  cross_reactive_cluster<-subset(cluster_result, clusterID %in% shared_cluster_IDs 
                                 & Epitope == target_epitope)
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
# Extract TCRs for LLLDRLNQL, KLWAQCVQL, YLQPRTFLL, ALWEIQQVV, YLFDESGEFKL, LLYDANYFL
LLL<-unique(ImmuneCODE[grep("LLLDRLNQL",ImmuneCODE$Epitope,fixed = T),])
KLW<-unique(ImmuneCODE[grep("KLWAQCVQL",ImmuneCODE$Epitope,fixed = T),])
YLQ<-unique(ImmuneCODE[grep("YLQPRTFLL",ImmuneCODE$Epitope,fixed = T),])
ALW<-unique(ImmuneCODE[grep("ALWEIQQVV",ImmuneCODE$Epitope,fixed = T),])
YLF<-unique(ImmuneCODE[grep("YLFDESGEFKL",ImmuneCODE$Epitope,fixed = T),]) #0 entries
LLY<-unique(ImmuneCODE[grep("LLYDANYFL",ImmuneCODE$Epitope,fixed = T),])
#### Cross reference by exact match ####
MART1_M1 <-subset(MART1,MART1$amino_acid%in% M1$amino_acid)
MART1_LLL<-subset(MART1,MART1$amino_acid%in% LLL$amino_acid)
SARS_M1 <-subset(ImmuneCODE, amino_acid %in% M1$amino_acid)
LLL_M1 <- subset(LLL,amino_acid %in% M1$amino_acid)
## Output peptide list ##
Peptide<-as.data.frame(table(SARS_M1$Epitope))
#### Cross reference by Similarity Match Using iSMART####
#LLL$Epitope<-"LLLDRLNQL"
MART1$Epitope<-"MART1"
M1$Epitope<-"M1"
All<-unique(rbind(select(LLL,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(KLW,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(YLQ,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(ALW,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(LLY,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(YLF,c("amino_acid","v_resolved","j_resolved","Epitope")),# 0 entries
           select(M1,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(MART1,c("amino_acid","v_resolved","j_resolved","Epitope"))))
write.table(
  All,
  file = paste0("Data/Output/iSMART", "/", "M1_MART1_TSCAN", "_iSMART.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names=TRUE
)
# After running iSMART, find the clusters containing the identified exact match TCRs. YLF is not recorded in ImmunCODE and thus not included in following analyses
cluster_result <- read.table("Data/Output/iSMART/M1_MART1_TSCAN_iSMART.tsv_clusteredCDR3s_7.5.txt",
                             sep = "\t",
                             header = T)
colnames(cluster_result)<-c("amino_acid","v_resolved","j_resolved","Epitope","clusterID")
MART1_Similar_M1<-find_cross_reactive(MART1,'M1',cluster_result)
LLL_similar_M1<-find_cross_reactive(LLL,'M1',cluster_result) # 121
KLW_similar_M1<-find_cross_reactive(KLW,'M1',cluster_result) # 67
YLQ_similar_M1<-find_cross_reactive(YLQ,'M1',cluster_result) # 162
ALW_similar_M1<-find_cross_reactive(ALW,'M1',cluster_result) # 0
LLY_similar_M1<-find_cross_reactive(LLY,'M1',cluster_result) # 330
# Summarize the numbers of similar TCRs and percentage
summary_epitopes<-c('LLL','KLW','YLQ','ALW','LLY')
query_combined <- list(LLL,KLW,YLQ,ALW,LLY)
similar_combined <- list(LLL_similar_M1,
                        KLW_similar_M1,
                        YLQ_similar_M1,
                        ALW_similar_M1,
                        LLY_similar_M1)
summary_similar<- data.frame(matrix(nrow= length(summary_epitopes), ncol = 4))
colnames(summary_similar)<-c("Epitope","ImmuneCODE entries","Similar M1 entries","Similar M1/ImmuneCODE")
for(i in 1:length(summary_epitopes)){
  summary_similar$Epitope[i]<-summary_epitopes[i]
  summary_similar$`ImmuneCODE entries`[i]<-nrow(query_combined[[i]])
  summary_similar$`Similar M1 entries`[i]<-nrow(similar_combined[[i]])
  summary_similar$`Similar M1/ImmuneCODE`[i]<-nrow(similar_combined[[i]])/nrow(query_combined[[i]])

}
write.table(
  summary_similar,
  file = paste0("Data/Output/iSMART", "/", "M1_MART1_TSCAN", "_iSMART_Summary.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names=TRUE
)

#### Cross reference by Similarity Match Using GIANA####
# GIANA and iSMART share the same input format
# After running GIANA, find the clusters containing the identified exact match TCRs. YLF is not recorded in ImmunCODE and thus not included in following analyses
cluster_result <- read.table("Data/Output/GIANA/M1_MART1_TSCAN_GIANA--RotationEncodingBL62.txt",
                             sep = "\t",
                             header = F)
colnames(cluster_result)<-c("amino_acid","clusterID","v_resolved","j_resolved","Epitope")
MART1_Similar_M1<-find_cross_reactive(MART1,'M1',cluster_result)
LLL_similar_M1<-find_cross_reactive(LLL,'M1',cluster_result) 
KLW_similar_M1<-find_cross_reactive(KLW,'M1',cluster_result) 
YLQ_similar_M1<-find_cross_reactive(YLQ,'M1',cluster_result) 
ALW_similar_M1<-find_cross_reactive(ALW,'M1',cluster_result) 
LLY_similar_M1<-find_cross_reactive(LLY,'M1',cluster_result) 
# Summarize the numbers of similar TCRs and percentage
summary_epitopes<-c('LLL','KLW','YLQ','ALW','LLY')
query_combined <- list(LLL,KLW,YLQ,ALW,LLY)
similar_combined <- list(LLL_similar_M1,
                         KLW_similar_M1,
                         YLQ_similar_M1,
                         ALW_similar_M1,
                         LLY_similar_M1)
summary_similar<- data.frame(matrix(nrow= length(summary_epitopes), ncol = 4))
colnames(summary_similar)<-c("Epitope","ImmuneCODE entries","Similar M1 entries","Similar M1/ImmuneCODE")
for(i in 1:length(summary_epitopes)){
  summary_similar$Epitope[i]<-summary_epitopes[i]
  summary_similar$`ImmuneCODE entries`[i]<-nrow(query_combined[[i]])
  summary_similar$`Similar M1 entries`[i]<-nrow(similar_combined[[i]])
  summary_similar$`Similar M1/ImmuneCODE`[i]<-nrow(similar_combined[[i]])/nrow(query_combined[[i]])
  
}
write.table(
  summary_similar,
  file = paste0("Data/Output/GIANA", "/", "M1_MART1_TSCAN", "_GIANA_Summary.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names=TRUE
)

#####Peptide Similarity Search#######
# Peptide comparison between TSCAN epitopes and M1 sequence. This returns the distance between 2 peptides
# COVID TSCAN sequences:  LLLDRLNQL, KLWAQCVQL, YLQPRTFLL, ALWEIQQVV, YLFDESGEFKL, LLYDANYFL
pairwiseAlignment(
  AAString('GILGFVFTL'),
  subject = AAString('LLYDANYFL'),
  type = 'global',
  substitutionMatrix = "BLOSUM62",
  gapOpening = -6,
  gapExtension = 0,
  scoreOnly = F
)
