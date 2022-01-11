## This is the R script to perform cross-reference with paired TRA-TRB info for cloning TCR
## VDJdb contains paired TRA-TRB information, whereas ImmuneCODE does not.
#### Peptides ####
# MART1: ELAGIGILTV
# M1:   GILGFVFTL
# LLL:  LLLDRLNQL
# LLY:  LLYDANYFL
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
library(RColorBrewer)
library(viridis)
library(ape)
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
  df<-subset(df,df$Gene=="TRB"|df$Gene=="TRA") # Include TRA
  names(df)[which(names(df)=="CDR3")]<-"amino_acid"
  names(df)[which(names(df)=="V")]<-"v_resolved"
  names(df)[which(names(df)=="J")]<-"j_resolved"
  df<-select(df,c(1:13,17)) # Include the last score column
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
M1<-read.csv("Data/VDJdb_M1_TCR_Paired.tsv",sep = "\t")
M1<-unique(rearrange_VDJdb(M1,peptide = "GILGFVFTL"))
# Read ImmuneCODE data
ImmuneCODE_raw<-rbind(read.csv("Data/peptide-detail-ci.csv"),
                      read.csv("Data/peptide-detail-cii.csv"))
ImmuneCODE<-rearrange_ImmuneCODE(ImmuneCODE_raw)
LLL<-unique(ImmuneCODE[grep("LLLDRLNQL",ImmuneCODE$Epitope,fixed = T),])
LLY<-unique(ImmuneCODE[grep("LLYDANYFL",ImmuneCODE$Epitope,fixed = T),])

#### Cross reference by exact match of CDR3b ####
LLL_M1_complex.id <- M1$complex.id[which(M1$amino_acid %in% LLL$amino_acid)]
LLL_M1<- subset(M1, M1$complex.id %in% LLL_M1_complex.id)

LLY_M1_complex.id <- M1$complex.id[which(M1$amino_acid %in% LLY$amino_acid)]
LLY_M1<- subset(M1, M1$complex.id %in% LLY_M1_complex.id)

M1_SARS <-subset(M1, amino_acid %in% ImmuneCODE$amino_acid) # In terms of VDJdb entry
SARS_M1 <-subset(ImmuneCODE, amino_acid %in% M1$amino_acid) # In terms of ImmuneCODE entries
## Output epitope list with frequencies ##
Peptide<-as.data.frame(table(SARS_M1$Epitope))

#### Cross reference by Similarity Match ####
# Combine M1 and SARS-CoV 2 data together and perform clustering first
# LLL$Epitope <- "LLL"
# LLY$Epitope <- "LLY"
# M1$Epitope <- "M1"
All<-rbind(select(ImmuneCODE,c("amino_acid","v_resolved","j_resolved","Epitope")),
           select(M1,c("amino_acid","v_resolved","j_resolved","Epitope")))
# Output into iSMART folder. 
write.table(
  All,
  file = paste0("Data/Output/iSMART", "/", "all", "_iSMART.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names=TRUE
)
# Output into GIANA folder. GIANA shares the same file structure as iSMART, but faster.
write.table(
  All,
  file = paste0("Data/Output/GIANA", "/", "all", "_GIANA.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names=TRUE
)

## Read cluster results
cluster_result <- read.table("Data/Output/GIANA/all_GIANA--RotationEncodingBL62.txt", sep = "\t")
colnames(cluster_result)<-c("amino_acid","clusterID","v_resolved","j_resolved","Epitope")

# Find the clusters containing the identified exact match TCRs
LLL_cross_reactive_cluster_ID <- subset(cluster_result, amino_acid %in% subset(LLL_M1, Gene == "TRB")$amino_acid)$clusterID
LLL_cross_reactive_cluster <-unique(subset(cluster_result, clusterID %in% LLL_cross_reactive_cluster_ID))
LLL_cross_reactive_similar<- LLL_cross_reactive_cluster[grep("LLYDANYFL",LLL_cross_reactive_cluster$Epitope,fixed = T),]
M1_cross_reactive_similar<- LLL_cross_reactive_cluster[grep("GILGFVFTL",LLL_cross_reactive_cluster$Epitope,fixed = T),]

## Plot tree diagram for similar TCRs in LLL dataset
# Score matrix for TCRs in LLL_cross_reactive_similar
score_mat <-
  data.frame(matrix(
    0,
    nrow(LLL_cross_reactive_similar),
    nrow(LLL_cross_reactive_similar)
  ))
colnames(score_mat) <- LLL_cross_reactive_similar$amino_acid
rownames(score_mat) <- LLL_cross_reactive_similar$amino_acid
for (i in 1:nrow(score_mat)) {
  for (j in 1:nrow(score_mat)) {
    score_mat[i,j] <- pairwiseAlignment(
      AAString(LLL_cross_reactive_similar$amino_acid[i]),
      subject = AAString(LLL_cross_reactive_similar$amino_acid[j]),
      type = 'local',
      substitutionMatrix = "BLOSUM62",
      gapOpening = 6,
      gapExtension = 6,
      scoreOnly = T
    )
  }
}
d <- as.dist((2-scale(score_mat)))
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "average" )
clus5 = cutree(hc1, 5)
# Plot the obtained dendrogram
plot(as.phylo(hc1),cex = 1,tip.color = brewer.pal(5, 'Set1')[clus5])
