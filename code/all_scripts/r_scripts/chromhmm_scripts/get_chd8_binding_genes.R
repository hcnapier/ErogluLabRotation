# Find genes chd8 binds in E4 and E5
# Hailey Napier
# December 13, 2022

# Last modified December 13, 2022

# 0.1 Load packages ----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("dplyr")
install.packages("stringr")
install.packages("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)


# 0.2 Load data ----
overlap_KO <- read.table("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/chd8_chromhmm/chd8_overlap/overlap_output_bed_files/KO_overlap_chd8.bed", sep = "\t", header = F)
overlap_WT <- read.table("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/chd8_chromhmm/chd8_overlap/overlap_output_bed_files/WT_overlap_chd8.bed", sep = "\t", header = F)


# 1.0 Filter bed files for E4 and E5 to assign genes ----
## KO ----
overlap_KO_E4 <- overlap_KO %>%
  filter(V4 == "E4")
overlap_KO_E5 <- overlap_KO %>%
  filter(V4 == "E5")
## WT ----
overlap_WT_E5 <- overlap_WT %>%
  filter(V4 == "E5")
overlap_WT_E4 <- overlap_WT %>%
  filter(V4 == "E4")
## save .bed files ----
setwd("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/chd8_chromhmm/gene_id/state_overlap_files")
write.table(overlap_WT_E4, "overlap_WT_E4.bed", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(overlap_WT_E5, "overlap_WT_E5.bed", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(overlap_KO_E4, "overlap_KO_E4.bed", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(overlap_KO_E5, "overlap_KO_E5.bed", col.names = F, row.names = F, quote = F, sep = "\t")


# 2.0 Read in output files from bedtools intersect ----
setwd("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/chd8_chromhmm/gene_id/overlap_out/refseq")
KO_E4 <- read.table("KO_E4_refseq_chromhmm_chd8_bind.bed", sep = "\t", header = F)
KO_E5 <- read.table("KO_E5_refseq_chromhmm_chd8_bind.bed", sep = "\t", header = F)
WT_E4 <- read.table("WT_E4_refseq_chromhmm_chd8_bind.bed", sep = "\t", header = F)
WT_E5 <- read.table("WT_E5_refseq_chromhmm_chd8_bind.bed", sep = "\t", header = F)

KO_E4_gene_names <- get_gene_names(KO_E4)
KO_E5_gene_names <- get_gene_names(KO_E5)
WT_E4_gene_names <- get_gene_names(WT_E4)
WT_E5_gene_names <- get_gene_names(WT_E5)


# 3.0 Extract unique gene names ----
get_gene_names <- function(df){
  gene_names <- df[,9] %>%
    as.vector() %>%
    word(1, sep = "; ") %>%
    word(2, sep = " ") %>%
    unique()
  
  return(gene_names)
}

KO_E4_known_gene_names <- get_gene_names(KO_E4_known)
KO_E5_gene_names <- get_gene_names(KO_E5)
WT_E4_gene_names <- get_gene_names(WT_E4)
WT_E5_gene_names <- get_gene_names(WT_E5)


# 4.0 Get GO terms ----
## Cellular component ----
KO_E4_ego_cc <- enrichGO(gene = KO_E4_gene_names, 
                      keyType = "SYMBOL",
                      OrgDb = org.Mm.eg.db, 
                      ont = "CC",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

KO_E5_ego_cc <- enrichGO(gene = KO_E5_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "CC",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

WT_E4_ego_cc <- enrichGO(gene = WT_E4_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "CC",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

WT_E5_ego_cc <- enrichGO(gene = WT_E5_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "CC",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

## Biological process ----
KO_E4_ego_bp <- enrichGO(gene = KO_E4_gene_names, 
                      keyType = "SYMBOL",
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

KO_E5_ego_bp <- enrichGO(gene = KO_E5_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

WT_E4_ego_bp <- enrichGO(gene = WT_E4_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

WT_E5_ego_bp <- enrichGO(gene = WT_E5_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

## Molecular function ----
KO_E4_ego_mf <- enrichGO(gene = KO_E4_gene_names, 
                      keyType = "SYMBOL",
                      OrgDb = org.Mm.eg.db, 
                      ont = "MF",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

KO_E5_ego_mf <- enrichGO(gene = KO_E5_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "MF",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

WT_E4_ego_mf <- enrichGO(gene = WT_E4_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "MF",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

WT_E5_ego_mf <- enrichGO(gene = WT_E5_gene_names, 
                      keyType = "SYMBOL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "MF",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

## Dotplots ----
### WT E4 ----
dotplot(WT_E4_ego_bp) + ggtitle("WT E4 Biological Process")
dotplot(WT_E4_ego_cc) + ggtitle("WT E4 Cellular Component")
dotplot(WT_E4_ego_mf) + ggtitle("WT Molecular Function")
### KO E4 ----
dotplot(KO_E4_ego_bp) + ggtitle("Nlgn2 KO E4 Biological Process")
dotplot(KO_E4_ego_cc) + ggtitle("Nlgn2 KO E4 Cellular Component")
dotplot(KO_E4_ego_mf) + ggtitle("Nlgn2 KO E4 Molecular Function")
### WT E5 ----
dotplot(WT_E5_ego_bp) + ggtitle("WT E5 Biological Process")
dotplot(WT_E5_ego_cc) + ggtitle("WT E5 Cellular ")
dotplot(WT_E5_ego_mf)
### KO E5 ----
dotplot(KO_E5_ego_bp)
dotplot(KO_E5_ego_cc)
dotplot(KO_E5_ego_mf)

