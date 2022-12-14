# Find genes chd8 binds in E4 and E5 and get GO terms for those genes
# Hailey Napier
# December 13, 2022

# Last modified December 14, 2022

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
library(ggvenn)


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

KO_E4_gene_names <- get_gene_names(KO_E4)
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
dotplot(WT_E4_ego_bp) + ggtitle("E4 WT Biological Process") + theme(axis.text.y = element_text(size = 10))
dotplot(WT_E4_ego_cc) + ggtitle("E4 WT Cellular Component") + theme(axis.text.y = element_text(size = 10))
dotplot(WT_E4_ego_mf) + ggtitle("E4 WT Molecular Function") + theme(axis.text.y = element_text(size = 10))
### KO E4 ----
dotplot(KO_E4_ego_bp) + ggtitle("E4 Nlgn2 KO Biological Process") + theme(axis.text.y = element_text(size = 10))
dotplot(KO_E4_ego_cc) + ggtitle("E4 Nlgn2 KO Cellular Component") + theme(axis.text.y = element_text(size = 10))
dotplot(KO_E4_ego_mf) + ggtitle("E4 Nlgn2 KO Molecular Function") + theme(axis.text.y = element_text(size = 10))
### WT E5 ----
dotplot(WT_E5_ego_bp) + ggtitle("E5 WT Biological Process") + theme(axis.text.y = element_text(size = 10))
dotplot(WT_E5_ego_cc) + ggtitle("E5 WT Cellular Component") + theme(axis.text.y = element_text(size = 10))
dotplot(WT_E5_ego_mf) + ggtitle("E5 WT Molecular Function") + theme(axis.text.y = element_text(size = 10))
### KO E5 ----
dotplot(KO_E5_ego_bp) + ggtitle("E5 Nlgn2 KO Biological Process") + theme(axis.text.y = element_text(size = 10))
dotplot(KO_E5_ego_cc) + ggtitle("E5 Nlgn2 KO Cellular Component") + theme(axis.text.y = element_text(size = 10))
dotplot(KO_E5_ego_mf) + ggtitle("E5 Nlgn2 KO Molecular Function") + theme(axis.text.y = element_text(size = 10))

## Gene difference GO terms ----
### Genes that are in KO E4 but not WT E4 ----
ko_not_wt_e4_genes <- KO_E4_gene_names[!KO_E4_gene_names %in% WT_E4_gene_names]
#### Find GO terms ----
ko_not_wt_e4_ego_mf <- enrichGO(gene = ko_not_wt_e4_genes, 
                                keyType = "SYMBOL", 
                                OrgDb = org.Mm.eg.db, 
                                ont = "MF",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
ko_not_wt_e4_ego_bp <- enrichGO(gene = ko_not_wt_e4_genes, 
                                keyType = "SYMBOL", 
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
ko_not_wt_e4_ego_cc <- enrichGO(gene = ko_not_wt_e4_genes, 
                                keyType = "SYMBOL", 
                                OrgDb = org.Mm.eg.db, 
                                ont = "CC",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
test <- KO_E4_ego_bp@result$geneID[which(KO_E4_ego_bp$Description == "histone modification")]
test
#### Plot ----
dotplot(ko_not_wt_e4_ego_bp) + ggtitle("E4 Nlgn2 KO genes not in WT\nBiological Process") + theme(axis.text.y = element_text(size = 10))
dotplot(ko_not_wt_e4_ego_cc) + ggtitle("E4 Nlgn2 KO genes not in WT\nCellular Component") + theme(axis.text.y = element_text(size = 10))
dotplot(ko_not_wt_e4_ego_mf) + ggtitle("E4 Nlgn2 KO genes not in WT\nMolecular Function") + theme(axis.text.y = element_text(size = 10))

### Genes that are in WT E5 but not KO E5
wt_not_ko_e5_genes <- WT_E5_gene_names[!WT_E5_gene_names %in% KO_E5_gene_names]
#### Find GO terms ----
wt_not_ko_e5_ego_cc <- enrichGO(gene = wt_not_ko_e5_genes, 
                                keyType = "SYMBOL", 
                                OrgDb = org.Mm.eg.db, 
                                ont = "CC",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
wt_not_ko_e5_ego_bp <- enrichGO(gene = wt_not_ko_e5_genes, 
                                keyType = "SYMBOL", 
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
wt_not_ko_e5_ego_mf <- enrichGO(gene = wt_not_ko_e5_genes, 
                                keyType = "SYMBOL", 
                                OrgDb = org.Mm.eg.db, 
                                ont = "MF",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
#### Plot ----
dotplot(wt_not_ko_e5_ego_bp) + ggtitle("E5 WT genes not in  Nlgn2 KO\nBiological Process") + theme(axis.text.y = element_text(size = 10))
dotplot(wt_not_ko_e5_ego_cc) + ggtitle("E5 WT genes not in  Nlgn2 KO\nCellular Component") + theme(axis.text.y = element_text(size = 10))
dotplot(wt_not_ko_e5_ego_mf) + ggtitle("E5 WT genes not in  Nlgn2 KO\nMolecular Function") + theme(axis.text.y = element_text(size = 10))

### Venn diagrams to show gene similarities/differences ----
#### E4 ----
e4_genes_wt_ko <- list(WT_E4_gene_names, KO_E4_gene_names)
names(e4_genes_wt_ko) <- c("WT", "Nlgn2 KO")
ggvenn(e4_genes_wt_ko, c("WT", "Nlgn2 KO")) + ggtitle("E4 genes")
#### E5 ----
e5_genes_wt_ko <- list(WT_E5_gene_names, KO_E5_gene_names)
names(e5_genes_wt_ko) <- c("WT", "Nlgn2 KO")
ggvenn(e5_genes_wt_ko, c("WT", "Nlgn2 KO")) + ggtitle("E5 genes")
