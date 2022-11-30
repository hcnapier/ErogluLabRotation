# chd8 peak overlap with WT and KO chromatin states
# Hailey Napier
# November 29, 2022

# 0.0 Load packages ----
install.packages("dplyr")
install.packages("ggplot2")
install.packages("RColorBrewer")
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 0.1 Load data ----
overlap_KO <- read.table("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/nlgn2_chromhmm/chd8_overlap/overlap_output_bed_files/KO_overlap_chd8.bed", sep = "\t", header = F)
overlap_WT <- read.table("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/nlgn2_chromhmm/chd8_overlap/overlap_output_bed_files/WT_overlap_chd8.bed", sep = "\t", header = F)


# 1.0 Find percentage binding for each state -----
## WT ----
WT_binding <- data.frame(table(overlap_WT$V4)) %>%
  rename(state = Var1, 
         binding = Freq)

## KO ----
KO_binding <- data.frame(table(overlap_KO$V4)) %>%
  rename(state = Var1, 
         binding = Freq)

## Order states ----
states_ordered <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13", "E14", "E15")


# 2.0 Get percentage binding for each state ----
## WT ----
wt_sum_binding <- sum(WT_binding$binding)
WT_binding$percent_binding <- WT_binding$binding/wt_sum_binding*100
WT_binding$genotype <- "wt"

## KO ----
ko_sum_binding <- sum(KO_binding$binding)
KO_binding$percent_binding <- KO_binding$binding/ko_sum_binding*100
KO_binding$genotype <- "ko"

## Join KO and WT into one dataframe ----
WT_KO_binding <- bind_rows(KO_binding, WT_binding)
WT_KO_binding$state <- factor(WT_KO_binding$state, levels = states_ordered)
WT_KO_binding$genotype <- factor(WT_KO_binding$genotype, levels = c("wt", "ko"))


# 3.0 Plot percent binding ----
ggplot(WT_KO_binding) + 
  geom_col(aes(state, percent_binding, fill = genotype), 
           position = "dodge") +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal()

