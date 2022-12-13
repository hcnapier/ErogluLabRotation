# chd8 peak overlap with WT and KO chromatin states
# Hailey Napier
# November 29, 2022
# Last modified December 8,2022

# 0.0 Load packages ----
install.packages("dplyr")
install.packages("ggplot2")
install.packages("RColorBrewer")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# 0.1 Load data ----
overlap_KO <- read.table("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/chd8_chromhmm/chd8_overlap/overlap_output_bed_files/KO_overlap_chd8.bed", sep = "\t", header = F)
overlap_WT <- read.table("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/chd8_chromhmm/chd8_overlap/overlap_output_bed_files/WT_overlap_chd8.bed", sep = "\t", header = F)


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
WT_binding$Genotype <- "WT"

## KO ----
ko_sum_binding <- sum(KO_binding$binding)
KO_binding$percent_binding <- KO_binding$binding/ko_sum_binding*100
KO_binding$Genotype <- "Nlgn2 KO"

## Join KO and WT into one dataframe ----
WT_KO_binding <- bind_rows(KO_binding, WT_binding)
WT_KO_binding$state <- factor(WT_KO_binding$state, levels = states_ordered)
WT_KO_binding$Genotype <- factor(WT_KO_binding$Genotype, levels = c("WT", "Nlgn2 KO"))


# 3.0 Plot percent binding ----
ggplot(WT_KO_binding) + 
  geom_col(aes(state, percent_binding, fill = Genotype), 
           position = "dodge") +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Chromatin State", y = "Percent Binding") +
  theme_minimal()

# 4.0 Statistics, all states at once ----
## Contingency table ----
percent_binding_table <- xtabs(binding ~ state + Genotype, data = WT_KO_binding)
## Chi squared test excluding 14th state (since there's no binding there there aren't enough observations to have the x-squared work)
chisq.test(percent_binding_table[-14])
#x-squared = 49832, df = 28, p-value < 2.2e-16


# 5.0 Statistics, each state individually ----
## Function to run chi-squared test on an individual state
chisq_ind_states <- function(in_state){
  state_binding <- WT_KO_binding %>%
    filter(state == in_state)
  con_table <- xtabs(binding ~ Genotype, data = state_binding)
  state_chisq_test <- chisq.test(con_table)
  return(state_chisq_test)
}

## Make list of states
state_list <- paste("E", c(1:13, 15), sep = "")

## Make output dataframe
state_stats <- data.frame("state" = state_list, "x_squared" = rep(NA, 14), "p_value" = rep(NA, 14))

## Loop through states and run chisq_ind_states for each one, save X-squared value and p-value in dataframe
for(i in state_list){
  i_state <- i
  chisq <- chisq_ind_states(i_state)
  state_stats$x_squared[which(state_stats$state == i_state)] <- chisq[[1]]
  state_stats$p_value[which(state_stats$state == i_state)] <- chisq[[3]]
}


chisq_ind_states("E1")[[3]]


# 6.0 Filter bed files for E4 and E5 to assign genes ----
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



