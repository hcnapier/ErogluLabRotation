# Connexin 43 Image Ananlysis
# All brains
# Hailey Napier
# December 6, 2022

# 0.0 Setup ----
## Load packages ----
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tibble")
install.packages("RColorBrewer")
library(dplyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)

## Load colocalization data ----
setwd("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/image_analysis/analysis_output/Cx43/")
concat_20221201 <- read.csv("20221201_Cx43_colocalization_concat.csv")
concat_20221201$X.2 <- NULL
concat_20221201$X.1 <- NULL
concat_20221202 <- read.csv("20221202_Cx43_colocalization_concat.csv")
cx43_all <- bind_rows(concat_20221201, concat_20221202)
brain_genotype <- cx43_all %>%
  select(Brain) %>%
  unique()
brain_genotype$Genotype <- c("WT", "KO", "KO", "WT", "WT")
label_brain_genotype <- cx43_all %>%
  select(Label, Brain) %>%
  unique()
label_brain_genotype <- full_join(label_brain_genotype, brain_genotype, by = "Brain")
cx43_all <- full_join(cx43_all, label_brain_genotype)
cx43_all$Genotype <- factor(cx43_all$Genotype, levels = c("WT", "KO"))

# 1.0 Number of puncta ----
## Find number of puncta per image ----
num_puncta_label <- data.frame(table(cx43_all$Label))
colnames(num_puncta_label) <- c("Label", "NumberPuncta")
num_puncta_label <- left_join(num_puncta_label, label_brain_genotype, by = "Label")

## Find average number of puncta per brain ----
av_puncta_brain <- data.frame(tapply(num_puncta_label$NumberPuncta, num_puncta_label$Brain, mean))
av_puncta_brain <- rownames_to_column(av_puncta_brain, "Brain")
colnames(av_puncta_brain) <- c("Brain", "AvNumPuncta")
av_puncta_brain <- left_join(av_puncta_brain, brain_genotype)

av_puncta_summary <- av_puncta_brain %>%
  group_by(Genotype) %>%
  summarize(
    sd = sd(AvNumPuncta), 
    mean = mean(AvNumPuncta)) %>%
  data.frame()


## Plot ----
dodge <- position_dodge(width = 0.25)
level_order = c("WT", "KO")
ggplot(data = av_puncta_brain, aes(x = Genotype, y = AvNumPuncta, fill = Brain)) + 
  geom_dotplot(data = num_puncta_label, 
               aes(x = Genotype, y = NumberPuncta), 
               binaxis = "y", 
               stackdir = "center",
               dotsize = 0.75,
               alpha = 0.25) + 
  geom_errorbar(inherit.aes = FALSE, data = av_puncta_summary, 
                aes(x = Genotype, ymin = mean-sd, ymax = mean+sd), width = 0.25) +
  geom_dotplot(binaxis = "y", 
               stackdir = "center", 
               position = dodge) + 
  labs(y="Average Number of Puncta Per Image") + 
  scale_x_discrete(limits = level_order) +
  theme_minimal()

# 2.0 Area of Puncta ----
## Find average area per image ----
av_area_label <- data.frame(tapply(cx43_all$Area, cx43_all$Label, mean))
av_area_label <- rownames_to_column(av_area_label, "Label")
colnames(av_area_label) <- c("Label", "AvAreaPuncta")
av_area_label <- left_join(av_area_label, label_brain_genotype)

## Find average area per brain ----
av_area_brain <- data.frame(tapply(cx43_all$Area, cx43_all$Brain, mean))
av_area_brain <- rownames_to_column(av_area_brain, "Brain")
colnames(av_area_brain) <- c("Brain", "AvAreaPuncta")
av_area_brain <- left_join(av_area_brain, brain_genotype)

av_area_summary <- av_area_brain %>%
  group_by(Genotype) %>%
  summarize(
    sd = sd(AvAreaPuncta), 
    mean = mean(AvAreaPuncta)) %>%
  data.frame()

## Plot ----
# Not sure what the units here are
dodge <- position_dodge(width = 0.25)
ggplot(data = av_area_brain, aes(x = Genotype, y = AvAreaPuncta, fill = Brain)) + 
  geom_dotplot(data = av_area_label, 
               aes(x = Genotype, y = AvAreaPuncta), 
               binaxis = "y", 
               stackdir = "center",
               dotsize = 0.75,
               alpha = 0.25) + 
  geom_errorbar(inherit.aes = FALSE, data = av_area_summary, 
                aes(x = Genotype, ymin = mean-sd, ymax = mean+sd), width = 0.25) +
  geom_dotplot(binaxis = "y", 
               stackdir = "center") + 
  labs(y=bquote("Average Puncta Area Per Image " (pixels ^2))) + 
  scale_x_discrete(limits = level_order) +
  theme_minimal()


# 3.0 Puncta Integrated Density ----
## Find average integrated density per image ----
av_intden_label <- data.frame(tapply(cx43_all$IntDen, cx43_all$Label, mean))
av_intden_label <- rownames_to_column(av_intden_label, "Label")
colnames(av_intden_label) <- c("Label", "AvPunctaIntDen")
av_intden_label <- left_join(av_intden_label, label_brain_genotype)

## Find average area per brain ----
av_intden_brain <- data.frame(tapply(cx43_all$IntDen, cx43_all$Brain, mean))
av_intden_brain <- rownames_to_column(av_intden_brain, "Brain")
colnames(av_intden_brain) <- c("Brain", "AvPunctaIntDen")
av_intden_brain <- left_join(av_intden_brain, brain_genotype)

av_intden_summary <- av_intden_brain %>%
  group_by(Genotype) %>%
  summarize(
    sd = sd(AvPunctaIntDen), 
    mean = mean(AvPunctaIntDen)) %>%
  data.frame()

## Plot ----
dodge <- position_dodge(width = 0.25)
ggplot(data = av_intden_brain, aes(x = Genotype, y = AvPunctaIntDen, fill = Brain)) + 
  geom_dotplot(data = av_intden_label, 
               aes(x = Genotype, y = AvPunctaIntDen), 
               binaxis = "y", 
               stackdir = "center",
               dotsize = 0.75,
               alpha = 0.25) +
  geom_errorbar(inherit.aes = FALSE, data = av_intden_summary, 
                aes(x = Genotype, ymin = mean-sd, ymax = mean+sd), width = 0.25) +
  geom_dotplot(binaxis = "y", 
               stackdir = "center") +
  labs(y="Average Puncta Integrated Density Per Image (a.u.)") + 
  scale_x_discrete(limits = level_order) +
  theme_minimal()
