# ChromHMM binarization file setup
# Hailey Napier
# November 22, 2022

# 0.0 Install packages ----
install.packages("tidyr")
install.packages("tidyverse")
library(tidyr)
library(tidyverse)


# 1.0 Create binarization dataframes ----
bams <- c("H3K27ac_HET_female.bam", 
          "H3K27ac_HET_male.bam",
          "H3K27ac_KO_female.bam",
          "H3K27ac_KO_male.bam",
          "H3K27ac_WT_female.bam",
          "H3K27ac_WT_male.bam",
          "H3K27me3_HET_female.bam",
          "H3K27me3_HET_male.bam",
          "H3K27me3_KO_female.bam",
          "H3K27me3_KO_male.bam",
          "H3K27me3_WT_female.bam",
          "H3K27me3_WT_male.bam",
          "H3K36me3_HET_female.bam",
          "H3K36me3_HET_male.bam",
          "H3K36me3_KO_female.bam",
          "H3K36me3_KO_male.bam",
          "H3K36me3_WT_female.bam",
          "H3K36me3_WT_male.bam",
          "H3K4me1_HET_female.bam",
          "H3K4me1_HET_male.bam",
          "H3K4me1_KO_female.bam",
          "H3K4me1_KO_male.bam",
          "H3K4me1_WT_female.bam",
          "H3K4me1_WT_male.bam",
          "H3K4me3_HET_female.bam",
          "H3K4me3_HET_male.bam",
          "H3K4me3_KO_female.bam",
          "H3K4me3_KO_male.bam",
          "H3K4me3_WT_female.bam",
          "H3K4me3_WT_male.bam",
          "H3K9me3_HET_female.bam",
          "H3K9me3_HET_male.bam",
          "H3K9me3_KO_female.bam",
          "H3K9me3_KO_male.bam",
          "H3K9me3_WT_female.bam",
          "H3K9me3_WT_male.bam")
bam_controls <- c("H3Total_HET_female.bam",
                  "H3Total_HET_male.bam",
                  "H3Total_KO_female.bam",
                  "H3Total_KO_male.bam",
                  "H3Total_WT_female.bam",
                  "H3Total_WT_male.bam")

bam_binarization <- data.frame(bams) 
bam_binarization <- bam_binarization %>% separate(bams, c("mark", "cell_type"), sep = "_", remove = F, convert = F)

# Change order of columns
bam_binarization <- bam_binarization[, c(3,2,1)]

# Add a new column for controls
bam_binarization <- bam_binarization %>%
  mutate(
    control = case_when(
      str_detect(bams, "HET") & str_detect(bams, "female") ~ bam_controls[1],
      str_detect(bams, "HET") & str_detect(bams, "male") ~ bam_controls[2],
      str_detect(bams, "KO") & str_detect(bams, "female") ~ bam_controls[3],
      str_detect(bams, "KO") & str_detect(bams, "male") ~ bam_controls[4], 
      str_detect(bams, "WT") & str_detect(bams, "female") ~ bam_controls[5], 
      str_detect(bams, "WT") & str_detect(bams, "male") ~ bam_controls[6]))

# Make binarization file for concatenated version
concat_bam_binarization <- bam_binarization

# 2.0 Save binarization files ----
# cat file
write.table(concat_bam_binarization, "/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/nlgn2_chromhmm/bam_binarization/concatenated/nlgn2_concatenated_bam_binarization.tsv", col.names = F, row.names = F, quote = F, sep = "\t")

