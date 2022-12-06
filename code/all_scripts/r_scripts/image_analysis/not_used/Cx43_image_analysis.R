# Cx43 Image Analysis
# Hailey Napier
# December 1, 2022

# 0.0 Load packages ----
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)

# 0.1 Load data ----
# Images that couldn't be analyzed are excluded
## 34872-3 ----
setwd("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/image_analysis/analysis_output/Cx43/34872-3")
cx43_34872.3_20221117_2 <- read.csv("Colocalization_20221117_34872-3_cx43-488_mcherry-568_dapi_2.csv")
#cx43_34872.3_20221117_2$Image <- "20221117_2"
cx43_34872.3_20221117_3 <- read.csv("Colocalization_20221117_34872-3_cx43-488_mcherry-568_dapi_3.csv")
#cx43_34872.3_20221117_3$Image <- "20221117_3"
cx43_34872.3_20221117_4 <- read.csv("Colocalization_20221117_34872-3_cx43-488_mcherry-568_dapi_4.csv")
#cx43_34872.3_20221117_4$Image <- "20221117_4"
cx43_34872.3_20221117_5 <- read.csv("Colocalization_20221117_34872-3_cx43-488_mcherry-568_dapi_5.csv")
#cx43_34872.3_20221117_5$Image <- "20221117_5"
cx43_34872.3_20221117_6 <- read.csv("Colocalization_20221117_34872-3_cx43-488_mcherry-568_dapi_6.csv")
#cx43_34872.3_20221117_6$Image <- "20221117_6"
cx43_34872.3 <- bind_rows(cx43_34872.3_20221117_2, cx43_34872.3_20221117_3, cx43_34872.3_20221117_4, cx43_34872.3_20221117_5, cx43_34872.3_20221117_6)
cx43_34872.3$Brain <- "34872-3"

## 36092-1 ----
setwd("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/image_analysis/analysis_output/Cx43/36092-1")
cx43_36092.1_20221118_3 <- read.csv("Colocalization_20221118_36092-1_cx43-488_mcherry-568_dapi_3.csv")
#cx43_36092.1_20221118_3$Image <- "20221118_3"
cx43_36092.1_20221118_4 <- read.csv("Colocalization_20221118_36092-1_cx43-488_mcherry-568_dapi_4.csv")
#cx43_36092.1_20221118_4$Image <- "20221118_4"
cx43_36092.1 <- bind_rows(cx43_36092.1_20221118_3, cx43_36092.1_20221118_4)
cx43_36092.1$Brain <- "36092-1"

cx43_data_1 <- bind_rows(cx43_36092.1, cx43_34872.3)
write.csv(cx43_data_1,"../20221201_Cx43_colocalization_concat.csv")

# 1.0 Make combined dataframes for each brain ----
## Join dataframes ----
cx43_34872.3 <- bind_rows(cx43_34872.3_20221117_2, cx43_34872.3_20221117_3, cx43_34872.3_20221117_4, cx43_34872.3_20221117_5, cx43_34872.3_20221117_6)
cx43_36092.1 <- bind_rows(cx43_36092.1_20221118_3, cx43_36092.1_20221118_4)


# 2.0 Plot 36092-1 ----
## Area ----
ggplot(data = cx43_36092.1, aes(x = Image, y = Area)) +
  geom_boxplot() + 
  labs(title = "Connexin 43\nPuncta Area\n36092-1") +
  theme_minimal()

## Integrated density ----
ggplot(data = cx43_36092.1, aes(x = Image, y = IntDen)) +
  geom_boxplot() + 
  labs(title = "Connexin 43\nIntegrated Density\n36092-1") +
  theme_minimal()


# 3.0 Plot 34872-3 ----
## Area ----
ggplot(data = cx43_34872.3, aes(x = Image, y = Area)) +
  geom_boxplot() + 
  labs(title = "Connexin 43\nPuncta Area\n36092-1") +
  theme_minimal()

## Integrated density
ggplot(data = cx43_34872.3, aes(x = Image, y = IntDen)) +
  geom_boxplot() + 
  labs(title = "Connexin 43\nIntegrated Density\n36092-1") +
  theme_minimal()

