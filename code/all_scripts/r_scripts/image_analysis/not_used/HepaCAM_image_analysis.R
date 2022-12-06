# HepaCAM Image Analysis
# Hailey Napier
# December 1, 2022

# 0.0 Load packages ----
install.packages("ggplot2")
library(ggplot2)

# 0.1 Load data ----
## 36092-1 ----
setwd("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/image_analysis/analysis_output/HepaCAM/36092-1")
hepa_36092.1_20221118_1 <- read.csv("Colocalization_20221118_63x_36092-1_hepacam-488_mcherry-594_dapi.csv")
hepa_36092.1_20221118_1$Image <- "20221118_1"
hepa_36092.1_20221121_1 <- read.csv("Colocalization_20221121_36092-1_63x_hepacam-488_mcherry-594_dapi.csv")
hepa_36092.1_20221121_1$Image <- "20221121_1"
hepa_36092.1_20221121_2 <- read.csv("Colocalization_20221121_63x_36092-1_hepacam-488_mcherry-594_dapi_2.csv")
hepa_36092.1_20221121_2$Image <- "20221121_2"
hepa_36092.1_20221121_4 <- read.csv("Colocalization_20221121_63x_36092-1_hepacam-488_mcherry-594_dapi_4.csv")
hepa_36092.1_20221121_4$Image <- "20221121_4"
hepa_36092.1_20221121_3 <- read.csv("Colocalization_20221121_6xx_36092-1_hepacam-488_mcherry-594_dapi_3.csv")
hepa_36092.1_20221121_3$Image <- "20221121_3"
## 34872-3 ----
setwd("/Users/haileynapier/Work/Rotations/Rotation2_Eroglu/image_analysis/analysis_output/HepaCAM/34872-3")
hepa_34872.3_20221118_1 <- read.csv("Colocalization_20221118_63x_34872-3_hepacam-488_mcherry-594_dapi.csv")
hepa_34872.3_20221118_1$Image <- "20221118_1"
hepa_34872.3_20221121_1 <- read.csv("Colocalization_20221121_63x_34872-3_hepacam-488_mcherry-594_dapi.csv")
hepa_34872.3_20221121_1$Image <- "20221121_1"
hepa_34872.3_20221121_2 <- read.csv("Colocalization_20221121_63x_34892-3_hepacam-488_mcherry-594_dapi_2.csv")
hepa_34872.3_20221121_2$Image <- "20221121_2"


# 1.0 Make combined dataframes for each brain ----
## Join dataframes ----
hepa_36092.1 <- bind_rows(hepa_36092.1_20221118_1, hepa_36092.1_20221121_1, hepa_36092.1_20221121_2, hepa_36092.1_20221121_3, hepa_36092.1_20221121_4)
hepa_34872.3 <- bind_rows(hepa_34872.3_20221118_1, hepa_34872.3_20221121_1, hepa_34872.3_20221121_2)



# 2.0 Plot 36092-1 ----
## Area ----
ggplot(data = hepa_36092.1, aes(x = Image, y = Area)) +
  geom_boxplot() + 
  labs(title = "HepaCAM\nPuncta Area\n36092-1") +
  theme_minimal()

## Integrated Density ----
ggplot(data = hepa_36092.1, aes(x = Image, y = IntDen)) +
  geom_boxplot() + 
  labs(title = "HepaCAM\nIntegrated Density\n36092-1") + 
  theme_minimal()


# 3.0 Plot 34872.3 ----
## Area ----
ggplot(data = hepa_34872.3, aes(x = Image, y = Area)) + 
  geom_boxplot() +
  labs(title = "HepaCAM\nPuncta Area\n34872-3") + 
  theme_minimal()

## Integrated Density ----
ggplot(data = hepa_34872.3, aes(x = Image, y = IntDen)) + 
  geom_boxplot() + 
  labs(title = "HepaCAM\nIntegrated Density\n34872-3") + 
  theme_minimal()
