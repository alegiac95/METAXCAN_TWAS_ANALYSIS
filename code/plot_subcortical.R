library(ggseg)
library(ggplot2)
library(tidyverse)
library(readxl)
library(ggsegExtra)
library(ggseg3d)

data_file <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/ENIGMA/ENIGMA_structural_organised.xlsx"
subcortical_data <- readxl::read_xlsx(data_file, sheet = 3)

area <- c('amygdala', 
          'caudate', 
          'hippocampus', 
          'pallidum', 
          'putamen', 
          'thalamus proper')

some_data <- data.frame(
  region = area,
  cohen = subcortical_data$SCZ[2:7],
  #hemi = rep("left", 6),
  #side = rep('coronal', 6),
  stringsAsFactors = FALSE)

ggplot() +
  geom_brain(atlas = aseg,
        data = some_data,
        mapping = aes(fill = cohen),
        #hemi = "left",
        #side = "coronal",
        colour = 'black') +
  scale_fill_gradient2(name = "Cohen's d", midpoint = 0, 
                       low = "navyblue", mid = "snow", high = "firebrick",
                       limits = c(-.5,.5),
                       na.value = "grey95") +
  theme_brain() +
  theme(legend.position = "left",
        legend.title.position = "left",
        legend.title = element_text(hjust = .5, angle = 90),
        legend.text = element_text(hjust = .5),
        axis.text.x=element_blank(),  
        axis.ticks.x=element_blank(),      
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(linetype = "solid", fill = "grey15"))

