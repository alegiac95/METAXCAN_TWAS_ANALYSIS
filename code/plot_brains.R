library(ggseg) # subcortical data plot
library(ggplot2) # general plotting
library(tidyverse) # data science general
library(readxl) # read excel files
library(ggsegExtra) # additional ggseg stuff
library(RColorBrewer) # nice colors
library(fsbrain) # cortical plotting
library(magick) # image manipulation



# ---- 
# Helper functions and path set
zscore <- function(column){
  m <- mean(column)
  s <- sd(column)
  z <- ((column - m) / s)
  return(z)
}

setwd("/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/")
base_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/"
figures_path <- "figures/"


template_subjects_dir <- fsbrain::fsaverage.path()   
atlas <- 'aparc'
template_subject <- 'fsaverage'

# color palettes
colFn_imaging_data3d <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, name="RdBu")))
colFn_imaging_data <- rev(RColorBrewer::brewer.pal(11, name="RdBu"))
colFn_gene_data3d <- colorRampPalette(RColorBrewer::brewer.pal(11, name="PiYG"))
colFn_gene_data <- RColorBrewer::brewer.pal(11, name="PiYG")
# 3D plot options
makecmap_options <-list('range'=c(-3,3),
                        'colFn'=colFn_imaging_data3d)
makecmap_options2 <-list('range'=c(-2,2),
                        'colFn'=colFn_gene_data3d)
# ---- 
# Structural imaging



data_file <- "data/ENIGMA/ENIGMA_structural_organised.xlsx"

subcortical_data <- readxl::read_xlsx(data_file, sheet = 3)|>
  dplyr::mutate(Structure = stringr::str_to_lower(Structure)) |>
  dplyr::mutate(Structure = stringr::str_replace(Structure, 
                                                 "thalamusproper",
                                                 "thalamus"))

cortical_data <- readxl::read_xlsx(data_file, sheet = 2) |>
  dplyr::mutate(Structure = stringr::str_replace_all(
    stringr::str_replace_all(Structure,"M_",""),"_thickavg", ""))

tdata <- dplyr::bind_rows(cortical_data, subcortical_data) |>
  dplyr::mutate_at(vars(ADHD, ASD, AN, BD, MDD, OCD, SCZ), zscore) 

area <- c('accumbens area',
          'amygdala', 
          'caudate', 
          'hippocampus', 
          'pallidum', 
          'putamen', 
          'thalamus')

sub_data <- tdata |>
  dplyr::filter(Structure %in% area)
cort_data <- tdata |>
  dplyr::filter(!Structure %in% area) 

diseases <- c("ADHD", "ASD", "AN", "BD", "MDD", "OCD", "SCZ")

for (dis in diseases){
  
  save_dis_folder <- paste0(figures_path, dis)
  if (!dir.exists(save_dis_folder)) {
    # If the folder does not exist, create it
    dir.create(save_dis_folder)
    cat("Folder created:", save_dis_folder, "\n")
  }
  
  sub_save <- paste0(base_path,save_dis_folder,"/",dis,"_subcortical_enigma.png")
  cort_save <- paste0(base_path, save_dis_folder,"/",dis,"_cortical_enigma.png")
  
  # Subcortical PLot 
  some_data <- data.frame(
    region = sub_data$Structure,
    cohen = sub_data[[dis]],
    #hemi = rep("left", 6),
    #side = rep('coronal', 6),
    stringsAsFactors = FALSE) 
  
  
  
  
  sub_plot <- ggplot() +
    geom_brain(atlas = aseg2,
          data = some_data,
          #mapping = aes(fill = cohen),
          hemi = "right",
          colour = 'black') +
    #scale_fill_gradientn(colors = colFn_imaging_data, breaks = c(-3, 0, 3), limits = c(-3,3), oob = scales::squish) +
    theme_brain() +
    theme(legend.position = "bottom",
          legend.title.position = "right",
          legend.title = element_text(hjust = .5, angle = 90),
          legend.text = element_text(hjust = .5),
          axis.text.x=element_blank(),  
          axis.ticks.x=element_blank(),      
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_rect(linetype = "solid", fill = "white"))
  
  
  ggsave(sub_save, plot = sub_plot)
  
  # Cortical plot
  
  # For the left hemi, we manually set data values for some regions.
  lh_region_value_list <- as.list(cort_data[[dis]])
  names(lh_region_value_list) <- cort_data$Structure 
  
  
  
  
  cortical_structure <- vis.region.values.on.subject(template_subjects_dir, 
                                                     template_subject, atlas,
                                                     surface = 'pial',
                                                     lh_region_value_list, 
                                                     rh_region_value_list = NULL, 
                                                     makecmap_options = makecmap_options,
                                                     draw_colorbar = TRUE)
  
  
  img <- export(cortical_structure, colorbar_legend="Z score", 
               draw_colorbar = 'horizontal',
               view_angles = get.view.angle.names(angle_set = 'lh'),
               output_img = cort_save)
  
  rgl::close3d(dev = rgl::rgl.dev.list())
  
  
  cortical_img <- magick::image_read(cort_save)
  sub_img <- magick::image_read(sub_save) |>
    magick::image_scale("670x")
  
  img <- c(cortical_img, sub_img)
  
  stacked_img <- image_append(img, stack = TRUE)
  
  image_write(stacked_img, 
              path = paste0(save_dis_folder,"/", dis, "_enigma.png"), 
              format = "png")
}

# ----------
# GENETIC TPRS PLOTS
weight_scheme <- c("", "_pos", "_neg")
thresholds <- c(10, 5, 1)

for (dis in diseases){
  
  save_dis_folder <- paste0(figures_path, dis)
  if (!dir.exists(save_dis_folder)) {
    # If the folder does not exist, create it
    dir.create(save_dis_folder)
    cat("Folder created:", save_dis_folder, "\n")
  }
  
  for (weight in weight_scheme){
    if (weight == "") {
      save_n <- "wavg"
    } else {
      save_n <- weight
    }
    save_weigh_folder <- paste0(save_dis_folder, "/", save_n)
    if (!dir.exists(save_weigh_folder)) {
      # If the folder does not exist, create it
      dir.create(save_weigh_folder)
      cat("Folder created:", save_weigh_folder, "\n")
    }
    
    for (thr in thresholds){
      
      gdata_file <- paste0(base_path, "results/", dis, "/",
                           dis, "_TPRS", weight, "_thr_", thr, ".tsv")
      
      gene_data <- readr::read_delim(file = gdata_file)|>
        dplyr::filter(hemisphere == "L") |>
        dplyr::filter(structure %in% c("cortex", "subcortex/brainstem")) |>
        dplyr::mutate(Z = zscore(weighted_avg))
      
      cort_gene <- gene_data |>
        dplyr::filter(structure == "cortex")
      
      sub_gene <- gene_data |>
        dplyr::filter(structure == "subcortex/brainstem") |>
        dplyr::mutate(label= stringr::str_replace(label, 
                                                    "thalamusproper",
                                                    "thalamus proper")) |>
        dplyr::filter(label %in% area)
        
      
      # save paths
      sub_save <- paste0(base_path,save_weigh_folder,"/",dis,
                         "_TPRS",weight,"_subcortical_thr_",thr ,".png")
      cort_save <- paste0(base_path,save_weigh_folder,"/",dis,
                          "_TPRS",weight,"_cortical_thr_",thr ,".png")
      
      # Subcortical PLot 
      some_data <- data.frame(
        region = sub_gene$label,
        cohen = sub_gene$Z,
        #hemi = rep("left", 6),
        #side = rep('coronal', 6),
        stringsAsFactors = FALSE) 
      
      
      
      
      sub_plot <- ggplot() +
        geom_brain(atlas = aseg,
                   data = some_data,
                   mapping = aes(fill = cohen),
                   hemi = "left",
                   #side = "coronal",
                   colour = 'black') +
        scale_fill_gradientn(colors = colFn_gene_data, breaks = c(-2, 0, 2),
                             limits = c(-2,2), oob = scales::squish) +
        theme_brain() +
        theme(legend.position = "none",
              legend.title.position = "left",
              legend.title = element_text(hjust = .5, angle = 90),
              legend.text = element_text(hjust = .5),
              axis.text.x=element_blank(),  
              axis.ticks.x=element_blank(),      
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.background = element_rect(linetype = "solid", fill = "white"))
      
      
      ggsave(sub_save, plot = sub_plot)
      
      # Cortical
      lh_region_value_list <- as.list(cort_gene$Z)
      names(lh_region_value_list) <- cort_gene$label 
      
      
      
      
      cortical_structure <- vis.region.values.on.subject(template_subjects_dir, 
                                                         template_subject, atlas,
                                                         surface = 'pial',
                                                         lh_region_value_list, 
                                                         rh_region_value_list = NULL, 
                                                         makecmap_options = makecmap_options2,
                                                         draw_colorbar = TRUE)
      
      img <- export(cortical_structure, colorbar_legend="Z score", 
                    draw_colorbar = 'horizontal',
                    view_angles = get.view.angle.names(angle_set = 'lh'),
                    output_img = cort_save)
      
      rgl::close3d(dev = rgl::rgl.dev.list())
      
      cortical_img <- magick::image_read(cort_save)
      sub_img <- magick::image_read(sub_save) |>
        magick::image_scale("670x")
      
      img <- c(cortical_img, sub_img)
      
      stacked_img <- image_append(img, stack = TRUE)
      
      image_write(stacked_img, 
                  path = paste0(save_weigh_folder,"/", dis, "_TPRS",
                                weight,"_thr_",thr,".png"), 
                  format = "png")
      
    }
  }
}
  

  
  
  