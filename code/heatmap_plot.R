library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)
library(readxl)
library(patchwork)
library(ComplexHeatmap)

twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"
data_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/"

figures_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/figures/"

adhd_corr_cort <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_cortical_abs_correlations_results.tsv")) |>
  mutate(disease = "ADHD")
adhd_corr_subcort <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_subcortical_abs_correlations_results.tsv"))|>
  mutate(disease = "ADHD")
asd_corr_cort <- readr::read_tsv(paste0(twas_path,"ASD/ASD_cortical_abs_correlations_results.tsv"))|>
  mutate(disease = "ASD")
asd_corr_subcort <- readr::read_tsv(paste0(twas_path,"ASD/ASD_subcortical_abs_correlations_results.tsv"))|>
  mutate(disease = "ASD")
an_corr_cort <- readr::read_tsv(paste0(twas_path,"AN/AN_cortical_abs_correlations_results.tsv"))|>
  mutate(disease = "AN")
an_corr_subcort <- readr::read_tsv(paste0(twas_path,"AN/AN_subcortical_abs_correlations_results.tsv"))|>
  mutate(disease = "AN")
bd_corr_cort <- readr::read_tsv(paste0(twas_path,"BD/BD_cortical_abs_correlations_results.tsv"))|>
  mutate(disease = "BD")
bd_corr_subcort <- readr::read_tsv(paste0(twas_path,"BD/BD_subcortical_abs_correlations_results.tsv"))|>
  mutate(disease = "BD")
mdd_corr_cort <- readr::read_tsv(paste0(twas_path,"MDD/MDD_cortical_abs_correlations_results.tsv"))|>
  mutate(disease = "MDD")
mdd_corr_subcort <- readr::read_tsv(paste0(twas_path,"MDD/MDD_subcortical_abs_correlations_results.tsv"))|>
  mutate(disease = "MDD")
ocd_corr_cort <- readr::read_tsv(paste0(twas_path,"OCD/OCD_cortical_abs_correlations_results.tsv"))|>
  mutate(disease = "OCD")
ocd_corr_subcort <- readr::read_tsv(paste0(twas_path,"OCD/OCD_subcortical_abs_correlations_results.tsv"))|>
  mutate(disease = "OCD")
scz_corr_cort <- readr::read_tsv(paste0(twas_path,"SCZ/SCZ_cortical_abs_correlations_results.tsv"))|>
  mutate(disease = "SCZ")
scz_corr_subcort <- readr::read_tsv(paste0(twas_path,"SCZ/SCZ_subcortical_abs_correlations_results.tsv"))|>
  mutate(disease = "SCZ")

cortical_df <- dplyr::bind_rows(adhd_corr_cort, 
                                  asd_corr_cort, 
                                  an_corr_cort, 
                                  bd_corr_cort,
                                  mdd_corr_cort,
                                  ocd_corr_cort,
                                  scz_corr_cort) |>
  mutate(disease = as.factor(disease))

subcortical_df <- dplyr::bind_rows(adhd_corr_subcort, 
                                asd_corr_subcort, 
                                an_corr_subcort, 
                                bd_corr_subcort,
                                mdd_corr_subcort,
                                ocd_corr_subcort,
                                scz_corr_subcort) |>
  mutate(disease = as.factor(disease))


cort_mat <- as.matrix(cortical_df |>
  dplyr::select(threshold, disease, corr) |>
  tidyr::spread(key = disease, value = corr) |>
  tibble::column_to_rownames(var = "threshold"))

cort_annotation_matrix <- as.matrix(cortical_df |>
  dplyr::select(threshold, disease, p) |>
  tidyr::spread(key = disease, value = p)|>
  tibble::column_to_rownames(var = "threshold"))


cort_heat <- ComplexHeatmap::Heatmap(cort_mat,
                        name = "Corr",
                        row_title = "Threshold",
                        column_title = "Brain Disorder",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        row_names_side = "left",
                        rect_gp = gpar(col = "white", lwd = 2),
                        heatmap_width = unit(13, "cm"), 
                        heatmap_height = unit(7, "cm"),
                        col = circlize::colorRamp2(c(-1,0,1), c("blue", "#EEEEEE", "red")),
                        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(cort_annotation_matrix[i, j] < 0.05)
                            grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 15, fontface = "bold"))
                        })
cort_heat 

subcort_mat <- subcortical_df |>
  dplyr::select(threshold, disease, corr) |>
  tidyr::spread(key = disease, value = corr) |>
  tibble::column_to_rownames(var = "threshold")

sub_annotation_matrix <- as.matrix(subcortical_df %>%
                                 select(threshold, disease, p) %>%
                                 spread(key = disease, value = p) %>%
                                 column_to_rownames(var = "threshold"))

subcort_heat <- ComplexHeatmap::Heatmap(as.matrix(subcort_mat),
                                     name = "Corr",
                                     row_title = "Threshold",
                                     column_title = "Brain Disorder",
                                     cluster_rows = FALSE,
                                     cluster_columns = FALSE,
                                     row_names_side = "left",
                                     rect_gp = gpar(col = "white", lwd = 2),
                                     heatmap_width = unit(13, "cm"), 
                                     heatmap_height = unit(7, "cm"),
                                     col = circlize::colorRamp2(c(-1,0,1), c("blue", "#EEEEEE", "red")),
                                     column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                                     row_title_gp = gpar(fontsize = 15, fontface = "bold"),
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                       if(sub_annotation_matrix[i, j] < 0.05)
                                         grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 15, fontface = "bold"))
                                     })
subcort_heat 

heat_list <- cort_heat %v% subcort_heat
png(paste0(figures_path,"Correlations_abs.png"),width=14.5,height=16,units="cm",res=300)
heat_list
dev.off()

