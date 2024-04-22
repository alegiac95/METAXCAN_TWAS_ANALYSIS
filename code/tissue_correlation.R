rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

library(tidyverse)
library(ComplexHeatmap)


twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"
figures_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/figures/"

diseases <- c("ADHD", "ASD", "AN", "BD", "MDD", "OCD", "SCZ")


for (disease in diseases){
  
  twas_path <- paste0(twas_path, disease)
  
  predicxcan_files <- list.files(path = twas_path,
                                 pattern = "en_PGC_.*_en_Brain_.*\\.csv$",
                                 full.names = TRUE)
  
  dfs <- list()
  
  for (pred_file in predicxcan_files){
    
    # Extract file name
    f_name <- stringr::str_replace(pred_file, twas_path, "")|> 
      stringr::str_replace(paste0("/en_PGC_", disease, "_en_Brain_"), "") |> 
      stringr::str_replace(".csv", "")
    
    # Extract the data
    data <- readr::read_csv(pred_file) |>
      tidyr::drop_na(pvalue) |>
      dplyr::select(gene_name) |>
      dplyr::arrange(gene_name)
    
    # Save the data in the list
    dfs[[f_name]] <- data
  }
  
  ribiosUtils::jaccardIndex(dfs[[1]]$gene_name, dfs[[2]]$gene_name)
  
  jaccard_matrix <- ribiosUtils::listOverlapCoefficient(dfs)
  breaks <- seq(0, 1, length.out = 15)
  
  lab_names <- stringr::str_replace_all(
    base::rownames(jaccard_matrix), "_basal_ganglia", "") |> 
    stringr::str_replace_all("_", "\n") |> 
    str_remove("BA9|BA24")
  
  p<-ComplexHeatmap::pheatmap(jaccard_matrix,
                           breaks = breaks,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           labels_row = lab_names,
                           labels_col = lab_names,
                           color = rev(RColorBrewer::brewer.pal(9, name = "GnBu")),
                           display_numbers = TRUE,
                           number_color = "black",
                           cellwidth = 40, cellheight = 40,
                           name = "Gene\nOverlap",
                           main = paste0(
                             disease, 
                             " PrediXcan Overlap of Genes Across GTEx Tissues"),
                          )
  dev.new()
  png(file=paste0(figures_path, disease, "_gtex_similarity.png"),
      width=8.5,height=8.5,units="in",res = 310)
  draw(p)
  dev.off()
  twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"

}


  