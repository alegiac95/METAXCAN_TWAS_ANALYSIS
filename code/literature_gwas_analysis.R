rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

library(tidyverse)
library(readxl)

# ----- Helper functions
top_gene_frac <- function(df1, thr){
  df1 <- df1 |>
    dplyr::arrange(p.value) |>
    dplyr::top_n(base::round(dplyr::n() * thr), wt = -p.value) |>
    dplyr::select(gene_name,logFC , p.value) 
  
  return(df1)
}

# Calculate weighted average
calculate_weighted_avg <- function(twas_df, gene_data_df) {
  weighted_avg <- gene_data_df |>
    tidyr::pivot_longer(cols = -label, 
                        names_to = "gene_name", 
                        values_to = "expression") |>
    dplyr::left_join(twas_df, by = "gene_name") |>
    dplyr::select(gene_name, label, expression, logFC) |> 
    tidyr::drop_na(logFC) |>
    group_by(label) |>
    summarise(weighted_avg = sum(expression * logFC) / sum(logFC))
}


#-----

twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/GWAS/mmc2(1).xlsx"
figures_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/figures/"
save_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"
ahba_data <- readr::read_csv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHBA_data_no_norm.csv") # data from the allen human brain atlas
ahpa_data <- readr::read_tsv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHPA_mrna_brain.tsv") # data from the allen human protein atlas


dk_labels <- readr::read_csv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/DK-files/DesikanKilliany_atlas.csv") |>
  dplyr::select(id, label, hemisphere, structure)

# Find intersection of the two datasets in the gene names
brain_genes <- colnames(ahba_data)
common_genes <- intersect(brain_genes, ahpa_data$Gene)


ahba_data_clean <- ahba_data |> 
  dplyr::select(dplyr::one_of(common_genes))
ahba_data_clean$label <- ahba_data$label

diseases <- c("ADHD", "ASD", "AN", "BD", "MDD", "OCD", "SCZ")

gene_data <- readxl::read_xlsx(twas_path, sheet = 2) |>
  dplyr::select(gene_name, logFC, p.value, Dx) |>
  dplyr::mutate(Dx = stringr::str_to_upper(Dx)) |>
  dplyr::mutate(Dx = if_else(Dx == "BP", "BD", Dx)) |>
  dplyr::mutate(Dx = as.factor(Dx)) |>
  dplyr::filter(Dx %in% diseases)


gene_data |>
  ggplot(aes(x = logFC, y = -log10(p.value), colour = -log10(p.value))) +
  facet_wrap(vars(Dx)) +
  scale_colour_viridis_c(option = "B") +
  geom_point() + 
  theme_classic()

# Calculate transcriptomic risk score
asd_genes <- gene_data |>
  dplyr::filter(Dx == "ASD") |>
  dplyr::select(gene_name, logFC, p.value) |>
  dplyr::filter(gene_name %in% common_genes) |>
  dplyr::arrange(p.value)

bd_genes <- gene_data |>
  dplyr::filter(Dx == "BD") |>
  dplyr::select(gene_name, logFC, p.value) |>
  dplyr::filter(gene_name %in% common_genes) |>
  dplyr::arrange(p.value)

mdd_genes <- gene_data |>
  dplyr::filter(Dx == "MDD") |>
  dplyr::select(gene_name, logFC, p.value) |>
  dplyr::filter(gene_name %in% common_genes) |>
  dplyr::arrange(p.value)

scz_genes <- gene_data |>
  dplyr::filter(Dx == "SCZ") |>
  dplyr::select(gene_name, logFC, p.value) |>
  dplyr::filter(gene_name %in% common_genes) |>
  dplyr::arrange(p.value)

thresholds <- c(0.1, 0.05, 0.01)


for (thr in thresholds) {
  # ----
  # ASD
  top_asd <- top_gene_frac(asd_genes, thr = thr)
  readr::write_tsv(top_asd, paste0(save_path,"ASD/ASD_literature_top_genes_thr_",thr*100,".tsv"))
  # weighted average
  gw_asd <- calculate_weighted_avg(top_asd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_asd, paste0(save_path,"ASD/ASD_literature_TPRS_thr_",thr*100,".tsv"))
  #-----
  # BD
  top_bd <- top_gene_frac(bd_genes, thr = thr)
  readr::write_tsv(top_bd, paste0(save_path,"BD/BD_literature_top_genes_thr_",thr*100,".tsv"))
  # weighted average
  gw_bd <- calculate_weighted_avg(top_bd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_bd, paste0(save_path,"BD/BD_literature_TPRS_thr_",thr*100,".tsv"))
  #------
  # MDD
  top_mdd <- top_gene_frac(mdd_genes, thr = thr)
  readr::write_tsv(top_mdd, paste0(save_path,"MDD/MDD_literature_top_genes_thr_",thr*100,".tsv"))
  # weighted average
  gw_mdd <- calculate_weighted_avg(top_mdd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_asd, paste0(save_path,"MDD/MDD_literature_TPRS_thr_",thr*100,".tsv"))
  #--------
  # SCZ
  top_scz <- top_gene_frac(scz_genes, thr = thr)
  readr::write_tsv(top_scz, paste0(save_path,"SCZ/SCZ_literature_top_genes_thr_",thr*100,".tsv"))
  # weighted average
  gw_scz <- calculate_weighted_avg(top_scz, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_scz, paste0(save_path,"SCZ/SCZ_literature_TPRS_thr_",thr*100,".tsv"))
}
