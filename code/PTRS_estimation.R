# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects

library(tidyverse)

# Load the Gene data of the brain
ahba_data <- readr::read_csv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHBA_data_no_norm.csv") # data from the allen human brain atlas
ahpa_data <- readr::read_tsv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHPA_mrna_brain.tsv") # data from the allen human protein atlas


dk_labels <- readr::read_csv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/DK-files/DesikanKilliany_atlas.csv") |>
  dplyr::select(id, label, hemisphere, structure)
# Path for TWAS files
twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"

# Find intersection of the two datasets in the gene names
brain_genes <- colnames(ahba_data)
common_genes <- intersect(brain_genes, ahpa_data$Gene)


ahba_data_clean <- ahba_data |> 
  dplyr::select(dplyr::one_of(common_genes))
ahba_data_clean$label <- ahba_data$label


rm(brain_genes)


# Filter the data to only include the common genes, i.e., that are in the AHBA and are also mostly expressed in the brain.


adhd_twas <- readr::read_tsv(paste0(twas_path,"ADHD/PGC_ADHD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

asd_twas <- readr::read_tsv(paste0(twas_path,"ASD/PGC_ASD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

an_twas <- readr::read_tsv(paste0(twas_path,"AN/PGC_AN_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

bd_twas <- readr::read_tsv(paste0(twas_path,"BD/PGC_BD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |> 
  tidyr::drop_na(pvalue)

mdd_twas <- readr::read_tsv(paste0(twas_path,"MDD/PGC_MDD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

ocd_twas <- readr::read_tsv(paste0(twas_path,"OCD/PGC_OCD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

scz_twas <- readr::read_tsv(paste0(twas_path,"SCZ/PGC_SCZ_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)


## ------------ TOP PERCENT OF GENES PTRS  ---------------------------
# Here for each disease the transcriptomic polygenic risk scores are calculated,
# using a weighte average approach, where the weights is given by the Z-score of 
# the top percent of genes for each disease, while the score is from the AHBA.

thresholds <- c(0.1, 0.05, 0.01, 0.005, 0.001) # thresholds to be used for the 
# analysis. This are the top 10, 5, 1, 0.5 and 0.1 percent of genes.

# Calculate weighted average
calculate_weighted_avg <- function(twas_df, gene_data_df) {
  weighted_avg <- gene_data_df |>
    tidyr::pivot_longer(cols = -label, 
                        names_to = "gene_name", 
                        values_to = "expression") |>
    dplyr::left_join(twas_df, by = "gene_name") |>
    dplyr::select(gene_name, label, expression, z_mean) |> 
    tidyr::drop_na(z_mean) |>
    group_by(label) |>
    summarise(weighted_avg = sum(expression * z_mean) / sum(z_mean))
}

calculate_weighted_avg_abs <- function(twas_df, gene_data_df) {
  weighted_avg <- gene_data_df |>
    tidyr::pivot_longer(cols = -label, 
                        names_to = "gene_name", 
                        values_to = "expression") |>
    dplyr::left_join(twas_df, by = "gene_name") |>
    dplyr::select(gene_name, label, expression, z_mean) |> 
    tidyr::drop_na(z_mean) |>
    group_by(label) |>
    summarise(weighted_avg = sum(expression * abs(z_mean)) / sum(abs(z_mean)))
}

top_gene_frac <- function(df1, thr){
  df1 <- df1 |>
    dplyr::arrange(pvalue) |>
    dplyr::top_n(base::round(dplyr::n() * thr), wt = -pvalue) |>
    dplyr::select(gene_name, z_mean, pvalue) 
  
  return(df1)
}


top_positive_gene <- function(df1, thr){
  df1 <- df1 |>
    dplyr::select(gene_name, z_mean, pvalue) |>
    dplyr::mutate(upregulated = dplyr::if_else(z_mean > 0, "UP", "DOWN")) |>
    dplyr::filter(upregulated == "UP") |>
    dplyr::arrange(pvalue) |>
    dplyr::filter(row_number() <= n() * thr)
  
  return(df1)
}

top_negative_gene <- function(df1, thr){
  df1 <- df1 |>
    dplyr::select(gene_name, z_mean, pvalue) |>
    dplyr::mutate(upregulated = dplyr::if_else(z_mean > 0, "UP", "DOWN")) |>
    dplyr::filter(upregulated == "DOWN") |>
    dplyr::arrange(pvalue) |>
    dplyr::filter(row_number() <= n() * thr)
  
  return(df1)
}


for (thr in thresholds) {
  # ----
  # ADHD
  top_adhd <- top_gene_frac(adhd_twas, thr = thr)
  readr::write_tsv(top_adhd, paste0(twas_path,"ADHD/ADHD_top_genes_thr_",thr*100,".tsv"))
  # weighted average
  gw_adhd <- calculate_weighted_avg(top_adhd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_adhd, paste0(twas_path,"ADHD/ADHD_TPRS_thr_",thr*100,".tsv"))
  # weighted average absolute value
  gw_adhd_abs <- calculate_weighted_avg_abs(top_adhd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_adhd_abs, paste0(twas_path,"ADHD/ADHD_TPRS_abs_thr_",thr*100,".tsv"))
  # top positive
  top_pos_adhd <- top_positive_gene(adhd_twas, thr = thr)
  readr::write_tsv(top_pos_adhd, paste0(twas_path,"ADHD/ADHD_top_pos_genes_thr_",thr*100,".tsv"))
  
  gw_adhd_abs <- calculate_weighted_avg_abs(top_pos_adhd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_adhd_abs, paste0(twas_path,"ADHD/ADHD_TPRS_pos_thr_",thr*100,".tsv"))
  
  # top negative
  top_neg_adhd <- top_negative_gene(adhd_twas, thr = thr)
  readr::write_tsv(top_neg_adhd, paste0(twas_path,"ADHD/ADHD_top_neg_genes_thr_",thr*100,".tsv"))
  
  gw_adhd_abs <- calculate_weighted_avg_abs(top_neg_adhd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_adhd_abs, paste0(twas_path,"ADHD/ADHD_TPRS_neg_thr_",thr*100,".tsv"))
  rm(top_adhd, gw_adhd, gw_adhd_abs, top_pos_adhd, top_neg_adhd)
  # ----
  # ASD
  top_asd <- top_gene_frac(asd_twas, thr = thr)
  readr::write_tsv(top_asd, paste0(twas_path,"ASD/ASD_top_genes_thr_",thr*100,".tsv"))
  gw_asd <- calculate_weighted_avg(top_asd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_asd, paste0(twas_path,"ASD/ASD_TPRS_thr_",thr*100,".tsv"))
  gw_asd_abs <- calculate_weighted_avg_abs(top_asd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_asd_abs, paste0(twas_path,"ASD/ASD_TPRS_abs_thr_",thr*100,".tsv"))
  # top positive
  top_pos_asd <- top_positive_gene(asd_twas, thr = thr)
  readr::write_tsv(top_pos_asd, paste0(twas_path,"ASD/ASD_top_pos_genes_thr_",thr*100,".tsv"))
  
  gw_asd_abs <- calculate_weighted_avg_abs(top_pos_asd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_asd_abs, paste0(twas_path,"ASD/ASD_TPRS_pos_thr_",thr*100,".tsv"))
  
  # top negative
  top_neg_asd <- top_negative_gene(asd_twas, thr = thr)
  readr::write_tsv(top_neg_asd, paste0(twas_path,"ASD/ASD_top_neg_genes_thr_",thr*100,".tsv"))
  
  gw_asd_abs <- calculate_weighted_avg_abs(top_neg_asd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_asd_abs, paste0(twas_path,"ASD/ASD_TPRS_neg_thr_",thr*100,".tsv"))
  rm(top_asd, gw_asd, gw_asd_abs, top_pos_asd, top_neg_asd)
  
  # ----
  # AN
  top_an<- top_gene_frac(an_twas, thr = thr)
  readr::write_tsv(top_an, paste0(twas_path,"AN/AN_top_genes_thr_",thr*100,".tsv"))
  gw_an <- calculate_weighted_avg(top_an, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_an, paste0(twas_path,"AN/AN_TPRS_thr_",thr*100,".tsv"))
  gw_an_abs <- calculate_weighted_avg_abs(top_an, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_an_abs, paste0(twas_path,"AN/AN_TPRS_abs_thr_",thr*100,".tsv"))
  
  # top positive
  top_pos_an <- top_positive_gene(an_twas, thr = thr)
  readr::write_tsv(top_pos_an, paste0(twas_path,"AN/AN_top_pos_genes_thr_",thr*100,".tsv"))
  
  gw_an_abs <- calculate_weighted_avg_abs(top_pos_an, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_an_abs, paste0(twas_path,"AN/AN_TPRS_pos_thr_",thr*100,".tsv"))
  
  # top negative
  top_neg_an <- top_negative_gene(an_twas, thr = thr)
  readr::write_tsv(top_neg_an, paste0(twas_path,"AN/AN_top_neg_genes_thr_",thr*100,".tsv"))
  
  gw_an_abs <- calculate_weighted_avg_abs(top_neg_an, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_an_abs, paste0(twas_path,"AN/AN_TPRS_neg_thr_",thr*100,".tsv"))
  rm(top_an, gw_an, gw_an_abs, top_pos_an, top_neg_an)
  
  # ----
  # BD
  top_bd <- top_gene_frac(bd_twas, thr = thr)
  readr::write_tsv(top_bd, paste0(twas_path,"BD/BD_top_genes_thr_",thr*100,".tsv"))
  gw_bd <- calculate_weighted_avg(top_bd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_bd, paste0(twas_path,"BD/BD_TPRS_thr_",thr*100,".tsv"))
  gw_bd_abs <- calculate_weighted_avg_abs(top_bd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_bd_abs, paste0(twas_path,"BD/BD_TPRS_abs_thr_",thr*100,".tsv"))
  
  # top positive
  top_pos_bd <- top_positive_gene(bd_twas, thr = thr)
  readr::write_tsv(top_pos_bd, paste0(twas_path,"BD/BD_top_pos_genes_thr_",thr*100,".tsv"))
  
  gw_bd_abs <- calculate_weighted_avg_abs(top_pos_bd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_bd_abs, paste0(twas_path,"BD/BD_TPRS_pos_thr_",thr*100,".tsv"))
  
  # top negative
  top_neg_bd <- top_negative_gene(bd_twas, thr = thr)
  readr::write_tsv(top_neg_bd, paste0(twas_path,"BD/BD_top_neg_genes_thr_",thr*100,".tsv"))
  
  gw_bd_abs <- calculate_weighted_avg_abs(top_neg_bd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_bd_abs, paste0(twas_path,"BD/BD_TPRS_neg_thr_",thr*100,".tsv"))
  rm(top_bd, gw_bd, gw_bd_abs, top_pos_bd, top_neg_bd)
  
  # ----
  # MDD
  top_mdd <- top_gene_frac(mdd_twas, thr = thr)
  readr::write_tsv(top_mdd, paste0(twas_path,"MDD/MDD_top_genes_thr_",thr*100,".tsv"))
  gw_mdd <- calculate_weighted_avg(top_mdd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_mdd, paste0(twas_path,"MDD/MDD_TPRS_thr_",thr*100,".tsv"))
  gw_mdd_abs <- calculate_weighted_avg_abs(top_mdd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_mdd_abs, paste0(twas_path,"MDD/MDD_TPRS_abs_thr_",thr*100,".tsv"))
  
  # top positive
  top_pos_mdd <- top_positive_gene(mdd_twas, thr = thr)
  readr::write_tsv(top_pos_mdd, paste0(twas_path,"MDD/MDD_top_pos_genes_thr_",thr*100,".tsv"))
  
  gw_mdd_abs <- calculate_weighted_avg_abs(top_pos_mdd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_mdd_abs, paste0(twas_path,"MDD/MDD_TPRS_pos_thr_",thr*100,".tsv"))
  
  # top negative
  top_neg_mdd <- top_negative_gene(mdd_twas, thr = thr)
  readr::write_tsv(top_neg_mdd, paste0(twas_path,"MDD/MDD_top_neg_genes_thr_",thr*100,".tsv"))
  
  gw_mdd_abs <- calculate_weighted_avg_abs(top_neg_mdd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_mdd_abs, paste0(twas_path,"MDD/MDD_TPRS_neg_thr_",thr*100,".tsv"))
  rm(top_mdd, gw_mdd, gw_mdd_abs, top_pos_mdd, top_neg_mdd)
  
  # ----
  # OCD
  top_ocd <- top_gene_frac(ocd_twas, thr = thr)
  readr::write_tsv(top_ocd, paste0(twas_path,"OCD/OCD_top_genes_thr_",thr*100,".tsv"))
  gw_ocd <- calculate_weighted_avg(top_ocd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_ocd, paste0(twas_path,"OCD/OCD_TPRS_thr_",thr*100,".tsv"))
  gw_ocd_abs <- calculate_weighted_avg_abs(top_ocd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_ocd_abs, paste0(twas_path,"OCD/OCD_TPRS_abs_thr_",thr*100,".tsv"))
  
  # top positive
  top_pos_ocd <- top_positive_gene(ocd_twas, thr = thr)
  readr::write_tsv(top_pos_ocd, paste0(twas_path,"OCD/OCD_top_pos_genes_thr_",thr*100,".tsv"))
  
  gw_ocd_abs <- calculate_weighted_avg_abs(top_pos_ocd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_ocd_abs, paste0(twas_path,"OCD/OCD_TPRS_pos_thr_",thr*100,".tsv"))
  
  # top negative
  top_neg_ocd <- top_negative_gene(ocd_twas, thr = thr)
  readr::write_tsv(top_neg_ocd, paste0(twas_path,"OCD/OCD_top_neg_genes_thr_",thr*100,".tsv"))
  
  gw_ocd_abs <- calculate_weighted_avg_abs(top_neg_ocd, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_ocd_abs, paste0(twas_path,"OCD/OCD_TPRS_neg_thr_",thr*100,".tsv"))
  rm(top_ocd, gw_ocd, gw_ocd_abs, top_pos_ocd, top_neg_ocd)
  
  # ----
  # SCZ
  top_scz <- top_gene_frac(scz_twas, thr = thr)
  readr::write_tsv(top_scz, paste0(twas_path,"SCZ/SCZ_top_genes_thr_",thr*100,".tsv"))
  gw_scz <- calculate_weighted_avg(top_scz, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_scz, paste0(twas_path,"SCZ/SCZ_TPRS_thr_",thr*100,".tsv"))
  gw_scz_abs <- calculate_weighted_avg_abs(top_scz, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_scz_abs, paste0(twas_path,"SCZ/SCZ_TPRS_abs_thr_",thr*100,".tsv"))
  
  # top positive
  top_pos_scz <- top_positive_gene(scz_twas, thr = thr)
  readr::write_tsv(top_pos_scz, paste0(twas_path,"SCZ/SCZ_top_pos_genes_thr_",thr*100,".tsv"))
  
  gw_scz_abs <- calculate_weighted_avg_abs(top_pos_scz, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_scz_abs, paste0(twas_path,"SCZ/SCZ_TPRS_pos_thr_",thr*100,".tsv"))
  
  # top negative
  top_neg_scz <- top_negative_gene(scz_twas, thr = thr)
  readr::write_tsv(top_neg_scz, paste0(twas_path,"SCZ/SCZ_top_neg_genes_thr_",thr*100,".tsv"))
  
  gw_scz_abs <- calculate_weighted_avg_abs(top_neg_scz, ahba_data_clean)|>
    dplyr::left_join(dk_labels, by = join_by(label == id), keep = FALSE) |>
    dplyr::select(label.y, hemisphere, structure, everything()) |>
    dplyr::select(-label) |>
    dplyr::rename(label = label.y)
  readr::write_tsv(gw_scz_abs, paste0(twas_path,"SCZ/SCZ_TPRS_neg_thr_",thr*100,".tsv"))
  rm(top_scz, gw_scz, gw_scz_abs, top_pos_scz, top_neg_scz)
}






