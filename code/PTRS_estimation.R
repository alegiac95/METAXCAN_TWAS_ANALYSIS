library(tidyverse)
library(ggplot2)
library(ggsankey)

# Load the Gene data of the brain
ahba_data <- readr::read_csv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHBA_data_no_norm.csv") # data from the allen human brain atlas
ahpa_data <- readr::read_tsv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHPA_mrna_brain.tsv") # data from the allen human protein atlas

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

top_gene_frac <- function(df1, thr){
  df1 <- df1 |>
    dplyr::arrange(pvalue) |>
    dplyr::top_n(base::round(dplyr::n() * thr), wt = -pvalue) |>
    dplyr::select(gene_name, z_mean, pvalue) 
  
  return(df1)
}

for (thr in thresholds) {
  # ADHD
  top_adhd <- top_gene_frac(adhd_twas, thr = thr)
  readr::write_tsv(top_adhd, paste0(twas_path,"ADHD/ADHD_top_genes_thr_",thr*100,".tsv"))
  gw_adhd <- calculate_weighted_avg(top_adhd, ahba_data_clean)
  readr::write_tsv(gw_adhd, paste0(twas_path,"ADHD/ADHD_TPRS_thr_",thr*100,".tsv"))
  # ASD
  top_asd <- top_gene_frac(asd_twas, thr = thr)
  readr::write_tsv(top_asd, paste0(twas_path,"ASD/ASD_top_genes_thr_",thr*100,".tsv"))
  gw_asd <- calculate_weighted_avg(top_adhd, ahba_data_clean)
  readr::write_tsv(gw_asd, paste0(twas_path,"ASD/ASD_TPRS_thr_",thr*100,".tsv"))
  # AN
  top_an<- top_gene_frac(an_twas, thr = thr)
  readr::write_tsv(top_an, paste0(twas_path,"AN/AN_top_genes_thr_",thr*100,".tsv"))
  gw_an <- calculate_weighted_avg(top_an, ahba_data_clean)
  readr::write_tsv(gw_an, paste0(twas_path,"AN/AN_TPRS_thr_",thr*100,".tsv"))
  # BD
  top_bd <- top_gene_frac(bd_twas, thr = thr)
  readr::write_tsv(top_bd, paste0(twas_path,"BD/BD_top_genes_thr_",thr*100,".tsv"))
  gw_bd <- calculate_weighted_avg(top_bd, ahba_data_clean)
  readr::write_tsv(gw_bd, paste0(twas_path,"BD/BD_TPRS_thr_",thr*100,".tsv"))
  # MDD
  top_mdd <- top_gene_frac(mdd_twas, thr = thr)
  readr::write_tsv(top_mdd, paste0(twas_path,"MDD/MDD_top_genes_thr_",thr*100,".tsv"))
  gw_mdd <- calculate_weighted_avg(top_mdd, ahba_data_clean)
  readr::write_tsv(gw_mdd, paste0(twas_path,"MDD/MDD_TPRS_thr_",thr*100,".tsv"))
  # OCD
  top_ocd <- top_gene_frac(ocd_twas, thr = thr)
  readr::write_tsv(top_ocd, paste0(twas_path,"OCD/OCD_top_genes_thr_",thr*100,".tsv"))
  gw_ocd <- calculate_weighted_avg(top_ocd, ahba_data_clean)
  readr::write_tsv(gw_ocd, paste0(twas_path,"OCD/OCD_TPRS_thr_",thr*100,".tsv"))
  # SCZ
  top_scz <- top_gene_frac(scz_twas, thr = thr)
  readr::write_tsv(top_scz, paste0(twas_path,"SCZ/SCZ_top_genes_thr_",thr*100,".tsv"))
  gw_scz <- calculate_weighted_avg(top_scz, ahba_data_clean)
  readr::write_tsv(gw_scz, paste0(twas_path,"SCZ/SCZ_TPRS_thr_",thr*100,".tsv"))
}




rm(top_adhd, top_asd, top_an, top_bd, top_mdd, top_ocd, top_scz,
   gw_adhd, gw_asd, gw_an, gw_bd, gw_mdd, gw_ocd, gw_scz)
