
# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

library(ggplot2)
library(ggsankey)


twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"
ahba_data <- readr::read_csv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHBA_data_no_norm.csv") # data from the allen human brain atlas
ahpa_data <- readr::read_tsv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHPA_mrna_brain.tsv") # data from the allen human protein atlas
figures_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/figures/"
data_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/"

brain_gmt <- fgsea::gmtPathways(paste0(data_path, "Geneset/BrainGMTv2_wGO_HumanOrthologs.gmt.txt"))

# Load genes df

adhd_10 <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_top_genes_thr_10.tsv")) |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))
asd_10 <- readr::read_tsv(paste0(twas_path,"ASD/ASD_top_genes_thr_10.tsv")) |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))
an_10 <- readr::read_tsv(paste0(twas_path,"AN/AN_top_genes_thr_10.tsv")) |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))
bd_10 <- readr::read_tsv(paste0(twas_path,"BD/BD_top_genes_thr_10.tsv")) |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))
mdd_10 <- readr::read_tsv(paste0(twas_path,"MDD/MDD_top_genes_thr_10.tsv")) |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))
ocd_10 <- readr::read_tsv(paste0(twas_path,"OCD/OCD_top_genes_thr_10.tsv")) |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))
scz_10 <- readr::read_tsv(paste0(twas_path,"SCZ/SCZ_top_genes_thr_10.tsv")) |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))
a <- a |>
  dplyr::mutate(rank = sign(z_mean) * (-log10(pvalue))) |>
  dplyr::arrange(desc(rank))

gene_list <- a$rank
names(gene_list) <- a$gene_name


test_result <- clusterProfiler::gseGO(geneList = gene_list,
                                      ont = "ALL",
                                      keyType = "SYMBOL",
                                      minGSSize = 3, 
                                      maxGSSize = 800, 
                                      pvalueCutoff = 0.05, 
                                      scoreType = "neg",
                                      verbose = TRUE, 
                                      OrgDb = 'org.Hs.eg.db', 
                                      pAdjustMethod = "none")


dotplot(test_result, showCategory=15, split=".sign", size = 9) + facet_grid(.~.sign) + theme(text = element_text(size = 9))

t_res <- fgsea::fgsea(pathways = brain_gmt, stats = gene_list, minSize = 5, maxSize = 500)


# ------------------------------
# Heatmap with all TWAS results

adhd <- adhd_twas |>
  dplyr::select(gene_name, z_mean, pvalue) |>
  dplyr:: mutate(disease = "ADHD")

asd <- asd_twas |>
  dplyr::select(gene_name, z_mean, pvalue) |>
  dplyr::mutate(disease = "ASD")

an <- an_twas |>
  dplyr::select(gene_name, z_mean, pvalue) |>
  dplyr::mutate(disease = "AN")

bd <- bd_twas |>
  dplyr::select(gene_name, z_mean, pvalue) |>
  dplyr::mutate(disease = "BD")

mdd <- mdd_twas |>
  dplyr::select(gene_name, z_mean, pvalue) |>
  dplyr::mutate(disease = "MDD")

ocd <- ocd_twas |>
  dplyr::select(gene_name, z_mean, pvalue) |>
  dplyr::mutate(disease = "OCD")

scz <- scz_twas |>
  dplyr::select(gene_name, z_mean, pvalue) |>
  dplyr::mutate(disease = "SCZ")


adhd_twas |>
  dplyr::select(gene_name, pvalue, z_mean) |>
  dplyr::mutate(regulation = if_else(z_mean > 0 , "UP", "DOWN")) |>
  ggplot(aes(x=abs(z_mean), y = -log10(pvalue))) +
  facet_grid(. ~ regulation) +
  geom_point(alpha = .7, aes(color=-log10(pvalue))) +
  scale_color_viridis_c(option = "G") + theme_classic()
  
