## ------------ GENE ANALYSIS OF GENE INTERSECTION ---------------------------
# In this part the analysis will focus on the results from the TWAS analysis. 
# The analysis to be done will be to try to make sense of the genes involved in
# the neuropsychiatric diseases and see whether common genes shared by diseases
# (if any) can help in the understanding of the heterogeneity of symptoms.

# Overlap of genes between disorders (ADHD, AN, ASD, BD, MDD, OCD, SCZ) 
# Top 10% only
library(tidyverse)
library(ggplot2)
library(ggsankey)
library(rlang)

twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"

adhd <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_top_genes_thr_10.tsv"))
an <- readr::read_tsv(paste0(twas_path,"AN/AN_top_genes_thr_10.tsv"))
asd <- readr::read_tsv(paste0(twas_path,"ASD/ASD_top_genes_thr_10.tsv"))
bd <- readr::read_tsv(paste0(twas_path,"BD/BD_top_genes_thr_10.tsv"))
mdd <- readr::read_tsv(paste0(twas_path,"MDD/MDD_top_genes_thr_10.tsv"))
ocd <- readr::read_tsv(paste0(twas_path,"OCD/OCD_top_genes_thr_10.tsv"))
scz <- readr::read_tsv(paste0(twas_path,"SCZ/SCZ_top_genes_thr_10.tsv"))


