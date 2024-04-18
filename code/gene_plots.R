library(tidyverse)
library(ggplot2)
library(ggsankey)
library(viridis)
library(ggrepel)

theme_set(theme_light() 
          +theme(
            axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
            axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
            plot.title = element_text(hjust = 0.5, face = 'bold')
          ))
twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"
ahba_data <- readr::read_csv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHBA_data_no_norm.csv") # data from the allen human brain atlas
ahpa_data <- readr::read_tsv("/Users/alessiogiacomel/Library/CloudStorage/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/AHBA/AHPA_mrna_brain.tsv") # data from the allen human protein atlas
figures_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/figures/"


# Find intersection of the two datasets in the gene names
brain_genes <- colnames(ahba_data)
common_genes <- intersect(brain_genes, ahpa_data$Gene)


adhd_twas <- readr::read_tsv(paste0(twas_path,"ADHD/PGC_ADHD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)


adhd_10 <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_top_genes_thr_10.tsv"))
adhd_05 <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_top_genes_thr_0.5.tsv"))
yint <- min(-log10(adhd_10$pvalue))


adhd_twas |>
  mutate(log10 = -log10(pvalue)) |>
  ggplot(aes(x=z_mean, y=log10)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey50", alpha = .5) +
  geom_hline(yintercept = 1.364, color = "#440154FF", linetype = 'dashed') +
  geom_point(alpha = .7, color = "grey", size =2) +
  geom_point(data = adhd_10, aes(x = z_mean, y = -log10(pvalue), color = -log10(pvalue)), size =3) +
  scale_color_viridis_c(name = "Top Genes",option="H", 
                        breaks = c(min(-log10(adhd_10$pvalue)), 
                                   max(-log10(adhd_10$pvalue))),
                        labels = c("10%", "0.1%"),
                        guide = guide_colorbar(position = "right",
    theme = theme(legend.ticks = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 12)))) +
  geom_text_repel(data = adhd_05, 
                  aes(x  = z_mean, 
                      y = -log10(pvalue), 
                      label = gene_name),
                  box.padding = .6, direction = "b") +
  annotate("text", x = 3.5, y = yint + 0.1, label = "Top 10%", 
           hjust = -0.1, vjust = 0, color = "#440154FF", size = 5, fontface = "bold") +
  xlim(-4.9, 5.1) +
  xlab("Z score") +
  ylab(expression(-log[10](p[value]))) +
  ggtitle("ADHD TWAS Genes")

ggsave(paste0(figures_path,"adhd_twas.png"), width = 7, height = 6, device = "png", dpi = "retina")

asd_twas <- readr::read_tsv(paste0(twas_path,"ASD/PGC_ASD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)


asd_10 <- readr::read_tsv(paste0(twas_path,"ASD/ASD_top_genes_thr_10.tsv"))
asd_05 <- readr::read_tsv(paste0(twas_path,"ASD/ASD_top_genes_thr_0.5.tsv"))
yint <- min(-log10(asd_10$pvalue))

asd_twas |>
  mutate(log10 = -log10(pvalue)) |>
  ggplot(aes(x=z_mean, y=log10)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey50", alpha = .5) +
  geom_hline(yintercept = yint, color = "#440154FF", linetype = 'dashed') +
  geom_point(alpha = .7, color = "grey", size = 2) +
  geom_point(data = asd_10, aes(x = z_mean, y = -log10(pvalue), color = -log10(pvalue)), size = 3) +
  scale_color_viridis_c(name = "Top Genes",option="H", 
                        breaks = c(min(-log10(asd_10$pvalue)), 
                                   max(-log10(asd_10$pvalue))),
                        labels = c("10%", "0.1%"),
                        guide = guide_colorbar(position = "right",
                                               theme = theme(legend.ticks = element_blank(),
                                                             legend.title = element_blank(),
                                                             legend.text = element_text(size = 12)))) +
  geom_text_repel(data = asd_05, 
                  aes(x  =z_mean, 
                      y = -log10(pvalue), 
                      label = gene_name),
                  box.padding = .6, direction = "b", max.overlaps = 20) +
  annotate("text", x = 3.5, y = yint + 0.1, label = "Top 10%", 
           hjust = -0.1, vjust = 0, color = "#440154FF", size = 5, fontface = "bold") +
  xlim(-6, 6) +
  xlab("Z score") +
  ylab(expression(-log[10](p[value]))) +
  ggtitle("ASD TWAS Genes")

ggsave(paste0(figures_path,"asd_twas.png"), width = 7, height = 6, device = "png", dpi = "retina")

an_twas <- readr::read_tsv(paste0(twas_path,"AN/PGC_AN_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

an_10 <- readr::read_tsv(paste0(twas_path,"AN/AN_top_genes_thr_10.tsv"))
an_05 <- readr::read_tsv(paste0(twas_path,"AN/AN_top_genes_thr_0.5.tsv"))
yint <- min(-log10(an_10$pvalue))

an_twas |>
  mutate(log10 = -log10(pvalue)) |>
  ggplot(aes(x=z_mean, y=log10)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey50", alpha = .5) +
  geom_hline(yintercept = yint, color = "#440154FF", linetype = 'dashed') +
  geom_point(alpha = .7, color = "grey", size = 2) +
  geom_point(data = an_10, aes(x = z_mean, y = -log10(pvalue), color = -log10(pvalue)), size = 3) +
  scale_color_viridis_c(name = "Top Genes",option="H", 
                        breaks = c(min(-log10(an_10$pvalue)), 
                                   max(-log10(an_10$pvalue))),
                        labels = c("Top 10%", "Top 0.1%"),
                        guide = guide_colorbar(position = "right",
                                               theme = theme(legend.ticks = element_blank(),
                                                             legend.title = element_blank(),
                                                             legend.text = element_text(size = 12)))) +
  geom_text_repel(data = an_05, 
                  aes(x  =z_mean, 
                      y = -log10(pvalue), 
                      label = gene_name),
                  box.padding = .6, direction = "b", max.overlaps = 20) +
  annotate("text", x = 3.5, y = yint + 0.1, label = "Top 10%", 
           hjust = -0.1, vjust = 0, color = "#440154FF", size = 5, fontface = "bold") +
  xlim(-6, 6) +
  xlab("Z score") +
  ylab(expression(-log[10](p[value]))) +
  ggtitle("AN TWAS Genes")

ggsave(paste0(figures_path,"an_twas.png"), width = 7, height = 6, device = "png", dpi = "retina")


bd_twas <- readr::read_tsv(paste0(twas_path,"BD/PGC_BD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |> 
  tidyr::drop_na(pvalue)

bd_10 <- readr::read_tsv(paste0(twas_path,"BD/BD_top_genes_thr_10.tsv"))
bd_05 <- readr::read_tsv(paste0(twas_path,"BD/BD_top_genes_thr_0.5.tsv"))
yint <- min(-log10(bd_10$pvalue))

bd_twas |>
  mutate(log10 = -log10(pvalue)) |>
  ggplot(aes(x=z_mean, y=log10)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey50", alpha = .5) +
  geom_hline(yintercept = yint, color = "#440154FF", linetype = 'dashed') +
  geom_point(alpha = .7, color = "grey", size =2) +
  geom_point(data = bd_10, aes(x = z_mean, y = -log10(pvalue), color = -log10(pvalue)), size = 3) +
  scale_color_viridis_c(name = "Top Genes",option="H", 
                        breaks = c(min(-log10(bd_10$pvalue)), 
                                   max(-log10(bd_10$pvalue))),
                        labels = c("10%", "0.1%"),
                        guide = guide_colorbar(position = "right",
                                               theme = theme(legend.ticks = element_blank(),
                                                             legend.title = element_blank(),
                                                             legend.text = element_text(size = 12)))) +
  geom_text_repel(data = bd_05, 
                  aes(x  =z_mean, 
                      y = -log10(pvalue), 
                      label = gene_name),
                  box.padding = .6, direction = "b", max.overlaps = 20) +
  annotate("text", x = 3.5, y = 1.45, label = "Top 10%", 
           hjust = -0.1, vjust = 0, color = "#440154FF", size = 5, fontface = "bold") +
  xlim(-6, 6) +
  xlab("Z score") +
  ylab(expression(-log[10](p[value]))) +
  ggtitle("BD TWAS Genes")

ggsave(paste0(figures_path,"bd_twas.png"), width = 7, height = 6, device = "png", dpi = "retina")


mdd_twas <- readr::read_tsv(paste0(twas_path,"MDD/PGC_MDD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

mdd_10 <- readr::read_tsv(paste0(twas_path,"MDD/MDD_top_genes_thr_10.tsv"))
mdd_05 <- readr::read_tsv(paste0(twas_path,"MDD/MDD_top_genes_thr_0.5.tsv"))
yint <- min(-log10(mdd_10$pvalue))

mdd_twas |>
  mutate(log10 = -log10(pvalue)) |>
  ggplot(aes(x=z_mean, y=log10)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey50", alpha = .5) +
  geom_hline(yintercept = yint, color = "#440154FF", linetype = 'dashed') +
  geom_point(alpha = .7, color = "grey", size =2) +
  geom_point(data = mdd_10, aes(x = z_mean, y = -log10(pvalue), color = -log10(pvalue)), size = 3) +
  scale_color_viridis_c(name = "Top Genes",option="H", 
                        breaks = c(min(-log10(mdd_10$pvalue)), 
                                   max(-log10(mdd_10$pvalue))),
                        labels = c("10%", "0.1%"),
                        guide = guide_colorbar(position = "right",
                                               theme = theme(legend.ticks = element_blank(),
                                                             legend.title = element_blank(),
                                                             legend.text = element_text(size = 12)))) +
  geom_text_repel(data = mdd_05, 
                  aes(x  =z_mean, 
                      y = -log10(pvalue), 
                      label = gene_name),
                  box.padding = .6, direction = "b", max.overlaps = 20) +
  annotate("text", x = 3.5, y = yint + 0.1, label = "Top 10%", 
           hjust = -0.1, vjust = 0, color = "#440154FF", size = 5, fontface = "bold") +
  xlim(-6, 6) +
  xlab("Z score") +
  ylab(expression(-log[10](p[value]))) +
  ggtitle("MDD TWAS Genes")

ggsave(paste0(figures_path,"mdd_twas.png"), width = 7, height = 6, device = "png", dpi = "retina")


ocd_twas <- readr::read_tsv(paste0(twas_path,"OCD/PGC_OCD_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

ocd_10 <- readr::read_tsv(paste0(twas_path,"OCD/OCD_top_genes_thr_10.tsv"))
ocd_05 <- readr::read_tsv(paste0(twas_path,"OCD/OCD_top_genes_thr_0.5.tsv"))
yint <- min(-log10(ocd_10$pvalue))

ocd_twas |>
  mutate(log10 = -log10(pvalue)) |>
  ggplot(aes(x=z_mean, y=log10)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey50", alpha = .5) +
  geom_hline(yintercept = yint, color = "#440154FF", linetype = 'dashed') +
  geom_point(alpha = .7, color = "grey", size = 2) +
  geom_point(data = ocd_10, aes(x = z_mean, y = -log10(pvalue), color = -log10(pvalue)), size = 3) +
  scale_color_viridis_c(name = "Top Genes",option="H", 
                        breaks = c(min(-log10(ocd_10$pvalue)), 
                                   max(-log10(ocd_10$pvalue))),
                        labels = c("10%", "0.1%"),
                        guide = guide_colorbar(position = "right",
                                               theme = theme(legend.ticks = element_blank(),
                                                             legend.title = element_blank(),
                                                             legend.text = element_text(size = 12)))) +
  geom_text_repel(data = ocd_05, 
                  aes(x  =z_mean, 
                      y = -log10(pvalue), 
                      label = gene_name),
                  box.padding = .6, direction = "b", max.overlaps = 20) +
  annotate("text", x = 3.5, y = 1.45, label = "Top 10%", 
           hjust = -0.1, vjust = 0, color = "#440154FF", size = 5, fontface = "bold") +
  xlim(-6, 6) +
  xlab("Z score") +
  ylab(expression(-log[10](p[value]))) +
  ggtitle("OCD TWAS Genes")

ggsave(paste0(figures_path,"ocd_twas.png"), width = 7, height = 6, device = "png", dpi = "retina")

scz_twas <- readr::read_tsv(paste0(twas_path,"SCZ/PGC_SCZ_SMetaXcan_en.tsv")) |>
  dplyr::filter(gene_name %in% common_genes) |>
  tidyr::drop_na(pvalue)

scz_10 <- readr::read_tsv(paste0(twas_path,"SCZ/SCZ_top_genes_thr_10.tsv"))
scz_05 <- readr::read_tsv(paste0(twas_path,"SCZ/SCZ_top_genes_thr_0.5.tsv"))
yint <- min(-log10(scz_10$pvalue))

scz_twas |>
  mutate(log10 = -log10(pvalue)) |>
  ggplot(aes(x=z_mean, y=log10)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey50", alpha = .5) +
  geom_hline(yintercept = yint, color = "#440154FF", linetype = 'dashed') +
  geom_point(alpha = .7, color = "grey", size =2) +
  geom_point(data = scz_10, aes(x = z_mean, y = -log10(pvalue), color = -log10(pvalue)), size =3) +
  scale_color_viridis_c(name = "Top Genes",option="H", 
                        breaks = c(min(-log10(scz_10$pvalue)), 
                                   max(-log10(scz_10$pvalue))),
                        labels = c("10%", "0.1%"),
                        guide = guide_colorbar(position = "right",
                                               theme = theme(legend.ticks = element_blank(),
                                                             legend.text = element_text(size = 12)))) +
  geom_text_repel(data = scz_05, 
                  aes(x = z_mean, 
                      y = -log10(pvalue), 
                      label = gene_name),
                  box.padding = .6, direction = "b", max.overlaps = 20) +
  annotate("text", x = 3.5, y = yint + 0.1, label = "Top 10%", 
           hjust = -0.1, vjust = 0, color = "#440154FF", size = 5, fontface = "bold") +
  xlim(-8, 8) +
  xlab("Z score") +
  ylab(expression(-log[10](p[value]))) +
  ggtitle("SCZ TWAS Genes")

ggsave(paste0(figures_path,"scz_twas.png"), width = 7, height = 6, device = "png", dpi = "retina")

