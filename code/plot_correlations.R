library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)
library(readxl)
library(patchwork)


theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# Plot the correlations in ggplot from the spin-permutation analysis

twas_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/"
data_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data/"

figures_path <- "/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/figures/"

enigma_data_cortical <- readxl::read_xlsx(paste0(data_path, "ENIGMA/ENIGMA_structural_organised.xlsx"),sheet = 2) |>
  mutate(Structure = str_replace(Structure, "^M_", "")) |>
  mutate(Structure = str_replace(Structure, "_thickavg$", ""))

enigma_data_subcortical <- readxl::read_xlsx(paste0(data_path, "ENIGMA/ENIGMA_structural_organised.xlsx"),sheet = 3) |>
  mutate(Structure = str_replace(Structure, "^M_", "")) |>
  mutate(Structure = str_to_lower(str_replace(Structure, "_thickavg$", "")))


# ADHD

all_corr_cort <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_cortical_correlations_results.tsv"))
all_corr_subcort <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_subcortical_correlations_results.tsv"))

# --------------------
# 10
thr <- 10
gene_data <- readr::read_tsv(paste0(twas_path, "ADHD/ADHD_TPRS_thr_",thr,".tsv"))



# Cortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_cortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_cort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_cortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = -.08, y = 47, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = -.075, y = 46.4, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()


hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = 0.36, y = 1.5, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = 0.38, y = 1.34, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("Cortical ADHD Correlation with Top ", thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))
ggsave(paste0(figures_path, "ADHD_", thr, "_cort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)

# Subcortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_subcortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_subcort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_subcortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = -.12, y = 49, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = -.116, y = 48.4, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()


hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = 0.65, y = 1.4, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = 0.67, y = 1.3, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()
hiplot

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("ADHD Subcortical Correlation with Top ",thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))

ggsave(paste0(figures_path, "ADHD_", thr, "_subcort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)

#----------
# 5
thr <- 5
gene_data <- readr::read_tsv(paste0(twas_path, "ADHD/ADHD_TPRS_thr_",thr,".tsv"))



# Cortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_cortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_cort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_cortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = -.08, y = 18.9, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = -.074, y = 18.65, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()
splot

hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = 0.28, y = 1.9, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = 0.30, y = 1.75, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()
hiplot

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("ADHD Cortical Correlation with Top ", thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))
ggsave(paste0(figures_path, "ADHD_", thr, "_cort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)

# Subcortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_subcortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_subcort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_subcortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = -.12, y = 18.5, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = -.117, y = 18.2, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()
splot

hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = 0.25, y = 1.5, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = 0.27, y = 1.38, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()
hiplot

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("ADHD Subcortical Correlation with Top ",thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))
ggsave(paste0(figures_path, "ADHD_", thr, "_subcort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)

                  
#----------
# 1

thr <- 1
gene_data <- readr::read_tsv(paste0(twas_path, "ADHD/ADHD_TPRS_thr_",thr,".tsv"))



# Cortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_cortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_cort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_cortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = -.08, y = 12, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = -.075, y = 11.7, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()
splot

hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = 0.26, y = 1.8, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = 0.28, y = 1.64, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()
hiplot

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("Cortical ADHD Correlation with Top ", thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))
ggsave(paste0(figures_path, "ADHD_", thr, "_cort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)

# Subcortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_subcortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_subcort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_subcortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = -.12, y = 11.5, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = -.116, y = 11, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()
splot

hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = -0.2, y = 1.4, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = -0.17, y = 1.3, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()
hiplot

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("ADHD Subcortical Correlation with Top ",thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))

ggsave(paste0(figures_path, "ADHD_", thr, "_subcort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)

#----------
# 0.5
thr <- 0.5
gene_data <- readr::read_tsv(paste0(twas_path, "ADHD/ADHD_TPRS_thr_",thr,".tsv"))



# Cortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_cortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_cort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_cortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = .08, y = 7.35, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = .08, y = 7.3, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()
splot

hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = -0.28, y = 1.85, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = -0.27, y = 1.75, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()
hiplot

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("Cortical ADHD Correlation with Top ", thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))
ggsave(paste0(figures_path, "ADHD_", thr, "_cort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)

# Subcortical data
corr_data <- readr::read_tsv(paste0(twas_path,"ADHD/ADHD_subcortical_corr_", thr,".tsv"), col_names = FALSE)
colnames(corr_data) <- c("Correlation")

corr_line <- all_corr_subcort |>
  filter(threshold == thr) |>
  select(corr, p)

data_scatter <- enigma_data_subcortical |>
  select(Structure, ADHD) |>
  left_join(gene_data, by = join_by(Structure == label)) |>
  filter(hemisphere == "L") |>
  select(ADHD, weighted_avg)

splot <- data_scatter |>
  ggplot(aes(x = ADHD, y = weighted_avg)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "#A8221C") +
  annotate(geom = "text", x = -.12, y = 7.4, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4)+
  annotate(geom = "text", x = -.116, y = 7.2, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  xlab("Cohen's D") +
  ylab("TPRS") +
  theme_Publication()
splot


hiplot <- corr_data |>
  ggplot(aes(x = Correlation)) +
  geom_histogram(aes(y =after_stat(density)),bins= 40, binwidth = .03, fill = "#A8221C",color="#e9ecef", alpha = .7) +
  geom_vline(xintercept = corr_line$corr, linetype = 2, lwd = 1) +
  annotate(geom = "text", x = 0.65, y = 1.4, label = paste0("r = ", round(corr_line$corr, digits = 3)), size = 4) + 
  annotate(geom = "text", x = 0.67, y = 1.3, label =  paste("p[spin] ==", corr_line$p), parse = T, size = 4) +
  ylab("Density") +
  xlab("Null Correlations") +
  theme_Publication()
hiplot

adhd_thr <- free(splot) / free(hiplot) +
  plot_annotation(title = paste0("ADHD Subcortical Correlation with Top ",thr,"% of TWAS Genes"),
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = .5)))

ggsave(paste0(figures_path, "ADHD_", thr, "_subcort_corr.png"), plot = adhd_thr, dpi = 300, width = 6, height = 8)
