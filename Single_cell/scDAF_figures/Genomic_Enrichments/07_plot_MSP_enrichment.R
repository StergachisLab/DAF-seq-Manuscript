library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(glue)
library(data.table)
library(tools)
library(ggrepel)
library(scales)


my_ggsave <- function(file, ...){
  file = glue(file)
  print(file)
  ext = file_ext(file)
  file_without_ext = tools::file_path_sans_ext(file)
  tbldir = paste0(dirname(file), "/Tables/")
  tblfile = paste0(tbldir, basename(file_without_ext), ".tbl.gz")
  ggsave(glue("{file}"), bg='transparent', ...)
  cmd = glue("cp {file} tmp.{ext}")
  fwrite(last_plot()$data, tblfile, sep="\t")
  print(cmd)
  system(cmd)
}

# TSS enrichment scores ----------------------------------------------
enrichments <- read_csv('ALL_2kb_tss_counts_FIRE_peaks_Stranded.csv')

FS_only <- enrichments %>% filter(assay == "Fiber-seq")
ATAC_only <- enrichments %>% filter(assay == "ATAC-seq")
DAF_only <- enrichments %>% filter(assay == "scDAF-seq")

DAF_median <- DAF_only %>% group_by(position) %>% summarize(med = median(enrichment))
DAF_median$sample <- "Median"

FS_line_tss_enrichment_FIRE_stranded <- ggplot(FS_only, aes(x = position, y = enrichment, group=sample, color=sample)) +
  geom_line() +
  geom_line(data = DAF_median, aes(x = position, y = med), linetype="dashed", color = "black") +
  scale_y_continuous(breaks = seq(0,30, by = 10), limits=c(0,30)) +
  theme_bw() +
  xlab("Position relative to TSSs") + ylab("TSS Enrichment Score") + ggtitle("Fiber-seq")
my_ggsave("figures/tss_enrichment_FS.pdf", FS_line_tss_enrichment_FIRE_stranded, width = 8, height = 5)

ATAC_line_tss_enrichment_FIRE_stranded <- ggplot(ATAC_only, aes(x = position, y = enrichment, group=sample, color=sample)) +
  geom_line() +
  geom_line(data = DAF_median, aes(x = position, y = med), linetype="dashed", color = "black") +
  scale_y_continuous(breaks = seq(0,30, by = 10), limits=c(0,30)) +
  theme_bw() +
  xlab("Position relative to TSSs") + ylab("TSS Enrichment Score") + ggtitle("ATAC-seq")
my_ggsave("figures/tss_enrichment_ATAC.pdf", ATAC_line_tss_enrichment_FIRE_stranded, width = 8, height = 5)

DAF_line_tss_enrichment_FIRE_stranded <- ggplot(DAF_only, aes(x = position, y = enrichment, group=sample, color=sample)) +
  geom_line() +
  geom_line(data = DAF_median, aes(x = position, y = med), linetype="dashed", color = "black") +
  scale_y_continuous(breaks = seq(0,30, by = 10), limits=c(0,30)) +
  theme_bw() +
  xlab("Position relative to TSSs") + ylab("TSS Enrichment Score") + ggtitle("scDAF-seq")
my_ggsave("figures/tss_enrichment_DAF.pdf", DAF_line_tss_enrichment_FIRE_stranded, width = 6, height = 5)


# correlation between scDAF-seq TSS enrichment and deamination rate
da_rates <- read_csv('../deamination_rate_by_cell.csv')
daf_max <- DAF_only %>% group_by(sample) %>% summarize(m = max(enrichment)) %>% rename(Cell = sample) # Max enrichment score
tss_da <- merge(da_rates, daf_max)

tss_da_corr <- ggscatter(tss_da, "prop_da_both_strands", "m", add = "reg.line") +
 xlab("Percentage of bases deaminated") + ylab("TSS Enrichment Score") +
 scale_x_continuous(labels = scales::percent, limits=c(0.15,0.3)) +
 ylim(0,20) +
 stat_cor(label.x = 0.2, label.y = 17) +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())
my_ggsave('figures/tss_enrichment_by_da_correlation.pdf', tss_da_corr, width = 4, height=4)


# Fiber-seq % m6A
fs_m6a <- read_tsv('tss_enrichment_intermediate/FS_m6A_prop.tsv', col_names = c("sample","prop"))
fs_max <- FS_only %>% group_by(sample) %>% summarize(m = max(enrichment)) # Max enrichment score
tss_fs <- merge(fs_m6a, fs_max)

tss_FS_corr <- ggscatter(tss_fs, "prop", "m", add = "reg.line") +
 xlab("Percentage of bases methylated (m6A)") + ylab("TSS Enrichment Score") +
 scale_x_continuous(labels = scales::percent, limits=c(0,0.1)) +
 ylim(0,30) +
 stat_cor(label.x = 0, label.y = 5) +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())
my_ggsave('figures/tss_enrichment_by_FS_m6A_correlation.pdf', tss_FS_corr, width = 4, height=4)

