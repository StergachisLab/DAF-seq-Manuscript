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
window_len = 2000
flank_len = 100

DAF_counts_stranded <- read_csv("DAF_2kb_tss_counts_FIRE_peaks_Stranded.csv")
DAF_flanking_mean_df_stranded <- DAF_counts_stranded %>% filter(position < (-window_len+flank_len) | position > (window_len-flank_len))
DAF_flanking_mean_stranded <- mean(DAF_flanking_mean_df_stranded$count)
DAF_counts_stranded$enrichment <- DAF_counts_stranded$count / DAF_flanking_mean_stranded

FS_counts_stranded <- read_csv("FS_2kb_tss_counts_FIRE_peaks_Stranded.csv")
FS_flanking_mean_df_stranded <- FS_counts_stranded %>% filter(position < (-window_len+flank_len) | position > (window_len-flank_len))
FS_flanking_mean_stranded <- mean(FS_flanking_mean_df_stranded$count)
FS_counts_stranded$enrichment <- FS_counts_stranded$count / FS_flanking_mean_stranded

DAF_line_tss_enrichment_FIRE_stranded <- ggplot(DAF_counts_stranded, aes(x = position, y = enrichment)) +
  geom_line() +
  scale_y_continuous(breaks = seq(0,12, by = 3), limits=c(0,12)) +
  theme_bw() +
  xlab("Position relative to TSSs") + ylab("TSS Enrichment Score") + ggtitle("scDAF-seq (all cells)")

FS_line_tss_enrichment_FIRE_stranded <- ggplot(FS_counts_stranded, aes(x = position, y = enrichment)) +
  geom_line() +
  scale_y_continuous(breaks = seq(0,30, by = 10), limits=c(0,30)) +
  theme_bw() +
  xlab("Position relative to TSSs") + ylab("TSS Enrichment Score") + ggtitle("Fiber-seq")

tss_grid <- plot_grid(DAF_line_tss_enrichment_FIRE_stranded, FS_line_tss_enrichment_FIRE_stranded, nrow=1)
ggsave("tss_enrichment_FIRE_stranded_DAF_vs_FS.pdf", tss_grid, width = 12, height = 6)


# FIRE enrichment scores -------------------------------------------------------------
DAF_counts_FIRE <- read_csv("DAF_MSPs_2kb_counts_binned_FIRE_peaks.csv")
DAF_counts_FIRE$bin_prop <- DAF_counts_FIRE$bin * 0.1
DAF_counts_FIRE$bin <- as.factor(DAF_counts_FIRE$bin)
DAF_flanking_mean_df_FIRE <- DAF_counts_FIRE %>% filter(position < (-window_len+flank_len) | position > (window_len-flank_len)) %>% group_by(bin) %>% summarize(avg = mean(count))
DAF_counts_FIRE$enrichment <- 0
for (i in 1:nrow(DAF_counts_FIRE)){
  row_bin <- DAF_counts_FIRE[i,]$bin
  DAF_counts_FIRE[i,]$enrichment <- DAF_counts_FIRE[i,]$count / DAF_flanking_mean_df_FIRE[DAF_flanking_mean_df_FIRE$bin == row_bin,]$avg
}

FS_counts_FIRE <- read_csv("FS_MSPs_2kb_counts_binned_FIRE_peaks.csv")
FS_counts_FIRE$bin_prop <- FS_counts_FIRE$bin * 0.1
FS_counts_FIRE$bin <- as.factor(FS_counts_FIRE$bin)
FS_flanking_mean_df_FIRE <- FS_counts_FIRE %>% filter(position < (-window_len+flank_len) | position > (window_len-flank_len)) %>% group_by(bin) %>% summarize(avg = mean(count))
FS_counts_FIRE$enrichment <- 0
for (i in 1:nrow(FS_counts_FIRE)){
  row_bin <- FS_counts_FIRE[i,]$bin
  FS_counts_FIRE[i,]$enrichment <- FS_counts_FIRE[i,]$count / FS_flanking_mean_df_FIRE[FS_flanking_mean_df_FIRE$bin == row_bin,]$avg
}


labels=seq(.1,.9,.1)
DAF_line_FIRE_enrichment <- ggplot(DAF_counts_FIRE, aes(x = position, y = enrichment, color = bin_prop-0.5, fill=NULL, group=bin_prop)) +
  geom_line() +
  scale_colour_steps2("", breaks=labels-0.5, labels=labels, low = muted("blue"), mid = "white", high = muted("red")) +
  scale_y_continuous(breaks = seq(0,15, by = 5), limits=c(0,15)) +
  theme_bw() +
  xlab("Position relative to FIRE Peak") + ylab("MSP Enrichment Score") + ggtitle("scDAF-seq (all cells)")

FS_line_FIRE_enrichment <- ggplot(FS_counts_FIRE, aes(x = position, y = enrichment, color = bin_prop-0.5, fill=NULL, group=bin_prop)) +
  geom_line() +
  scale_colour_steps2("", breaks=labels-0.5, labels=labels, low = muted("blue"), mid = "white", high = muted("red")) +
  scale_y_continuous(breaks = seq(0,45, by = 15), limits=c(0,45)) +
  theme_bw() +
  xlab("Position relative to FIRE Peak") + ylab("MSP Enrichment Score") + ggtitle("Fiber-seq")

FIRE_grid <- plot_grid(DAF_line_FIRE_enrichment, FS_line_FIRE_enrichment, nrow=1)
ggsave("FIRE_enrichment_DAF_vs_FS.pdf", FIRE_grid, width = 12, height = 6)



# MSP enrichment in Repeat Masker classes ------------------------------------------------------------------------
DAF_MSP_enrichments <- read_csv("rm_grouped_MSP_counts_DAF.csv")
FS_MSP_enrichments <- read_csv("rm_grouped_MSP_counts_FS.csv")

DAF_MSP_enrichments$enrichment <- DAF_MSP_enrichments$n_msp / DAF_MSP_enrichments$rm_bases
DAF_MSP_enrichments$scaled <- DAF_MSP_enrichments$enrichment / (DAF_MSP_enrichments %>% filter(repeat_class == "ALL"))$enrichment
FS_MSP_enrichments$enrichment <- FS_MSP_enrichments$n_msp / FS_MSP_enrichments$rm_bases
FS_MSP_enrichments$scaled <- FS_MSP_enrichments$enrichment / (FS_MSP_enrichments %>% filter(repeat_class == "ALL"))$enrichment


# Merged data
DAF_MSP_enrichments$assay <- "DAF-seq"
FS_MSP_enrichments$assay <- "Fiber-seq"

Merged_MSP_enrichment <- rbind(DAF_MSP_enrichments, FS_MSP_enrichments)

# filter out "?" repeat classes
Merged_MSP_enrichment_filt <- Merged_MSP_enrichment[!grepl("\\?", Merged_MSP_enrichment$repeat_class), ]
Merged_MSP_enrichment_filt$repeat_class <- factor(Merged_MSP_enrichment_filt$repeat_class, levels = c("DNA","LINE","LTR","RC","RNA","Retroposon","SINE","Satellite","Simple_repeat","Unknown","rRNA","scRNA","snRNA","tRNA","non_RM","FIRE_peaks"))

msp_Merged_rm <- ggbarplot(Merged_MSP_enrichment_filt %>% filter(repeat_class != "ALL"), "repeat_class", "scaled", orientation = "horiz", fill = "assay", color = "assay", position = position_dodge(0.9), palette = c("#00AFBB", "#E7B800")) +
                geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
                xlab("Repeat Class") + ylab("Enrichment over All Regions") +
                theme_bw(15)

my_ggsave("MSP_enrichment_RM_MERGED.pdf", msp_Merged_rm, height = 6, width = 10)

