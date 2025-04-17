library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(cowplot) # for theme_minimal_hgrid()
library(colorspace) # for darken()
library(ggpubr)
library(gplots)
library(ggside)
library(RColorBrewer)
library(glue)
library(data.table)
library(tools)
library(tidyr)
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


Colon <- read_tsv('Colon_SLC39A4_PS00682_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample.tsv')
Colon$tissue <- 'Colon_SLC39A4'
GM12878 <- read_tsv('GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample.tsv')
GM12878$tissue <- 'GM12878_SLC39A4'
Heart <- read_tsv('Heart_SLC39A4_PS00681_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample.tsv')
Heart$tissue <- 'Heart_SLC39A4'
Liver <- read_tsv('Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample.tsv')
Liver$tissue <- 'Liver_SLC39A4'
NAPA <- read_tsv('PS00626_map-pb_corrected_realigned_NAPA_region_DAF_DeDuplicated_stats_Downsample.tsv')
NAPA$tissue <- 'GM12878_NAPA'
WASF1 <- read_tsv('PS00626_map-pb_corrected_realigned_WASF1_region_DAF_DeDuplicated_stats_Downsample.tsv')
WASF1$tissue <- 'GM12878_WASF1'
UBA1 <- read_tsv('PS00685_GM12878_UBA1_map-pb_corrected_realigned_UBA1_region_DAF_DeDuplicated_stats_Downsample.tsv')
UBA1$tissue <- 'GM12878_UBA1'
COLO <- read_tsv('PS00719_COLO_Region2_map-pb_corrected_realigned_COLO_region_DAF_DeDuplicated_stats_Downsample.tsv')
COLO$tissue <- 'COLO829BLT'

merged_all <- rbind(Colon, GM12878, Heart, Liver, NAPA, WASF1, UBA1, COLO)


perc_by_tissue <- ggbarplot(merged_all, x = "tissue", y = "percent_dup", orientation = "horiz",
    fill = "tissue", facet.by = "total_reads",
    position = position_dodge(0.9)) +
    scale_y_continuous(labels = scales::comma) +
    theme(legend.position="none") +
    xlab("Sample") + ylab('Percent Duplicates') +
    scale_fill_brewer(palette="Set1")
    
my_ggsave('../figures/dup_perc_by_tissue.pdf', perc_by_tissue)


max_dup_by_tissue <- ggbarplot(merged_all, x = "tissue", y = "max_copies", orientation = "horiz",
    fill = "tissue", facet.by = "total_reads",
    position = position_dodge(0.9)) +
    scale_y_continuous(labels = scales::comma) +
    theme(legend.position="none") +
    xlab("Sample") + ylab('Max Read Count') +
    scale_fill_brewer(palette="Set1")

my_ggsave('../figures/max_dup_by_tissue.pdf', max_dup_by_tissue)



# Final Bar plots (25k reads, all reads, number of reads per sample)

# Dup stats from ALL reads
Colon_ALL <- read_tsv('Colon_SLC39A4_PS00682_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats.tsv')
Colon_ALL$tissue <- 'Colon_SLC39A4'
GM12878_ALL <- read_tsv('GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats.tsv')
GM12878_ALL$tissue <- 'GM12878_SLC39A4'
Heart_ALL <- read_tsv('Heart_SLC39A4_PS00681_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats.tsv')
Heart_ALL$tissue <- 'Heart_SLC39A4'
Liver_ALL <- read_tsv('Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats.tsv')
Liver_ALL$tissue <- 'Liver_SLC39A4'
NAPA_ALL <- read_tsv('PS00626_map-pb_corrected_realigned_NAPA_region_DAF_DeDuplicated_stats.tsv')
NAPA_ALL$tissue <- 'GM12878_NAPA'
WASF1_ALL <- read_tsv('PS00626_map-pb_corrected_realigned_WASF1_region_DAF_DeDuplicated_stats.tsv')
WASF1_ALL$tissue <- 'GM12878_WASF1'
UBA1_ALL <- read_tsv('PS00685_GM12878_UBA1_map-pb_corrected_realigned_UBA1_region_DAF_DeDuplicated_stats.tsv')
UBA1_ALL$tissue <- 'GM12878_UBA1'
COLO_ALL <- read_tsv('PS00719_COLO_Region2_map-pb_corrected_realigned_COLO_region_DAF_DeDuplicated_stats.tsv')
COLO_ALL$tissue <- 'COLO829BLT'

merged_all2 <- rbind(Colon_ALL, GM12878_ALL, Heart_ALL, Liver_ALL, NAPA_ALL, WASF1_ALL, UBA1_ALL, COLO_ALL)

# 25,000 read downsample
down_25k <- ggbarplot(merged_all %>% filter(total_reads == 25000), x = "tissue", y = "percent_dup", orientation = "horiz",
    fill = "tissue",
    position = position_dodge(0.9)) +
    ylim(0,30) +
    theme(legend.position="none") +
    xlab("Sample") + ylab('Percent Duplicates') +
    scale_fill_brewer(palette="Set1")

# All read percent duplicates
perc_dub_all <- ggbarplot(merged_all2, x = "tissue", y = "percent_dup", orientation = "horiz",
    fill = "tissue",
    position = position_dodge(0.9)) +
    ylim(0,30) +
    theme(legend.position="none") +
    xlab("Sample") + ylab('Percent Duplicates') +
    scale_fill_brewer(palette="Set1")

# Total read count by sample
n_reads <- ggbarplot(merged_all2, x = "tissue", y = "total_reads", orientation = "horiz",
    fill = "tissue",
    position = position_dodge(0.9)) +
    ylim(0,1500000) +
    theme(legend.position="none") +
    xlab("Sample") + ylab('Total Reads') +
    scale_fill_brewer(palette="Set1")


my_ggsave('../figures/percent_dup_25k.pdf', down_25k)
my_ggsave('../figures/percent_dup_all.pdf', perc_dub_all)
my_ggsave('../figures/read_count.pdf', n_reads)


# line plot of sequencing depth on X-axis and duplication rate on y-axis 

Colon <- read_tsv('Colon_SLC39A4_PS00682_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
Colon$tissue <- 'Colon_SLC39A4'
GM12878 <- read_tsv('GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
GM12878$tissue <- 'GM12878_SLC39A4'
Heart <- read_tsv('Heart_SLC39A4_PS00681_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
Heart$tissue <- 'Heart_SLC39A4'
Liver <- read_tsv('Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
Liver$tissue <- 'Liver_SLC39A4'
NAPA <- read_tsv('PS00626_map-pb_corrected_realigned_NAPA_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
NAPA$tissue <- 'GM12878_NAPA'
WASF1 <- read_tsv('PS00626_map-pb_corrected_realigned_WASF1_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
WASF1$tissue <- 'GM12878_WASF1'
UBA1 <- read_tsv('PS00685_GM12878_UBA1_map-pb_corrected_realigned_UBA1_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
UBA1$tissue <- 'GM12878_UBA1'
COLO <- read_tsv('PS00719_COLO_Region2_map-pb_corrected_realigned_COLO_region_DAF_DeDuplicated_stats_Downsample_highRange.tsv')
COLO$tissue <- 'COLO829BLT'

merged_all <- rbind(Colon, GM12878, Heart, Liver, NAPA, WASF1, UBA1, COLO)


# merged_all$prop <- merged_all$percent_dup / 100

downsample_line <- ggline(merged_all %>% filter(total_reads >= 4000), "total_reads", "percent_dup",
   linetype = "tissue",
   color = "tissue") +
   scale_y_continuous(trans='log2') +
   scale_color_brewer(palette="Set1")

my_ggsave('../figures/percent_dup_downsample_line_log2.pdf', downsample_line, width = 9, height = 8, units = "in")


downsample_line <- ggline(merged_all %>% filter(total_reads >= 4000), "total_reads", "percent_dup",
   linetype = "tissue",
   color = "tissue") +
   scale_color_brewer(palette="Set1")

my_ggsave('../figures/percent_dup_downsample_line.pdf', downsample_line, width = 9, height = 8, units = "in")

