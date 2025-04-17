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

# TSS % MSP overlaps ----------------------------------------------------------------------------

tss_df <- read_tsv('tss_uba1_masSeq.bed', col_names=c('chr','start','end','PB_ID','score','strand')) # TSS coordinates from MAS-seq

ft_names=c('chr','start','end','coverage','fire_coverage','score','nuc_coverage','msp_coverage','m6a_coverage')
h1_ft <- read_tsv('pileups/ft_pileup_H1_uba1_masSeq_TSS.bed', col_names=ft_names)
h1_ft$prop <- h1_ft$msp_coverage / h1_ft$coverage
h1_ft$hap <- "H1"
h2_ft <- read_tsv('pileups/ft_pileup_H2_uba1_masSeq_TSS.bed', col_names=ft_names)
h2_ft$prop <- h2_ft$msp_coverage / h2_ft$coverage
h2_ft$hap <- "H2"

both_haps <- rbind(h1_ft, h2_ft)

# tss_df$TSS_ID <- c(1:4)
both_haps$TSS_ID <- as.factor(both_haps$start)

tss_perc_acc <- ggbarplot(both_haps, "TSS_ID", "prop",
  fill = "hap", color = "hap", palette = "Paired",
  position = position_dodge(0.9)) +
  xlab("TSS Coordinate") +
  ylab("Percent Actuation") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic(20)

my_ggsave('figures/tss_perc_acc_UBA1.pdf', tss_perc_acc, width = 8, height = 6, units = "in")



# Canonical promoter TF % occupancy ----------------------------------------------------------------------------
# bedtools intersect -wa -a UBA1_GM12878_footprinted_motifs_H1_sorted.bed -b canonical_uba1_promoter.bed > canonical_uba1_promoter_footprinted_motifs_H1.bed
footprint_stats = read_csv('footprinting/footprint_stats.csv')
filt_motifs <- read_tsv('footprinting/canonical_uba1_promoter_footprinted_motifs_H1.bed', col_names=c('chr','start','end','coverage','name','prop_H1','strand'))
filt_coords <- paste0('chrX:', filt_motifs$start, '-', filt_motifs$end)
h1_filt <- footprint_stats %>% filter(motif %in% filt_coords)
h1_filt[h1_filt == 'CTCF_Mod2'] <- 'CTCF_L'

TF_occ_long <- gather(h1_filt, hap, prop, c(H1_prop,H2_prop), factor_key=TRUE)
levels(TF_occ_long$hap) <- c('H1', 'H2')
name_split <- strsplit(TF_occ_long$name, split=".", fixed=TRUE)
TF_occ_long$name_short <- sapply(name_split, tail, 1)

# TF_occ_long %>% select(motif, name, name_short, prop, hap)

# Sort TFs by coordinate
coord_split <- strsplit(TF_occ_long$motif, split='[:-]+', fixed=FALSE)
starts <- c()
for (i in 1:length(coord_split)){
  starts <- c(starts, coord_split[[i]][2])
}
TF_occ_long$start <- as.numeric(starts)
TF_occ_long <- TF_occ_long[order(TF_occ_long$start),]

TF_perc_occ <- ggbarplot(TF_occ_long, "name_short", "prop",
  fill = "hap", color = "hap", palette = "Paired",
  position = position_dodge(0.9)) +
  xlab("TF Motif") +
  ylab("Percent Occupancy") +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
  theme_classic(20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

my_ggsave('figures/TF_perc_occ_UBA1_canonical.pdf', TF_perc_occ, width = 12, height = 6, units = "in")


h1_regions <- read_tsv('footprinting/canonical_uba1_promoter_footprinted_regions_H1.bed', col_names=c('chr','start','end','region','prop','strand'))
h2_regions <- read_tsv('footprinting/canonical_uba1_promoter_footprinted_regions_H2.bed', col_names=c('chr','start','end','region','prop','strand'))
h1_regions$hap <- 'H1'
h2_regions$hap <- 'H2'

reg_both_haps <- rbind(h1_regions,h2_regions)

region_perc_occ <- ggbarplot(reg_both_haps, "start", "prop",
  fill = "hap", color = "hap", palette = "Paired",
  position = position_dodge(0.9)) +
  xlab("TF Region") +
  ylab("Percent Occupancy") +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
  theme_classic(20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

my_ggsave('figures/region_perc_occ_UBA1_canonical.pdf', region_perc_occ, width = 12, height = 6, units = "in")


