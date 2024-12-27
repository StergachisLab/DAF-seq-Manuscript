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

# Pileups generated in pileups/ft_extract_stats_UBA1.sh
# bedtools intersect -a ft_pileup_GM12878_UBA1_H1.bed -b tss_uba1_masSeq.bed > ft_pileup_H1_uba1_masSeq_TSS.bed --> ran in bash
# bedtools intersect -a ft_pileup_GM12878_UBA1_H2.bed -b tss_uba1_masSeq.bed > ft_pileup_H2_uba1_masSeq_TSS.bed --> ran in bash

ft_names=c('chr','start','end','coverage','fire_coverage','score','nuc_coverage','msp_coverage','m6a_coverage')
h1_ft <- read_tsv('ft_pileup_H1_uba1_masSeq_TSS.bed', col_names=ft_names)
h1_ft$prop <- h1_ft$msp_coverage / h1_ft$coverage
h1_ft$hap <- "H1"
h2_ft <- read_tsv('ft_pileup_H2_uba1_masSeq_TSS.bed', col_names=ft_names)
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


# Canonical promoter TF % occupancy BY REGION ----------------------------------------------------------------------------
# bedtools intersect -wa -a region_footprint_props_UBA1_H1.bed -b canonical_uba1_promoter.bed > canonical_uba1_promoter_footprinted_regions_H1.bed
# bedtools intersect -wa -a region_footprint_props_UBA1_H2.bed -b canonical_uba1_promoter.bed > canonical_uba1_promoter_footprinted_regions_H2.bed


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



# h1z = pd.read_csv('zmw_footprint_regions_UBA1_H1.csv')
# h2z = pd.read_csv('zmw_footprint_regions_UBA1_H2.csv')

# Counter(h1z['31'])
# Counter(h2z['31'])
# len(h1z[h1z['31'].isin([0,1])])
# len(h2z[h2z['31'].isin([0,1])])
# sum(Counter(h1z['31']).values())
# sum(Counter(h2z['31']).values())
# len(h1z[h1z['31'].isin([0,1])])/sum(Counter(h1z['31']).values()) # 0.5366057902389635
# len(h2z[h2z['31'].isin([0,1])])/sum(Counter(h2z['31']).values()) # 0.46560854600849433


# Counter(h1z['32'])
# Counter(h2z['32'])
# len(h1z[h1z['32'].isin([0,1])])
# len(h2z[h2z['32'].isin([0,1])])
# sum(Counter(h1z['32']).values())
# sum(Counter(h2z['32']).values())
# len(h1z[h1z['32'].isin([0,1])])/sum(Counter(h1z['32']).values()) # 0.6622385980404845
# len(h2z[h2z['32'].isin([0,1])])/sum(Counter(h2z['32']).values()) # 0.6296147154242993




# # Codependency and co-occupancy tracks for TSS MSP overlaps ----------------------------------------------------------------------------
# # ALL data is from H1!!!
# scores <- read_csv('codep_scores_pairs_H1.csv') # data to use for codep and co-occ heatmaps!!!
# scores <- scores[order(scores$score),]
# scores$rank <- as.integer(rownames(scores))
# scores$name <- str_c(scores$reg1,":",scores$reg2)
# # nfilt = 3
# # scores_filt <- scores %>% filter(rank <= nfilt | rank > (nrow(scores)-nfilt))

# score_point <- ggplot(scores, aes(x=reorder(rank, score), y=score)) +
#     geom_point(size = 1.5, color = "red") +
#     xlab('Peak Pair') +
#     ylab('Codependency Score') +
#     theme_bw(20) +
#     theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x = element_blank()) +
#     ylim(-1,1)

# my_ggsave('figures/codep_score_point_UBA1_H1.pdf', score_point, width = 10, height = 5)






