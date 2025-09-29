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
library(ggbeeswarm)


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


stat_df <- read_tsv("collapsed_read_stats.tsv")
stat_df$Gb_seq <- stat_df$bases / 1000000000
stat_df$longKb <- stat_df$longest / 1000
stat_df$UL_Mb <- stat_df$ultra_long_bases / 1000000
stat_df$cov_prop_all <- stat_df$ALL / stat_df$mappable_ref_bp
stat_df$cov_prop_H1 <- stat_df$H1 / stat_df$mappable_ref_bp
stat_df$cov_prop_H2 <- stat_df$H2 / stat_df$mappable_ref_bp

# Swarmplot --- Sequenced bases by Cell
stat_df$swarm_factor <- 'ALL'
stat_df$swarm_factor <- factor(stat_df$swarm_factor)

# Subset example cells for Figure 6
include_cells <- c("PS00756","PS00757","PS00758","PS00867")
stat_df_subset <- stat_df %>% filter(Cell %in% include_cells)


# Barplot --- Sequenced bases by Cell
pSeqBases <- ggbarplot(stat_df_subset, "Cell", "Gb_seq",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out", lab.nb.digits = 1) +
 ylim(0,150) +
 ylab("HiFi bases (Gb)") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/sequenced_bases_barplot.pdf', pSeqBases, width = 8, height=4)

# Barplot --- N50 by Cell
pN50 <- ggbarplot(stat_df_subset, "Cell", "N50",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out") +
 ylim(0,40000) +
 ylab("Collapsed Read N50") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/N50_barplot.pdf', pN50, width = 8, height=4)

# Barplot --- Number of Ultra Long reads per by Cell
pUltraLongCount <- ggbarplot(stat_df_subset, "Cell", "num_100Kb",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out") +
 ylim(0,6000) +
 ylab("Ultra Long (100Kb) Read Count") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/ultraLong_Count_consensus_barplot.pdf', pUltraLongCount, width = 8, height=4)

pUltraLongCount_swarm <- ggplot(stat_df, aes(x=swarm_factor, y=num_100Kb, size = 5)) +
  geom_beeswarm(cex = 7) +
  ylim(0,6000) +
  geom_label_repel(data = stat_df_subset, aes(label=Cell), box.padding = 1,  segment.color = "red", nudge_y = 1) +
  ylab("Ultra Long (100Kb) Read Count") +
  theme_bw(10) +
  theme(legend.position="none")

my_ggsave('figures/ultraLong_Count_consensus_swarm.pdf', pUltraLongCount_swarm)


# Barplot --- Bases in Ultra Long reads per by Cell
pUltraLongBases <- ggbarplot(stat_df_subset, "Cell", "UL_Mb",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out", lab.nb.digits = 1) +
 ylim(0,800) +
 ylab("Mb of sequence in Ultra Long (100Kb)") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/ultraLong_nBases_consensus_barplot.pdf', pUltraLongBases, width = 8, height=4)

pUltraLongBases_swarm <- ggplot(stat_df, aes(x=swarm_factor, y=UL_Mb, size = 5)) +
  geom_beeswarm(cex = 7) +
  ylim(0,800) +
  geom_label_repel(data = stat_df_subset, aes(label=Cell), box.padding = 1,  segment.color = "red", nudge_y = 1) +
  ylab("Mb of sequence in Ultra Long (100Kb)") +
  theme_bw(10) +
  theme(legend.position="none")

my_ggsave('figures/ultraLong_nBases_consensus_swarm.pdf', pUltraLongBases_swarm)


# Barplot --- % bp of mappable genome covered (HG38, low cov regions) --- ALL, H1, H2
cov_df_all <- stat_df %>% select(Cell, cov_prop_all) %>% rename("cov" = cov_prop_all)
cov_df_all$hap <- "All"
cov_df_H1 <- stat_df %>% select(Cell, cov_prop_H1) %>% rename("cov" = cov_prop_H1)
cov_df_H1$hap <- "H1"
cov_df_H2 <- stat_df %>% select(Cell, cov_prop_H2) %>% rename("cov" = cov_prop_H2)
cov_df_H2$hap <- "H2"
cov_merged <- rbind(cov_df_all, cov_df_H1, cov_df_H2)
cov_merged$Gb_coverged <- (cov_merged$cov*stat_df$mappable_ref_bp)/1000000000

cov_merged_subset <- cov_merged %>% filter(Cell %in% include_cells)

# percent covered
pCoverage <- ggbarplot(cov_merged_subset, "hap", "cov",
 fill = "Cell", color = "Cell",
 position = position_dodge(0.9)) +
 scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
 ylab("Mappable hg38 - Percent Coverage") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/percent_coverage_barplot.pdf', pCoverage)

# Gb covered
pCoverageGb <- ggbarplot(cov_merged_subset, "Cell", "Gb_coverged",
 fill = "hap", color = "hap",
 position = position_dodge(0.9)) +
 ylab("Mappable hg38 - Covered Bases (Gb)") +
 theme_bw(10) +
 ylim(0,3) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/bases_coverage_barplot.pdf', pCoverageGb)


# Beeswarm --------------------
pCoverage_swarm <- ggplot(cov_merged, aes(x=hap, y=cov, colour=hap, size = 10)) +
  geom_beeswarm(cex = 3) +
  scale_color_manual(values=c("#CC79A7","#00AFBB", "#FC4E07")) +
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  geom_label_repel(data = cov_merged_subset, aes(label=Cell), box.padding = 1,  segment.color = "black") +
  ylab("Mappable hg38 - Percent Coverage") +
  theme_bw(10) +
  theme(legend.position="none")

my_ggsave('figures/percent_coverage_swarmplot.pdf', pCoverage_swarm)


pSeqBases_swarm <- ggplot(stat_df, aes(x=swarm_factor, y=Gb_seq, size = 20), color = Cell) +
  geom_beeswarm(cex = 3.5, color = "steelblue") +
  geom_label_repel(data = stat_df %>% filter(Cell %in% include_cells), aes(label=Cell), box.padding = 2,  segment.color = "black", nudge_x = 0.1, nudge_y = -1) +
  ylim(0,150) +
  ylab("HiFi bases (Gb)") +
  theme_bw(15) +
  theme(legend.position="none")

my_ggsave('figures/pSeqBases_swarm.pdf', pSeqBases_swarm)


# Sequencing depth of the single cell vs % of genome covered
depth_vs_coverage_corr <- ggscatter(stat_df, x = "Gb_seq", y = "cov_prop_all") +
            xlim(0,150) +
            scale_y_continuous(labels = scales::percent, limits=c(0,1))
my_ggsave('figures/depth_vs_coverage_corr.pdf', depth_vs_coverage_corr)

dvc_all <- stat_df %>% select(Cell, cov_prop_all, Gb_seq) %>% rename("prop" = cov_prop_all)
dvc_all$type <- "All"
dvc_H1 <- stat_df %>% select(Cell, cov_prop_H1, Gb_seq) %>% rename("prop" = cov_prop_H1)
dvc_H1$type <- "H1"
dvc_H2 <- stat_df %>% select(Cell, cov_prop_H2, Gb_seq) %>% rename("prop" = cov_prop_H2)
dvc_H2$type <- "H2"
dvc_df <- rbind(dvc_all, dvc_H1, dvc_H2)

# By hap
depth_vs_coverage_corr_by_hap <- ggscatter(dvc_df, x = "Gb_seq", y = "prop", color = "type", palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
            xlim(0,150) +
            scale_y_continuous(labels = scales::percent, limits=c(0,1))
my_ggsave('figures/depth_vs_coverage_corr_by_hap.pdf', depth_vs_coverage_corr_by_hap)

# By ALL vs Phased
dvc_phased_all <- stat_df %>% select(Cell, ALL, Gb_seq, mappable_ref_bp) %>% rename("bp" = ALL)
dvc_phased_all$prop <- dvc_phased_all$bp / dvc_phased_all$mappable_ref_bp
dvc_phased_all$type <- "All"
dvc_phased_H1 <- stat_df %>% select(Cell, H1, Gb_seq, mappable_ref_bp) %>% rename("bp" = H1)
dvc_phased_H2 <- stat_df %>% select(Cell, H2, Gb_seq, mappable_ref_bp) %>% rename("bp" = H2)
all_phased <- data.frame(Cell = dvc_phased_H1$Cell, bp = dvc_phased_H1$bp + dvc_phased_H2$bp, Gb_seq = dvc_phased_H1$Gb_seq, mappable_ref_bp = dvc_phased_H1$mappable_ref_bp)
all_phased$prop <- all_phased$bp / (all_phased$mappable_ref_bp*2)
all_phased$type <- "Phased"
dvc_phased_merged <- rbind(dvc_phased_all, all_phased)

depth_vs_coverage_corr_phased <- ggscatter(dvc_phased_merged, x = "Gb_seq", y = "prop", color = "type", palette = c("#00AFBB", "#FC4E07")) +
            xlim(0,150) +
            scale_y_continuous(labels = scales::percent, limits=c(0,1))
my_ggsave('figures/depth_vs_coverage_corr_phased.pdf', depth_vs_coverage_corr_phased)


# SAVE COVERAGE STATS AS A TABLE
write.csv(cov_merged, "figures/coverage_proportion_by_cell.csv", row.names = FALSE)


# Consensus read phasing percentages
perc_phased <- c()
for (i in 1:nrow(stat_df)){
  phasing <- read.table(paste0("../phasing/haplotagged/", stat_df$Cell[i] ,"_consensus_BothStrands_HG38_corrected.haplotag_list.txt.gz"))
  n_phased <- nrow(phasing %>% filter(V2 != "none"))
  perc_phased <- c(perc_phased, (n_phased/nrow(phasing))*100)
}
perc_phased_df <- data.frame(Cell = stat_df$Cell, percent_phased = perc_phased)
write_tsv(perc_phased_df, "figures/percent_phased_consensus.tsv")


# N50 of raw reads
raw_read_lengths <- read_tsv('raw_read_lengths.tsv.gz', col_names=c('length','Cell'))

calculate_N50 <- function(lengths){
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  cumulative_lengths <- cumsum(sorted_lengths)
  half_total_length <- sum(lengths) / 2
  N50 <- sorted_lengths[cumulative_lengths >= half_total_length][1]
  return(N50)
}
N50_all_cells <- calculate_N50(raw_read_lengths$length)

cell_names <- unique(raw_read_lengths$Cell)
n50s <- c()
for (i in 1:length(cell_names)){
  temp_df <- raw_read_lengths %>% filter(Cell == cell_names[i])
  n50s <- c(n50s, calculate_N50(temp_df$length))
}

# Save Raw Read N50 data
raw_N50_df <- data.frame(N50_all = N50_all_cells)
for (i in 1:length(cell_names)){
  cname <- paste0(cell_names[i],"_N50")
  raw_N50_df[[cname]] <- n50s[i]
}
write.csv(raw_N50_df, "figures/raw_read_N50.csv", row.names = FALSE, quote=FALSE)


consensus_lengths <- read_tsv('consensus_read_lengths.tsv', col_names=c('length','Cell'))

# Subset example cells for Figure 6
raw_read_lengths_subset <- raw_read_lengths %>% filter(Cell %in% include_cells)
consensus_lengths_subset <- consensus_lengths %>% filter(Cell %in% include_cells)


# Read length Density Plots weighted by total bp
raw_read_length_density_weight <- ggplot(raw_read_lengths_subset, aes(x=length, color=Cell)) +
                            geom_density(alpha=0.6, aes(weight=length)) +
                            scale_x_continuous(trans = 'log2', breaks = scales::breaks_log(2, n = 7), limits = c(1, 524288)) +
                            scale_color_manual(values = c("#00A14B","#F26522","#00AEEF","#7F3F98")) +
                            theme_bw(10)
my_ggsave('figures/raw_read_length_density_weighted.pdf', raw_read_length_density_weight, width = 8, height=4)

collapsed_read_length_density_weight <- ggplot(consensus_lengths_subset, aes(x=length, color=Cell)) +
                            geom_density(alpha=0.6, aes(weight=length)) +
                            scale_x_continuous(trans = 'log2', breaks = scales::breaks_log(2, n = 6), limits = c(1, 524288)) +
                            scale_color_manual(values = c("#00A14B","#F26522","#00AEEF","#7F3F98")) +
                            theme_bw(10)
my_ggsave('figures/collapsed_read_length_density_weighted.pdf', collapsed_read_length_density_weight, width = 8, height=4)



# Karyoplots -----------------------------------------------------------------------------------------------------------

# Karyoplot coverage --- chr2 uncollapsed MAPPABLE reads
library(karyoploteR)
library(GenomicAlignments)


# COLLAPSED Reads filtered by MAPPABLE regoins AND High Coverage Regions

# Raw clipped reads FILTERED ---------------------------------------------
raw_clipped_BothStrands_FILT_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00756_clipped_FILTERED.bed.gz', col_names = c('chrom','start','end','name'))
gr_PS00756_clipped_FILT_BothStrands <- makeGRangesFromDataFrame(raw_clipped_BothStrands_FILT_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL,
                                       seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)

# Raw Reads
png("figures/PS00756_raw_reads_BothStrandsMerged_clipped_FILT_karyo_chr13.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr13")
kp <- kpPlotCoverage(kp, data=gr_PS00756_clipped_FILT_BothStrands, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()
# Save AXIS LABELS as a PDF
clip_filt_max <- kp$latest.plot$computed.values$max.coverage
pdf("figures/PS00756_raw_reads_BothStrandsMerged_clipped_FILT_karyo_chr13_AXIS.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr13")
kpAxis(kp, ymax=clip_filt_max, data.panel = 1, numticks=6)
dev.off()

pdf("figures/PS00756_raw_reads_BothStrandsMerged_clipped_FILT_karyo_chr13.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr13")
kp <- kpPlotCoverage(kp, data=gr_PS00756_clipped_FILT_BothStrands, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()


collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00756_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt6 <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL,
                                       seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long6 <- granges_collapsed_filt6[width(granges_collapsed_filt6) >= 10000]

png("figures/PS00756_BothStrandsMerged_Mappable_10Kb_karyo_chr13.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr13")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long6, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()
# Save AXIS LABELS as a PDF
collapsed_filt_max <- kp$latest.plot$computed.values$max.coverage
pdf("figures/PS00756_BothStrandsMerged_Mappable_10Kb_karyo_chr13_AXIS.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr13")
kpAxis(kp, ymax=collapsed_filt_max, data.panel = 1, numticks=6)
dev.off()


pdf("figures/PS00756_BothStrandsMerged_Mappable_10Kb_karyo_chr13.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr13")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long6, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()


# Extended Data Figure - Karyoplots for all chromosomes and cells ----------------------------------------------
png("figures/PS00756_BothStrandsMerged_Mappable_10Kb_karyo_AllChroms.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2)
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long6, data.panel = 1)
dev.off()

# chr1 from each Cell
collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00718_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00718_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00756_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00756_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00757_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00757_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00758_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00758_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00867_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00867_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00868_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00868_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00869_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00869_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=5)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00870_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00870_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00871_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00871_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00872_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00872_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=6)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00873_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00873_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=5)
dev.off()

collapsed_filt_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00874_no_high_STRAND_coverage_subtract.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt <- makeGRangesFromDataFrame(collapsed_filt_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL, seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_long_collapsed_long <- granges_collapsed_filt[width(granges_collapsed_filt) >= 10000]
png("figures/PS00874_BothStrandsMerged_Mappable_10Kb_karyo_chr1.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr1")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_long, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=5)
dev.off()




# MSP Actuation --------------------------------------------------------------------

# Actuation at FIRE peaks binned by Fiber-seq actuation at different MSP lengths
msp_150 <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP150.tsv')

# Filter each peakset for peaks with width +/- 25bp to each MSP cutoff (ex. 175-225bp wide for MSP 200 cutoff)
msp_150$peak_width <- msp_150$end - msp_150$start
msp_150 <- msp_150 %>% filter(peak_width >= 125, peak_width <= 175)

# for each bin: count the number of haps x peaks in that bin covered & the total accessible haps
PS00718_bins <- msp_150 %>% filter(PS00718_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00718_nCov), nMSP = sum(PS00718_nMSP))
PS00756_bins <- msp_150 %>% filter(PS00756_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00756_nCov), nMSP = sum(PS00756_nMSP))
PS00757_bins <- msp_150 %>% filter(PS00757_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00757_nCov), nMSP = sum(PS00757_nMSP))
PS00758_bins <- msp_150 %>% filter(PS00758_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00758_nCov), nMSP = sum(PS00758_nMSP))
msp150_bins <- data.frame(bin = PS00718_bins$bin)
msp150_bins$nCov <- PS00718_bins$nCov + PS00756_bins$nCov + PS00757_bins$nCov + PS00758_bins$nCov
msp150_bins$nMSP <- PS00718_bins$nMSP + PS00756_bins$nMSP + PS00757_bins$nCov + PS00758_bins$nMSP
msp150_bins$prop <- msp150_bins$nMSP / msp150_bins$nCov


# Autocorrelation -------------------------------------------------------------------------------------------------------------------------

acf_merged <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(acf_merged) = c('acf','lag','corr','cell')
for (cn in cell_names){
  acf_df <- read_tsv(paste0('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/',cn,'_ft_acf.tsv'), col_names=c('acf','lag','corr'))
  acf_df$cell <- cn
  acf_merged <- rbind(acf_merged, acf_df)
}

acf_merged$change <- c(0, (acf_merged$corr[1:(length(acf_merged$corr)-1)] * acf_merged$corr[2:length(acf_merged$corr)]) <= 0)

acf1 <- ggplot(acf_merged %>% filter(lag > 25), aes(x=lag, y=corr, color=cell)) +
  geom_hline(aes(yintercept=0), color="darkblue", linewidth=1, linetype="dashed") + 
  geom_vline(data=NULL, aes(xintercept=147), linetype="dashed", color="black", alpha=0.5) +
  geom_line() +
  # scale_color_brewer(palette = "Set3") +
  geom_label_repel(data = . %>% filter(change == 1, (lag > 130 & lag < 200)),
    aes(y=0, x = lag, label=lag),
    min.segment.length = 0, # draw all line segments
    nudge_x=-5,
    nudge_y=0.01,
    show.legend = FALSE)+
  theme_minimal_grid() + 
  scale_y_continuous("Autocorrelation between deamination events") + 
  scale_x_continuous("Lag between deamination events") + 
  guides(color = guide_legend(override.aes = list(size = 2, shape="")))+
  theme(legend.position = "top", legend.text=element_text(size=4))
my_ggsave('figures/auto_corr_line_allCells.pdf', acf1)

# Deamination rate per cell
da_rates <- read_csv('deamination_rate_by_cell.csv')
da_rates$rate_Kb <- (da_rates$DA_bp/da_rates$Tot_bp)*1000

pDaRate_prop_select <- ggbarplot(da_rates %>% filter(Cell %in% include_cells), "Cell", "prop_da_both_strands", fill = "Cell", color = FALSE) +
 ylab("Deaminations per Kb") +
 scale_y_continuous(labels = scales::percent, limits=c(0,0.4)) +
 scale_fill_manual(values = c("#00A14B","#F26522","#00AEEF","#7F3F98")) +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())
my_ggsave('figures/deamination_prop_select_cell_barplot.pdf', pDaRate_prop_select, width = 4, height=4)

da_rates$x <- "ALL"
pDaRate_prop_violin <- ggviolin(da_rates, "x", "prop_da_both_strands", fill = "#E7B800", add = "boxplot", add.params = list(fill = "white")) +
 scale_y_continuous(labels = scales::percent, limits=c(0,0.3)) +
 ylab("Deaminations per Kb") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())
my_ggsave('figures/deamination_prop_by_cell_violin.pdf', pDaRate_prop_violin, width = 4, height=4)

write.csv(da_rates, "figures/deaminations_per_1Kb_by_cell.csv", row.names = FALSE)


# Correlate mapping % by deamination rate
perc_mapping <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/data/percent_mapped_aggregated.tsv')
da_rates$Percent_Mapped <- perc_mapping$Percent_Mapped[match(da_rates$Cell, perc_mapping$Sample)]
da_rates$exact_prop_da <- da_rates$DA_bp/da_rates$Tot_bp

da_corr <- ggscatter(da_rates, x = "exact_prop_da", y = "Percent_Mapped", add = "reg.line") +
            xlim(0.15,0.30) + ylim(99.9,100) +
            stat_cor(label.x = 0.25, label.y = 99.9)
my_ggsave('figures/percent_mapping_by_da_corr.pdf', da_corr)




# FIRE Plots -------------------------------------------------------------------------------------------------------------------------

# Barplot --- Number of FIRE peaks covered by Cell --- ALL, H1, H2
stat_df$fire_cov_prop_all <- stat_df$FIRE_cov_ALL / stat_df$mappable_peaks
stat_df$fire_cov_prop_H1 <- stat_df$FIRE_cov_H1 / stat_df$mappable_peaks
stat_df$fire_cov_prop_H2 <- stat_df$FIRE_cov_H2 / stat_df$mappable_peaks
stat_df$fire_cov_prop_both <- stat_df$FIRE_cov_BOTH_haps / stat_df$mappable_peaks

fire_cov_df_all <- stat_df %>% select(Cell, fire_cov_prop_all) %>% rename("fire_cov" = fire_cov_prop_all)
fire_cov_df_all$hap <- "All"
fire_cov_df_H1 <- stat_df %>% select(Cell, fire_cov_prop_H1) %>% rename("fire_cov" = fire_cov_prop_H1)
fire_cov_df_H1$hap <- "H1"
fire_cov_df_H2 <- stat_df %>% select(Cell, fire_cov_prop_H2) %>% rename("fire_cov" = fire_cov_prop_H2)
fire_cov_df_H2$hap <- "H2"
fire_cov_df_both <- stat_df %>% select(Cell, fire_cov_prop_both) %>% rename("fire_cov" = fire_cov_prop_both)
fire_cov_df_both$hap <- "Both"
fire_cov_merged <- rbind(fire_cov_df_all, fire_cov_df_H1, fire_cov_df_H2, fire_cov_df_both)
fire_cov_merged$fire_cov_count <- fire_cov_merged$fire_cov * stat_df$mappable_peaks[1]

fire_cov_merged$hap <- factor(fire_cov_merged$hap, levels=c('All','H1','H2','Both'))

# SAVE FIRE STATS AS A TABLE
write.csv(fire_cov_merged, "figures/FIRE_coverage_by_cell.csv", row.names = FALSE, quote=FALSE)


# Barplot --- Number of FIRE peaks ACTUATED by Cell --- ALL, H1, H2
fire_df <- read_tsv('scDAF_FIRE_actuation_MSP150_withProximal.tsv') # promoter proximal info added with "sc_actuation_comparison.py"

fire_act <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(fire_act) = c('Cell','hap','count')
for (cn in cell_names){
  temp_act <- data.frame(Cell = c(cn,cn,cn))
  temp_act$hap <- c("All","H1","H2")
  temp_act$count <- c(table(fire_df[[paste0(cn,'_ALL')]])['1'], table(fire_df[[paste0(cn,'_H1')]])['1'], table(fire_df[[paste0(cn,'_H2')]])['1'])
  fire_act <- rbind(fire_act, temp_act)
}

# SAVE FIRE ACTUATION AS A TABLE
write.csv(fire_act, "figures/FIRE_actuation_by_cell.csv", row.names = FALSE,quote=FALSE)



# Actuation heterogeneity between cells --------------------------------------------------------------------------------------------------
jaccard_dist_df <- read_tsv('jaccard_distances_actuation_by_cell.tsv') # computed in "sc_actuation_comparison.py"

jaccard_dist_df$Opposite <- FALSE
jaccard_dist_df[grepl("Opposite_Hap", jaccard_dist_df$Group), ]$Opposite <- TRUE
jaccard_dist_df <- jaccard_dist_df %>% mutate(Group2 = gsub("_Opposite_Hap", "", Group))

jaccard_dist_df$WithinCell <- FALSE
jaccard_dist_df[grepl("Within_Cell", jaccard_dist_df$Group), ]$WithinCell <- TRUE
jaccard_dist_df <- jaccard_dist_df %>% mutate(Group3 = gsub("_Within_Cell", "", Group2))

jaccard_dist_df$Type <- "Same Hap"
jaccard_dist_df[jaccard_dist_df$Opposite == TRUE,]$Type <- "Opposite Hap"
jaccard_dist_df[jaccard_dist_df$WithinCell == TRUE,]$Type <- "Within Cell"

jaccard_dist_df$Type <- factor(jaccard_dist_df$Type, levels=c("Within Cell", "Same Hap", "Opposite Hap"))


jaccard_violin <- ggviolin(jaccard_dist_df, "Group3", "Jaccard_Dist", color = "Type", fill = "Type",
                        palette = c("#00AFBB", "#E7B800", "#FC4E07"), position = position_dodge(1)) +
                        ylim(0.2,1) +
                        stat_compare_means(aes(group = Type, label = after_stat(p.format)), method = "t.test")+
                        theme_bw()
my_ggsave('figures/jaccard_dist_actuation_violin.pdf', jaccard_violin, width = 8, height = 5)

# T-tests
c1 <- jaccard_dist_df %>% filter(Type == "Within Cell", Group3 == "Overall")
c2 <- jaccard_dist_df %>% filter(Type == "Same Hap", Group3 == "Overall")
t.test(c1$Jaccard_Dist, c2$Jaccard_Dist)

c1 <- jaccard_dist_df %>% filter(Type == "Within Cell", Group3 == "Distal")
c2 <- jaccard_dist_df %>% filter(Type == "Same Hap", Group3 == "Distal")
t.test(c1$Jaccard_Dist, c2$Jaccard_Dist)

c1 <- jaccard_dist_df %>% filter(Type == "Within Cell", Group3 == "Proximal")
c2 <- jaccard_dist_df %>% filter(Type == "Same Hap", Group3 == "Proximal")
t.test(c1$Jaccard_Dist, c2$Jaccard_Dist)


my_comparisons <- list(c("Proximal", "Distal"), c("Proximal", "Overall"), c("Distal", "Overall"))
jaccard_violin_no_Opp <- ggviolin(jaccard_dist_df %>% filter(Opposite == FALSE, WithinCell == FALSE), "Group", "Jaccard_Dist", fill = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),) +
                        ylim(0.2,1) +
                        stat_compare_means(aes(label = after_stat(p.format)), comparisons = my_comparisons, method = "t.test", label.y = c(0.75,0.80,0.85))+
                        theme_bw()
my_ggsave('figures/jaccard_dist_actuation_violin_noOpposite.pdf', jaccard_violin_no_Opp, width = 8, height = 5)




# % actuation --> promoter-proximal vs promoter-distal
prop_act_proximal <- read_tsv('prop_actuation_proximal_vs_distal_by_cell.tsv')

pActuationByProximal_violin <- ggviolin(prop_act_proximal, "is_proximal", "prop", fill = "is_proximal",
                          palette = c("#00AFBB", "#E7B800"), position = position_dodge(1), add = "boxplot", add.params = list(fill = "white")) +
                          scale_y_continuous(labels = scales::percent, limits=c(0.2,1.0)) +
                          stat_compare_means(aes(group = is_proximal, label = after_stat(p.format)), method = "t.test") +
                          theme_bw() +
                          ylab("scDAF-seq Percent Actuation") + xlab("Promoter-distal vs. Promoter-proximal FIRE Peaks")
my_ggsave('figures/actuation_perc_promoter_proximal_distal_violin.pdf', pActuationByProximal_violin)







# Fiber-seq % actuation by proximal
FS_prop_act <- fire_df %>% group_by(is_proximal)  %>% summarise(prop = mean(prop_acc))

pFSActuationByProximal <- ggboxplot(fire_df, "is_proximal", "prop_acc", fill = "is_proximal", palette = c("#26AAE1", "#FCB041")) +
 scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
 stat_compare_means(aes(group = is_proximal, label = after_stat(p.format)), method = "t.test", label.y = c(0.02)) +
 ylab("Fiber-seq Percent Actuation") + xlab("Promoter-distal vs. Promoter-proximal FIRE Peaks") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())
my_ggsave('figures/actuation_perc_promoter_proximal_distal_Fiber-seq.pdf', pFSActuationByProximal)


# Fiber-seq Haplotype Selective peaks -----------------------------------------------------------------------------
fire_df$ID <- paste0(fire_df$chrom,":",fire_df$start,"-",fire_df$end)
volcano_df <- read_tsv('/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/HG002/hap1-vs-hap2/hap1-vs-hap2-volcano.tbl.gz')
volcano_df$ID <- paste0(volcano_df$'#chrom',":",volcano_df$peak_start,"-",volcano_df$peak_end)

# identify "discordant" peaks within single cells (actuated on H1 but not H2 or visa versa)
discordant_plots <- list()
density_plots <- list()
for (cn in cell_names){
  h1_name <- paste0(cn,"_H1")
  h2_name <- paste0(cn,"_H2")
  discordant_df <- fire_df %>% filter((.data[[h1_name]] == 1 & .data[[h2_name]] == 0) | (.data[[h1_name]] == 0 & .data[[h2_name]] == 1)) %>% select('chrom','start','end','prop_acc','bin',h1_name,h2_name,'ID')
  cov_df <- fire_df %>% filter(.data[[h1_name]] >= 0, .data[[h2_name]] >= 0) # require coverage on both haplotypes within the cell !!!
  filt_volcano <- volcano_df %>% filter(ID %in% fire_df$ID, ID %in% cov_df$ID)
  filt_volcano <- filt_volcano %>% mutate(discordant = ifelse(ID %in% discordant_df$ID, 'discordant', 'same'))
  discordant_plots[[cn]] <- ggplot(data = filt_volcano, aes(x = diff, y = -log10(p_value), col = discordant)) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_x_continuous(limits=c(-1,1)) +
    scale_y_continuous(limits=c(0,17)) +
    theme_classic(8) +
    theme(legend.position = "none") +
    xlab("Difference between paternal and maternal accessibility") +
    ggtitle(cn) +
    scale_color_manual(values = c("grey","#BE1E2D"), # to set the colours of our variable
                      labels = c("Consistent", "Discordant"))
  # Density plot of of the X-axis of the volcano plot for: 1) ALL HG002 peaks, 2) Discordant peaks
  dens_df <- rbind(data.frame(group='Same', diff=(filt_volcano %>% filter(!ID %in% discordant_df$ID))$diff), data.frame(group='Discordant', diff=(filt_volcano %>% filter(ID %in% discordant_df$ID))$diff))
  dens_df$group <- factor(dens_df$group, levels=c('Same','Discordant'))
  density_plots[[cn]] <- ggplot(dens_df, aes(x=diff, color=group)) +
                              geom_density(alpha=0.6) +
                              geom_vline(xintercept = c(-0.5, 0, 0.5), col = "gray", linetype = 'dashed') +
                              ylim(0,3.5) +
                              scale_x_continuous(limits=c(-1,1)) +
                              xlab("Difference between paternal and maternal accessibility") +
                              ggtitle(cn) +
                              scale_color_manual(values = c("grey","#BE1E2D")) +
                              theme_classic(8) +
                              theme(legend.position = "none")
}

disc_grid <- plot_grid(plotlist = discordant_plots)
ggsave("figures/discordant_peaks_volcano_hg002.pdf", disc_grid, width = 12, height = 8)

density_grid <- plot_grid(plotlist = density_plots)
ggsave("figures/discordant_peaks_density_hapDiff.pdf", density_grid, width = 12, height = 8)


