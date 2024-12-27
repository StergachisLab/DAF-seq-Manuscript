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


# Barplot --- Sequenced bases by Cell
stat_df$Gb_seq <- stat_df$bases / 1000000000

pSeqBases <- ggbarplot(stat_df, "Cell", "Gb_seq",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out", lab.nb.digits = 1) +
 ylab("HiFi bases (Gb)") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/sequenced_bases_barplot.pdf', pSeqBases, width = 4, height=4)

# Barplot --- N50 by Cell
pN50 <- ggbarplot(stat_df, "Cell", "N50",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out") +
 ylim(0,40000) +
 ylab("Collapsed Read N50") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/N50_barplot.pdf', pN50, width = 4, height=4)

# Barplot --- Longest collapse by Cell
stat_df$longKb <- stat_df$longest / 1000


# Barplot --- Number of Ultra Long reads per by Cell
pUltraLongCount <- ggbarplot(stat_df, "Cell", "num_100Kb",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out") +
#  ylim(0,500) +
 ylab("Ultra Long (100Kb) Read Count") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/ultraLong_Count_consensus_barplot.pdf', pUltraLongCount, width = 4, height=4)

# Barplot --- Bases in Ultra Long reads per by Cell

stat_df$UL_Mb <- stat_df$ultra_long_bases / 1000000

pUltraLongBases <- ggbarplot(stat_df, "Cell", "UL_Mb",
 fill = "steelblue", color = "steelblue",
 label = TRUE, lab.pos = "out", lab.nb.digits = 1) +
#  ylim(0,700) +
 ylab("Mb of sequence in Ultra Long (100Kb)") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/ultraLong_nBases_consensus_barplot.pdf', pUltraLongBases, width = 4, height=4)



# Barplot --- % bp of mappable genome covered (HG38, low cov regions) --- ALL, H1, H2
stat_df$cov_prop_all <- stat_df$ALL / stat_df$mappable_ref_bp
stat_df$cov_prop_H1 <- stat_df$H1 / stat_df$mappable_ref_bp
stat_df$cov_prop_H2 <- stat_df$H2 / stat_df$mappable_ref_bp

cov_df_all <- stat_df %>% select(Cell, cov_prop_all) %>% rename("cov" = cov_prop_all)
cov_df_all$hap <- "All"
cov_df_H1 <- stat_df %>% select(Cell, cov_prop_H1) %>% rename("cov" = cov_prop_H1)
cov_df_H1$hap <- "H1"
cov_df_H2 <- stat_df %>% select(Cell, cov_prop_H2) %>% rename("cov" = cov_prop_H2)
cov_df_H2$hap <- "H2"
cov_merged <- rbind(cov_df_all, cov_df_H1, cov_df_H2)

# percent covered
pCoverage <- ggbarplot(cov_merged, "Cell", "cov",
 fill = "hap", color = "hap",
 position = position_dodge(0.9)) +
 scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
 ylab("Mappable hg38 - Percent Coverage") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/percent_coverage_barplot.pdf', pCoverage)

# Gb covered
cov_merged$Gb_coverged <- (cov_merged$cov*stat_df$mappable_ref_bp)/1000000000

pCoverageGb <- ggbarplot(cov_merged, "Cell", "Gb_coverged",
 fill = "hap", color = "hap",
 position = position_dodge(0.9)) +
 ylab("Mappable hg38 - Covered Bases (Gb)") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/bases_coverage_barplot.pdf', pCoverageGb)

# SAVE COVERAGE STATS AS A TABLE
write.csv(cov_merged, "figures/coverage_proportion_by_cell.csv", row.names = FALSE)


# N50 of raw reads
raw_read1 <- read_tsv('PS00718_raw_read_lengths.tsv.gz', col_names=c('length','Cell'))
raw_read2 <- read_tsv('PS00756_raw_read_lengths.tsv.gz', col_names=c('length','Cell'))
raw_read3 <- read_tsv('PS00757_raw_read_lengths.tsv.gz', col_names=c('length','Cell'))
raw_read4 <- read_tsv('PS00758_raw_read_lengths.tsv.gz', col_names=c('length','Cell'))
raw_read_len <- rbind(raw_read1,raw_read2,raw_read3,raw_read4)

calculate_N50 <- function(lengths){
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  cumulative_lengths <- cumsum(sorted_lengths)
  half_total_length <- sum(lengths) / 2
  N50 <- sorted_lengths[cumulative_lengths >= half_total_length][1]
  return(N50)
}
N50_all_reads <- calculate_N50(raw_read_len$length)

PS00718_N50 <- calculate_N50(raw_read1$length)
PS00756_N50 <- calculate_N50(raw_read2$length)
PS00757_N50 <- calculate_N50(raw_read3$length)
PS00758_N50 <- calculate_N50(raw_read4$length)

# Save Raw Read N50 data
raw_N50_df <- data.frame(N50_all = N50_all_reads, PS00718_N50 = PS00718_N50, PS00756_N50 = PS00756_N50, PS00757_N50 = PS00757_N50, PS00758_N50 = PS00758_N50)
write.csv(raw_N50_df, "figures/raw_read_N50.csv", row.names = FALSE, quote=FALSE)


# density plot of raw read lengths colored by cell
raw_read_length_density <- ggplot(raw_read_len, aes(x=length, color=Cell)) +
                            geom_density(alpha=0.6) +
                            scale_x_continuous(trans = 'log2', breaks = scales::breaks_log(2, n = 7), limits = c(1, 524288)) +
                            scale_color_manual(values = c("#00A14B","#F26522","#00AEEF","#7F3F98")) +
                            theme_bw(10)
my_ggsave('figures/raw_read_length_density.pdf', raw_read_length_density)

# density plot of COLLAPSED read lengths colored by cell
consensus_lengths <- read_tsv('consensus_read_lengths.tsv', col_names=c('length','Cell'))

collapsed_read_length_density <- ggplot(consensus_lengths, aes(x=length, color=Cell)) +
                            geom_density(alpha=0.6) +
                            scale_x_continuous(trans = 'log2', breaks = scales::breaks_log(2, n = 6), limits = c(1, 524288)) +
                            scale_color_manual(values = c("#00A14B","#F26522","#00AEEF","#7F3F98",'red')) +
                            theme_bw(10)
my_ggsave('figures/collapsed_read_length_density.pdf', collapsed_read_length_density)


# Read length Density Plots weighted by total bp
raw_read_length_density_weight <- ggplot(raw_read_len, aes(x=length, color=Cell)) +
                            geom_density(alpha=0.6, aes(weight=length)) +
                            scale_x_continuous(trans = 'log2', breaks = scales::breaks_log(2, n = 7), limits = c(1, 524288)) +
                            scale_color_manual(values = c("#00A14B","#F26522","#00AEEF","#7F3F98")) +
                            theme_bw(10)
my_ggsave('figures/raw_read_length_density_weighted.pdf', raw_read_length_density_weight)

collapsed_read_length_density_weight <- ggplot(consensus_lengths, aes(x=length, color=Cell)) +
                            geom_density(alpha=0.6, aes(weight=length)) +
                            scale_x_continuous(trans = 'log2', breaks = scales::breaks_log(2, n = 6), limits = c(1, 524288)) +
                            scale_color_manual(values = c("#00A14B","#F26522","#00AEEF","#7F3F98",'red')) +
                            theme_bw(10)
my_ggsave('figures/collapsed_read_length_density_weighted.pdf', collapsed_read_length_density_weight)



# Yak switch rates -----------------------------------------------------------------------
yak_df <- read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/yak/yak_error_rate_summary.csv')

yak_switch <- yak_df %>% filter(yak_type == 'W')

yak_switch_p <- ggbarplot(yak_switch, "sample", "rate",
                  fill = "read_type", color = "read_type", palette = "Paired",
                  facet.by = "strand", position = position_dodge(0.9))

ggsave('figures/yak_switch_rate_barplot.png', yak_switch_p)



# Karyoplots -----------------------------------------------------------------------------------------------------------

# Karyoplot coverage --- chr2 uncollapsed MAPPABLE reads
library(karyoploteR)
library(GenomicAlignments)

#
# Regions for Coverage were precomputed in 03_mappable_coverage_precompute.sh !!!!
#

auto_chroms <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')


# Raw clipped reads ---------------------------------------------
raw_clipped_CT_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00758_clipped_CT.bed.gz', col_names = c('chrom','start','end','name'))
raw_clipped_GA_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00758_clipped_GA.bed.gz', col_names = c('chrom','start','end','name'))
gr_PS00758_clipped_CT <- makeGRangesFromDataFrame(raw_clipped_CT_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL,
                                       seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
gr_PS00758_clipped_GA <- makeGRangesFromDataFrame(raw_clipped_GA_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL,
                                       seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)

# Raw clipped reads FILTERED ---------------------------------------------
raw_clipped_BothStrands_FILT_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00758_clipped_FILTERED.bed.gz', col_names = c('chrom','start','end','name'))
gr_PS00758_clipped_FILT_BothStrands <- makeGRangesFromDataFrame(raw_clipped_BothStrands_FILT_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL,
                                       seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)



# COLLAPSED Reads filtered by MAPPABLE regoins AND High Coverage Regions
collapsed_filt_CT_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00758_CT_WIDE_high_cov_FILTERED_consensus_reads.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
collapsed_filt_GA_BED <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/PS00758_GA_WIDE_high_cov_FILTERED_consensus_reads.bed.gz', col_names = c('chrom','start','end','name','score','strand'))
granges_collapsed_filt_CT <- makeGRangesFromDataFrame(collapsed_filt_CT_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL,
                                       seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)
granges_collapsed_filt_GA <- makeGRangesFromDataFrame(collapsed_filt_GA_BED, keep.extra.columns=FALSE, ignore.strand = TRUE, seqinfo=NULL,
                                       seqnames.field=c("chrom"),start.field="start",end.field=c("end"), starts.in.df.are.0based=TRUE)

# COLLAPSED Read UNFILTERED
collapsed_ALL_CT_BAM <- '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/haplotagged/PS00758_consensus_CT_HG38_corrected.haplotagged.bam'
gAlign_collapsed_ALL_CT <- readGAlignments(collapsed_ALL_CT_BAM, use.names = TRUE)
granges_collapsed_ALL_CT <- granges(gAlign_collapsed_ALL_CT)
collapsed_ALL_GA_BAM <- '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/haplotagged/PS00758_consensus_GA_HG38_corrected.haplotagged.bam'
gAlign_collapsed_ALL_GA <- readGAlignments(collapsed_ALL_GA_BAM, use.names = TRUE)
granges_collapsed_ALL_GA <- granges(gAlign_collapsed_ALL_GA)


# COMPARE Raw Coverage by Strand -------------------
# png("figures/PS00758_raw_reads_clipped_bothStrands_ALL_karyo_chr2.png", width = 3200, height = 3200, res = 600)
pdf("figures/PS00758_raw_reads_clipped_bothStrands_ALL_karyo_chr2.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
kp <- kpPlotCoverage(kp, data=gr_PS00758_clipped_CT, col="#26AAE1", data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1)
kp <- kpPlotCoverage(kp, data=gr_PS00758_clipped_GA, col="#FCB041", data.panel = 2) # GREEN
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 2)
dev.off()

# COMPARE Collapsed Unfiltered by Strand -------------------
# png("figures/PS00758_raw_reads_clipped_bothStrands_ALL_karyo_chr2.png", width = 3200, height = 3200, res = 600)
pdf("figures/PS00758_Collapsed_UNFILTERED_clipped_bothStrands_ALL_karyo_chr2.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
kp <- kpPlotCoverage(kp, data=granges_collapsed_ALL_CT, col="#26AAE1", data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1)
kp <- kpPlotCoverage(kp, data=granges_collapsed_ALL_GA, col="#FCB041", data.panel = 2) # GREEN
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 2)
dev.off()


# COMPARE Coverage Raw vs Collapsed > 10 Kb -------------------
gr_long_collapsed_CT <- granges_collapsed_filt_CT[width(granges_collapsed_filt_CT) >= 10000]
gr_long_collapsed_GA <- granges_collapsed_filt_GA[width(granges_collapsed_filt_GA) >= 10000]

# png("figures/PS00758_raw_reads_CT_clipped_ALL_vs_Mappable_10Kb_karyo.png", width = 3200, height = 3200, res = 600)
pdf("figures/PS00758_raw_reads_CT_clipped_ALL_vs_Mappable_10Kb_karyo_chr2.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
kp <- kpPlotCoverage(kp, data=gr_PS00758_clipped_CT, col="#26AAE1", data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=4)
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_CT, col="#26AAE1", data.panel = 2)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 2, numticks=4)
dev.off()

# find regions above the cutoff
highCov_CT <- kp$latest.plot$computed.values$coverage[elementMetadata(kp$latest.plot$computed.values$coverage)$coverage >= 3]
width(highCov_CT)


# COMPARE Coverage Raw vs Collapsed > 10 Kb -- MERGED STRANDS !!!!!! -------------------
# gr_PS00758_clipped_BothStrands <- c(gr_PS00758_clipped_CT,gr_PS00758_clipped_GA)
gr_long_collapsed_BothStrands <- c(gr_long_collapsed_CT,gr_long_collapsed_GA)

# FIGURE 6 PLOT! -----------------------------
# SAVE ONE AT A Time
png("figures/PS00758_raw_reads_BothStrandsMerged_clipped_FILT_karyo_chr2.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
kp <- kpPlotCoverage(kp, data=gr_PS00758_clipped_FILT_BothStrands, data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1, numticks=5)
dev.off()
# Save AXIS LABELS as a PDF
clip_filt_max <- kp$latest.plot$computed.values$max.coverage
pdf("figures/PS00758_raw_reads_BothStrandsMerged_clipped_FILT_karyo_chr2_AXIS.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
kpAxis(kp, ymax=clip_filt_max, data.panel = 1, numticks=5)
dev.off()

png("figures/PS00758_raw_reads_BothStrandsMerged_Mappable_10Kb_karyo_chr2.png", width = 6400, height = 6400, res = 1200)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
kp <- kpPlotCoverage(kp, data=gr_long_collapsed_BothStrands, data.panel = 1, r0=0, r1=0.875)
kpAxis(kp, ymax=8, data.panel = 1, numticks=5)
dev.off()
# Save AXIS LABELS as a PDF
pdf("figures/PS00758_raw_reads_BothStrandsMerged_Mappable_10Kb_karyo_chr2_AXIS.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
kpAxis(kp, ymax=8, data.panel = 1, numticks=5)
dev.off()



# # ----------------------------------------------------------------------

# # COMPARE Collapsed > 10 Kb --- Top vs Bottom -------------------
# # png("figures/PS00758_raw_reads_CT_vs_GA_Mappable_10Kb_karyo.png", width = 3200, height = 3200, res = 600)
# pdf("figures/PS00758_raw_reads_CT_vs_GA_Mappable_10Kb_karyo_chr2.pdf")
# kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = "chr2")
# kp <- kpPlotCoverage(kp, data=gr_long_collapsed_CT, col="#26AAE1", data.panel = 1, r0=0, r1=0.75)
# kpAxis(kp, ymax=4, data.panel = 1, numticks=5)
# kp <- kpPlotCoverage(kp, data=gr_long_collapsed_GA, col="#FCB041", data.panel = 2)
# kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 2, numticks=5)
# dev.off()

# # find regions above the cutoff
# highCov_GA <- kp$latest.plot$computed.values$coverage[elementMetadata(kp$latest.plot$computed.values$coverage)$coverage >= 3]
# width(highCov_GA)


# # COMPARE Collapsed > 10 Kb --- Top ONLY ALL Chromosomes -------------------
# # png("figures/PS00758_raw_reads_CT_AllChroms_Mappable_10Kb_karyo.png", width = 4800, height = 9600, res = 900)
# pdf("figures/PS00758_raw_reads_CT_AllChroms_Mappable_10Kb_karyo.pdf")
# kp <- plotKaryotype(genome="hg38", chromosomes = auto_chroms)
# kp <- kpPlotCoverage(kp, data=gr_long_collapsed_CT, col="#26AAE1", data.panel = 1)
# dev.off()

# # Look at RAW coverage over all chromosomes ------------------
# # png("figures/PS00758_raw_reads_clipped_bothStrands_AllChroms_karyo.png", width = 4800, height = 9600, res = 900)
# pdf("figures/PS00758_raw_reads_clipped_bothStrands_AllChroms_karyo.pdf")
# gr_PS00758_clipped_BothStrands <- c(gr_PS00758_clipped_CT, gr_PS00758_clipped_GA)
# kp <- plotKaryotype(genome="hg38", chromosomes = auto_chroms)
# kp <- kpPlotCoverage(kp, data=gr_PS00758_clipped_BothStrands)
# dev.off() # Close the device




# MSP Actuation --------------------------------------------------------------------

# Actuation at FIRE peaks binned by Fiber-seq actuation at different MSP lengths

msp_150 <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP150.tsv')
msp_200 <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP200.tsv')
msp_300 <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP300.tsv')

# Filter each peakset for peaks with width +/- 25bp to each MSP cutoff (ex. 175-225bp wide for MSP 200 cutoff)
msp_150$peak_width <- msp_150$end - msp_150$start
msp_150 <- msp_150 %>% filter(peak_width >= 125, peak_width <= 175)
msp_200$peak_width <- msp_200$end - msp_200$start
msp_200 <- msp_200 %>% filter(peak_width >= 175, peak_width <= 225)
msp_300$peak_width <- msp_300$end - msp_300$start
msp_300 <- msp_300 %>% filter(peak_width >= 275, peak_width <= 325)

# # Coverage on BOTH haplotypes
# msp_150_both <- msp_150 %>% select(chrom:bin)
# msp_150 %>% filter(PS00718_nCov == 2) %>% group_by(bin) %>% select(PS00718_nMSP) %>% sum()


# for each bin: count the number of haps x peaks in that bin covered & the total accessible haps
PS00718_bins <- msp_150 %>% filter(PS00718_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00718_nCov), nMSP = sum(PS00718_nMSP))
PS00756_bins <- msp_150 %>% filter(PS00756_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00756_nCov), nMSP = sum(PS00756_nMSP))
PS00757_bins <- msp_150 %>% filter(PS00757_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00757_nCov), nMSP = sum(PS00757_nMSP))
PS00758_bins <- msp_150 %>% filter(PS00758_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00758_nCov), nMSP = sum(PS00758_nMSP))
msp150_bins <- data.frame(bin = PS00718_bins$bin)
msp150_bins$nCov <- PS00718_bins$nCov + PS00756_bins$nCov + PS00757_bins$nCov + PS00758_bins$nCov
msp150_bins$nMSP <- PS00718_bins$nMSP + PS00756_bins$nMSP + PS00757_bins$nCov + PS00758_bins$nMSP
msp150_bins$prop <- msp150_bins$nMSP / msp150_bins$nCov

PS00718_bins <- msp_200 %>% filter(PS00718_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00718_nCov), nMSP = sum(PS00718_nMSP))
PS00756_bins <- msp_200 %>% filter(PS00756_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00756_nCov), nMSP = sum(PS00756_nMSP))
PS00757_bins <- msp_200 %>% filter(PS00757_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00757_nCov), nMSP = sum(PS00757_nMSP))
PS00758_bins <- msp_200 %>% filter(PS00758_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00758_nCov), nMSP = sum(PS00758_nMSP))
msp200_bins <- data.frame(bin = PS00718_bins$bin)
msp200_bins$nCov <- PS00718_bins$nCov + PS00756_bins$nCov + PS00757_bins$nCov + PS00758_bins$nCov
msp200_bins$nMSP <- PS00718_bins$nMSP + PS00756_bins$nMSP + PS00757_bins$nCov + PS00758_bins$nMSP
msp200_bins$prop <- msp200_bins$nMSP / msp200_bins$nCov

PS00718_bins <- msp_300 %>% filter(PS00718_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00718_nCov), nMSP = sum(PS00718_nMSP))
PS00756_bins <- msp_300 %>% filter(PS00756_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00756_nCov), nMSP = sum(PS00756_nMSP))
PS00757_bins <- msp_300 %>% filter(PS00757_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00757_nCov), nMSP = sum(PS00757_nMSP))
PS00758_bins <- msp_300 %>% filter(PS00758_nCov >= 1) %>% group_by(bin) %>% summarise(nCov = sum(PS00758_nCov), nMSP = sum(PS00758_nMSP))
msp300_bins <- data.frame(bin = PS00718_bins$bin)
msp300_bins$nCov <- PS00718_bins$nCov + PS00756_bins$nCov + PS00757_bins$nCov + PS00758_bins$nCov
msp300_bins$nMSP <- PS00718_bins$nMSP + PS00756_bins$nMSP + PS00757_bins$nCov + PS00758_bins$nMSP
msp300_bins$prop <- msp300_bins$nMSP / msp300_bins$nCov


bin_names <- c('10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%')

msp150_barplot <- ggbarplot(msp150_bins, "bin", "prop",
 fill = "bin", color = "bin",
 position = position_dodge(0.9)) +
 scale_x_discrete(labels= bin_names) +
 scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
 ylab("Average Actuation") +
 theme_bw(5) +
 theme(panel.grid.major.x = element_blank(), legend.position="none") +
 ggtitle("FIRE peak bins - 150 bp MSP")

msp200_barplot <- ggbarplot(msp200_bins, "bin", "prop",
 fill = "bin", color = "bin",
 position = position_dodge(0.9)) +
 scale_x_discrete(labels= bin_names) +
 scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
 ylab("Average Actuation") +
 theme_bw(5) +
 theme(panel.grid.major.x = element_blank(), legend.position="none") +
 ggtitle("FIRE peak bins - 200 bp MSP")

 msp300_barplot <- ggbarplot(msp300_bins, "bin", "prop",
 fill = "bin", color = "bin",
 position = position_dodge(0.9)) +
 scale_x_discrete(labels= bin_names) +
 scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
 ylab("Average Actuation") +
 theme_bw(5) +
 theme(panel.grid.major.x = element_blank(), legend.position="none") +
 ggtitle("FIRE peak bins - 300 bp MSP")


msp_bins_grid <- plot_grid(msp150_barplot, msp200_barplot, msp300_barplot, nrow = 1)

ggsave('figures/actuation_bins_by_MSP_Length_barplot.pdf', msp_bins_grid, width = 8.5, height = 3)


# Autocorrelation -------------------------------------------------------------------------------------------------------------------------

PS00718_acf <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/PS00718_ft_acf.tsv', col_names=c('acf','lag','corr'))
PS00756_acf <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/PS00756_ft_acf.tsv', col_names=c('acf','lag','corr'))
PS00757_acf <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/PS00757_ft_acf.tsv', col_names=c('acf','lag','corr'))
PS00758_acf <- read_tsv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/PS00758_ft_acf.tsv', col_names=c('acf','lag','corr'))

PS00718_acf$cell <- 'PS00718'
PS00756_acf$cell <- 'PS00756'
PS00757_acf$cell <- 'PS00757'
PS00758_acf$cell <- 'PS00758'

acf_merged <- rbind(PS00718_acf,PS00756_acf,PS00757_acf,PS00758_acf)

acf_merged$change <- c(0, (acf_merged$corr[1:(length(acf_merged$corr)-1)] * acf_merged$corr[2:length(acf_merged$corr)]) <= 0)

acf1 <- ggplot(acf_merged %>% filter(lag > 25), aes(x=lag, y=corr, color=cell)) +
  geom_hline(aes(yintercept=0), color="darkblue", linewidth=1, linetype="dashed") + 
  geom_vline(data=NULL, aes(xintercept=147), linetype="dashed", color="black", alpha=0.5) +
  geom_line() +
  scale_color_manual(values = c("PS00718" = "#00A14B", "PS00756" = "#F26522", "PS00757" = "#00AEEF", "PS00758" = "#7F3F98")) +
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

pDaRateKb <- ggbarplot(da_rates, "Cell", "rate_Kb",
 fill = "steelblue", color = "steelblue") +
 ylab("Deaminations per Kb") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/deaminations_per_1Kb_by_cell_barplot.pdf', pDaRateKb, width = 4, height=4)

write.csv(da_rates, "figures/deaminations_per_1Kb_by_cell.csv", row.names = FALSE)



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

pFireCoverage <- ggbarplot(fire_cov_merged, "Cell", "fire_cov",
 fill = "hap", color = "hap",
 position = position_dodge(0.9)) +
 scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
 ylab("Mappable hg38 - Percent FIRE Peaks Covered") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/percent_fire_coverage_barplot.pdf', pFireCoverage, width = 6, height=4)

pFireCount <- ggbarplot(fire_cov_merged, "Cell", "fire_cov_count",
 fill = "hap", color = "hap",
 position = position_dodge(0.9)) +
 ylab("Mappable hg38 - Number of FIRE Peaks Covered") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/fire_coverage_count_barplot.pdf', pFireCount, width = 6, height=4)

# SAVE FIRE STATS AS A TABLE
write.csv(fire_cov_merged, "figures/FIRE_coverage_by_cell.csv", row.names = FALSE, quote=FALSE)


# Barplot --- Number of FIRE peaks ACTUATED by Cell --- ALL, H1, H2

# 150 bp MSP Length !!!!!!!!!!!!!!!!!!
fire_df <- read_tsv('scDAF_FIRE_actuation_MSP150_withProximal.tsv') # promoter proximal info added with "sc_actuation_comparison.py"

cellA <- data.frame(Cell = c("PS00718","PS00718","PS00718"))
cellA$hap <- c("All","H1","H2")
cellA$count <- c(table(fire_df$PS00718_ALL)['1'], table(fire_df$PS00718_H1)['1'], table(fire_df$PS00718_H2)['1'])
cellB <- data.frame(Cell = c("PS00756","PS00756","PS00756"))
cellB$hap <- c("All","H1","H2")
cellB$count <- c(table(fire_df$PS00756_ALL)['1'], table(fire_df$PS00756_H1)['1'], table(fire_df$PS00756_H2)['1'])
cellC <- data.frame(Cell = c("PS00757","PS00757","PS00757"))
cellC$hap <- c("All","H1","H2")
cellC$count <- c(table(fire_df$PS00757_ALL)['1'], table(fire_df$PS00757_H1)['1'], table(fire_df$PS00757_H2)['1'])
cellD <- data.frame(Cell = c("PS00758","PS00758","PS00758"))
cellD$hap <- c("All","H1","H2")
cellD$count <- c(table(fire_df$PS00758_ALL)['1'], table(fire_df$PS00758_H1)['1'], table(fire_df$PS00758_H2)['1'])

fire_act <- rbind(cellA, cellB, cellC, cellD)

pActFireCount <- ggbarplot(fire_act, "Cell", "count",
 fill = "hap", color = "hap",
 position = position_dodge(0.9)) +
 ylab("Mappable hg38 - Number of Actuated FIRE Peaks") + 
 theme_bw(10) +
 ylim(0,80000) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/actuated_fire_count_barplot.pdf', pActFireCount, width = 4.5, height=4)

# SAVE FIRE ACTUATION AS A TABLE
write.csv(fire_act, "figures/FIRE_actuation_by_cell.csv", row.names = FALSE,quote=FALSE)


# Actuation heterogeneity between cells --------------------------------------------------------------------------------------------------
jaccard_dist_df <- read_tsv('jaccard_distances_actuation_by_cell.tsv') # computed in "sc_actuation_comparison.py"

jaccard_dotplot <- ggdotplot(jaccard_dist_df, "Proximal", "Jaccard_Dist", fill = "Proximal", color = "Proximal",
                        palette = c("#26AAE1", "#FCB041",'red')) +
                        ylim(0.4,0.8)
my_ggsave('figures/jaccard_dist_actuation_dotplot.pdf', jaccard_dotplot)

# % actuation --> promoter-proximal vs promoter-distal
prop_act_proximal <- read_tsv('prop_actuation_proximal_vs_distal_by_cell.tsv')

pActuationByProximal <- ggbarplot(prop_act_proximal, "Cell", "prop",
 fill = "is_proximal", color = "is_proximal", palette = c("#26AAE1", "#FCB041"), position = position_dodge(0.9)) +
 scale_y_continuous(labels = scales::percent) +
 ylab("Mean Percent Actuation") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/actuation_perc_promoter_proximal_distal.pdf', pActuationByProximal)

# Fiber-seq % actuation by proximal

FS_prop_act <- fire_df %>% group_by(is_proximal)  %>% summarise(prop = mean(prop_acc))

pFSActuationByProximal <- ggbarplot(FS_prop_act, "is_proximal", "prop",
 fill = "is_proximal", color = "is_proximal", palette = c("#26AAE1", "#FCB041"), position = position_dodge(0.9)) +
 scale_y_continuous(labels = scales::percent) +
 ylab("Mean Percent Actuation") +
 theme_bw(10) +
 theme(panel.grid.major.x = element_blank())

my_ggsave('figures/actuation_perc_promoter_proximal_distal_Fiber-seq.pdf', pFSActuationByProximal)


# Fiber-seq Haplotype Selective peaks -----------------------------------------------------------------------------
fire_df$ID <- paste0(fire_df$chrom,":",fire_df$start,"-",fire_df$end)
volcano_df <- read_tsv('/mmfs1/gscratch/stergachislab/mvollger/projects/FIREv2.0/results/HG002/hap1-vs-hap2/hap1-vs-hap2-volcano.tbl.gz')
volcano_df$ID <- paste0(volcano_df$'#chrom',":",volcano_df$peak_start,"-",volcano_df$peak_end)


# identify "discordant" peaks within single cells (actuated on H1 but not H2 or visa versa)

# PS00718
discordant_PS00718 <- fire_df %>% filter((PS00718_H1==1 & PS00718_H2 == 0) | (PS00718_H1==0 & PS00718_H2 == 1)) %>% select('chrom','start','end','prop_acc','bin','PS00718_H1','PS00718_H2','ID')
cov_PS00718 <- fire_df %>% filter(PS00718_H1>=0, PS00718_H2>=0) # require coverage on both haplotypes within the cell !!!
filt_volcano_PS00718 <- volcano_df %>% filter(ID %in% fire_df$ID, ID %in% cov_PS00718$ID)
filt_volcano_PS00718 <- filt_volcano_PS00718 %>% mutate(discordant = ifelse(ID %in% discordant_PS00718$ID, 'discordant', 'same'))
volcano_hg002_PS00718 <- ggplot(data = filt_volcano_PS00718, aes(x = diff, y = -log10(p_value), col = discordant)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_x_continuous(limits=c(-1,1)) +
  scale_y_continuous(limits=c(0,17)) +
  theme_classic(10) +
  theme(legend.position = "none") +
  xlab("Difference between paternal and maternal accessibility") +
  ggtitle("PS00718") +
  scale_color_manual(values = c("grey","#BE1E2D"), # to set the colours of our variable
                     labels = c("Consistent", "Discordant"))

# PS00756
discordant_PS00756 <- fire_df %>% filter((PS00756_H1==1 & PS00756_H2 == 0) | (PS00756_H1==0 & PS00756_H2 == 1)) %>% select('chrom','start','end','prop_acc','bin','PS00756_H1','PS00756_H2','ID')
cov_PS00756 <- fire_df %>% filter(PS00756_H1>=0, PS00756_H2>=0) # require coverage on both haplotypes within the cell !!!
filt_volcano_PS00756 <- volcano_df %>% filter(ID %in% fire_df$ID, ID %in% cov_PS00756$ID)
filt_volcano_PS00756 <- filt_volcano_PS00756 %>% mutate(discordant = ifelse(ID %in% discordant_PS00756$ID, 'discordant', 'same'))
volcano_hg002_PS00756 <- ggplot(data = filt_volcano_PS00756, aes(x = diff, y = -log10(p_value), col = discordant)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_x_continuous(limits=c(-1,1)) +
  scale_y_continuous(limits=c(0,17)) +
  theme_classic(10) +
  theme(legend.position = "none") +
  xlab("Difference between paternal and maternal accessibility") +
  ggtitle("PS00756") +
  scale_color_manual(values = c("grey","#BE1E2D"), # to set the colours of our variable
                     labels = c("Consistent", "Discordant"))

# PS00757
discordant_PS00757 <- fire_df %>% filter((PS00757_H1==1 & PS00757_H2 == 0) | (PS00757_H1==0 & PS00757_H2 == 1)) %>% select('chrom','start','end','prop_acc','bin','PS00757_H1','PS00757_H2','ID')
cov_PS00757 <- fire_df %>% filter(PS00757_H1>=0, PS00757_H2>=0) # require coverage on both haplotypes within the cell !!!
filt_volcano_PS00757 <- volcano_df %>% filter(ID %in% fire_df$ID, ID %in% cov_PS00757$ID)
filt_volcano_PS00757 <- filt_volcano_PS00757 %>% mutate(discordant = ifelse(ID %in% discordant_PS00757$ID, 'discordant', 'same'))
volcano_hg002_PS00757 <- ggplot(data = filt_volcano_PS00757, aes(x = diff, y = -log10(p_value), col = discordant)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_x_continuous(limits=c(-1,1)) +
  scale_y_continuous(limits=c(0,17)) +
  theme_classic(10) +
  theme(legend.position = "none") +
  xlab("Difference between paternal and maternal accessibility") +
  ggtitle("PS00757") +
  scale_color_manual(values = c("grey","#BE1E2D"), # to set the colours of our variable
                     labels = c("Consistent", "Discordant"))

# PS00758
discordant_PS00758 <- fire_df %>% filter((PS00758_H1==1 & PS00758_H2 == 0) | (PS00758_H1==0 & PS00758_H2 == 1)) %>% select('chrom','start','end','prop_acc','bin','PS00758_H1','PS00758_H2','ID')
cov_PS00758 <- fire_df %>% filter(PS00758_H1>=0, PS00758_H2>=0) # require coverage on both haplotypes within the cell !!!
filt_volcano_PS00758 <- volcano_df %>% filter(ID %in% fire_df$ID, ID %in% cov_PS00758$ID)
filt_volcano_PS00758 <- filt_volcano_PS00758 %>% mutate(discordant = ifelse(ID %in% discordant_PS00758$ID, 'discordant', 'same'))
volcano_hg002_PS00758 <- ggplot(data = filt_volcano_PS00758, aes(x = diff, y = -log10(p_value), col = discordant)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_x_continuous(limits=c(-1,1)) +
  scale_y_continuous(limits=c(0,17)) +
  theme_classic(10) +
  theme(legend.position = "none") +
  xlab("Difference between paternal and maternal accessibility") +
  ggtitle("PS00758") +
  scale_color_manual(values = c("grey","#BE1E2D"), # to set the colours of our variable
                     labels = c("Consistent", "Discordant"))

my_ggsave('figures/volcano_hg002_PS00718.pdf', volcano_hg002_PS00718)
my_ggsave('figures/volcano_hg002_PS00756.pdf', volcano_hg002_PS00756)
my_ggsave('figures/volcano_hg002_PS00757.pdf', volcano_hg002_PS00757)
my_ggsave('figures/volcano_hg002_PS00758.pdf', volcano_hg002_PS00758)


# Density plot of of the X-axis of the volcano plot for: 1) ALL HG002 peaks, 2) Discordant peaks
dens_df_PS00718 <- rbind(data.frame(group='Same', diff=(filt_volcano_PS00718 %>% filter(!ID %in% discordant_PS00718$ID))$diff), data.frame(group='Discordant', diff=(filt_volcano_PS00718 %>% filter(ID %in% discordant_PS00718$ID))$diff))
dens_df_PS00718$group <- factor(dens_df_PS00718$group, levels=c('Same','Discordant'))
density_p_PS00718 <- ggplot(dens_df_PS00718, aes(x=diff, color=group)) +
                            geom_density(alpha=0.6) +
                            geom_vline(xintercept = c(-0.5, 0, 0.5), col = "gray", linetype = 'dashed') +
                            ylim(0,3.5) +
                            scale_x_continuous(limits=c(-1,1)) +
                            xlab("Difference between paternal and maternal accessibility") +
                            ggtitle("PS00718") +
                            scale_color_manual(values = c("grey","#BE1E2D")) +
                            theme_classic(10) +
                            theme(legend.position = "none")

dens_df_PS00756 <- rbind(data.frame(group='Same', diff=(filt_volcano_PS00756 %>% filter(!ID %in% discordant_PS00756$ID))$diff), data.frame(group='Discordant', diff=(filt_volcano_PS00756 %>% filter(ID %in% discordant_PS00756$ID))$diff))
dens_df_PS00756$group <- factor(dens_df_PS00756$group, levels=c('Same','Discordant'))
density_p_PS00756 <- ggplot(dens_df_PS00756, aes(x=diff, color=group)) +
                            geom_density(alpha=0.6) +
                            geom_vline(xintercept = c(-0.5, 0, 0.5), col = "gray", linetype = 'dashed') +
                            ylim(0,3.5) +
                            scale_x_continuous(limits=c(-1,1)) +
                            xlab("Difference between paternal and maternal accessibility") +
                            ggtitle("PS00756") +
                            scale_color_manual(values = c("grey","#BE1E2D")) +
                            theme_classic(10) +
                            theme(legend.position = "none")

dens_df_PS00757 <- rbind(data.frame(group='Same', diff=(filt_volcano_PS00757 %>% filter(!ID %in% discordant_PS00757$ID))$diff), data.frame(group='Discordant', diff=(filt_volcano_PS00757 %>% filter(ID %in% discordant_PS00757$ID))$diff))
dens_df_PS00757$group <- factor(dens_df_PS00757$group, levels=c('Same','Discordant'))
density_p_PS00757 <- ggplot(dens_df_PS00757, aes(x=diff, color=group)) +
                            geom_density(alpha=0.6) +
                            geom_vline(xintercept = c(-0.5, 0, 0.5), col = "gray", linetype = 'dashed') +
                            ylim(0,3.5) +
                            scale_x_continuous(limits=c(-1,1)) +
                            xlab("Difference between paternal and maternal accessibility") +
                            ggtitle("PS00757") +
                            scale_color_manual(values = c("grey","#BE1E2D")) +
                            theme_classic(10) +
                            theme(legend.position = "none")

dens_df_PS00758 <- rbind(data.frame(group='Same', diff=(filt_volcano_PS00758 %>% filter(!ID %in% discordant_PS00758$ID))$diff), data.frame(group='Discordant', diff=(filt_volcano_PS00758 %>% filter(ID %in% discordant_PS00758$ID))$diff))
dens_df_PS00758$group <- factor(dens_df_PS00758$group, levels=c('Same','Discordant'))
density_p_PS00758 <- ggplot(dens_df_PS00758, aes(x=diff, color=group)) +
                            geom_density(alpha=0.6) +
                            geom_vline(xintercept = c(-0.5, 0, 0.5), col = "gray", linetype = 'dashed') +
                            ylim(0,3.5) +
                            scale_x_continuous(limits=c(-1,1)) +
                            xlab("Difference between paternal and maternal accessibility") +
                            ggtitle("PS00758") +
                            scale_color_manual(values = c("grey","#BE1E2D")) +
                            theme_classic(10) +
                            theme(legend.position = "none")


my_ggsave('figures/density_hapDiff_PS00718.pdf', density_p_PS00718)
my_ggsave('figures/density_hapDiff_PS00756.pdf', density_p_PS00756)
my_ggsave('figures/density_hapDiff_PS00757.pdf', density_p_PS00757)
my_ggsave('figures/density_hapDiff_PS00758.pdf', density_p_PS00758)




# # Chi-square contingency tables --------------------------------------------------------------------------------------------------

# # for each cell, for each peak covered by both haps, track Neither, H1 only, H2 only, or both
# neither = 0
# h1_only = 0
# h2_only = 0
# both = 0

# neither = neither + nrow(fire_df %>% filter(PS00718_H1 == 0, PS00718_H2 == 0))
# h1_only = h1_only + nrow(fire_df %>% filter(PS00718_H1 == 1, PS00718_H2 == 0))
# h2_only = h2_only + nrow(fire_df %>% filter(PS00718_H1 == 0, PS00718_H2 == 1))
# both = both + nrow(fire_df %>% filter(PS00718_H1 == 1, PS00718_H2 == 1))

# neither = neither + nrow(fire_df %>% filter(PS00756_H1 == 0, PS00756_H2 == 0))
# h1_only = h1_only + nrow(fire_df %>% filter(PS00756_H1 == 1, PS00756_H2 == 0))
# h2_only = h2_only + nrow(fire_df %>% filter(PS00756_H1 == 0, PS00756_H2 == 1))
# both = both + nrow(fire_df %>% filter(PS00756_H1 == 1, PS00756_H2 == 1))

# neither = neither + nrow(fire_df %>% filter(PS00757_H1 == 0, PS00757_H2 == 0))
# h1_only = h1_only + nrow(fire_df %>% filter(PS00757_H1 == 1, PS00757_H2 == 0))
# h2_only = h2_only + nrow(fire_df %>% filter(PS00757_H1 == 0, PS00757_H2 == 1))
# both = both + nrow(fire_df %>% filter(PS00757_H1 == 1, PS00757_H2 == 1))

# neither = neither + nrow(fire_df %>% filter(PS00758_H1 == 0, PS00758_H2 == 0))
# h1_only = h1_only + nrow(fire_df %>% filter(PS00758_H1 == 1, PS00758_H2 == 0))
# h2_only = h2_only + nrow(fire_df %>% filter(PS00758_H1 == 0, PS00758_H2 == 1))
# both = both + nrow(fire_df %>% filter(PS00758_H1 == 1, PS00758_H2 == 1))

# # (quick sanity check)
# A_table <- table((fire_df %>% filter(PS00718_nCov == 2) %>% select(PS00718_nCov, PS00718_nMSP))$PS00718_nMSP)


# hap_act <- matrix(c(neither, h1_only, h2_only, both), ncol=2, byrow=FALSE)
# colnames(hap_act) <- c('Closed','Actuated')
# rownames(hap_act) <- c('Closed','Actuated')
# hap_act <- as.table(hap_act)

# chi <- chisq.test(hap_act)

# chi$expected

# chi$residuals

#         Pearson's Chi-squared test with Yates' continuity correction

# data:  hap_act
# X-squared = 23133, df = 1, p-value < 2.2e-16

# > hap_act
#          Closed Actuated
# Closed   115246    42813
# Actuated  43326    57342
# > chi$expected

# chi$residuals
#            Closed Actuated
# Closed   96873.27 61185.73
# Actuated 61698.73 38969.27
#             Closed  Actuated
# Closed    59.02984 -74.27600
# Actuated -73.96657  93.07058



# # Bar plot by bin --> expected prop actuation on both haps --------------------------------------------------------------------------------------------------
# # expected both = prop^2, one = (prop*(1-prop))*2, neither = (1-prop)^2
# # vs ACTUAL fraction of neither, one, both

# # Group FIRE peaks by FS actuation bins --> calc total scDAF % actuation in each bin
# fire_df <- read_tsv('scDAF_FIRE_actuation_MSP150_withProximal.tsv') # promoter proximal info added with "sc_actuation_comparison.py"
# bin_names <- c('10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%')

# scDAF_bin_act <- fire_df %>% group_by(bin) %>% summarise(nCov = sum(PS00718_H1>=0,PS00718_H2>=0,PS00756_H1>=0,PS00756_H2>=0,PS00757_H1>=0,PS00757_H2>=0,PS00758_H1>=0,PS00758_H2>=0),
#                                         nMSP = sum(PS00718_H1>=1,PS00718_H2>=1,PS00756_H1>=1,PS00756_H2>=1,PS00757_H1>=1,PS00757_H2>=1,PS00758_H1>=1,PS00758_H2>=1))

# scDAF_bin_act$prop_act <- scDAF_bin_act$nMSP / scDAF_bin_act$nCov

# scDAF_bin_act$exp_neither <- (1-scDAF_bin_act$prop_act)**2
# scDAF_bin_act$exp_H1 <- scDAF_bin_act$prop_act * (1-scDAF_bin_act$prop_act)
# scDAF_bin_act$exp_H2 <- scDAF_bin_act$prop_act * (1-scDAF_bin_act$prop_act)
# scDAF_bin_act$exp_both <- scDAF_bin_act$prop_act**2

# scDAF_bin_act %>% select(exp_neither:exp_both) %>% rowSums()

# act_df <- read_csv('hap_actuation_bins.csv')

# actuation_by_bin <- act_df %>% group_by(bin) %>% summarise(neither=sum(neither), H1_only=sum(H1_only), H2_only=sum(H2_only), both=sum(both))

# scDAF_bin_act$neither <- actuation_by_bin$neither
# scDAF_bin_act$H1_only <- actuation_by_bin$H1_only
# scDAF_bin_act$H2_only <- actuation_by_bin$H2_only
# scDAF_bin_act$both <- actuation_by_bin$both

# scDAF_bin_act$total <- scDAF_bin_act %>% select(neither, H1_only, H2_only, both) %>% rowSums()
# scDAF_bin_act$pneither <- scDAF_bin_act$neither / scDAF_bin_act$total
# scDAF_bin_act$pH1_only <- scDAF_bin_act$H1_only / scDAF_bin_act$total
# scDAF_bin_act$pH2_only <- scDAF_bin_act$H2_only / scDAF_bin_act$total
# scDAF_bin_act$pboth <- scDAF_bin_act$both / scDAF_bin_act$total

# # create single, long dataframe for plotting

# hap_act_long <- scDAF_bin_act %>% select(bin,exp_neither)
# hap_act_long <- hap_act_long %>% rename(prop = exp_neither)
# hap_act_long$type <- 'Expected'
# hap_act_long$group <- 'Neither'

# temp <- scDAF_bin_act %>% select(bin,exp_H1)
# temp <- temp %>% rename(prop = exp_H1)
# temp$type <- 'Expected'
# temp$group <- 'H1'
# hap_act_long <- rbind(hap_act_long, temp)

# temp <- scDAF_bin_act %>% select(bin,exp_H2)
# temp <- temp %>% rename(prop = exp_H2)
# temp$type <- 'Expected'
# temp$group <- 'H2'
# hap_act_long <- rbind(hap_act_long, temp)

# temp <- scDAF_bin_act %>% select(bin,exp_both)
# temp <- temp %>% rename(prop = exp_both)
# temp$type <- 'Expected'
# temp$group <- 'Both'
# hap_act_long <- rbind(hap_act_long, temp)

# temp <- scDAF_bin_act %>% select(bin,pneither)
# temp <- temp %>% rename(prop = pneither)
# temp$type <- 'Observed'
# temp$group <- 'Neither'
# hap_act_long <- rbind(hap_act_long, temp)

# temp <- scDAF_bin_act %>% select(bin,pH1_only)
# temp <- temp %>% rename(prop = pH1_only)
# temp$type <- 'Observed'
# temp$group <- 'H1'
# hap_act_long <- rbind(hap_act_long, temp)

# temp <- scDAF_bin_act %>% select(bin,pH2_only)
# temp <- temp %>% rename(prop = pH2_only)
# temp$type <- 'Observed'
# temp$group <- 'H2'
# hap_act_long <- rbind(hap_act_long, temp)

# temp <- scDAF_bin_act %>% select(bin,pboth)
# temp <- temp %>% rename(prop = pboth)
# temp$type <- 'Observed'
# temp$group <- 'Both'
# hap_act_long <- rbind(hap_act_long, temp)

# scDAF_bin_act %>% select(exp_neither, pneither)
# scDAF_bin_act %>% select(exp_H1, pH1_only)
# scDAF_bin_act %>% select(exp_H2, pH2_only)
# scDAF_bin_act %>% select(exp_both, pboth)

# hap_act_long$bin <- factor(hap_act_long$bin)

# plot_act_obs_exp_bins <- ggplot(hap_act_long, aes(x=bin, y=prop, fill=group)) +
#                           geom_bar(stat = "identity") +
#                           facet_wrap(~type) +
#                           theme_bw()

# ggsave('figures/hap_actuation_by_bins.png', plot_act_obs_exp_bins, width=10, height=6)


