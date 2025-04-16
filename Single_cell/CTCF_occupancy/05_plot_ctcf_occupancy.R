library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
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


da_rate <- read_csv('da_rate_by_prop_bound.csv')

da_corr <- ggscatter(da_rate, x = "prop_da_both_strands", y = "prop_bound", add = "reg.line") +
            xlim(0.15,0.30) + ylim(0.25, 0.75) +
            stat_cor(label.x = 0.25, label.y = 0.5)
my_ggsave('figures/da_ctcf_bound_corr.pdf', da_corr)



# ChIA-PET scDAF CTCF co-occupnacy
chia_co_occ <- read_csv('ctcf_ChIA-PET_co-occupancy.csv')

no_na <- chia_co_occ %>% filter(prop_co_occ != 'NA')
no_na$looped <- FALSE
no_na[no_na$pet_score > 0,]$looped <- TRUE
no_na$dist <- 1

unlooped <- no_na %>% filter(looped == FALSE)
unlooped$mot_strands <- ''
unlooped[unlooped$strand_1 == '+' & unlooped$strand_2 == '+',]$mot_strands <- '++'
unlooped[unlooped$strand_1 == '+' & unlooped$strand_2 == '-',]$mot_strands <- '+-'
unlooped[unlooped$strand_1 == '-' & unlooped$strand_2 == '+',]$mot_strands <- '-+'
unlooped[unlooped$strand_1 == '-' & unlooped$strand_2 == '-',]$mot_strands <- '--'
unlooped$pet_bin = 10

unlooped_co_occ_motStrands <- ggboxplot(unlooped, "dist", "prop_co_occ", color = "mot_strands")
my_ggsave('figures/unlooped_co_occ_motStrands.pdf', unlooped_co_occ_motStrands)

# t-test of genome-wide co-occupancy vs within loop anchor regions
in_loops <- no_na %>% filter(pet_score > 0)
t.test(in_loops$prop_co_occ, no_na$prop_co_occ)



ft_codes <- read_csv('ctcf_footprint_codes.csv')

# percent aaccesible (both covered)
example_fibers1 <- ft_codes %>% filter(chrom == "chr1", start == 1122246, code >= 0) %>% select(cell, hap)
example_fibers2 <- ft_codes %>% filter(chrom == "chr1", start == 248807309, code >= 0) %>% select(cell, hap)

example_fibers1$hap_strand <- paste0(example_fibers1$cell, "_", example_fibers1$hap)
example_fibers2$hap_strand <- paste0(example_fibers2$cell, "_", example_fibers2$hap)
both_cov <- intersect(example_fibers1$hap_strand, example_fibers2$hap_strand) # 10 hap_strands

# percent co-occupied
example_fibers1 <- ft_codes %>% filter(chrom == "chr1", start == 1122246, code >= 1) %>% select(cell, hap)
example_fibers2 <- ft_codes %>% filter(chrom == "chr1", start == 248807309, code >= 1) %>% select(cell, hap)
example_fibers1$hap_strand <- paste0(example_fibers1$cell, "_", example_fibers1$hap)
example_fibers2$hap_strand <- paste0(example_fibers2$cell, "_", example_fibers2$hap)
both_acc <- intersect(example_fibers1$hap_strand, example_fibers2$hap_strand) # 6 hap_strands


# COMPARE BY LOOP REGIONS --------------------------------------------------------------------------------------------------------------------------------

vibrant_palette_four <- c("#CC79A7","#00AFBB", "#E7B800", "#FC4E07")

# co-occ as >= 1 co-occ by strand
loop_co_occ_any <- read_tsv('loop_coocc_by_strand_any.tsv')

loop_co_occ_any_uniq <- loop_co_occ_any %>% select(loop_ID, pet_score) %>% unique()
loop_co_occ_any_uniq_ranked <- loop_co_occ_any_uniq[order(-loop_co_occ_any_uniq$pet_score),]
loop_co_occ_any_uniq_ranked$rank <- 1:nrow(loop_co_occ_any_uniq_ranked)
loop_co_occ_any_ranked <- loop_co_occ_any[order(-loop_co_occ_any$pet_score),]
loop_co_occ_any_ranked$rank <- 0
for (i in 1:nrow(loop_co_occ_any_ranked)){
  loop_co_occ_any_ranked[i,]$rank <- loop_co_occ_any_uniq_ranked[loop_co_occ_any_uniq_ranked$loop_ID == loop_co_occ_any_ranked[i,]$loop_ID,]$rank
}

loop_co_occ_any_ranked$pet_bin <- 1
bin_cutoffs <- c(20, 100, 250, 500, 750, 1000, 2000, 4000)
bin_numbers <- 2:9
for (i in bin_numbers){
  loop_co_occ_any_ranked[loop_co_occ_any_ranked$rank > bin_cutoffs[i-1], ]$pet_bin <- i
}


# compare to NULL unlooped anchor regions as controls
null_anchors <- read_tsv('loop_coocc_by_strand_any_SHUFFLED.tsv')
null_anchors$looped <- FALSE
null_anchors$pet_bin <- 10

loop_co_occ_any_ranked$looped <- TRUE
loop_co_occ_any_ranked_plus_CTRL <- loop_co_occ_any_ranked %>% select(pet_bin, prop_co_occ, mot_strands, looped)
loop_co_occ_any_ranked_plus_CTRL <- rbind(loop_co_occ_any_ranked_plus_CTRL, null_anchors %>% select(pet_bin, prop_co_occ, mot_strands, looped))

loop_co_occ_pet_bin_any_motStrands_wCTRL <- ggboxplot(loop_co_occ_any_ranked_plus_CTRL %>% filter(mot_strands == "+-"), "pet_bin", "prop_co_occ", fill = "looped", palette = vibrant_palette_four) +
                                              scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
                                              xlab("ChIA-PET Score Bin") + ylab("Percentage of Loop Anchors Co-occupied") + ggtitle("CTCF +/- Orientation Co-occupancy within loop anchors") +
                                              theme_bw()

my_ggsave('figures/loop_co_occ_pet_bin_any_motStrands_wCTRL.pdf', loop_co_occ_pet_bin_any_motStrands_wCTRL)


# Boxplot of co-occ by bin --> only show +-, separate boxplot of bin1 for each strand combo (maybe different bin size to overcome small n)

# Use larger bins for top 500 PET Scores
loop_co_occ_any_large_bin <- loop_co_occ_any_ranked %>% filter(rank <= 500)
loop_co_occ_any_large_bin$mot_strands <- factor(loop_co_occ_any_large_bin$mot_strands, levels = c("+-","--","++","-+"))

my_comparisons <- list( c("--","++"), c("+-", "--"), c("+-","++"), c("+-","-+") )

loop_co_occ_any_large_bin_by_strand_boxplot <- ggboxplot(loop_co_occ_any_large_bin, "mot_strands", "prop_co_occ", fill = "mot_strands", palette = vibrant_palette_four) +
                                              scale_y_continuous(labels = scales::percent, limits=c(0,1.5)) +
                                              stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(1, 1.1, 1.2, 1.3)) +
                                              xlab("CTCF Orientation") + ylab("Percentage of Loop Anchors Co-occupied") + ggtitle("Top 500 ChIA-PET Scores (Bins 1-4)") +
                                              theme_bw()
my_ggsave('figures/loop_co_occ_any_large_bin_by_strand_boxplot.pdf', loop_co_occ_any_large_bin_by_strand_boxplot)


# Composition of bins
bin_composition <- loop_co_occ_any_ranked_plus_CTRL %>% group_by(pet_bin, mot_strands) %>%
                                                        summarise(Count = n(), .groups = "drop") %>%
                                                        group_by(pet_bin) %>%
                                                        mutate(bin_total = sum(Count))
bin_composition$bin_prop <- bin_composition$Count / bin_composition$bin_total
bin_composition$mot_strands <- factor(bin_composition$mot_strands, levels = c("-+","++","--","+-"))

bin_composition_barplot <- ggbarplot(bin_composition, "pet_bin", "bin_prop", fill = "mot_strands", color = "mot_strands", palette = rev(vibrant_palette_four)) +
                                scale_y_continuous(labels = scales::percent) +
                                xlab("ChIA-PET Score Bin") + ylab("Strand orientation by bin") +
                                theme_bw()
my_ggsave('figures/bin_composition_barplot.pdf', bin_composition_barplot)


# label bins on ranked ChIA-PET scores
loop_co_occ_any_uniq_ranked$pet_bin <- 1
for (i in bin_numbers){
  loop_co_occ_any_uniq_ranked[loop_co_occ_any_uniq_ranked$rank > bin_cutoffs[i-1], ]$pet_bin <- i
}
loop_co_occ_any_uniq_ranked$pet_bin <- factor(loop_co_occ_any_uniq_ranked$pet_bin)

pet_rank_Boxplot_by_loop_any_markedBins <- ggboxplot(loop_co_occ_any_uniq_ranked, "pet_bin", "pet_score", fill = "pet_bin", palette = "Set1") +
  xlab("ChIA-PET Score Bin") + ylab("ChIA-PET Score") +
  theme_bw()
my_ggsave("figures/pet_rank_Boxplot_by_loop_any_markedBins.pdf", pet_rank_Boxplot_by_loop_any_markedBins)


