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

# Lineplots of binned codep ---------------------------------------------------------

# Same haplotype
avg_scores <- read_csv('bin_codep_DAF_exp_avg_score.csv')
all_scores_df <- read_csv('codep_DAF_exp_df.csv')

# Opposite Haplotype codependency
avg_scores_Opp <- read_csv('bin_codep_DAF_exp_avg_score_OPPOSITE_HAP.csv')
all_scores_Opp_df <- read_csv('codep_DAF_exp_OPPOSITE_HAP_df.csv')


# 1) Codep Same Hap vs. Opposite Hap
same_vs_opp_df1 <- all_scores_df
same_vs_opp_df1$hap <- "Same"
same_vs_opp_df1$dist_bin_log <- round(log2(same_vs_opp_df1$dist))

same_vs_opp_df2 <- all_scores_Opp_df
same_vs_opp_df2$hap <- "Opposite"
same_vs_opp_df2$dist_bin_log <- round(log2(same_vs_opp_df2$dist))
same_vs_opp_df <- rbind(same_vs_opp_df1, same_vs_opp_df2)


# log2 binning --------------------------------------------------------------------------------------------
same_vs_opp_df$bin_log2 <- floor(log2(same_vs_opp_df$dist))
same_vs_opp_df[same_vs_opp_df$bin_log2 < 8, ]$bin_log2 <- 8 # < 512 bp in same bin

same_bin2_avg_scores_log <- same_vs_opp_df %>% filter(hap == "Same") %>% group_by(bin_log2) %>% summarize(bin_mean = mean(score))
opp_bin2_avg_scores_log <- same_vs_opp_df %>% filter(hap == "Opposite") %>% group_by(bin_log2) %>% summarize(bin_mean = mean(score))
same_bin2_avg_scores_log$hap <- "Same"
opp_bin2_avg_scores_log$hap <- "Opposite"
same_opp_bin_idx_mean_log2 <- rbind(same_bin2_avg_scores_log, opp_bin2_avg_scores_log)

same_vs_opp_avg_codep_line_log2 <- ggline(same_opp_bin_idx_mean_log2, "bin_log2", "bin_mean",
   color = "hap", palette = c("#00AFBB", "#E7B800")) +
   ylim(0,0.04) +
   scale_x_discrete(labels = c("<512","> 1 kb","> 5 kb", "> 16 kb", "> 66 kb", "> 262 kb", "> 1 Mb", "> 4 Mb", "> 17 Mb", "> 67 Mb", "> 134 Mb"), breaks = c(8,10,12,14,16,18,20,22,24,26,27)) +
   xlab("Distance bins (bp between peaks)") + ylab("Mean binned codependency scores") +
   theme_bw(10) +
   theme(legend.position = "none")


same_bin2_avg_scores_log$diff <- same_bin2_avg_scores_log$bin_mean - opp_bin2_avg_scores_log$bin_mean
same_vs_opp_codep_DIFF_line_log2 <- ggline(same_bin2_avg_scores_log, "bin_log2", "diff") +
   ylim(-0.005,0.03) +
   scale_x_discrete(labels = c("<512","> 1 kb","> 5 kb", "> 16 kb", "> 66 kb", "> 262 kb", "> 1 Mb", "> 4 Mb", "> 17 Mb", "> 67 Mb", "> 134 Mb"), breaks = c(8,10,12,14,16,18,20,22,24,26,27)) +
   xlab("Distance bins (bp between peaks)") + ylab("Mean binned codependency scores (Same - Opposite Haplotype)") +
   theme_bw(10)


codep_grid_log2 <- plot_grid(same_vs_opp_avg_codep_line_log2, same_vs_opp_codep_DIFF_line_log2, ncol = 1)
ggsave("figures/codep_grid_log2.pdf", codep_grid_log2, height = 8, width = 10)


# T-test between haplotypes for each distance bin
bin_names <- sort(unique(same_vs_opp_df$bin_log2))

p_values <- c()
t_test_n <- c()
test_stats <- c()
t_deg_freedon <- c()
CI_1 <- c()
CI_2 <- c()
for (b in bin_names){
   same_vals <- same_vs_opp_df %>% filter(bin_log2 == b, hap == "Same") %>% select(score)
   opp_vals <- same_vs_opp_df %>% filter(bin_log2 == b, hap == "Opposite") %>% select(score)
   t <- t.test(same_vals, opp_vals, alternative = "greater")
   n <- nrow(same_vals) + nrow(opp_vals)
   t_stat <- t$statistic
   tdf <- t$parameter
   CI1 <- t$conf.int[1]
   CI2 <- t$conf.int[2]
   p_values <- c(p_values, t$p.value)
   t_test_n <- c(t_test_n, n)
   test_stats <- c(test_stats, t_stat)
   t_deg_freedon <- c(t_deg_freedon, tdf)
   CI_1 <- c(CI_1, CI1)
   CI_2 <- c(CI_2, CI2)
}

codep_p_values <- data.frame(bin = bin_names, p_values = p_values, n = t_test_n, t_statistic = test_stats, degrees_of_freedom = t_deg_freedon, confidence_interval_95_lower = CI_1, confidence_interval_95_upper = CI_2)
write_csv(codep_p_values, "msp_codependency_by_hap_t-test.csv")


# Compare between different haplotypes & cells to generate a NULL codependency dataset --------------------------------------------------------------------------
diffCells_scores_df <- read_csv('codep_DAF_exp_DIFF_CELLS_df.csv.gz')

same_vs_opp_df3 <- diffCells_scores_df
same_vs_opp_df3$hap <- "DiffCells"
same_vs_opp_df3$dist_bin_log <- round(log2(same_vs_opp_df3$dist))
same_vs_opp_vs_diffCells_df <- rbind(same_vs_opp_df1, same_vs_opp_df2, same_vs_opp_df3)


same_vs_opp_vs_diffCells_df$bin_log2 <- floor(log2(same_vs_opp_vs_diffCells_df$dist))
same_vs_opp_vs_diffCells_df[same_vs_opp_vs_diffCells_df$bin_log2 < 8, ]$bin_log2 <- 8 # < 512 bp in same bin

same_bin2_avg_scores_log <- same_vs_opp_vs_diffCells_df %>% filter(hap == "Same") %>% group_by(bin_log2) %>% summarize(bin_mean = mean(score))
opp_bin2_avg_scores_log <- same_vs_opp_vs_diffCells_df %>% filter(hap == "Opposite") %>% group_by(bin_log2) %>% summarize(bin_mean = mean(score))
diffCells_bin2_avg_scores_log <- same_vs_opp_vs_diffCells_df %>% filter(hap == "DiffCells") %>% group_by(bin_log2) %>% summarize(bin_mean = mean(score))
same_bin2_avg_scores_log$hap <- "Same Fiber"
opp_bin2_avg_scores_log$hap <- "Opposite Hap"
diffCells_bin2_avg_scores_log$hap <- "Different Cells"
same_opp_diffCells_bin_idx_mean_log2 <- rbind(same_bin2_avg_scores_log, opp_bin2_avg_scores_log, diffCells_bin2_avg_scores_log)


same_vs_opp_vs_diffCells_avg_codep_line_log2 <- ggline(same_opp_diffCells_bin_idx_mean_log2, "bin_log2", "bin_mean",
   color = "hap", palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
   ylim(-0.01,0.04) +
   scale_x_discrete(labels = c("<512","> 1 kb","> 5 kb", "> 16 kb", "> 66 kb", "> 262 kb", "> 1 Mb", "> 4 Mb", "> 17 Mb", "> 67 Mb", "> 134 Mb"), breaks = c(8,10,12,14,16,18,20,22,24,26,27)) +
   xlab("Distance bins (bp between peaks)") + ylab("Mean binned codependency scores\n") +
   theme_bw(12) +
   theme(legend.position = "top")

same_bin2_avg_scores_log$diff_same_diffCells <- same_bin2_avg_scores_log$bin_mean - diffCells_bin2_avg_scores_log$bin_mean
same_vs_diffCells_codep_DIFF_line_log2 <- ggline(same_bin2_avg_scores_log, "bin_log2", "diff_same_diffCells") +
   ylim(0,0.04) +
   scale_x_discrete(labels = c("<512","> 1 kb","> 5 kb", "> 16 kb", "> 66 kb", "> 262 kb", "> 1 Mb", "> 4 Mb", "> 17 Mb", "> 67 Mb", "> 134 Mb"), breaks = c(8,10,12,14,16,18,20,22,24,26,27)) +
   xlab("Distance bins (bp between peaks)") + ylab("Mean binned codependency scores\n(Same Hap & Cell - Different Cells)") +
   theme_bw(12)

opp_bin2_avg_scores_log$diff_opp_diffCells <- opp_bin2_avg_scores_log$bin_mean - diffCells_bin2_avg_scores_log$bin_mean
opp_vs_diffCells_codep_DIFF_line_log2 <- ggline(opp_bin2_avg_scores_log, "bin_log2", "diff_opp_diffCells") +
   ylim(0,0.04) +
   scale_x_discrete(labels = c("<512","> 1 kb","> 5 kb", "> 16 kb", "> 66 kb", "> 262 kb", "> 1 Mb", "> 4 Mb", "> 17 Mb", "> 67 Mb", "> 134 Mb"), breaks = c(8,10,12,14,16,18,20,22,24,26,27)) +
   xlab("Distance bins (bp between peaks)") + ylab("Mean binned codependency scores\n(Opposite Haplotype within Cell - Different Cells)") +
   theme_bw(12)


codep_grid_log2_with_DiffCells <- plot_grid(same_vs_opp_vs_diffCells_avg_codep_line_log2, same_vs_diffCells_codep_DIFF_line_log2, opp_vs_diffCells_codep_DIFF_line_log2, ncol = 1)
ggsave("figures/codep_grid_log2_with_DiffCells.pdf", codep_grid_log2_with_DiffCells, height = 12, width = 10)

