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


codep_data <- read_csv("napa_TF_footprinting_codep_chi_square.csv")
codep_data$chi2_stat_log10 <- log10(codep_data$chi2_stat)

# color by Bonferroni corrected p-value significance of <= 0.01
codep_data$bon_sig <- FALSE
codep_data[codep_data$adj_p_value <= 0.01,]$bon_sig <- TRUE
codep_data$highlight <- FALSE
codep_data[codep_data$adj_p_value == 0,]$highlight <- TRUE
codep_data$pair_name <- paste0(codep_data$reg1,"-",codep_data$reg2)


score_chiStat_sp <- ggscatter(codep_data, x = "score", y = "chi2_stat_log10", color = "bon_sig", palette = c("blue","red")) +
                        xlab("Codependency Score") + ylab("Chi-square Test Statistic (log10)") + xlim(-0.1,0.6) +
                        geom_label_repel(data = subset(codep_data, highlight == TRUE), box.padding = 0.4,
                            aes(label = pair_name))

my_ggsave("figures/score_by_chiStat.pdf", score_chiStat_sp)


