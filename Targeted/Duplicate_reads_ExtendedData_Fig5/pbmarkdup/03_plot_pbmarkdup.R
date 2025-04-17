library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(glue)
library(data.table)
library(tools)


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

dup_rates <- read_csv('duplicate_rates_pbmarkdub.csv')

dup_rates$sample <- factor(dup_rates$sample, levels = rev(c("COLO","UBA1","WASF1","napa","liver","heart","GM12878","colon")))

perc_by_tissue <- ggbarplot(dup_rates, x = "sample", y = "dup_rate", orientation = "horiz", fill = "sample", position = position_dodge(0.9)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
    theme(legend.position="none") +
    xlab("Sample") + ylab('Percent Duplicates') +
    scale_fill_brewer(palette="Set1")
    
my_ggsave('percent_dup_by_tissue_pbmarkdup.pdf', perc_by_tissue)

