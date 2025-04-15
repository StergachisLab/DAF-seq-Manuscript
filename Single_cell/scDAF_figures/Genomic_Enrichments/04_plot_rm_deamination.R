library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(cowplot)
library(colorspace)
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

grouped_rm <- read_csv("rm_grouped_deamination_counts.csv")

grouped_rm$total <- (grouped_rm$deam + grouped_rm$unmod)
grouped_rm$prop <- grouped_rm$deam / grouped_rm$total

# filter out "?" repeat classes
grouped_rm_filt <- grouped_rm[!grepl("\\?", grouped_rm$repeat_class), ]

prop_bar <- ggbarplot(grouped_rm_filt, "repeat_class", "prop", fill = "repeat_class", ylab = FALSE) +
                ylim(0,0.5) +
                theme(legend.position = "none") +
                rotate()

total_bar <- ggbarplot(grouped_rm_filt, "repeat_class", "total", fill = "repeat_class", ylab = FALSE) +
                theme(legend.position = "none") +
                rotate()

total_bar_log <- ggbarplot(grouped_rm_filt, "repeat_class", "total", fill = "repeat_class", ylab = FALSE) +
                theme(legend.position = "none") +
                yscale("log10", .format = TRUE) +
                rotate()

p_grid <- plot_grid(prop_bar, total_bar, total_bar_log, nrow=1)

ggsave("rm_group_grid.pdf", p_grid, height=6, width = 15)
