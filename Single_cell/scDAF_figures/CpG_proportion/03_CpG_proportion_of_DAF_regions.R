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


cpg_daf_props <- read_csv('cpg_prop_of_c_in_DAF_regions.csv')

cpg_bar <- ggbarplot(cpg_daf_props, "category", "prop",
  fill = "CpG", color = "CpG", palette = "Paired",
  label = TRUE, lab.nb.digits = 2,
  position = position_dodge(0.9)) +
  xlab('')

my_ggsave('cpg_prop_of_c_in_DAF_regions.pdf', cpg_bar)

