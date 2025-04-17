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


run_dens <- read_tsv("run_of_3_da_occurences_all.tsv", col_names=c('chrom','pos','cov'))

# 47515062 is only on the top strand
run_dens <- run_dens %>% filter(pos > 47515062)

plot <- ggplot(run_dens, aes(x = pos+1, y = cov)) +
            geom_line() +
            theme_bw() +
            xlab("NAPA Promoter") + ylab("Count of 3 straight non-deaminated bases")

ggsave("run_of_3_lineplot.png", plot, width=8, height=6)


# highlight 
footprint_regions <- read_tsv("merged_ft_on_both_strands.bed", col_names=c('chrom','start','end'))


