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
library(scales)

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


deamination_df <- read_csv('deamination_summary.csv')

# prop deamination (merged bases)
deamin_prop_bar <- ggbarplot(deamination_df, "Region", "prop_da",
  fill = "group", palette = "Paired",
  label = TRUE, lab.pos = "out", lab.nb.digits = 3,
  position = position_dodge(0.9)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.5))

my_ggsave('figures/CpG_deamination_prop_barplot.pdf', deamin_prop_bar)


# Quantify control CpG methylation
NAPA_pos <- read_tsv('CpG_tools_out/PS00724_NAPA_CpG.combined.bed', col_names=c('chrom','start','end','mod_score','hap','coverage','mod_count','unmod_count','avg_mod_score','avg_unmod_score'))
UBA1_pos <- read_tsv('CpG_tools_out/PS00725_UBA1_CpG.combined.bed', col_names=c('chrom','start','end','mod_score','hap','coverage','mod_count','unmod_count','avg_mod_score','avg_unmod_score'))
NAPA_neg <- read_tsv('CpG_tools_out/PS00726_NAPA_NegCtrl.combined.bed', col_names=c('chrom','start','end','mod_score','hap','coverage','mod_count','unmod_count','avg_mod_score','avg_unmod_score'))
UBA1_neg <- read_tsv('CpG_tools_out/PS00727_UBA1_NegCtrl.combined.bed', col_names=c('chrom','start','end','mod_score','hap','coverage','mod_count','unmod_count','avg_mod_score','avg_unmod_score'))

df1 <- data.frame(score=NAPA_pos$mod_score)
df1$region <- 'NAPA'
df1$type <- 'Pos'
df2 <- data.frame(score=UBA1_pos$mod_score)
df2$region <- 'UBA1'
df2$type <- 'Pos'
df3 <- data.frame(score=NAPA_neg$mod_score)
df3$region <- 'NAPA'
df3$type <- 'Neg'
df4 <- data.frame(score=UBA1_neg$mod_score)
df4$region <- 'UBA1'
df4$type <- 'Neg'

cpg_cntl_df <- rbind(df1,df2,df3,df4)

cpg_meth_violins <- ggviolin(cpg_cntl_df, "region", "score",
                      fill = "type", add = "boxplot")

my_ggsave('figures/control_CpG_meth_violins.pdf', cpg_meth_violins)

