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

FS_props <- read_csv("saturation_Fiber-seq_NAPA_m6a_ref_coord_props.csv")
FS_props$assay <- "Fiber-seq"

# 4 uM for 10 minutes
DAF_props <- read_tsv("/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/qc_plots/prop_da_titration_NAPA.tsv")
ft_pos <- read_csv('napa_footprint_positions.txt', col_names = c('position'))
NAPA_prom_DAF <- DAF_props %>% select(position, PS00626) %>% filter(position >= 47515063, position <= 47515660, !position %in% ft_pos$position) %>% rename(ref_coord = position, m6a_prop = PS00626)
NAPA_prom_DAF$sample <- "Targeted DAF-seq (GM12878)"
NAPA_prom_DAF$assay <- "DAF-seq"

merged <- rbind(FS_props, NAPA_prom_DAF)
merged$sample <- factor(merged$sample, levels=c("Targeted DAF-seq (GM12878)","GM12878","HG002","CHM13","PS00338_COLO829BL_1","PS00356_COLO829BL_2","COLO_T_2_PS00_418_451_488",
                                                "PS00272","PS00321","PS00327","PS00381","PS00382","PS00383","PS00384","PS30743","ST001-lung","ST001-liver"))

merged_violin <- ggboxplot(merged, "sample", "m6a_prop", fill = "assay") +
                        ylim(0,1) +
                        scale_y_continuous(labels = percent) +
                        theme_bw(8) +
                        theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                        ylab("Percentage bases modified") +
                        xlab("") +
                        ggtitle("Modification Percentage within NAPA Promoter (Outside footprint regions)")

my_ggsave('figures/saturation_NAPA_DAF_vs_FS.pdf', merged_violin, width = 8.5, height = 5)


# DAF-seq vs. Fiber-seq stat test
D <- merged %>% filter(sample == "Targeted DAF-seq (GM12878)")
F <- merged %>% filter(sample == "GM12878")

t.test(D$m6a_prop, F$m6a_prop, alternative = "greater")

# merged %>% group_by(sample) %>% summarize(mean = mean(m6a_prop), median = median(m6a_prop))


