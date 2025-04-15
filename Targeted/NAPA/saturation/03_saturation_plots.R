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


# Condition Order: ['position', '4 uM 10 min', '1 uM 10 min', '0.25 uM 10 min', '4 uM 20 min', '1 uM 20 min', '0.25 uM 20 min']
# "PS00626","PS00627","PS00628","PS00629","PS00630","PS00631"

non_footprint_napa <- read_tsv("non_ft_da_props_NAPA.tsv")
ctcf_mod2_napa <- read_tsv("ctcf_mod2_da_props_NAPA.tsv")

wasf1_peak1 <- read_tsv("non_ft_da_props_WASF1_peak1.tsv")
wasf1_peak2 <- read_tsv("non_ft_da_props_WASF1_peak2.tsv")
wasf1_nuc <- read_tsv("non_ft_da_props_WASF1_nuc.tsv")


# convert wide to long
non_footprint_napa_long <- non_footprint_napa %>% pivot_longer(cols = `PS00626`:`PS00631`, names_to = "treatment", values_to = "proportion")
non_footprint_napa_long$treatment <- as.factor(non_footprint_napa_long$treatment)
ctcf_mod2_napa <- ctcf_mod2_napa %>% pivot_longer(cols = `PS00626`:`PS00631`, names_to = "treatment", values_to = "proportion")
ctcf_mod2_napa$treatment <- as.factor(ctcf_mod2_napa$treatment)

wasf1_peak1 <- wasf1_peak1 %>% pivot_longer(cols = `PS00626`:`PS00631`, names_to = "treatment", values_to = "proportion")
wasf1_peak1$treatment <- as.factor(wasf1_peak1$treatment)
wasf1_peak2 <- wasf1_peak2 %>% pivot_longer(cols = `PS00626`:`PS00631`, names_to = "treatment", values_to = "proportion")
wasf1_peak2$treatment <- as.factor(wasf1_peak2$treatment)
wasf1_nuc <- wasf1_nuc %>% pivot_longer(cols = `PS00626`:`PS00631`, names_to = "treatment", values_to = "proportion")
wasf1_nuc$treatment <- as.factor(wasf1_nuc$treatment)

limits <- c("PS00629","PS00626","PS00630","PS00627","PS00631","PS00628")
labels <- c("4 uM\n20 min","4 uM\n10 min","1 uM\n20 min","1 uM\n10 min","0.25 uM\n20 min","0.25 uM\n10 min")


#  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Line plot with DddA conc on the X and % deamination on the Y. One line per position. Separate by between footprints, CTCF, well-positioned nucleosome, random non-FIRE

# add median line for each group
non_footprint_napa_long_medians <- non_footprint_napa_long %>% group_by(treatment) %>% summarize(proportion = median(proportion))
non_footprint_napa_long_medians$position <- "median"

line_non_ft_napa <- ggplot(non_footprint_napa_long, aes(x = treatment, y = proportion, group = position)) +
  geom_line() +
  geom_point() +
  ylim(0,1) +
  theme_bw() +
  scale_x_discrete(limits = limits, label = labels) +
  xlab("DddA Treatment") + ylab("Proportion Deaminated") + ggtitle("NAPA Promoter - bases outside TF footprint regions")
line_non_ft_napa <- line_non_ft_napa + geom_line(data = non_footprint_napa_long_medians, aes(x = treatment, y = proportion), color = "blue", linewidth = 2)

my_ggsave("figures/lineplot_non_ft_da_by_treatment_napa.pdf", line_non_ft_napa)


# add median line for each group
ctcf_mod2_napa_medians <- ctcf_mod2_napa %>% group_by(treatment) %>% summarize(proportion = median(proportion))
ctcf_mod2_napa_medians$position <- "median"

line_ctcf_napa <- ggplot(ctcf_mod2_napa, aes(x = treatment, y = proportion, group = position)) +
  geom_line() +
  geom_point() +
  ylim(0,1) +
  theme_bw() +
  scale_x_discrete(limits = limits, label = labels) +
  xlab("DddA Treatment") + ylab("Proportion Deaminated") + ggtitle("NAPA Promoter - bases within footprint 11 (CTCF)")
line_ctcf_napa <- line_ctcf_napa + geom_line(data = ctcf_mod2_napa_medians, aes(x = treatment, y = proportion), color = "blue", linewidth = 1)

my_ggsave("figures/lineplot_ctcf_mod2_da_by_treatment_napa.pdf", line_ctcf_napa)


# add median line for each group
wasf1_p1_long_medians <- wasf1_peak1 %>% group_by(treatment) %>% summarize(proportion = median(proportion))
wasf1_p1_long_medians$position <- "median"

line_non_ft_wasf1_peak1 <- ggplot(wasf1_peak1, aes(x = treatment, y = proportion, group = position)) +
  geom_line() +
  geom_point() +
  ylim(0,1) +
  theme_bw() +
  scale_x_discrete(limits = limits, label = labels) +
  xlab("DddA Treatment") + ylab("Proportion Deaminated") + ggtitle("WASF1 Promoter")
line_non_ft_wasf1_peak1 <- line_non_ft_wasf1_peak1 + geom_line(data = wasf1_p1_long_medians, aes(x = treatment, y = proportion), color = "blue", linewidth = 2)

my_ggsave("figures/lineplot_non_ft_da_by_treatment_wasf1_peak1.pdf", line_non_ft_wasf1_peak1)


# add median line for each group
wasf1_p2_long_medians <- wasf1_peak2 %>% group_by(treatment) %>% summarize(proportion = median(proportion))
wasf1_p2_long_medians$position <- "median"

line_non_ft_wasf1_peak2 <- ggplot(wasf1_peak2, aes(x = treatment, y = proportion, group = position)) +
  geom_line() +
  geom_point() +
  ylim(0,1) +
  theme_bw() +
  scale_x_discrete(limits = limits, label = labels) +
  xlab("DddA Treatment") + ylab("Proportion Deaminated") + ggtitle("CD40 Promoter")
line_non_ft_wasf1_peak2 <- line_non_ft_wasf1_peak2 + geom_line(data = wasf1_p2_long_medians, aes(x = treatment, y = proportion), color = "blue", linewidth = 2)

my_ggsave("figures/lineplot_non_ft_da_by_treatment_wasf1_peak2.pdf", line_non_ft_wasf1_peak2)


# add median line for each group
wasf1_nuc_medians <- wasf1_nuc %>% group_by(treatment) %>% summarize(proportion = median(proportion))
wasf1_nuc_medians$position <- "median"

line_non_ft_wasf1_nuc <- ggplot(wasf1_nuc, aes(x = treatment, y = proportion, group = position)) +
  geom_line() +
  geom_point() +
  ylim(0,1) +
  theme_bw() +
  scale_x_discrete(limits = limits, label = labels) +
  xlab("DddA Treatment") + ylab("Proportion Deaminated") + ggtitle("WASF1 well-positioned nucleosome core (center 20 bp)")
line_non_ft_wasf1_nuc <- line_non_ft_wasf1_nuc + geom_line(data = wasf1_nuc_medians, aes(x = treatment, y = proportion), color = "blue", linewidth = 1)

my_ggsave("figures/lineplot_non_ft_da_by_treatment_wasf1_nuc.pdf", line_non_ft_wasf1_nuc)


# Violin of % deamination by DddA conc
violin_non_ft_napa_TC <- ggviolin(non_footprint_napa_long, "treatment", "proportion", fill = "TC",
   add.params = list(fill = "white"), position = position_dodge(1)) +
   scale_x_discrete(limits = limits, label = labels) +
   xlab("DddA Treatment") + ylab("Proportion Deaminated") + ggtitle("NAPA Promoter - bases outside TF footprint regions")

my_ggsave("figures/violin_by_TC_non_ft_da_by_treatment_napa.pdf", violin_non_ft_napa_TC)

violin_non_ft_napa_cpg <- ggviolin(non_footprint_napa_long, "treatment", "proportion", fill = "CpG",
   add.params = list(fill = "white"), position = position_dodge(1)) +
   scale_x_discrete(limits = limits, label = labels) +
   xlab("DddA Treatment") + ylab("Proportion Deaminated") + ggtitle("NAPA Promoter - bases outside TF footprint regions")

my_ggsave("figures/violin_by_CpG_non_ft_da_by_treatment_napa.pdf", violin_non_ft_napa_cpg)


# highlight TC positions
line_non_ft_napa_TC <- ggplot(non_footprint_napa_long, aes(x = treatment, y = proportion, group = position)) +
  geom_line(color = "gray", alpha = 0.5) +
  geom_point(color = "gray") +
  ylim(0,1) +
  theme_bw() +
  scale_x_discrete(limits = limits, label = labels) +
  xlab("DddA Treatment") + ylab("Proportion Deamianted") + ggtitle("NAPA Promoter - bases outside TF footprint regions (TC in red)")
line_non_ft_napa_TC <- line_non_ft_napa_TC +
                        geom_line(data = non_footprint_napa_long %>% filter(TC == "TC"), aes(x = treatment, y = proportion), color = "red") +
                        geom_point(data = non_footprint_napa_long %>% filter(TC == "TC"), color = "red")
ggsave("figures/lineplot_TC_highlight_non_ft_da_by_treatment_napa.pdf", line_non_ft_napa_TC)


# highlight CpG positions
line_non_ft_napa_cpg <- ggplot(non_footprint_napa_long, aes(x = treatment, y = proportion, group = position)) +
  geom_line(color = "gray", alpha = 0.5) +
  geom_point(color = "gray") +
  ylim(0,1) +
  theme_bw() +
  scale_x_discrete(limits = limits, label = labels) +
  xlab("DddA Treatment") + ylab("Proportion Deamianted") + ggtitle("NAPA Promoter - bases outside TF footprint regions (CpG in blue)")
line_non_ft_napa_cpg <- line_non_ft_napa_cpg +
                        geom_line(data = non_footprint_napa_long %>% filter(CpG == "CpG"), aes(x = treatment, y = proportion), color = "blue") +
                        geom_point(data = non_footprint_napa_long %>% filter(TC == "CpG"), color = "blue")
ggsave("figures/lineplot_CpG_highlight_non_ft_da_by_treatment_napa.pdf", line_non_ft_napa_cpg)


# T-tests of deamination frequencies for TC and CpG ---------------------------------

# 4 uM for 10 min
TC_high <- non_footprint_napa_long %>% filter(TC == "TC", treatment == "PS00626")
non_TC_high <- non_footprint_napa_long %>% filter(TC != "TC", treatment == "PS00626")
t.test(TC_high$proportion, non_TC_high$proportion)

TC_low <- non_footprint_napa_long %>% filter(TC == "TC", treatment == "PS00628")
non_TC_low <- non_footprint_napa_long %>% filter(TC != "TC", treatment == "PS00628")
t.test(TC_low$proportion, non_TC_low$proportion)

CpG_high <- non_footprint_napa_long %>% filter(CpG == "CpG", treatment == "PS00626")
non_CpG_high <- non_footprint_napa_long %>% filter(CpG != "CpG", treatment == "PS00626")
t.test(CpG_high$proportion, non_CpG_high$proportion)

CpG_low <- non_footprint_napa_long %>% filter(CpG == "CpG", treatment == "PS00628")
non_CpG_low <- non_footprint_napa_long %>% filter(CpG != "CpG", treatment == "PS00628")
t.test(CpG_low$proportion, non_CpG_low$proportion)


# t-test of 4 uM 10 vs 20 minutes -------------------------------------------------

# CTCF footprint

# 4 uM
ctcf_ten <- ctcf_mod2_napa %>% filter(treatment == "PS00626")
ctcf_twenty <- ctcf_mod2_napa %>% filter(treatment == "PS00629")
t.test(ctcf_ten$proportion, ctcf_twenty$proportion)

# 1 uM
ctcf_ten <- ctcf_mod2_napa %>% filter(treatment == "PS00627")
ctcf_twenty <- ctcf_mod2_napa %>% filter(treatment == "PS00630")
t.test(ctcf_ten$proportion, ctcf_twenty$proportion)

# 0.25 uM
ctcf_ten <- ctcf_mod2_napa %>% filter(treatment == "PS00628")
ctcf_twenty <- ctcf_mod2_napa %>% filter(treatment == "PS00631")
t.test(ctcf_ten$proportion, ctcf_twenty$proportion)


# 4 uM vs 1 uM
ctcf_high <- ctcf_mod2_napa %>% filter(treatment == "PS00626")
ctcf_low <- ctcf_mod2_napa %>% filter(treatment == "PS00627")
t.test(ctcf_high$proportion, ctcf_low$proportion)

# 4 uM vs 0.25 uM
ctcf_high <- ctcf_mod2_napa %>% filter(treatment == "PS00626")
ctcf_low <- ctcf_mod2_napa %>% filter(treatment == "PS00628")
t.test(ctcf_high$proportion, ctcf_low$proportion)

# 1 uM vs 0.25 uM
ctcf_high <- ctcf_mod2_napa %>% filter(treatment == "PS00627")
ctcf_low <- ctcf_mod2_napa %>% filter(treatment == "PS00628")
t.test(ctcf_high$proportion, ctcf_low$proportion)

# Condition Order: ['position', '4 uM 10 min', '1 uM 10 min', '0.25 uM 10 min', '4 uM 20 min', '1 uM 20 min', '0.25 uM 20 min']
# "PS00626","PS00627","PS00628","PS00629","PS00630","PS00631"




# well-positioned nucleosome
nuc_ten <- wasf1_nuc %>% filter(treatment == "PS00626")
nuc_twenty <- wasf1_nuc %>% filter(treatment == "PS00629")
t.test(nuc_ten$proportion, nuc_twenty$proportion)

# well-positioned nucleosome
napa_ten <- non_footprint_napa_long %>% filter(treatment == "PS00626")
napa_twenty <- non_footprint_napa_long %>% filter(treatment == "PS00629")
t.test(napa_ten$proportion, napa_twenty$proportion)


