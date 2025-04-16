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



# Upset plot of footprint combos -------------
library(UpSetR)

z_df <- read_csv('zmw_footprint_regions.csv')
z_df[z_df == -1] <- 0 # replace -1 with 0 for upset
z_df <- as.data.frame(z_df)
z_df$"11" <- NULL # removing CTCF for upset plot


z_small <- z_df %>% select('1','2','3')
upset_plot_small <- upset(z_small, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
      )

# can't use ggsave for upset plot
pdf(file='figures/upset_footprint_combos_small.pdf')
upset_plot_small
dev.off()


# average the RAW of each edge in each graph ------------------
df_raw <- read_tsv("conditional_codependency_graphs_raw_1-2-3.tsv")
df_raw <- df_raw %>% filter(peak_omitted != "baseline")
df_raw$peak_omitted <- as.factor(df_raw$peak_omitted)

bar_raw <- ggplot(data=df_raw, aes(x=reorder(peak_omitted, -codep_ratio), y=codep_ratio, fill=peak_omitted)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_brewer(palette="Dark2") +
  xlab("Peak Omitted from Graph") +
  ylab("Codependency Ratio") +
  theme(legend.position="none")
my_ggsave("figures/conditional_codep_raw_bar_1-2-3.pdf", bar_raw)
