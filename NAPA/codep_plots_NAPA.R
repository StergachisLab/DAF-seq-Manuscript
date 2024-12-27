library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)


scores <- read_csv('codep_scores_pairs.csv')
scores <- scores[order(scores$score),]
scores$rank <- as.integer(rownames(scores))
scores$name <- str_c(scores$reg1,":",scores$reg2)
nfilt = 3
scores_filt <- scores %>% filter(rank <= nfilt | rank > (nrow(scores)-nfilt))

score_point <- ggplot(scores, aes(x=reorder(rank, score), y=score)) +
    geom_point(size = 1.5, color = "red") +
    xlab('Peak Pair') +
    ylab('Codependency Score') +
    theme_bw(20) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x = element_blank()) +
    ylim(-0.2, 0.6) +
    geom_label_repel(data = scores_filt, aes(label = name), size = 3)

ggsave('figures/codep_score_point.png', score_point, width = 10, height = 5)


# Upset plot of footprint combos -------------
library(UpSetR)

z_df <- read_csv('zmw_footprint_regions.csv')
z_df[z_df == -1] <- 0 # replace -1 with 0 for upset
z_df <- as.data.frame(z_df)
z_df$"11" <- NULL # removing CTCF for upset plot

upset_plot <- upset(z_df, 
      nintersects = 40, 
      nsets = 15, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
      )

# can't use ggsave for upset plot
pdf(file='figures/upset_footprint_combos.pdf')
upset_plot
dev.off()


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
ggsave("figures/conditional_codep_raw_bar_1-2-3.png", bar_raw)
ggsave("figures/conditional_codep_raw_bar_1-2-3.pdf", bar_raw)

df_raw <- read_tsv("conditional_codependency_graphs_raw_1-2-3-4-5-6-7-8-9-10-11.tsv")
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
#   scale_fill_brewer(palette="Dark2") +
  xlab("Peak Omitted from Graph") +
  ylab("Codependency Ratio") +
  theme(legend.position="none")
ggsave("figures/conditional_codep_raw_bar_all.png", bar_raw)
ggsave("figures/conditional_codep_raw_bar_all.pdf", bar_raw)

