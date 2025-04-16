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

# Cluster composition by tissue ---------------------------------------------------
cluster_counts <- read_csv("cluster_tissue_counts_GA_FINAL.csv")
liver_prop <- as.data.frame(cluster_counts$cluster)
names(liver_prop) <- c('cluster')
liver_prop$prop <- cluster_counts$num_Liver / (cluster_counts$num_Liver+cluster_counts$num_GM12878)
liver_prop$tissue <- 'liver'
gm_prop <- as.data.frame(cluster_counts$cluster)
names(gm_prop) <- c('cluster')
gm_prop$prop <- cluster_counts$num_GM12878 / (cluster_counts$num_Liver+cluster_counts$num_GM12878)
gm_prop$tissue <- 'GM12878'
cluster_merged_GA <- rbind(liver_prop, gm_prop)

clust_bar_GA <- ggbarplot(cluster_merged_GA, "cluster", "prop",
  fill = "tissue", color = "tissue", palette = "Paired",
  label = FALSE, lab.col = "white", lab.pos = "in") +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1))

my_ggsave('figures/cluster_composition_GA.pdf', clust_bar_GA)

# By Haplotype
H1_counts <- read_csv("cluster_tissue_counts_GA_FINAL_H1.csv")
H2_counts <- read_csv("cluster_tissue_counts_GA_FINAL_H2.csv")
merged_hap_counts <- merge(H1_counts, H2_counts, by = "cluster")
merged_hap_counts$total <- merged_hap_counts %>% select(Liver_H1:GM12878_H2) %>% rowSums()
merged_hap_counts$H1L_prop <- merged_hap_counts$Liver_H1 / merged_hap_counts$total
merged_hap_counts$H2L_prop <- merged_hap_counts$Liver_H2 / merged_hap_counts$total
merged_hap_counts$H1G_prop <- merged_hap_counts$GM12878_H1 / merged_hap_counts$total
merged_hap_counts$H2G_prop <- merged_hap_counts$GM12878_H2 / merged_hap_counts$total

H1L <- as.data.frame(merged_hap_counts$cluster)
names(H1L) <- c('cluster')
H1L$prop <- merged_hap_counts$H1L_prop
H1L$tissue <- 'Liver_H1'
H2L <- as.data.frame(merged_hap_counts$cluster)
names(H2L) <- c('cluster')
H2L$prop <- merged_hap_counts$H2L_prop
H2L$tissue <- 'Liver_H2'
H1G <- as.data.frame(merged_hap_counts$cluster)
names(H1G) <- c('cluster')
H1G$prop <- merged_hap_counts$H1G_prop
H1G$tissue <- 'GM12878_H1'
H2G <- as.data.frame(merged_hap_counts$cluster)
names(H2G) <- c('cluster')
H2G$prop <- merged_hap_counts$H2G_prop
H2G$tissue <- 'GM12878_H2'
cluster_hap_merged_GA <- rbind(H1L, H2L, H1G, H2G)

clust_hap_bar_GA <- ggbarplot(cluster_hap_merged_GA, "cluster", "prop",
  fill = "tissue", color = "tissue",
  label = FALSE, lab.col = "white", lab.pos = "in") +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1))

my_ggsave('figures/cluster_composition_with_hap_GA.pdf', clust_hap_bar_GA)




# plot raw counts in each cluster by tissue & by hap

H1L <- as.data.frame(merged_hap_counts$cluster)
names(H1L) <- c('cluster')
H1L$count <- merged_hap_counts$Liver_H1
H1L$tissue <- 'Liver_H1'
H2L <- as.data.frame(merged_hap_counts$cluster)
names(H2L) <- c('cluster')
H2L$count <- merged_hap_counts$Liver_H2
H2L$tissue <- 'Liver_H2'
H1G <- as.data.frame(merged_hap_counts$cluster)
names(H1G) <- c('cluster')
H1G$count <- merged_hap_counts$GM12878_H1
H1G$tissue <- 'GM12878_H1'
H2G <- as.data.frame(merged_hap_counts$cluster)
names(H2G) <- c('cluster')
H2G$count <- merged_hap_counts$GM12878_H2
H2G$tissue <- 'GM12878_H2'
cluster_raw_counts_hap_merged_GA <- rbind(H1L, H2L, H1G, H2G)

clust_hap_raw_count_bar_GA <- ggbarplot(cluster_raw_counts_hap_merged_GA, "cluster", "count",
  fill = "tissue", color = "tissue",
  label = FALSE, lab.col = "white", lab.pos = "in")

my_ggsave('figures/cluster_counts_with_hap_GA.pdf', clust_hap_raw_count_bar_GA)


# Sub Clustering Cluster 6 --------------------------------------------

# By Haplotype
H1_counts <- read_csv("cluster6_tissue_counts_GA_FINAL_H1.csv")
H2_counts <- read_csv("cluster6_tissue_counts_GA_FINAL_H2.csv")
merged_hap_counts_sub <- merge(H1_counts, H2_counts, by = "cluster")
merged_hap_counts_sub$total <- merged_hap_counts_sub %>% select(Liver_H1:GM12878_H2) %>% rowSums()
merged_hap_counts_sub$H1L_prop <- merged_hap_counts_sub$Liver_H1 / merged_hap_counts_sub$total
merged_hap_counts_sub$H2L_prop <- merged_hap_counts_sub$Liver_H2 / merged_hap_counts_sub$total
merged_hap_counts_sub$H1G_prop <- merged_hap_counts_sub$GM12878_H1 / merged_hap_counts_sub$total
merged_hap_counts_sub$H2G_prop <- merged_hap_counts_sub$GM12878_H2 / merged_hap_counts_sub$total

H1L <- as.data.frame(merged_hap_counts_sub$cluster)
names(H1L) <- c('cluster')
H1L$prop <- merged_hap_counts_sub$H1L_prop
H1L$tissue <- 'Liver_H1'
H2L <- as.data.frame(merged_hap_counts_sub$cluster)
names(H2L) <- c('cluster')
H2L$prop <- merged_hap_counts_sub$H2L_prop
H2L$tissue <- 'Liver_H2'
H1G <- as.data.frame(merged_hap_counts_sub$cluster)
names(H1G) <- c('cluster')
H1G$prop <- merged_hap_counts_sub$H1G_prop
H1G$tissue <- 'GM12878_H1'
H2G <- as.data.frame(merged_hap_counts_sub$cluster)
names(H2G) <- c('cluster')
H2G$prop <- merged_hap_counts_sub$H2G_prop
H2G$tissue <- 'GM12878_H2'
cluster_hap_merged_GA <- rbind(H1L, H2L, H1G, H2G)

sub_clust_hap_bar_GA <- ggbarplot(cluster_hap_merged_GA, "cluster", "prop",
  fill = "tissue", color = "tissue",
  label = FALSE, lab.col = "white", lab.pos = "in") +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1))

my_ggsave('figures/cluster6_composition_with_hap_GA.pdf', sub_clust_hap_bar_GA)




# Sub Cluster 6 Liver Conditional Codependency ----------------------------------------------------------------------------
df_raw <- read_tsv("codependency_promoter/conditional_codependency_graphs_raw_clust6_Liver_H1_1-2-3-4-5.tsv")
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
ggsave("figures/conditional_codep_raw_bar_Clust6_Liver_H1_1-2-3-4-5.png", bar_raw)

df_raw <- read_tsv("codependency_promoter/conditional_codependency_graphs_raw_clust6_Liver_H2_1-2-3-4-5.tsv")
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
  xlab("Peak Omitted from Graph") +
  ylab("Codependency Ratio") +
  theme(legend.position="none")
ggsave("figures/conditional_codep_raw_bar_Clust6_Liver_H2_1-2-3-4-5.png", bar_raw)


# module usage plots -------------------------------------------

z_df_h1 <- read.csv('codependency_promoter/clust6_msp_df_Liver_H1.csv', row.names=1)
z_df_h2 <- read.csv('codependency_promoter/clust6_msp_df_Liver_H2.csv', row.names=1)

z_df_h1$rsum <- z_df_h1 %>% select('X1':'X5') %>% rowSums()
z_df_h2$rsum <- z_df_h2 %>% select('X1':'X5') %>% rowSums()

# save # of module actuated to csv
h1_table <- as.data.frame(table(z_df_h1$rsum))
h2_table <- as.data.frame(table(z_df_h2$rsum))
colnames(h1_table) <- c('n_acc','Count')
colnames(h2_table) <- c('n_acc','Count')
write_csv(h1_table, 'cluster6_numModulesBound_H1.csv')
write_csv(h2_table, 'cluster6_numModulesBound_H2.csv')


# H1
new_df_h1 <- data.frame(n_acc = numeric(), module = numeric(), count = numeric())
for (i in 1:5){
  sub_df <- z_df_h1 %>% filter(rsum == i) %>% select('X1':'X5')
  for (k in 1:ncol(sub_df)){
    count <- sum(sub_df[,k])
    new_df_h1 <- rbind(new_df_h1, c(i,k,count))
  }
}
names(new_df_h1) <- c('n_acc','module','count')
new_df_h1$n_acc <- factor(new_df_h1$n_acc)
new_df_h1$module <- factor(new_df_h1$module)

# H2
new_df_h2 <- data.frame(n_acc = numeric(), module = numeric(), count = numeric())
for (i in 1:5){
  sub_df <- z_df_h2 %>% filter(rsum == i) %>% select('X1':'X5')
  for (k in 1:ncol(sub_df)){
    count <- sum(sub_df[,k])
    new_df_h2 <- rbind(new_df_h2, c(i,k,count))
  }
}
names(new_df_h2) <- c('n_acc','module','count')
new_df_h2$n_acc <- factor(new_df_h2$n_acc)
new_df_h2$module <- factor(new_df_h2$module)


module_bar_h1 <- ggplot(data=new_df_h1, aes(x=module, y=count, fill=module)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_brewer(palette="Dark2") +
  xlab("Module") +
  ylab("Count") +
  ggtitle("Module occupancy by N modules bound -- Hap 1 (C)") +
  theme(legend.position="none") +
  facet_wrap(~n_acc, nrow=1)

module_bar_h2 <- ggplot(data=new_df_h2, aes(x=module, y=count, fill=module)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_brewer(palette="Dark2") +
  xlab("Module") +
  ylab("Count") +
  ggtitle("Module occupancy by N modules bound -- Hap 2 (T)") +
  theme(legend.position="none") +
  facet_wrap(~n_acc, nrow=1)


my_ggsave("figures/module_usage_bar_H1.pdf", module_bar_h1)
my_ggsave("figures/module_usage_bar_H2.pdf", module_bar_h2)


# stacked bar plot where modules add up to 100% -------------------------------------------
gsums <- c()
for (i in 1:nrow(new_df_h1)){
  temp <- new_df_h1 %>% filter(n_acc == new_df_h1[i,]$n_acc)
  gsums <- c(gsums, sum(temp$count))
}
new_df_h1$group_sum <- gsums
new_df_h1$prop_of_group <- new_df_h1$count / new_df_h1$group_sum

gsums <- c()
for (i in 1:nrow(new_df_h2)){
  temp <- new_df_h2 %>% filter(n_acc == new_df_h2[i,]$n_acc)
  gsums <- c(gsums, sum(temp$count))
}
new_df_h2$group_sum <- gsums
new_df_h2$prop_of_group <- new_df_h2$count / new_df_h2$group_sum

module_bar_h1 <- ggbarplot(new_df_h1, "n_acc", "prop_of_group",
  fill = "module", color = "module") +
  ggtitle("Module occupancy by N modules bound -- Hap 1 (C)")
my_ggsave("figures/module_usage_bar_stacked_H1.pdf", module_bar_h1)

module_bar_h2 <- ggbarplot(new_df_h2, "n_acc", "prop_of_group",
  fill = "module", color = "module") +
  ggtitle("Module occupancy by N modules bound -- Hap 2 (T)")
my_ggsave("figures/module_usage_bar_stacked_H2.pdf", module_bar_h2)


# 2 modules being on, break this out in terms of which two modules that is
# stacked bar plot of 1:2, 1:3, 1:4, 1:5, 2:3, 2:4, 2:5, 3:4, 3:5, 4:5

mod2_h1_on <- z_df_h1 %>% filter(rsum == 2)
mod2_h1_on$first <- 0
mod2_h1_on$second <- 0
for (i in 1:nrow(mod2_h1_on)){
  acc <- c()
  for (k in 1:5){
    if (mod2_h1_on[i,k] == 1){
      acc <- c(acc,k)
    }
  }
  mod2_h1_on[i,]$first <- acc[1]
  mod2_h1_on[i,]$second <- acc[2]
}
mod2_h1_on$combo <- paste0(mod2_h1_on$first,':',mod2_h1_on$second)

mod3_h1_on <- z_df_h1 %>% filter(rsum == 3)
mod3_h1_on$first <- 0
mod3_h1_on$second <- 0
mod3_h1_on$third <- 0
for (i in 1:nrow(mod3_h1_on)){
  acc <- c()
  for (k in 1:5){
    if (mod3_h1_on[i,k] == 1){
      acc <- c(acc,k)
    }
  }
  mod3_h1_on[i,]$first <- acc[1]
  mod3_h1_on[i,]$second <- acc[2]
  mod3_h1_on[i,]$third <- acc[3]
}
mod3_h1_on$combo <- paste0(mod3_h1_on$first,':',mod3_h1_on$second,':',mod3_h1_on$third)

mod4_h1_on <- z_df_h1 %>% filter(rsum == 4)
mod4_h1_on$first <- 0
mod4_h1_on$second <- 0
mod4_h1_on$third <- 0
mod4_h1_on$fourth <- 0
for (i in 1:nrow(mod4_h1_on)){
  acc <- c()
  for (k in 1:5){
    if (mod4_h1_on[i,k] == 1){
      acc <- c(acc,k)
    }
  }
  mod4_h1_on[i,]$first <- acc[1]
  mod4_h1_on[i,]$second <- acc[2]
  mod4_h1_on[i,]$third <- acc[3]
  mod4_h1_on[i,]$fourth <- acc[4]
}
mod4_h1_on$combo <- paste0(mod4_h1_on$first,':',mod4_h1_on$second,':',mod4_h1_on$third,':',mod4_h1_on$fourth)

# H2
mod2_h2_on <- z_df_h2 %>% filter(rsum == 2)
mod2_h2_on$first <- 0
mod2_h2_on$second <- 0
for (i in 1:nrow(mod2_h2_on)){
  acc <- c()
  for (k in 1:5){
    if (mod2_h2_on[i,k] == 1){
      acc <- c(acc,k)
    }
  }
  mod2_h2_on[i,]$first <- acc[1]
  mod2_h2_on[i,]$second <- acc[2]
}
mod2_h2_on$combo <- paste0(mod2_h2_on$first,':',mod2_h2_on$second)

mod3_h2_on <- z_df_h2 %>% filter(rsum == 3)
mod3_h2_on$first <- 0
mod3_h2_on$second <- 0
mod3_h2_on$third <- 0
for (i in 1:nrow(mod3_h2_on)){
  acc <- c()
  for (k in 1:5){
    if (mod3_h2_on[i,k] == 1){
      acc <- c(acc,k)
    }
  }
  mod3_h2_on[i,]$first <- acc[1]
  mod3_h2_on[i,]$second <- acc[2]
  mod3_h2_on[i,]$third <- acc[3]
}
mod3_h2_on$combo <- paste0(mod3_h2_on$first,':',mod3_h2_on$second,':',mod3_h2_on$third)

mod4_h2_on <- z_df_h2 %>% filter(rsum == 4)
mod4_h2_on$first <- 0
mod4_h2_on$second <- 0
mod4_h2_on$third <- 0
mod4_h2_on$fourth <- 0
for (i in 1:nrow(mod4_h2_on)){
  acc <- c()
  for (k in 1:5){
    if (mod4_h2_on[i,k] == 1){
      acc <- c(acc,k)
    }
  }
  mod4_h2_on[i,]$first <- acc[1]
  mod4_h2_on[i,]$second <- acc[2]
  mod4_h2_on[i,]$third <- acc[3]
  mod4_h2_on[i,]$fourth <- acc[4]
}
mod4_h2_on$combo <- paste0(mod4_h2_on$first,':',mod4_h2_on$second,':',mod4_h2_on$third,':',mod4_h2_on$fourth)


mod2_on_h1_groups <- data.frame(table(mod2_h1_on$combo))
names(mod2_on_h1_groups) <- c('group','count')
mod2_on_h1_groups$prop <- mod2_on_h1_groups$count / sum(mod2_on_h1_groups$count)
mod2_on_h1_groups$n_acc <- 2

mod3_on_h1_groups <- data.frame(table(mod3_h1_on$combo))
names(mod3_on_h1_groups) <- c('group','count')
mod3_on_h1_groups$prop <- mod3_on_h1_groups$count / sum(mod3_on_h1_groups$count)
mod3_on_h1_groups$n_acc <- 3

mod4_on_h1_groups <- data.frame(table(mod4_h1_on$combo))
names(mod4_on_h1_groups) <- c('group','count')
mod4_on_h1_groups$prop <- mod4_on_h1_groups$count / sum(mod4_on_h1_groups$count)
mod4_on_h1_groups$n_acc <- 4

mod3_on_h1_groups_merged <- rbind(mod2_on_h1_groups, mod3_on_h1_groups, mod4_on_h1_groups)

mod2_on_h2_groups <- data.frame(table(mod2_h2_on$combo))
names(mod2_on_h2_groups) <- c('group','count')
mod2_on_h2_groups$prop <- mod2_on_h2_groups$count / sum(mod2_on_h2_groups$count)
mod2_on_h2_groups$n_acc <- 2

mod3_on_h2_groups <- data.frame(table(mod3_h2_on$combo))
names(mod3_on_h2_groups) <- c('group','count')
mod3_on_h2_groups$prop <- mod3_on_h2_groups$count / sum(mod3_on_h2_groups$count)
mod3_on_h2_groups$n_acc <- 3

mod4_on_h2_groups <- data.frame(table(mod4_h2_on$combo))
names(mod4_on_h2_groups) <- c('group','count')
mod4_on_h2_groups$prop <- mod4_on_h2_groups$count / sum(mod4_on_h2_groups$count)
mod4_on_h2_groups$n_acc <- 4

mod3_on_h2_groups_merged <- rbind(mod2_on_h2_groups, mod3_on_h2_groups, mod4_on_h2_groups)

mod2_3_on_bar_h1 <- ggbarplot(mod3_on_h1_groups_merged, "n_acc", "prop",
  fill = "group", color = "group") +
  ggtitle("Module occupancy by N modules bound -- Hap 1 (C)")

mod2_3_on_bar_h2 <- ggbarplot(mod3_on_h2_groups_merged, "n_acc", "prop",
  fill = "group", color = "group") +
  ggtitle("Module occupancy by N modules bound -- Hap 2 (T)")


my_ggsave("figures/two_three_modules_on_bar_h1.pdf", mod2_3_on_bar_h1)
my_ggsave("figures/two_three_modules_on_bar_h2.pdf", mod2_3_on_bar_h2)



mod3_on_h1_groups_merged$hap <- "H1"
mod3_on_h2_groups_merged$hap <- "H2"
mod23_merged_all <- rbind(mod3_on_h1_groups_merged, mod3_on_h2_groups_merged)


mod2_3_on_bar_ALL <- ggbarplot(mod23_merged_all, "n_acc", "prop",
  fill = "group", color = "group", facet.by = "hap") +
  ggtitle("Module occupancy by N modules bound")


my_ggsave("figures/two_three_modules_on_bar_bothHaps.pdf", mod2_3_on_bar_ALL)


# by num mods
mod2_only_on_bar_ALL <- ggbarplot(mod23_merged_all %>% filter(n_acc == 2), "hap", "prop",
  fill = "group", color = "group") +
  ggtitle("Module occupancy by N modules bound")

my_ggsave("figures/two_modules_on_bar_bothHaps.pdf", mod2_only_on_bar_ALL)


mod3_only_on_bar_ALL <- ggbarplot(mod23_merged_all %>% filter(n_acc == 3), "hap", "prop",
  fill = "group", color = "group") +
  ggtitle("Module occupancy by N modules bound")

my_ggsave("figures/three_modules_on_bar_bothHaps.pdf", mod3_only_on_bar_ALL)


mod4_only_on_bar_ALL <- ggbarplot(mod23_merged_all %>% filter(n_acc == 4), "hap", "prop",
  fill = "group", color = "group") +
  ggtitle("Module occupancy by N modules bound")

my_ggsave("figures/four_modules_on_bar_bothHaps.pdf", mod4_only_on_bar_ALL)



# Hypothesis test of Cluster 6 (Cluster 1 in the figure) composition for Liver. Null is H1 == H2, one-sided ------------------------------------------------------
tot_clustered_H1 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H1") %>% select(count) %>% sum()
tot_clustered_H2 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H2") %>% select(count) %>% sum()
clust_6_H1 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H1", cluster == 6) %>% select(count) %>% sum()
clust_6_H2 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H2", cluster == 6) %>% select(count) %>% sum()

# H1 successes, H1 failures, H2 successes, H2 failures
fisher.test(matrix(c(clust_6_H1, tot_clustered_H1-clust_6_H1, clust_6_H2, tot_clustered_H2-clust_6_H2), nrow=2), alternative="less") # p-value < 2.2e-16


# Hypothesis test of sub-cluster 6 (Cluster 1.7 in the figure) composition for Liver. Null is H1 == H2, one-sided ------------------------------------------------------
sub_tot_clustered_H1 <- merged_hap_counts_sub %>% select(Liver_H1) %>% sum()
sub_tot_clustered_H2 <- merged_hap_counts_sub %>% select(Liver_H2) %>% sum()
sub_clust_6_H1 <- merged_hap_counts_sub %>% filter(cluster == 6) %>% select(Liver_H1) %>% sum()
sub_clust_6_H2 <- merged_hap_counts_sub %>% filter(cluster == 6) %>% select(Liver_H2) %>% sum()

# One-sided Fisher's exact test
fisher.test(matrix(c(sub_clust_6_H1, sub_tot_clustered_H1-sub_clust_6_H1, sub_clust_6_H2, sub_tot_clustered_H2-sub_clust_6_H2), nrow=2), alternative="less") # p-value 0.001158




# Specifically, clusters 2, 3, and 7 (4, 5, and 6 in Figure 4) which consist of variably positioned nucleosome arrays, were selective to reads from lymphoblastoid cells. ------------------------
clust_2_Liv <- cluster_counts %>% filter(cluster == 2) %>% select(num_Liver) %>% sum()
clust_2_GM <- cluster_counts %>% filter(cluster == 2) %>% select(num_GM12878) %>% sum()
clust_3_Liv <- cluster_counts %>% filter(cluster == 3) %>% select(num_Liver) %>% sum()
clust_3_GM <- cluster_counts %>% filter(cluster == 3) %>% select(num_GM12878) %>% sum()
clust_7_Liv <- cluster_counts %>% filter(cluster == 7) %>% select(num_Liver) %>% sum()
clust_7_GM <- cluster_counts %>% filter(cluster == 7) %>% select(num_GM12878) %>% sum()

tot_clustered_Liv <- sum(cluster_counts$num_Liver)
tot_clustered_GM <- sum(cluster_counts$num_GM12878)

fisher.test(matrix(c(clust_2_Liv, tot_clustered_Liv-clust_2_Liv, clust_2_GM, tot_clustered_GM-clust_2_GM), nrow=2), alternative="less") # p-value < 2.2e-16
fisher.test(matrix(c(clust_3_Liv, tot_clustered_Liv-clust_3_Liv, clust_3_GM, tot_clustered_GM-clust_3_GM), nrow=2), alternative="less") # p-value < 2.2e-16
fisher.test(matrix(c(clust_7_Liv, tot_clustered_Liv-clust_7_Liv, clust_7_GM, tot_clustered_GM-clust_7_GM), nrow=2), alternative="less") # p-value < 2.2e-16



# In contrast, clusters 4 and 5 (2 and 3 in Figure 4), have a preference for liver reads from the rs2280838-C haplotype
clust_4_H1 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H1", cluster == 4) %>% select(count) %>% sum()
clust_4_H2 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H2", cluster == 4) %>% select(count) %>% sum()
clust_5_H1 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H1", cluster == 5) %>% select(count) %>% sum()
clust_5_H2 <- cluster_raw_counts_hap_merged_GA %>% filter(tissue == "Liver_H2", cluster == 5) %>% select(count) %>% sum()

fisher.test(matrix(c(clust_4_H1, tot_clustered_H1 - clust_4_H1, clust_4_H2, tot_clustered_H2 - clust_4_H2), nrow=2), alternative="greater") # p-value < 2.2e-16
fisher.test(matrix(c(clust_5_H1, tot_clustered_H1 - clust_5_H1, clust_5_H2, tot_clustered_H2 - clust_5_H2), nrow=2), alternative="greater") # p-value < 2.2e-16


# Notably, rs2280838-C had a higher proportion of fibers with 0 modules actuated compared to rs2280838-T, corresponding to a closed chromatin state which is consistent with reduced expression from the rs2280838-C allele.
no_mods_H1 <- nrow(z_df_h1 %>% filter(rsum == 0))
no_mods_H2 <- nrow(z_df_h2 %>% filter(rsum == 0))

fisher.test(matrix(c(no_mods_H1, nrow(z_df_h1) - no_mods_H1, no_mods_H2, nrow(z_df_h2) - no_mods_H2), nrow=2), alternative="greater") # p-value < 2.2e-16


# For the 1-module actuated state, module C was preferentially actuated in rs2280838-T fibers (59.7%)
one_mod_H1 <- z_df_h1 %>% filter(rsum == 1)
one_mod_H2 <- z_df_h2 %>% filter(rsum == 1)
one_mod_C_H1 <- one_mod_H1 %>% filter(X3 == 1)
one_mod_C_H2 <- one_mod_H2 %>% filter(X3 == 1)
one_mod_B_H1 <- one_mod_H1 %>% filter(X2 == 1)
one_mod_E_H1 <- one_mod_H2 %>% filter(X5 == 1)


one_mod_C_H1_n <- nrow(one_mod_C_H1)
one_mod_C_H2_n <- nrow(one_mod_C_H2)

fisher.test(matrix(c(one_mod_C_H1_n, nrow(z_df_h1) - one_mod_C_H1_n, one_mod_C_H2_n, nrow(z_df_h2) - one_mod_C_H2_n), nrow=2), alternative="less") # p-value = 2.292e-07


# whereas modules B, C, or E were similarly actuated in rs2280838-C fibers (Fig. 4e) (30%, 31%, 30%)
nrow(one_mod_B_H1)/nrow(one_mod_H1)
nrow(one_mod_C_H1)/nrow(one_mod_H1)
nrow(one_mod_E_H1)/nrow(one_mod_H1)

# In contrast, the specific combinations of modules actuated in the 2-and 3-module bound states were more similar between haplotypes and largely included actuation of module C (two module 69% & 87%, three module 95% & 91%)
two_mod_H1 <- z_df_h1 %>% filter(rsum == 2)
two_mod_H2 <- z_df_h2 %>% filter(rsum == 2)
three_mod_H1 <- z_df_h1 %>% filter(rsum == 3)
three_mod_H2 <- z_df_h2 %>% filter(rsum == 3)

nrow(two_mod_H1 %>% filter(X3 == 1)) / nrow(two_mod_H1)
nrow(two_mod_H2 %>% filter(X3 == 1)) / nrow(two_mod_H2)
nrow(three_mod_H1 %>% filter(X3 == 1)) / nrow(three_mod_H1)
nrow(three_mod_H2 %>% filter(X3 == 1)) / nrow(three_mod_H2)


