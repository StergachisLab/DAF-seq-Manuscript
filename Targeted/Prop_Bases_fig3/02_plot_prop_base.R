library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(cowplot)
library(colorspace)
library(ggpubr)
library(hexbin)
library(RColorBrewer)


n_bins = 30

samples <- c('NAPA','WASF1','UBA1','GM12878_SLC','Liver_SLC', 'Heart_SLC', 'Colon_SLC','COLO_Mix')

for (i in 1:length(samples))
{
    df_a = read_csv(paste0('prop_by_strand/', samples[i], '_prop_A_by_strand.csv'))
    df_c = read_csv(paste0('prop_by_strand/', samples[i], '_prop_C_by_strand.csv'))
    df_g = read_csv(paste0('prop_by_strand/', samples[i], '_prop_G_by_strand.csv'))
    df_t = read_csv(paste0('prop_by_strand/', samples[i], '_prop_T_by_strand.csv'))
    plot_a <- ggplot(df_a, aes(x=bottom_prop, y=top_prop) ) +
        geom_hex(bins = n_bins) +
        # scale_fill_distiller("", palette = "Spectral", trans="log10") +
        scale_fill_continuous(type = "viridis", trans = "log10") +
        theme_bw(10) +
        xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
    plot_c <- ggplot(df_c, aes(x=bottom_prop, y=top_prop) ) +
        geom_hex(bins = n_bins) +
        # scale_fill_distiller("", palette = "Spectral", trans="log10") +
        scale_fill_continuous(type = "viridis", trans = "log10") +
        theme_bw(10) +
        xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
    plot_g <- ggplot(df_g, aes(x=bottom_prop, y=top_prop) ) +
        geom_hex(bins = n_bins) +
        # scale_fill_distiller("", palette = "Spectral", trans="log10") +
        scale_fill_continuous(type = "viridis", trans = "log10") +
        theme_bw(10) +
        xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
    plot_t <- ggplot(df_t, aes(x=bottom_prop, y=top_prop) ) +
        geom_hex(bins = n_bins) +
        # scale_fill_distiller("", palette = "Spectral", trans="log10") +
        scale_fill_continuous(type = "viridis", trans = "log10") +
        theme_bw(10) +
        xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
    merged_plots <- plot_grid(plot_c, plot_g, plot_t, plot_a, labels = c('C','G','T','A'), label_size = 12)
    ggsave(paste0(samples[i], '_prop_by_base.png'), plot=merged_plots, width = 8, height = 6)
}


# MERGE ALL SAMPLES

# initialize DFs
df_a = read_csv(paste0('prop_by_strand/', samples[1], '_prop_A_by_strand.csv'))
df_c = read_csv(paste0('prop_by_strand/', samples[1], '_prop_C_by_strand.csv'))
df_g = read_csv(paste0('prop_by_strand/', samples[1], '_prop_G_by_strand.csv'))
df_t = read_csv(paste0('prop_by_strand/', samples[1], '_prop_T_by_strand.csv'))
df_a$sample <- samples[1]
df_c$sample <- samples[1]
df_g$sample <- samples[1]
df_t$sample <- samples[1]
for (i in 2:length(samples))
{
    temp_a = read_csv(paste0('prop_by_strand/', samples[i], '_prop_A_by_strand.csv'))
    temp_c = read_csv(paste0('prop_by_strand/', samples[i], '_prop_C_by_strand.csv'))
    temp_g = read_csv(paste0('prop_by_strand/', samples[i], '_prop_G_by_strand.csv'))
    temp_t = read_csv(paste0('prop_by_strand/', samples[i], '_prop_T_by_strand.csv'))
    temp_a$sample <- samples[i]
    temp_c$sample <- samples[i]
    temp_g$sample <- samples[i]
    temp_t$sample <- samples[i]
    df_a <- rbind(df_a, temp_a)
    df_c <- rbind(df_c, temp_c)
    df_g <- rbind(df_g, temp_g)
    df_t <- rbind(df_t, temp_t)
}

plot_a <- ggplot(df_a, aes(x=bottom_prop, y=top_prop) ) +
    geom_hex(bins = n_bins) +
    # scale_fill_distiller("", palette = "Spectral", trans="log10") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_bw(10) +
    xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
plot_c <- ggplot(df_c, aes(x=bottom_prop, y=top_prop) ) +
    geom_hex(bins = n_bins) +
    # scale_fill_distiller("", palette = "Spectral", trans="log10") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_bw(10) +
    xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
plot_g <- ggplot(df_g, aes(x=bottom_prop, y=top_prop) ) +
    geom_hex(bins = n_bins) +
    # scale_fill_distiller("", palette = "Spectral", trans="log10") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_bw(10) +
    xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
plot_t <- ggplot(df_t, aes(x=bottom_prop, y=top_prop) ) +
    geom_hex(bins = n_bins) +
    # scale_fill_distiller("", palette = "Spectral", trans="log10") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_bw(10) +
    xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.")
merged_plots <- plot_grid(plot_c, plot_g, plot_t, plot_a, labels = c('C','G','T','A'), label_size = 12)
ggsave('ALL_samples_prop_by_base.pdf', plot=merged_plots, width = 8, height = 6)


# Zoom in of COLO BL:T C->T variants (19,447,245 and 19,447,246)
df_c = read_csv('prop_by_strand/COLO_Mix_prop_C_by_strand.csv')
df_t = read_csv('prop_by_strand/COLO_Mix_prop_T_by_strand.csv')

var_pos <- c(19447245,19447246)

df_t$var <- FALSE
df_t[df_t$position %in% var_pos,]$var <- TRUE

plot_t <- ggplot(df_t, aes(x=bottom_prop, y=top_prop) ) +
    geom_point(aes(colour=factor(var), fill = factor(var)), size = 3) +
    scale_fill_manual(values=c("blue", "red")) +
    scale_colour_manual(values=c("blue", "red")) +
    xlim(0,0.05) +
    ylim(0, 0.2) +
    theme_bw(15) +
    xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.") +
    theme(legend.position="none")


df_c %>% filter(position %in% var_pos)
df_c$var <- FALSE
df_c[df_c$position %in% var_pos,]$var <- TRUE

plot_c <- ggplot(df_c, aes(x=bottom_prop, y=top_prop) ) +
    geom_point(aes(colour=factor(var), fill = factor(var)), size = 3) +
    scale_fill_manual(values=c("blue", "red")) +
    scale_colour_manual(values=c("blue", "red")) +
    xlim(0.95,1) +
    ylim(0.8, 1) +
    theme_bw(15) +
    xlab("Bottom Strand Prop.") + ylab("Top Strand Prop.") +
    theme(legend.position="none")


merged_plots <- plot_grid(plot_c, plot_t, labels = c('C','T'), label_size = 18)
ggsave('COLO_Mix_Zoom_prop_by_base.pdf', plot=merged_plots, width = 8, height = 4)


