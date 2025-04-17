#!/bin/bash

sort-bed UBA1_GM12878_footprinted_motifs_H1.bed > UBA1_GM12878_footprinted_motifs_H1_sorted.bed
bedmap --max-element UBA1_GM12878_footprinted_motifs_H1_sorted.bed | sort | uniq | sort-bed - > UBA1_H1_max_motifs.bed

# Merge motifs to find footprinted regions, not motifs
bedops --everything UBA1_GM12878_footprinted_motifs_H1_sorted.bed UBA1_GM12878_footprinted_motifs_H1_sorted.bed | sort-bed - | bedmap --ec --fraction-either 0.5 --echo-map-range - | sort-bed - | uniq > temp_merged.bed
cat temp_merged.bed | bedmap --ec --fraction-either 0.1 --echo-map-range - | sort-bed - | uniq | bedmap --ec --fraction-either 0.1 --echo-map-range - | sort-bed - | uniq > merged_UBA1_GM12878_footprinted_motifs_H1_sorted.bed
rm temp_merged.bed

