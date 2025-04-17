#!/bin/bash

bedmap --ec --echo --echo-map-id tss_uba1_masSeq.bed <( sort-bed ../pileups/GM12878_UBA1_PS00685_H1_large_msp_positions.bed ) > large_msps_H1_fire_map.txt
bedmap --ec --echo --echo-map-id ../tss_uba1_masSeq.bed <( sort-bed ../pileups/GM12878_UBA1_PS00685_H2_large_msp_positions.bed ) > large_msps_H2_fire_map.txt

samtools view ../pileups/GM12878_UBA1_PS00685_nuc_H1.bam | awk '{ print $1 }' > H1_GA_GM12878_UBA1_zmws.txt
