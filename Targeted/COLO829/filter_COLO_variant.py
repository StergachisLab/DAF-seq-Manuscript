""" Filter COLO region2 GA reads to only include T hap2 reads with BOTH C->T variants. """

import pysam
import csv

# HEART SAMPLE ---------------------------------------------------------------------------------------------------------
bam_name = 'PS00719_COLO_Region2_m84055_240822_221245_s3.hifi_reads.bc2069.map-pb_corrected_realigned_GA_ONLY.bam'
new_bam = 'hap2_T_only.bam'
bam = pysam.AlignmentFile(bam_name, "rb")

vars = [19447245,19447246]

# chr17:19,446,516-19,450,322
reg_chrom = 'chr17'
reg_start = 19446516
reg_end = 19450322

# write corrected reads to new BAM
corrected_bam = pysam.AlignmentFile(new_bam, "wb", template=bam)
for read in bam.fetch(reg_chrom, reg_start, reg_end):
        pass_count = 0
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        for pos in pair:
            if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                mcoord = pos[0]
                ref_coord = pos[1] + 1
                if ref_coord in vars:
                    if seq[mcoord] == 'T':
                        pass_count += 1
        if pass_count == 2:
            corrected_bam.write(read)
bam.close()
corrected_bam.close()
