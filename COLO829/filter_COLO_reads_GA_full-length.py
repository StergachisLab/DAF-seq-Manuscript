""" Filter COLO region2 reads to only include GA reads that fully span the region. """

import pysam
import csv

# HEART SAMPLE ---------------------------------------------------------------------------------------------------------
bam_name = 'PS00719_COLO_Region2_m84055_240822_221245_s3.hifi_reads.bc2069.map-pb_corrected_realigned.bam'
new_bam = 'PS00719_COLO_Region2_m84055_240822_221245_s3.hifi_reads.bc2069.map-pb_corrected_realigned_GA_ONLY.bam'
bam = pysam.AlignmentFile(bam_name, "rb")

# chr17:19,446,516-19,450,322
reg_chrom = 'chr17'
reg_start = 19446516
reg_end = 19450322
bufferlen = 50

# write corrected reads to new BAM
corrected_bam = pysam.AlignmentFile(new_bam, "wb", template=bam)
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False and read.reference_start <= (reg_start+bufferlen) and read.reference_end >= (reg_end-bufferlen):
        strand = read.get_tag('ST')
        if (strand == 'GA' and 'Y' not in read.seq):
            corrected_bam.write(read)
bam.close()
corrected_bam.close()
