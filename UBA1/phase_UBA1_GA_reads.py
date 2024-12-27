import pysam

bam_name = 'UBA1_region_GA.bam'
bam = pysam.AlignmentFile(bam_name, "rb")

# chrX:47,194,331 C/T, H1/H2, Xa/Xi

variant = 47194331
new_bam = 'UBA1_region_GA_phased.bam'
phased_bam = pysam.AlignmentFile(new_bam, "wb", template=bam)
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        hap = None
        for pos in pair:
            if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                ref_coord = pos[1] + 1
                if ref_coord == variant:
                    mcoord = pos[0]
                    basecall = seq[mcoord]
                    if basecall == 'C':
                        hap = 1
                    elif basecall == 'T':
                        hap = 2
        if hap != None:
            read.set_tag('HP', hap, 'i')
            phased_bam.write(read)
bam.close()
phased_bam.close()

