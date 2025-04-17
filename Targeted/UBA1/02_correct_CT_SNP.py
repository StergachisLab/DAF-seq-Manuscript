import pysam
import csv


# CORRECT C->T basecalls at chrX:47,194,331 C/T Het (change from Y to T on CT strands)

bam = 'PS00685_GM12878_UBA1_m84046_240802_231812_s1.hifi_reads.bc2093.map-pb_corrected_realigned.bam'


# chrX:47,190,549-47,194,939 UBA1
reg_chrom = 'chrX'
reg_start = 47190549
reg_end = 47194939

SNP = 47194331

# correct reads by haplotype and filter deaminations by hap & position
def correct_SNP(read_obj, strand, SNP):
    seq = read_obj.query_sequence
    pair = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
    new_seq = ''
    amb_codes = {'CT':'Y', 'GA':'R'}
    deam_pos = [] # mol coordinates of likely base changes
    for pos in pair:
        if pos[0] == None: # deletion, ignore
            pass
        elif pos[1] == None: # insertion, use seq base
            qi = pos[0]
            new_seq += seq[qi]
        else:
            qi = pos[0]
            ref_coord = pos[1] + 1
            if ref_coord == SNP:
                new_base = None
                if seq[qi] == 'Y':
                    new_base = 'T'
                else:
                    new_base = seq[qi]
                new_seq += new_base
                if new_base == amb_codes[strand]:
                    deam_pos.append(qi+1) # track DA positions in 1-indexed
            else:
                new_seq += seq[qi] # don't correct base
                if seq[qi] == amb_codes[strand]:
                    deam_pos.append(qi+1) # track DA positions in 1-indexed
    return(new_seq, deam_pos)


# write corrected reads to new BAM
new_bam = 'UBA1_PS00685_GM12878_haplotype_corrected.bam'
corrected_bam = pysam.AlignmentFile(new_bam, "wb", template=bam)
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        strand = read.get_tag('ST')
        if strand in ['CT','GA']:
            # WRITE NEW seq with added tags
            new_seq, deam_pos = correct_SNP(read, strand, SNP)
            quals = read.query_qualities
            read.query_sequence = new_seq
            read.seq = new_seq
            read.query_qualities = quals
            read.set_tag('DA', deam_pos)
            read.set_tag('FD', deam_pos[0], "i")
            read.set_tag('LD', deam_pos[-1], "i")
            read.set_tag('ST', strand)
        else:
            read.set_tag('DA', [0])
            read.set_tag('FD', 0, "i")
            read.set_tag('LD', 0, "i")
            read.set_tag('ST', strand)
        corrected_bam.write(read)
bam.close()
corrected_bam.close()






