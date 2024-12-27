import pysam
import csv
import sys
import os


bams = sys.argv[1:]
if len(bams) < 1:
    sys.exit(f"BAM file not provided! Exiting.....")
for bf in bams:
    if not os.path.exists(bf):
        sys.exit(f"BAM file {bf} can't be found! Exiting.....")

name = os.path.basename(bams[0]).split('_')[0]

cutoff = 0.90

def determine_da_strand_MD(read_obj, cutoff):
    # based on the proportion of C->T & G->A determine the strand acted upon by DddA
    # only counting single base substitutions for now
    seq = read_obj.query_sequence
    pair = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
    c = 0
    g = 0
    total = 0
    for pos in pair:
        if pos[0] == None or pos[1] == None: # indel, ignore
            pass
        else:
            qi = pos[0]
            ref_pos = pos[2].upper()
            if seq[qi] != ref_pos:
                total += 1
                change = ref_pos + seq[qi]
                if change == "CT":
                    c += 1
                elif change == 'GA':
                    g += 1
    if c+g == 0:
        return('none')
    elif c/(c+g) >= cutoff:
        return('CT')
    elif g/(c+g) >= cutoff:
        return('GA')
    else:
        return('undetermined')

CT_seqs = []
GA_seqs = []
for bam_file in bams:
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam.fetch():
        if read.is_secondary == False and read.is_supplementary == False:
            strand = determine_da_strand_MD(read, cutoff)
            seq = read.seq
            pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
            clipped_seq = ''
            first_base = None
            final_base = None
            for pos in pairs:
                if pos[1] and pos[1] >= read.reference_start and pos[1] <= read.reference_end:
                    if pos[0]:
                        clipped_seq += seq[pos[0]]
                        if first_base == None:
                            first_base = pos[1]
                        final_base = pos[1]
                    else: # deletion
                        clipped_seq += 'D'
                        # don't count indel bases as first or last ref positions !!!!!!!!!!!!!!!!!
                        if first_base == None:
                            first_base = pos[1]
                        final_base = pos[1]
                elif pos[1] == None: # insertion
                    clipped_seq += seq[pos[0]].lower() # lowercase for inserted bases
            # clip inserted bases from the ends of the reference sequence to aid in kmer alignment (ref coordinates, accurate first_base and final_base)
            clipped_seq = clipped_seq.strip('gtca')
            if strand == 'CT':
                CT_seqs.append([read.qname, read.reference_name, first_base+1, final_base+1, clipped_seq])
            elif strand == 'GA':
                GA_seqs.append([read.qname, read.reference_name, first_base+1, final_base+1, clipped_seq])

with open(f'clipped_seqs/{name}_clipped_CT.tsv','w') as fw:
    writer = csv.writer(fw, delimiter='\t')
    for f in CT_seqs:
        writer.writerow(f)
        
with open(f'clipped_seqs/{name}_clipped_GA.tsv','w') as fw:
    writer = csv.writer(fw, delimiter='\t')
    for f in GA_seqs:
        writer.writerow(f)
