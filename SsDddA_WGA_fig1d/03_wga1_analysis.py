import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import logomaker as lm


"""
For each processed BAM file: make a table of DA kmers and PWM.
"""

def rc(s):
    bases = {'A':'T','C':'G','G':'C','T':'A','Y':'R','R':'Y'}
    new = ''
    for i in range(len(s)):
        new = bases[s[i]] + new
    return(new)

def gc_prop(seq):
    gc = 0
    tot = 0
    for i in range(len(seq)):
        tot += 1
        if seq[i] in {'G','C'}:
            gc += 1
    prop = gc/tot
    return(prop)

def generate_all_kmers(k):
    # recursive https://thecodingbiologist.com/posts/Find-all-possible-kmers-in-python
    bases = ['A','C','G','T']
    if k < 1:
        return([])
    if k == 1:
        return(bases)
    sub_sequences = generate_all_kmers(k-1)
    sequences = []
    for seq in sub_sequences:
        for b in bases:
                 sequences.append(b + seq)
    return(sequences)

def kmer_pwm(kmers, k):
    """ Create position weight matrix of deamination events based on a list of kmers """
    counts  = [{'A':0,'C':0,'G':0,'T':0} for i in range(k)]
    for kmer,kval in kmers.items():
        for i in range(k):
            counts[i][kmer[i]] += kval
    totals = [sum(i.values()) for i in counts]
    props  = [{'A':counts[i]['A']/totals[i],'C':counts[i]['C']/totals[i],'G':counts[i]['G']/totals[i],'T':counts[i]['T']/totals[i]} for i in range(k)]
    return(props)


def da_kmer(read_obj, k):
    # track kmers of size k, centered on DA
    # if k%2 == 1:
    #     pass
    padding = int((k-1)/2)
    ref_kmers = []
    read_kmers = []
    seq = read_obj.query_sequence
    pair = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
    DA = read_obj.get_tag('DA')
    for pos in DA:
        pos -= 1 # change back to 0-indexed
        pos_st = pos-padding
        pos_end = pos+padding
        if pos_st >= 0 and pos_end < len(seq):
            ref_kmer = ''
            read_kmer = ''
            if seq[pos] == 'Y':
                for i in range(pos_st,pos_end+1):
                    if i == pos: # middle of kmer should be original base (C or G)
                        ref_kmer += 'C'
                        read_kmer += 'C'
                    elif seq[i] == 'Y':
                        ref_kmer += 'C'
                        read_kmer += 'T'
                    else:
                        if pair[i][2] == None:
                            ref_kmer += seq[i]
                        else:
                            ref_kmer += pair[i][2].upper()
                        read_kmer += seq[i]
                if 'N' not in ref_kmer:
                    ref_kmers.append(ref_kmer)
                if 'N' not in read_kmer:
                    read_kmers.append(read_kmer)
            elif seq[pos] == 'R':
                for i in range(pos_st,pos_end+1):
                    if i == pos: # middle of kmer should be original base (C or G)
                        ref_kmer += 'G'
                        read_kmer += 'G'
                    elif seq[i] == 'R':
                        ref_kmer += 'G'
                        read_kmer += 'A'
                    else:
                        if pair[i][2] == None:
                            ref_kmer += seq[i]
                        else:
                            ref_kmer += pair[i][2].upper()
                        read_kmer += seq[i]
                # convert kmers from G/A strands back to C/T orientation
                if 'N' not in ref_kmer:
                    ref_kmer = rc(ref_kmer)
                    ref_kmers.append(ref_kmer)
                if 'N' not in read_kmer:
                    read_kmer = rc(read_kmer)
                    read_kmers.append(read_kmer)
    return([ref_kmers,read_kmers])


wt_bam = 'PS00535_K562_WT.map-pb_corrected_realigned.bam'
d5_bam = 'PS00536_K562_5.map-pb_corrected_realigned.bam'
hg002_wt_bam = 'PS00537_HG002_WT.map-pb_corrected_realigned.bam'
hg002_noAmp_bam = 'PS00538_HG002_5_noWGA.map-pb_corrected_realigned.bam'

# generate every possible kmer --> track surrounding nt context based on Ref sequence and read sequence
k = 7
all_kmers = {seq for seq in generate_all_kmers(k) if seq[3] == 'C'}
ref_kmers_wt = {seq:0 for seq in all_kmers}
read_kmers_wt = {seq:0 for seq in all_kmers}
ref_kmers_d5 = {seq:0 for seq in all_kmers}
read_kmers_d5 = {seq:0 for seq in all_kmers}

ref_kmers_hg = {seq:0 for seq in all_kmers}
read_kmers_hg = {seq:0 for seq in all_kmers}
ref_kmers_noAmp = {seq:0 for seq in all_kmers}
read_kmers_noAmp = {seq:0 for seq in all_kmers}

bam = pysam.AlignmentFile(wt_bam, "rb")
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        refk,readk = da_kmer(read, 7)
        for kmer in refk:
            if 'N' not in kmer:
                ref_kmers_wt[kmer] += 1
        for kmer in readk:
            if 'N' not in kmer:
                read_kmers_wt[kmer] += 1

ref_wt_pwm = kmer_pwm(ref_kmers_wt, 7)
read_wt_pwm = kmer_pwm(read_kmers_wt, 7)


bam = pysam.AlignmentFile(d5_bam, "rb")
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        refk,readk = da_kmer(read, 7)
        for kmer in refk:
            if 'N' not in kmer:
                ref_kmers_d5[kmer] += 1
        for kmer in readk:
            if 'N' not in kmer:
                read_kmers_d5[kmer] += 1

ref_d5_pwm = kmer_pwm(ref_kmers_d5, 7)
read_d5_pwm = kmer_pwm(read_kmers_d5, 7)


# HG002 ---------------------------------------------
bam = pysam.AlignmentFile(hg002_wt_bam, "rb")
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        refk,readk = da_kmer(read, 7)
        for kmer in refk:
            if 'N' not in kmer:
                ref_kmers_hg[kmer] += 1
        for kmer in readk:
            if 'N' not in kmer:
                read_kmers_hg[kmer] += 1

ref_hg_pwm = kmer_pwm(ref_kmers_hg, 7)
read_hg_pwm = kmer_pwm(read_kmers_hg, 7)

bam = pysam.AlignmentFile(hg002_noAmp_bam, "rb")
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        refk,readk = da_kmer(read, 7)
        for kmer in refk:
            if 'N' not in kmer:
                ref_kmers_noAmp[kmer] += 1
        for kmer in readk:
            if 'N' not in kmer:
                read_kmers_noAmp[kmer] += 1

ref_noAmp_pwm = kmer_pwm(ref_kmers_noAmp, 7)
read_noAmp_pwm = kmer_pwm(read_kmers_noAmp, 7)



# Available logo fonts
    # lm.list_font_names()

pwm = [ref_wt_pwm, read_wt_pwm, ref_d5_pwm, read_d5_pwm, ref_hg_pwm, read_hg_pwm, ref_noAmp_pwm, read_noAmp_pwm]
pwm_names = ['ref_wt', 'read_wt', 'ref_d5', 'read_d5', 'ref_hg', 'read_hg', 'ref_noAmp', 'read_noAmp']
for i in range(len(pwm)):
    df = pd.DataFrame(pwm[i])
    logo = lm.Logo(df, font_name = 'Arial Rounded MT Bold')
    plt.savefig(f'{pwm_names[i]}_test.pdf')




# base_counts = {'A':0, 'C':0, 'G':0, 'T':0, 'Y':0, 'R':0}
# base_counts = {'A':0, 'C':0, 'G':0, 'T':0, 'Y':0, 'R':0}
# for bam_name in ['PS00536_K562_5.map-ont_corrected.bam']:
#     bam = pysam.AlignmentFile(bam_name, "rb")
#     for read in bam.fetch("chr22"):
#         if read.is_secondary == False and read.is_supplementary == False:
#             print()

# read.get_cigar_stats()
# read.cigarstring
# read.is_mapped
# read.seq
# read.query_sequence
# read.qlen # reference span

# 143265541
# chr22
# 33950255

# 53152966
# m84046_240224_152130_s4/53152966/ccs
# chr22:30,360,473-30,364,300
# chr2:74,587,397-74,591,957

# PS00576
# chr4:2,631,175-2,638,339 Intra read chimera?
# m84055_240423_222813_s1/28049886/ccs 
# 4224 bp
# 2187
# 1977






