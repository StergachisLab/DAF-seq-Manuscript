""" Use the basecall % from the CT and GA reads to identify true genomic variants rather than DddA-induced base changes. """

import pysam
import csv

bam_name = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/PS00626.m84046_240619_124816_s1.bc2072.ft.map-pb_corrected_realigned.bam'
bam = pysam.AlignmentFile(bam_name, "rb")

# chr19:47,514,458-47,519,061 NAPA
reg_chrom = 'chr19'
reg_start = 47514458
reg_end = 47519061

# basecalls by CT & GA strands
ct_pos = dict()
ga_pos = dict()

# ref coordinates are 1-based !!!!!!!!!!
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        strand = read.get_tag('ST')
        if strand == 'CT':
            for pos in pair:
                if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                    mcoord = pos[0]
                    ref_coord = pos[1] + 1
                    if ref_coord not in ct_pos.keys():
                        ct_pos[ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0}
                    if seq[mcoord] == 'Y':
                        ct_pos[ref_coord]['T'] += 1
                    else:
                        ct_pos[ref_coord][seq[mcoord]] += 1
        elif strand == 'GA':
            for pos in pair:
                if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                    mcoord = pos[0]
                    ref_coord = pos[1] + 1
                    if ref_coord not in ga_pos.keys():
                        ga_pos[ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0}
                    if seq[mcoord] == 'R':
                        ga_pos[ref_coord]['A'] += 1
                    else:
                        ga_pos[ref_coord][seq[mcoord]] += 1

shared_keys = set(ct_pos.keys()).intersection(set(ga_pos.keys()))

# identify het variants
strand_cutoff = 0.9
hets = []
daf_basecalls = dict()
base_type = None
for sk in shared_keys:
    if (ct_pos[sk]['G'] / sum(ct_pos[sk].values())) >= strand_cutoff:
        base_type = 'G'
    elif (ga_pos[sk]['C'] / sum(ga_pos[sk].values())) >= strand_cutoff:
        base_type = 'C'
    elif (ct_pos[sk]['A'] / sum(ct_pos[sk].values())) >= strand_cutoff:
        base_type = 'A'
    elif (ct_pos[sk]['T'] / sum(ct_pos[sk].values())) >= strand_cutoff:
        base_type = 'T'
    else:
        base_type = 'HET'
        hets.append(sk)
    daf_basecalls[sk] = base_type

# Count combos of specified HET variants in reads
het_ref = [47515000, 47516117, 47518383]
ct_hets = dict()
ga_hets = dict()
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        strand = read.get_tag('ST')
        if strand == 'CT':
            read_calls = ['N' for i in range(len(het_ref))]
            for pos in pair:
                if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                    mcoord = pos[0]
                    ref_coord = pos[1] + 1
                    if ref_coord in het_ref:
                        basecall = seq[mcoord]
                        if basecall == 'Y':
                            basecall = 'T'
                        read_calls[het_ref.index(ref_coord)] = basecall
            call_seq = ''.join(read_calls)
            if call_seq not in ct_hets.keys():
                ct_hets[call_seq] = set()
            ct_hets[call_seq].add(read.qname)
        elif strand == 'GA':
            read_calls = ['N' for i in range(len(het_ref))]
            for pos in pair:
                if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                    mcoord = pos[0]
                    ref_coord = pos[1] + 1
                    if ref_coord in het_ref:
                        basecall = seq[mcoord]
                        if basecall == 'R':
                            basecall = 'A'
                        read_calls[het_ref.index(ref_coord)] = basecall
            call_seq = ''.join(read_calls)
            if call_seq not in ga_hets.keys():
                ga_hets[call_seq] = set()
            ga_hets[call_seq].add(read.qname)

ct_hets = dict(sorted(ct_hets.items(), key=lambda item: len(item[1]), reverse=True))
ga_hets = dict(sorted(ga_hets.items(), key=lambda item: len(item[1]), reverse=True))


# [47515000, 47516117, 47518383]
# H1: TAC
# H2: GTG

# treating C/T on CT reads as C, G/A on GA reads as G
h1_zmw = set()
h2_zmw = set()
unphased = set()
for k,v in ct_hets.items():
    h1_score = 0
    h2_score = 0
    # var 1
    if k[0] in ['T']:
        h1_score += 1
    elif k[0] in ['G']:
        h2_score += 1
    # var 2
    if k[1] in ['A']:
        h1_score += 1
    elif k[1] in ['T']:
        h2_score += 1
    # var 3
    if k[2] in ['C','T']:
        h1_score += 1
    elif k[2] in ['G']:
        h2_score += 1
    if (h1_score > 0 and h2_score > 0) == False:
        if h1_score > 0:
            if len(h1_zmw) == 0:
                h1_zmw = v
            h1_zmw = h1_zmw.union(v)
        elif h2_score > 0:
            if len(h2_zmw) == 0:
                h2_zmw = v
            h2_zmw = h2_zmw.union(v)
        else:
            if len(unphased) == 0:
                unphased = v
            unphased = unphased.union(v)
    else:
        if len(unphased) == 0:
            unphased = v
        unphased = unphased.union(v)
for k,v in ga_hets.items():
    h1_score = 0
    h2_score = 0
    # var 1
    if k[0] in ['T']:
        h1_score += 1
    elif k[0] in ['G','A']:
        h2_score += 1
    # var 2
    if k[1] in ['A']:
        h1_score += 1
    elif k[1] in ['T']:
        h2_score += 1
    # var 3
    if k[2] in ['C']:
        h1_score += 1
    elif k[2] in ['G','A']:
        h2_score += 1
    if (h1_score > 0 and h2_score > 0) == False:
        if h1_score > 0:
            if len(h1_zmw) == 0:
                h1_zmw = v
            h1_zmw = h1_zmw.union(v)
        elif h2_score > 0:
            if len(h2_zmw) == 0:
                h2_zmw = v
            h2_zmw = h2_zmw.union(v)
        else:
            if len(unphased) == 0:
                unphased = v
            unphased = unphased.union(v)
    else:
        if len(unphased) == 0:
            unphased = v
        unphased = unphased.union(v)

# Identify diff from the hg38 reference
ref_bases = {pos:'' for pos in shared_keys}
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        for pos in pair:
            if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                ref_coord = pos[1] + 1
                if ref_coord in shared_keys:
                    ref_bases[ref_coord] = pos[2].upper()
    if '' not in ref_bases.values():
        break

# generate consensus H1 & H2 basecalls
hap_basecalls = dict()
for k,v in ref_bases.items():
    if v != daf_basecalls[k]:
        hap_basecalls[k] = ['','']
        if daf_basecalls[k] != 'HET':
            hap_basecalls[k] = [daf_basecalls[k],daf_basecalls[k]]

hap_counts = {pos:[{'A':0, 'C':0, 'G':0, 'T':0}, {'A':0, 'C':0, 'G':0, 'T':0}] for pos in hets}
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        for pos in pair:
            if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                mcoord = pos[0]
                ref_coord = pos[1] + 1
                if ref_coord in hap_counts.keys():
                    base = seq[mcoord]
                    if base == 'R':
                        base = 'A'
                    elif base == 'Y':
                        base = 'T'
                    if read.qname in h1_zmw:
                        hap_counts[ref_coord][0][base] += 1
                    elif read.qname in h2_zmw:
                        hap_counts[ref_coord][1][base] += 1
for k,v in hap_counts.items():
    for i in range(2):
        called = False
        for k2,v2 in v[i].items():
            if (v2 / sum(v[i].values())) > strand_cutoff:
                hap_basecalls[k][i] = k2
                called = True
        if not called:
            if (v[i]['C']+v[i]['T'])/sum(v[i].values()) > strand_cutoff: # large mix of C->T
                hap_basecalls[k][i] = 'C'
            elif (v[i]['G']+v[i]['A'])/sum(v[i].values()) > strand_cutoff: # large mix of G->A
                hap_basecalls[k][i] = 'G'
            else:
                hap_basecalls[k][i] = 'N'

# positions in hap_basecalls show either differences from the reference or dubious mixed calls (not 50/50 but not low %) to omit from deamination calling
# 47518160 ['A', 'A'] contains 7% T
# 47518977 ['C', 'C'] contains 13% T, end of the read


# correct reads by haplotype and filter deaminations by hap & position
def correct_by_hap(read_obj, strand, haplotype, hap_basecalls):
    """ Correct previously identified DddA induced deaminations at sites with non-reference and heterozygous positions
    that might result in false-positives and false-negatives.
     Correct the original sequence using ambiguity codes (C|T: Y, G|A: R) and output new DA-tag positions.
     Limit detection to the previously identified DddA strand info (either ct or ga).
     Output everything in FIBER coordinates, not reference!
    """
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
            if ref_coord in hap_basecalls.keys():
                new_base = None
                if hap_basecalls[ref_coord][0] == hap_basecalls[ref_coord][1] and hap_basecalls[ref_coord][0] == pos[2].upper(): # identified as dubious basecall. Omit Da calls at this position.
                    if seq[qi] == 'Y':
                        new_base = 'T'
                    elif seq[qi] == 'R':
                        new_base = 'A'
                else:
                    if haplotype == 1:
                        hap_base = hap_basecalls[ref_coord][0]
                    elif haplotype == 2:
                        hap_base = hap_basecalls[ref_coord][1]
                    if haplotype == 'UNK': # convert DA positions back to original base for unphased reads
                        if seq[qi] == 'Y':
                            new_base = 'T'
                        elif seq[qi] == 'R':
                            new_base = 'A'
                    else:
                        if hap_base == 'T' and strand == 'CT':
                            if seq[qi] == 'Y':
                                new_base = 'T'
                        elif hap_base == 'A' and strand == 'GA':
                            if seq[qi] == 'R':
                                new_base = 'A'
                        elif hap_base == 'G' and strand == 'GA':
                            if seq[qi] == 'A':
                                new_base = 'R'
                        elif hap_base == 'C' and strand == 'CT':
                            if seq[qi] == 'T':
                                new_base = 'Y'
                if new_base == None:
                    new_base = seq[qi] # don't correct base
                new_seq += new_base
                if new_base == amb_codes[strand]:
                    deam_pos.append(qi+1) # track DA positions in 1-indexed
            else:
                new_seq += seq[qi] # don't correct base
                if seq[qi] == amb_codes[strand]:
                    deam_pos.append(qi+1) # track DA positions in 1-indexed
    return(new_seq, deam_pos)


# write corrected reads to new BAM
new_bam = 'NAPA_PS00626_haplotype_corrected.bam'
corrected_bam = pysam.AlignmentFile(new_bam, "wb", template=bam)
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        strand = read.get_tag('ST')
        if strand in ['CT','GA']:
            if read.qname in h1_zmw:
                hap = 1
            elif read.qname in h2_zmw:
                hap = 2
            else:
                hap = 'UNK'
            # WRITE NEW seq with added tags
            new_seq, deam_pos = correct_by_hap(read, strand, hap, hap_basecalls)
            quals = read.query_qualities
            read.query_sequence = new_seq
            read.seq = new_seq
            read.query_qualities = quals
            # read.set_tags([('DA', deam_pos), ('FD', deam_pos[0], "i"), ('LD', deam_pos[-1], "i"), ('ST', strand)])
            read.set_tag('DA', deam_pos)
            read.set_tag('FD', deam_pos[0], "i")
            read.set_tag('LD', deam_pos[-1], "i")
            read.set_tag('ST', strand)
            if hap != 'UNK':
                read.set_tag('HP', hap, 'i')
        else:
            # read.set_tags([('DA', [0]), ('FD', 0, "i"), ('LD', 0, "i"), ('ST', strand)])
            read.set_tag('DA', [0])
            read.set_tag('FD', 0, "i")
            read.set_tag('LD', 0, "i")
            read.set_tag('ST', strand)
        corrected_bam.write(read)
bam.close()
corrected_bam.close()

# write Hap basecalls used for phasing and correcting DAs
with open('NAPA_PS00626_haplotype_corrected_hap_basecalls.tsv','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(['Position','H1','H2'])
    for k,v in hap_basecalls.items():
        line = [k, v[0], v[1]]
        writer.writerow(line)
