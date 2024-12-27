""" Use the basecall % from the CT and GA reads to identify true genomic variants rather than DddA-induced base changes. """

import pysam

bam_dir = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1'

in_bams = ['PS00627.m84046_240619_124816_s1.bc2073.ft.map-pb_corrected_realigned.bam',
           'PS00628.m84046_240619_124816_s1.bc2074.ft.map-pb_corrected_realigned.bam',
           'PS00629.m84046_240619_124816_s1.bc2075.ft.map-pb_corrected_realigned.bam',
           'PS00630.m84046_240619_124816_s1.bc2076.ft.map-pb_corrected_realigned.bam',
           'PS00631.m84046_240619_124816_s1.bc2077.ft.map-pb_corrected_realigned.bam']

# chr19:47,514,458-47,519,061 NAPA
reg_chrom = 'chr19'
reg_start = 47514458
reg_end = 47519061

# verified variant calls from PS00626
# positions in hap_basecalls show either differences from the reference or dubious mixed calls (not 50/50 but not low %) to omit from deamination calling
# 47518160 ['A', 'A'] contains 7% T
# 47518977 ['C', 'C'] contains 13% T, end of the read
hap_basecalls = {47514700: ['T', 'C'],
 47515000: ['T', 'G'],
 47515603: ['G', 'G'],
 47516117: ['A', 'T'],
 47517160: ['C', 'T'],
 47517464: ['T', 'C'],
 47518032: ['T', 'T'],
 47518383: ['C', 'G'],
 47518977: ['C', 'C']}

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

for b in in_bams:
    bam = pysam.AlignmentFile(f'{bam_dir}/{b}', "rb")
    new_bam = f"NAPA_{b.split('.')[0]}_haplotype_corrected.bam"
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
    # write corrected reads to new BAM
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
                if read.qname == 'm84046_240619_124816_s1/231738587/ccs':
                    break
                new_seq, deam_pos = correct_by_hap(read, strand, hap, hap_basecalls)
                quals = read.query_qualities
                read.query_sequence = new_seq
                read.seq = new_seq
                read.query_qualities = quals
                if len(deam_pos) < 1:
                    read.set_tag('DA', [0])
                    read.set_tag('FD', 0, "i")
                    read.set_tag('LD', 0, "i")
                else:
                    read.set_tag('DA', deam_pos)
                    read.set_tag('FD', deam_pos[0], "i")
                    read.set_tag('LD', deam_pos[-1], "i")
                read.set_tag('ST', strand)
                if hap != 'UNK':
                    read.set_tag('HP', hap, 'i')
            else:
                read.set_tag('DA', [0])
                read.set_tag('FD', 0, "i")
                read.set_tag('LD', 0, "i")
                read.set_tag('ST', strand)
            corrected_bam.write(read)
    bam.close()
    corrected_bam.close()

