import csv
import pysam
from glob import glob


# one DA dict for each chromosome
ct_dicts = {f'chr{i}':{} for i in range(1,23)}
ct_dicts['chrX'] = {}
ct_dicts['chrY'] = {}

ga_dicts = {f'chr{i}':{} for i in range(1,23)}
ga_dicts['chrX'] = {}
ga_dicts['chrY'] = {}


# 2) for each BAM (aggregate data), count basecalls at each REF C or G position ------------------------------------------------------------------------------------

min_seq_len = 4000
read_buff_len = 300 # don't count bases from positions X bp from the ends of the read to avoid double counting collapse failures.

bam_dir = '../collapse/consensus_bams'
bams = glob(f'{bam_dir}/PS*_consensus_BothStrands_HG38_corrected.bam')

for i in range(len(bams)):
    bam = pysam.AlignmentFile(f'{bams[i]}', "rb")
    for read in bam.fetch():
        seq = read.query_sequence
        if len(seq) >= min_seq_len:
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            strand = read.get_tag('ST')
            chrom = read.reference_name
            if chrom in ct_dicts.keys():
                for pos in pair:
                    if pos[0] == None or pos[1] == None: # indel, ignore
                        pass
                    else:
                        qi = pos[0]
                        ref_coord = pos[1] + 1
                        ref_base = pos[2].upper()
                        base = seq[qi]
                        if base != 'N' and qi > read_buff_len and qi < (len(seq)-read_buff_len):
                            if strand == 'CT':
                                if ref_base == 'C':
                                    if ref_coord not in ct_dicts[chrom].keys():
                                        ct_dicts[chrom][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                                    ct_dicts[chrom][ref_coord][base] += 1
                            elif strand == 'GA':
                                if ref_base == 'G':
                                    if ref_coord not in ga_dicts[chrom].keys():
                                        ga_dicts[chrom][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                                    ga_dicts[chrom][ref_coord][base] += 1
    print(f'{bams[i]} COMPLETED')
    bam.close()


# Output TSV of per position count of DA and TOTAL COV

# Columns: chrom, start, end, n_unmod, n_deam
with open('da_counts_by_pos_all_cells.tsv','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for chrom,data in ct_dicts.items():
        for pos,counts in data.items():
            line = [chrom, pos-1, pos, counts['C'], counts['Y']]
            writer.writerow(line)
    for chrom,data in ga_dicts.items():
        for pos,counts in data.items():
            line = [chrom, pos-1, pos, counts['G'], counts['R']]
            writer.writerow(line)

