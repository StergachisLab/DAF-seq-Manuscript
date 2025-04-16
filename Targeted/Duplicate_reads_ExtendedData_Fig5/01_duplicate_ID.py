import pysam
import argparse
import csv
import os
import statistics as st


"""
Identify PCR duplciates in DAF-seq reads. Requires ST tag (CT or GA strand) to be set. Reads must contain ambiguous R/Y bases in their sequence.

1) Identify PCR duplicates based on the deamination status at each position (C|T: Y, G|A: R).
2) Group duplicate reads.
3) Output list of deduped ZMWs to subset BAM (one ZMW per duplicate group).
4) Output duplication stats.
"""

# parse command line arguments
parser = argparse.ArgumentParser(description = "DddA BAM preprocessing", epilog = "")
parser.add_argument("-b", "--bam", required = True, metavar = '', help = "DddA aligned BAM to correct")
parser.add_argument("-n", "--name", required = False, metavar = '', help = "Name to append to outputs")
args = parser.parse_args()

# identify fastq files in dir
bam_name = args.bam
if args.name:
    name = args.name
else:
    name = os.path.basename(bam_name).split('.')[0]

bam = pysam.AlignmentFile(bam_name, "rb")


# use order of ambiguous bases (Y/R)
# negative value for reference base, positive value for deamination
hash_dict = dict()
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False:
        strand = read.get_tag('ST')
        # C --> T
        if strand == 'CT':
            seq = read.seq
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            da_hash = []
            for pos in pair:
                mcoord = pos[0]
                if mcoord != None and pos[1] != None:
                    if pos[2].upper() == 'C': # reference base must be a C to be counted
                        if seq[mcoord] == 'Y':
                            da_hash.append(pos[1])
                        else:
                            da_hash.append(-pos[1])
        # G --> A
        elif strand == 'GA':
            seq = read.seq
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            da_hash = []
            for pos in pair:
                mcoord = pos[0]
                if mcoord != None and pos[1] != None:
                    if pos[2].upper() == 'G': # reference base must be a G to be counted
                        if seq[mcoord] == 'R':
                            da_hash.append(pos[1])
                        else:
                            da_hash.append(-pos[1])
        if strand in ['CT','GA']:              
            key = tuple(da_hash)
            if key not in hash_dict.keys():
                hash_dict[key] = {'count':1, 'zmws':[read.qname]}
            else:
                hash_dict[key]['count'] += 1
                hash_dict[key]['zmws'].append(read.qname)

# write groups of duplicates
with open(f'{name}_DAF_duplicate_reads.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    header = ['duplicate_count','zmws']
    writer.writerow(header)
    for k,v in hash_dict.items():
        writer.writerow([v['count'], ';'.join(v['zmws'])])

# write DeDuped ZMWs
# the first ZMW in the group will be retained
with open(f'{name}_DAF_DeDuplicated_zmws.txt','w') as fw:
    for k,v in hash_dict.items():
        fw.write(f"{v['zmws'][0]}\n")

# write deduplication stats
with open(f'{name}_DAF_DeDuplicated_stats.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    header = ['total_reads','duplicate_reads','percent_dup', 'mean_copies','median_copies','max_copies']
    writer.writerow(header)
    total=0
    dup=0
    n_copies=[]
    for k,v in hash_dict.items():
        total += v['count']
        dup += (v['count']-1)
        n_copies.append(v['count'])
    writer.writerow([total, dup, (dup/total)*100, st.mean(n_copies), st.median(n_copies), max(n_copies)])
