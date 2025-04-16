import csv
from itertools import combinations
from collections import Counter
import re
import argparse
import os
from math import ceil
import sys


# parse command line arguments
parser = argparse.ArgumentParser(description = "Exctract scDAF-seq read sequences for collapse", epilog = "")
parser.add_argument("-r", "--reads", required = True, metavar = '', help = "Clipped reads to collapse (TSV)")
parser.add_argument("-c", "--chr", required = True, metavar = '', help = "Chromosome to extract from")
args = parser.parse_args()

# identify fastq files in dir
read_file = args.reads
chrom = args.chr
name = os.path.basename(read_file).replace('.tsv','')
if 'CT' in read_file:
    strand = 'CT'
elif 'GA' in read_file:
    strand = 'GA'
else:
    sys.exit(1)


# chromosome size for bins /mmfs1/gscratch/stergachislab/assemblies/simple-names/hg38.fa.fai
sizes={'chr1':248956422,'chr2':242193529,'chr3':198295559,'chr4':190214555,'chr5':181538259,'chr6':170805979,'chr7':159345973,'chr8':145138636,'chr9':138394717,'chr10':133797422,'chr11':135086622,'chr12':133275309,
'chr13':114364328,'chr14':107043718,'chr15':101991189,'chr16':90338345,'chr17':83257441,'chr18':80373285,'chr19':58617616,'chr20':64444167,'chr21':46709983,'chr22':50818468,'chrX':156040895,'chrY':57227415}

cutoff = 0.99 # hamming dist identity cutoff
prop_bin_id = 0.80 # proportion of bins with >= hamming cutoff to group together

# define bins for region
bin_len=150
slide=25
bins=list(range(1, sizes[chrom], slide))
min_bins=11 # minimum overlapping bins to group reads (11 bins = 400 bp)
min_cons_bins=7 # Allow shorter overlap of 7 bins (300 bp) for grouping consensus sequences

reads = dict()
with open(read_file) as fr:
    reader = csv.DictReader(fr, delimiter="\t", fieldnames=['zmw','chr','start','end','seq'])
    for row in reader:
        if row['chr'] == chrom:
            zmw = row['zmw']
            del row['zmw']
            row['start'] = int(row['start'])
            row['end'] = int(row['end'])
            reads[zmw] = row

# assign reads to bins
bin_reads = dict()
for r in reads:
    cbin = ceil(reads[r]['start']/slide)
    while reads[r]['start'] <= bins[cbin] and reads[r]['end'] >= (bins[cbin]+bin_len):
        if bins[cbin] not in bin_reads.keys():
            bin_reads[bins[cbin]] = {r}
        else:
            bin_reads[bins[cbin]].add(r)
        cbin += 1

def hamming_dist(s1, s2):
    if len(s1) != len(s2):
        return(len(s1)) # return high dist value
    dist = 0
    for i in range(len(s1)):
        if s1[i] != s2[i] and s1[i] != 'D' and s2[i] != 'D': # don't count deletion bases in hamming distance
            dist += 1
    return(dist)

# find graph nodes and edges
nodes = {}
for b in bin_reads.keys():
    for c in combinations(bin_reads[b], 2): # compute hamming dist for pairs of reads by bin
        r1,r2 = c
        r1_st = b - reads[r1]['start'] - 1
        r2_st = b - reads[r2]['start'] - 1
        # ignore ins bases for hamming distance
        r1_ref = re.sub('[a,c,g,t]', '', reads[r1]['seq'])
        r2_ref = re.sub('[a,c,g,t]', '', reads[r2]['seq'])
        r1_seq = r1_ref[r1_st:r1_st+bin_len]
        r2_seq = r2_ref[r2_st:r2_st+bin_len]
        # account for insertion bases, compare only ref coordinates
        dist = hamming_dist(r1_seq, r2_seq)
        if r1 not in nodes.keys():
            nodes[r1] = {}
        if r2 not in nodes.keys():
            nodes[r2] = {}
        if b not in nodes[r1].keys():
            nodes[r1][b] = {}
        if b not in nodes[r2].keys():
            nodes[r2][b] = {}
        nodes[r1][b][r2] = dist
        nodes[r2][b][r1] = dist

# find neighborhoods of fibers
neighbors = {z:{} for z in nodes.keys()}
for k,v in nodes.items():
    for k2,v2 in v.items():
        for k3,v3 in v2.items():
            prop_id = (1-(v3/bin_len))
            if k3 not in neighbors[k].keys():
                neighbors[k][k3] = [prop_id]
            else:
                neighbors[k][k3].append(prop_id)

# create read groupings
def get_group(zmw, neighbors, cutoff, prop_bin_id, min_bins):
    neighbor_props = neighbors[zmw]
    group=set()
    for k2,v2 in neighbor_props.items():
        p_cutoff = len([p for p in v2 if p >= cutoff])/len(v2)
        if p_cutoff >= prop_bin_id and len(v2) >= min_bins:
            group.add(k2)
    return(group)

groups = []
for k,v in neighbors.items():
    g2 = get_group(k, neighbors, cutoff, prop_bin_id, min_bins)
    g2.add(k)
    groups.append(g2)

def merge_sets(sets):
    merged = []
    while sets:
        first, *rest = sets
        first = set(first)
        lf = -1
        while len(first) > lf:
            lf = len(first)
            rest2 = []
            for r in rest:
                if first.intersection(r):
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2
        merged.append(first)
        sets = rest
    return merged

groups2 = merge_sets(groups)


final_groups = []
for orig_group in groups2:
    regrouped = set() # skip over fibers that have been re-grouped
    new_groups = []
    for fiber in orig_group: # start new group
        if fiber not in regrouped:
            new_set = {fiber}
            regrouped.add(fiber)
            for fiber2 in orig_group:
                if fiber2 not in regrouped:
                    # compare new fiber against all other fibers in the new set
                    agree = True
                    fib_overlap = False # fiber2 must overlap at least one fiber in the set to be added
                    for ns_fiber in new_set:
                        if ns_fiber in neighbors[fiber2].keys():
                            fib_overlap = True
                            props = neighbors[fiber2][ns_fiber]
                            p_cutoff = len([p for p in props if p >= cutoff])/len(props)
                            if p_cutoff < prop_bin_id: # disagreement of fiber2 with one of the fibers in the new set. Move on
                                agree = False
                                break
                    if agree == True and fib_overlap == True:
                        new_set.add(fiber2)
                        regrouped.add(fiber2)
            new_groups.append(new_set)
    # repeat merging step with new_groups
    niter = 2
    old_groups = new_groups
    for i in range(niter):
        old_groups = new_groups
        regrouped = set()
        new_groups = []
        for group in old_groups:
            if len(group.intersection(regrouped)) == 0:
                new_set = group
                for group2 in old_groups:
                    if len(group2.intersection(regrouped)) == 0:
                        agree = True
                        fib_overlap = False # groups must overlap by at least one fiber
                        # compare each fiber between the two groups
                        for fiber in group:
                            if agree == False:
                                break
                            for fiber2 in group2:
                                if fiber in neighbors[fiber2].keys():
                                    fib_overlap = True
                                    props = neighbors[fiber][fiber2]
                                    p_cutoff = len([p for p in props if p >= cutoff])/len(props)
                                    if p_cutoff < prop_bin_id: # disagreement of fiber2 with one of the fibers in the new set. Move on
                                        agree = False
                                        break
                        if agree == True and fib_overlap == True:
                            new_set = new_set.union(group2)
                            regrouped = regrouped.union(group, group2)
                new_groups.append(new_set)
                regrouped = regrouped.union(group)
    for ng in new_groups:
        final_groups.append(ng)

# add reads that shared no bins with other reads
for r in reads.keys():
    if r not in nodes.keys():
        final_groups.append({r})

# create merged, final consensus sequence of grouped reads
def consensus_seq(group, read_dict, base_cutoff, strand):
    if strand == 'CT':
        da = ['C','T']
    elif strand == 'GA':
        da = ['G','A']
    else:
        return(None)
    con_first = min([read_dict[r]['start'] for r in group])
    con_end = max([read_dict[r]['end'] for r in group])
    new_seq = ''
    i = con_first
    while i < (con_end+1):
        bases = []
        for read in group:
            st = read_dict[read]['start']
            end = read_dict[read]['end']
            if (st <= i and end >= i):
                pos_base = read_dict[read]['seq'][i-st]
                if pos_base.islower():
                    ins_bases = []
                    for ins_base in read_dict[read]['seq'][i-st:]:
                        if ins_base.islower():
                            ins_bases.append(ins_base)
                        else:
                            break
                    ins = ''.join(ins_bases)
                    bases.append(ins)
                    # remove insertion bases from read in dict!!!
                    # avoids offset issue with uppercase Ref bases further on in the read.
                    read_dict[read]['seq'] = read_dict[read]['seq'].replace(ins,'',1)
                else:
                    bases.append(pos_base)
        clower = [cl for cl in bases if cl.islower()]
        if len(clower) > 0:
            if len(clower) == len(bases): # ins in all reads, take the shortest
                new_seq += min(clower, key=len)
                # NO NEED TO INCREASE i for ins, ref does not include
            else: # if ins in some but not all reads, repeat loop with ins removed
                continue
        # if no ins, take base
        else:
            bcounts = Counter(bases)
            # disagreements ----------------------------
            if len(bcounts) > 1:
                if (bcounts.most_common(1)[0][1] / sum(bcounts.values())) >= base_cutoff: # take consensus base
                    new_base = bcounts.most_common(1)[0][0]
                else:
                    top2_b = [mc[0] for mc in bcounts.most_common(2)]
                    top2_c = [mc[1] for mc in bcounts.most_common(2)]
                    if da[0] in top2_b and da[1] in top2_b and (sum(top2_c)/sum(bcounts.values()) >= base_cutoff): # C/T --> take C, or G/A --> take G, could be amplicon deamination
                        new_base = da[0]
                    elif 'D' in bcounts.keys():
                        new_base = [b for b,c in bcounts.most_common(2) if b != 'D'][0]
                    elif 'N' in top2_b and (sum([c for b,c in bcounts.most_common(2) if b != 'N'])/sum([c for b,c in bcounts.items() if b != 'N']) >= base_cutoff): # N/other --> take other
                        new_base = [b for b,c in bcounts.most_common(2) if b != 'N'][0]
                    else: # mix --> take N
                        new_base = 'N'
            # single basecall ----------------------------
            else:
                new_base = bases[0]
            # add new base
            i += 1
            new_seq += new_base
    return(new_seq, con_first, con_end)


# Repeat collapse on the consensus sequences
# Where there is a disgreement between grouped consensus sequences take base from the consensus group with the highest read count WITHIN THAT BIN
out_reads = []
hash_zmw = []
cons_zmw_dict = dict()
for g in final_groups:
    cs,con_first,con_end = consensus_seq(g, reads, 0.5, strand)
    joined_hash = str(hash('-'.join(list(g))))
    if cs != None:
        out_reads.append([joined_hash, chrom, con_first, con_end, cs, joined_hash])
        hash_zmw.append([joined_hash, f'{chrom}:{con_first}-{con_end}', ';'.join(sorted(g))])
        cons_zmw_dict[joined_hash] = set(g)

# bin consensus seqs
cons_reads = dict()
for i in range(len(final_groups)):
    read = out_reads[i]
    cons_name = read[0]
    cons_reads[cons_name] = {'chr':read[1], 'start':read[2], 'end':read[3], 'seq':read[4], 'group_num':len(hash_zmw[i][2].split(';'))}

# assign consensus reads to bins
bin_cons_reads = dict()
for r in cons_reads:
    cbin = ceil(cons_reads[r]['start']/slide)
    while cons_reads[r]['start'] <= bins[cbin] and cons_reads[r]['end'] >= (bins[cbin]+bin_len):
        if bins[cbin] not in bin_cons_reads.keys():
            bin_cons_reads[bins[cbin]] = {r}
        else:
            bin_cons_reads[bins[cbin]].add(r)
        cbin += 1

# group overlapping consensus seqs 
cons_nodes = {}
for b in bin_cons_reads.keys():
    for c in combinations(bin_cons_reads[b], 2): # compute hamming dist for pairs of reads by bin
        r1,r2 = c
        r1_st = b - cons_reads[r1]['start'] - 1
        r2_st = b - cons_reads[r2]['start'] - 1
        # ignore ins bases for hamming distance
        r1_ref = re.sub('[a,c,g,t]', '', cons_reads[r1]['seq'])
        r2_ref = re.sub('[a,c,g,t]', '', cons_reads[r2]['seq'])
        r1_seq = r1_ref[r1_st:r1_st+bin_len]
        r2_seq = r2_ref[r2_st:r2_st+bin_len]
        # account for insertion bases, compare only ref coordinates
        dist = hamming_dist(r1_seq, r2_seq)
        if r1 not in cons_nodes.keys():
            cons_nodes[r1] = {}
        if r2 not in cons_nodes.keys():
            cons_nodes[r2] = {}
        if b not in cons_nodes[r1].keys():
            cons_nodes[r1][b] = {}
        if b not in cons_nodes[r2].keys():
            cons_nodes[r2][b] = {}
        cons_nodes[r1][b][r2] = dist
        cons_nodes[r2][b][r1] = dist

# find neighborhoods of fibers
cons_neighbors = {z:{} for z in cons_nodes.keys()}
for k,v in cons_nodes.items():
    for k2,v2 in v.items():
        for k3,v3 in v2.items():
            prop_id = (1-(v3/bin_len))
            if k3 not in cons_neighbors[k].keys():
                cons_neighbors[k][k3] = [prop_id]
            else:
                cons_neighbors[k][k3].append(prop_id)

cons_groups = []
for k,v in cons_neighbors.items():
    g2 = get_group(k, cons_neighbors, cutoff, prop_bin_id, min_cons_bins)
    g2.add(k)
    cons_groups.append(g2)
cons_groups2 = merge_sets(cons_groups)

# add consensus reads that shared no bins with other reads
for r in cons_reads.keys():
    if r not in cons_nodes.keys():
        cons_groups2.append({r})


# merge consensus seqs into final consensus seqs
# create merged, final consensus sequence of grouped reads
def consensus_merge(group, read_dict, strand, slide, cons_zmw_dict, bin_reads, bins):
    if strand == 'CT':
        da = ['C','T']
    elif strand == 'GA':
        da = ['G','A']
    else:
        return(None)
    con_first = min([read_dict[r]['start'] for r in group])
    con_end = max([read_dict[r]['end'] for r in group])
    new_seq = ''
    i = con_first
    while i < (con_end+1):
        bases = []
        read_names = []
        for read in group:
            st = read_dict[read]['start']
            end = read_dict[read]['end']
            if (st <= i and end >= i):
                pos_base = read_dict[read]['seq'][i-st]
                if pos_base.islower():
                    ins_bases = []
                    for ins_base in read_dict[read]['seq'][i-st:]:
                        if ins_base.islower():
                            ins_bases.append(ins_base)
                        else:
                            break
                    ins = ''.join(ins_bases)
                    bases.append(ins)
                    read_names.append(read)
                    # remove insertion bases from read in dict!!!
                    # avoids offset issue with uppercase Ref bases further on in the read.
                    read_dict[read]['seq'] = read_dict[read]['seq'].replace(ins,'',1)
                else:
                    bases.append(pos_base)
                    read_names.append(read)
        clower = [cl for cl in bases if cl.islower()]
        if len(clower) > 0:
            if len(clower) == len(bases): # ins in all reads, take the shortest
                new_seq += min(clower, key=len)
                # NO NEED TO INCREASE i for ins, ref does not include
            else: # if ins in some but not all reads, repeat loop with ins removed
                continue
        # if no ins, take base
        else:
            bcounts = Counter([b for b in bases if b.isupper()])
            # disagreements ----------------------------
            if len(bcounts) > 1:
                if 'D' in bcounts.keys():
                    new_base = [b for b,c in bcounts.most_common(2) if b != 'D'][0]
                else: # mix --> take the base from the consensus seq with the highest number of raw reads within that bin
                    cbin = ceil(i/slide)
                    if bins[cbin] in bin_reads.keys():
                        num_bin_zmw = [len(cons_zmw_dict[name].intersection(bin_reads[bins[cbin]])) for name in read_names]
                        max_zmw_idx = num_bin_zmw.index(max(num_bin_zmw))
                        while bases[max_zmw_idx].islower():
                            num_bin_zmw.pop(max_zmw_idx)
                            bases.pop(max_zmw_idx)
                            max_zmw_idx = num_bin_zmw.index(max(num_bin_zmw))
                        new_base = bases[max_zmw_idx]
                    else: # when cbin is not in bin keys take the most common base
                        new_base = list(bcounts.keys())[0]
            # single basecall ----------------------------
            else:
                new_base = bases[0]
            # add new base
            i += 1
            if new_base != 'D': # don't write ref deletion bases (already wrote ins bases)
                new_seq += new_base
    return(new_seq, con_first, con_end)


# output consensus reads in fasta format using a group ZMW hash as the read name
out_reads2 = []
hash_zmw2 = []
for g in cons_groups2:
    cs,con_first,con_end = consensus_merge(g, cons_reads, strand, slide, cons_zmw_dict, bin_reads, bins)
    joined_hash = hash('-'.join(list(g)))
    if cs != None:
        out_reads2.append([joined_hash, chrom, con_first, con_end, cs.upper(), joined_hash])
        hash_zmw2.append([joined_hash, f'{chrom}:{con_first}-{con_end}', ';'.join(sorted(g))])

out_fasta = f'consensus_seqs/{name}_{chrom}_CONSENSUS.fa'
with open(out_fasta, 'w') as fw:
    for con_read in out_reads2:
        header = f'>{con_read[0]}'
        fw.write(f'{header}\n')
        fw.write(f'{con_read[4]}\n')

out_hash_groups = f'consensus_seqs/{name}_{chrom}_HASH_ZMW_Groups.tsv'
with open(out_hash_groups, 'w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(['hash','consensus_coord','grouped_ZMWs'])
    for group in hash_zmw2:
        writer.writerow(group)

