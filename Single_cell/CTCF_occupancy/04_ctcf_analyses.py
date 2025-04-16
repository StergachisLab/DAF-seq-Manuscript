import pandas as pd
import numpy as np
import itertools
import gzip
from statistics import mean,median
import csv


ctcf_df = pd.read_csv('ctcf_footprint_codes.csv')

samples = ctcf_df['cell'].unique()


# Count the number of occupied CTCF sites per cell (is it different globally)?
for cell in ctcf_df.groupby('cell'):
    cname = cell[0]
    ft_values = cell[1]['code'].value_counts()
    n_span = sum(ft_values[1:])
    n_msp = sum(ft_values[2:])
    n_bound = ft_values[2]
    # print(cname, n_span, n_msp, n_bound, n_bound/n_span, n_bound/n_msp)


# Does per-cell occupancy correlate with % deamination within the cell?
da_rate = pd.read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/figures_daf/deamination_rate_by_cell.csv')

da_rate['prop_bound'] = 0
for cell in ctcf_df.groupby('cell'):
    cname = cell[0]
    ft_values = cell[1]['code'].value_counts()
    n_msp = sum(ft_values[2:])
    n_bound = ft_values[2]
    da_rate.loc[da_rate['Cell'] == cname, 'prop_bound'] = n_bound / n_msp
    
# Pearson correlation coefficient
correlation = da_rate['prop_da_both_strands'].corr(da_rate['prop_bound'])
# print("Correlation:", correlation)

da_rate.to_csv('da_rate_by_prop_bound.csv', index=False)


# Create dictionary of haplotype-strand motif calls
hap_strand = [f'{samp}_H1' for samp in samples] + [f'{samp}_H2' for samp in samples]
hs_calls = {hs:{} for hs in hap_strand}

for samp in samples:
    for hap in ['H1','H2']:
        hs = f'{samp}_{hap}'
        temp_df = ctcf_df[(ctcf_df['cell'] == samp) & (ctcf_df['hap'] == hap)]
        code_list = temp_df.groupby('motif_ID')['code'].apply(list)
        for idx,mc in enumerate(code_list):
            if max(mc) == -1:
                hs_code = -1
            else:
                hs_code = min(code for code in mc if code >= 0)
            hs_calls[hs][code_list.index[idx]] = hs_code



# CTCF ChHI-PET ------------------------------------------------------------------------------------------------------------------------

# Fiber-seq motif data
fs_data = pd.read_csv('FiberSeq_merged_ctcf_footprint.csv')
hg002_fs = fs_data[fs_data['sample'] == 'HG002']
fs_dict = hg002_fs[['motif_ID','chrom','start','end','strand','prop_bound','prop_act']].set_index(['motif_ID']).to_dict(orient='index')

ctcf_id_dict = dict()
motif_id_dict = dict()
ctcf_df['m_coord'] = ctcf_df['chrom'] + ':' + ctcf_df['start'].astype(str) + ':' + ctcf_df['end'].astype(str)
temp_df = ctcf_df[(ctcf_df['cell'] == samples[0]) & (ctcf_df['hap'] == 'H1') & (ctcf_df['strand'] == 'Top')]
for row in temp_df.iterrows():
    ctcf_id_dict[row[1]['m_coord']] = row[1]['motif_ID']
    motif_id_dict[row[1]['motif_ID']] = row[1]['m_coord']


# Parse ChIA-PET data
loops = pd.read_csv('ENCODE/GM12878_CTCF_ChIA-PET_Loops_ENCFF780PGS.bedpe.gz', sep="\t", names=['chrom1','start1','end1','chrom2','start2','end2','pet_score'])
loops['m_coord_1'] = loops['chrom1'] + ':' + loops['start1'].astype(str) + ':' + loops['end1'].astype(str)
loops['m_coord_2'] = loops['chrom2'] + ':' + loops['start2'].astype(str) + ':' + loops['end2'].astype(str)


# track motif IDs in each half of the ChIA-PET loops
reg1_dict = {}
with gzip.open('Loops_first_regions_map.txt.gz', 'rt') as fr:
    for row in fr:
        chia,mapped = row.strip().split('|')
        chia_id = ':'.join(chia.split('\t'))
        mapped = mapped.split(';')
        reg1_dict[chia_id] = {ctcf_id_dict[':'.join(mot.split('\t')[:3])] for mot in mapped}

reg2_dict = {}
with gzip.open('Loops_second_regions_map.txt.gz', 'rt') as fr:
    for row in fr:
        chia,mapped = row.strip().split('|')
        chia_id = ':'.join(chia.split('\t'))
        mapped = mapped.split(';')
        reg2_dict[chia_id] = {ctcf_id_dict[':'.join(mot.split('\t')[:3])] for mot in mapped}


# create dict of motif:motif pairs with ChIA-PET scores for fast lookups
m_combo_dict = dict()
for loop in loops.iterrows():
    reg1 = loop[1]['m_coord_1']
    reg2 = loop[1]['m_coord_2']
    if reg1 in reg1_dict.keys() and reg2 in reg2_dict.keys():
        reg_combos = itertools.product(reg1_dict[reg1], reg2_dict[reg2])
        for rc in reg_combos:
            m_combo_dict[rc] = loop[1]['pet_score']
        # add combos in the reverse orientation incase the motif combos are not ordered by genome coordinates
        reg_combos = itertools.product(reg2_dict[reg2], reg1_dict[reg1])
        for rc in reg_combos:
            m_combo_dict[rc] = loop[1]['pet_score']



# calculate pair-wise co-occupancy
min_act = 0.3

# minimum number of fibers to count the pair
min_motif_cov = 4

ctcf_df_filt = ctcf_df[(ctcf_df['code'] >= 0) & (ctcf_df['hap'].isin(['H1','H2']))]

tot_motifs = 0
combos = []
for ch in ctcf_df_filt['chrom'].unique():
    ch_motifs = ctcf_df_filt[ctcf_df_filt['chrom'] == ch]['motif_ID'].unique()
    ch_motifs = [mot for mot in ch_motifs if (fs_dict[mot]['prop_act'] >=  min_act)] # filter by percent actuation
    tot_motifs += len(ch_motifs)
    for combo in itertools.combinations(ch_motifs, 2):
        if combo[0] != combo[1]:
            combos.append(combo)


# initialize codep DF
codep_df_rows = []
for c in combos:
    m1,m2 = c
    p1 = fs_dict[m1]['prop_bound']
    p2 = fs_dict[m2]['prop_bound']
    counts = [0,0,0,0]
    for k,v in hs_calls.items():
        c1 = v[m1]
        c2 = v[m2]
        # codep for occupancy, fiber must have a MSP spanning each motif
        if c1 >= 1 and c2 >= 1:
            if c1 == 2 and c2 == 2: # both
                counts[3] += 1
            elif c1 == 2: # first only
                counts[1] += 1
            elif c2 == 2: # second only
                counts[2] += 1
            else: # neither
                counts[0] += 1
    new_row = [m1, fs_dict[m1]['chrom'], fs_dict[m1]['start'], fs_dict[m1]['end'], fs_dict[m1]['strand'], p1, m2, fs_dict[m2]['chrom'], fs_dict[m2]['start'], fs_dict[m2]['end'], fs_dict[m2]['strand'], p2] + counts
    if sum(counts) >= min_motif_cov:
        co_occ = counts[3] / sum(counts)
    else:
        co_occ = None
    new_row.append(co_occ)
    # ADD ChiA-PET scores
    score = 0 # There should not be multiple motif:motif interactions in ChIA-PET data
    if (m1,m2) in m_combo_dict.keys():
        score = m_combo_dict[(m1,m2)]
    elif (m2,m1) in m_combo_dict.keys():
        score = m_combo_dict[(m2,m1)]
    new_row.append(score)
    codep_df_rows.append(new_row)

chia_coocc_df = pd.DataFrame(codep_df_rows, columns=['motif_1','chrom_1','start_1','end_1','strand_1','prop_1','motif_2','chrom_2','start_2','end_2','strand_2','prop_2','neither','first','second','both','prop_co_occ','pet_score'])
chia_coocc_df.to_csv('ctcf_ChIA-PET_co-occupancy.csv', index=False)



# ------------------------------------------------------------------------------------------------------------------------

loops['ID'] = loops.index + 1

# output new DF with one row per loop / hap_strand with occ calls
loop_motif_rows = [] # loop_ID, hap_strand, reg1_calls, reg2_calls, reg1_motifs, reg2_motifs, reg1_strand, reg2_strand

moitf_coocc_strand_rows = []
loop_coocc_strand_rows = []
loop_num_cooc = []

for loop in loops.iterrows():
    reg1 = loop[1]['m_coord_1']
    reg2 = loop[1]['m_coord_2']
    if reg1 in reg1_dict.keys() and reg2 in reg2_dict.keys():
        reg1_motifs = sorted([m for m in reg1_dict[reg1] if fs_dict[m]['prop_act'] >=  min_act]) # filter motifs by FS actuation
        reg2_motifs = sorted([m for m in reg2_dict[reg2] if fs_dict[m]['prop_act'] >=  min_act]) # filter motifs by FS actuation
        reg1_strands = [fs_dict[m]['strand'] for m in reg1_motifs]
        reg2_strands = [fs_dict[m]['strand'] for m in reg2_motifs]
        # track loop counts by motif-strand combo
        mot_strands_tot = {'--':0, '-+':0, '+-':0, '++':0}
        mot_strands_occ = {'--':0, '-+':0, '+-':0, '++':0}
        any_strands_tot = {'--':set(), '-+':set(), '+-':set(), '++':set()}
        any_strands_occ = {'--':set(), '-+':set(), '+-':set(), '++':set()}
        # raw count of how many co-occ per hap_strand
        num_co_occ = {'--':[], '-+':[], '+-':[], '++':[]}
        for hs,calls in hs_calls.items():
            hs_num_co_occ = {}
            r1_calls = [calls[m] for m in reg1_motifs]
            r2_calls = [calls[m] for m in reg2_motifs]
            if -1 not in r1_calls and -1 not in r2_calls:
                loop_row = [loop[1]['ID'], hs, r1_calls, r2_calls, reg1_motifs, reg2_motifs, reg1_strands, reg2_strands, reg1, reg2]
                loop_motif_rows.append(loop_row)
                for i in range(len(reg1_motifs)):
                    for j in range(len(reg2_motifs)):
                        st_combo = reg1_strands[i] + reg2_strands[j]
                        if r1_calls[i] > 0 and r2_calls[j] > 0:
                            if st_combo not in hs_num_co_occ.keys():
                                hs_num_co_occ[st_combo] = 0
                            if r1_calls[i] == 2 and r2_calls[j] == 2:
                                mot_strands_occ[st_combo] += 1
                                any_strands_occ[st_combo].add(hs)
                                hs_num_co_occ[st_combo] += 1
                            mot_strands_tot[st_combo] += 1
                            any_strands_tot[st_combo].add(hs)
            for k,v in hs_num_co_occ.items():
                num_co_occ[k].append(v)
        for k,v in mot_strands_tot.items():
            if v >= min_motif_cov:
                moitf_coocc_strand_rows.append([loop[1]['ID'], loop[1]['pet_score'], k, v, mot_strands_occ[k], reg1, reg2])
        # add >= 1 motif co-occ data (prop of fibers with >= 1 co-occ per orientation)
        for k,v in any_strands_tot.items():
            if len(v) >= min_motif_cov:
                loop_coocc_strand_rows.append([loop[1]['ID'], loop[1]['pet_score'], k, len(v), len(any_strands_occ[k]), reg1, reg2])
        # add number of co-occ per strand combo
        for k,v in mot_strands_tot.items():
            if v > 0 and len(num_co_occ[k]) > 0:
                loop_num_cooc.append([loop[1]['ID'], loop[1]['pet_score'], k, len(num_co_occ[k]), mean(num_co_occ[k]), median(num_co_occ[k]), num_co_occ[k], reg1, reg2])


loop_occ_combos_df = pd.DataFrame(loop_motif_rows, columns=['loop_ID','hap_strand','reg1_calls','reg2_calls','reg1_motifs','reg2_motifs','reg1_strands','reg2_strands','reg1_coord','reg2_coord'])
loop_occ_combos_df.to_csv('loop_motif_calls_full_fiber.tsv', sep="\t", index = False)



# (2) for regions with >0 ctcf sites, compare unique anchor vs non-unique anchor (>1 ctcf site in region) AND (3) Combination of occupancy in ChIA-PET regions when >1 ctcf site within either region of pair

# output dataframes with one row per loop & strand combo: -----------------------

# reg1 & reg2 co-occ by strand
moitf_coocc_strand_df = pd.DataFrame(moitf_coocc_strand_rows, columns=['loop_ID','pet_score','mot_strands','total','num_co_occ','reg1_coord','reg2_coord'])
moitf_coocc_strand_df['prop_co_occ'] = moitf_coocc_strand_df['num_co_occ'] / moitf_coocc_strand_df['total']
moitf_coocc_strand_df.to_csv('loop_coocc_by_strand_exact.tsv', sep="\t", index = False)

# co-occ as >= 1 co-occ by strand
loop_coocc_strand_df = pd.DataFrame(loop_coocc_strand_rows, columns=['loop_ID','pet_score','mot_strands','total','num_co_occ', 'reg1_coord','reg2_coord'])
loop_coocc_strand_df['prop_co_occ'] = loop_coocc_strand_df['num_co_occ'] / loop_coocc_strand_df['total']
loop_coocc_strand_df.to_csv('loop_coocc_by_strand_any.tsv', sep="\t", index = False)


# absolute counts of number of co-occ motifs by strand orientation per hap_strand fiber
loop_abs_nCoOcc_strand_df = pd.DataFrame(loop_num_cooc, columns=['loop_ID','pet_score','mot_strands','total','mean_num_co_occ','median_num_co_occ','raw_co_occ_counts', 'reg1_coord','reg2_coord'])
loop_abs_nCoOcc_strand_df.to_csv('loop_abs_nCoOcc_strand_df.tsv', sep="\t", index = False)


# Repeat Loop >1 for shuffled loop anchor regions (NULL distribution) ---------------------------------------------------------------------

# track motif IDs in each half of the ChIA-PET loops
reg1_dict_shuff = {}
with gzip.open('Loops_first_regions_SHUFFLED_map.txt.gz', 'rt') as fr:
    for row in fr:
        chia,mapped = row.strip().split('|')
        chia_id = ':'.join(chia.split('\t'))
        mapped = mapped.split(';')
        reg1_dict_shuff[chia_id] = {ctcf_id_dict[':'.join(mot.split('\t')[:3])] for mot in mapped}

reg2_dict_shuff = {}
with gzip.open('Loops_second_regions_SHUFFLED_map.txt.gz', 'rt') as fr:
    for row in fr:
        chia,mapped = row.strip().split('|')
        chia_id = ':'.join(chia.split('\t'))
        mapped = mapped.split(';')
        reg2_dict_shuff[chia_id] = {ctcf_id_dict[':'.join(mot.split('\t')[:3])] for mot in mapped}


# create loop anchor pairs from combos of reg1 and reg2 regions
loops_shuff_rows = []
reg_combos = itertools.product(reg1_dict_shuff.keys(), reg2_dict_shuff.keys())
for rc in reg_combos:
    chrom1,start1,end1 = rc[0].split(':')
    chrom2,start2,end2 = rc[1].split(':')
    if chrom1 == chrom2 and int(end1) < int(start2):
        loops_shuff_rows.append([chrom1, start1, end1, chrom2, start2, end2, 0])

loops_shuff = pd.DataFrame(loops_shuff_rows, columns=['chrom1','start1','end1','chrom2','start2','end2','pet_score'])
loops_shuff['m_coord_1'] = loops_shuff['chrom1'] + ':' + loops_shuff['start1'].astype(str) + ':' + loops_shuff['end1'].astype(str)
loops_shuff['m_coord_2'] = loops_shuff['chrom2'] + ':' + loops_shuff['start2'].astype(str) + ':' + loops_shuff['end2'].astype(str)

loops_shuff['ID'] = loops_shuff.index + 1

# output new DF with one row per loop / hap_strand with occ calls
loop_coocc_strand_rows_shuff = []
for loop in loops_shuff.iterrows():
    reg1 = loop[1]['m_coord_1']
    reg2 = loop[1]['m_coord_2']
    if reg1 in reg1_dict_shuff.keys() and reg2 in reg2_dict_shuff.keys():
        reg1_motifs = sorted([m for m in reg1_dict_shuff[reg1] if fs_dict[m]['prop_act'] >=  min_act]) # filter motifs by FS actuation !!!!!!!!!!!!!!
        reg2_motifs = sorted([m for m in reg2_dict_shuff[reg2] if fs_dict[m]['prop_act'] >=  min_act]) # filter motifs by FS actuation !!!!!!!!!!!!!!
        reg1_strands = [fs_dict[m]['strand'] for m in reg1_motifs]
        reg2_strands = [fs_dict[m]['strand'] for m in reg2_motifs]
        # track loop counts by motif-strand combo
        any_strands_tot = {'--':set(), '-+':set(), '+-':set(), '++':set()}
        any_strands_occ = {'--':set(), '-+':set(), '+-':set(), '++':set()}
        for hs,calls in hs_calls.items():
            r1_calls = [calls[m] for m in reg1_motifs]
            r2_calls = [calls[m] for m in reg2_motifs]
            if -1 not in r1_calls and -1 not in r2_calls:
                for i in range(len(reg1_motifs)):
                    for j in range(len(reg2_motifs)):
                        st_combo = reg1_strands[i] + reg2_strands[j]
                        if r1_calls[i] > 0 and r2_calls[j] > 0:
                            if r1_calls[i] == 2 and r2_calls[j] == 2:
                                any_strands_occ[st_combo].add(hs)
                            any_strands_tot[st_combo].add(hs)
        # add >= 1 motif co-occ data (prop of fibers with >= 1 co-occ per orientation)
        for k,v in any_strands_tot.items():
            if len(v) >= min_motif_cov:
                loop_coocc_strand_rows_shuff.append([loop[1]['ID'], loop[1]['pet_score'], k, len(v), len(any_strands_occ[k]), reg1, reg2])

# output dataframes with one row per loop & strand combo: -----------------------

# co-occ as >= 1 co-occ by strand
loop_coocc_strand_SHUFF_df = pd.DataFrame(loop_coocc_strand_rows_shuff, columns=['loop_ID','pet_score','mot_strands','total','num_co_occ', 'reg1_coord','reg2_coord'])
loop_coocc_strand_SHUFF_df['prop_co_occ'] = loop_coocc_strand_SHUFF_df['num_co_occ'] / loop_coocc_strand_SHUFF_df['total']
loop_coocc_strand_SHUFF_df.to_csv('loop_coocc_by_strand_any_SHUFFLED.tsv', sep="\t", index = False)
