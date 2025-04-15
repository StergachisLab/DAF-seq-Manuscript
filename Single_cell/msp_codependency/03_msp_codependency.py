import pandas as pd
import numpy as np
from itertools import combinations,product



act_df = pd.read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/figures_daf/scDAF_FIRE_actuation_MSP150_withProximal.tsv', sep="\t")

samples = set()
for col in act_df.columns:
    if col.startswith('PS'):
        samples.add(col.split('_')[0])


# Create dictionary of haplotype-strand MSP calls
hap_strand = [f'{samp}_H1' for samp in samples] + [f'{samp}_H2' for samp in samples]

hs_calls = dict()
for samp in samples:
    for hap in ['H1','H2']:
        hs = f'{samp}_{hap}'
        hs_calls[hs] = dict(zip(act_df['ID'], act_df[hs]))


# return scaled (-1 to 1) pairwise codependency score for provided counts: [neither, first, second, both]
def codep_score_fsExp(counts, exp):
    total = sum(counts)
    p1 = (counts[1]+counts[3]) / total
    p2 = (counts[2]+counts[3]) / total
    obs = counts[3] / total
    score = (obs-exp)*4
    return(score)

def codep_score(counts):
    total = sum(counts)
    p1 = (counts[1]+counts[3]) / total
    p2 = (counts[2]+counts[3]) / total
    exp = p1 * p2
    obs = counts[3] / total
    score = (obs-exp)*4
    return(score)


# By chromosome P vs Q arms (codep of p-p, q-q, p-q, cross centromeres) ------------------------------------------------------------------------
chrom_arms = pd.read_csv('FIRE-peaks_HG002_noUnreliableCoverage.chrom_arm_map.bed.gz', sep="\t", names = ['peak_chrom','peak_start','peak_end','strand','prop_act','spacer','arm_chrom','arm_start','arm_end','arm','overlap_bp'])
chrom_arms['ID'] = chrom_arms['peak_chrom'] + ':' + chrom_arms['peak_start'].astype(str) + '-' + chrom_arms['peak_end'].astype(str)
act_df['arm'] = act_df['ID'].map(chrom_arms.set_index('ID')['arm'])

# MSP codep for A/B compartments ------------------------------------------------------------------------
A_comp = pd.read_csv('FIRE-peaks_HG002_noUnreliableCoverage.A_compartments.bed.gz', sep="\t", names = ['peak_chrom','peak_start','peak_end','strand','prop_act','spacer','compartment_chrom','compartment_start','compartment_end','overlap_bp'])
B_comp = pd.read_csv('FIRE-peaks_HG002_noUnreliableCoverage.B_compartments.bed.gz', sep="\t", names = ['peak_chrom','peak_start','peak_end','strand','prop_act','spacer','compartment_chrom','compartment_start','compartment_end','overlap_bp'])

# track identity of the overlapping compartment
A_comp['compartment_key'] = A_comp['compartment_chrom'] + ':' + A_comp['compartment_start'].astype(str) + ':' + A_comp['compartment_end'].astype(str)
B_comp['compartment_key'] = B_comp['compartment_chrom'] + ':' + B_comp['compartment_start'].astype(str) + ':' + B_comp['compartment_end'].astype(str)
A_comp_idx = A_comp[['compartment_chrom','compartment_start','compartment_end','compartment_key']].drop_duplicates()
A_comp_idx['compartment_id'] = range(1,len(A_comp_idx)+1)
B_comp_idx = B_comp[['compartment_chrom','compartment_start','compartment_end','compartment_key']].drop_duplicates()
B_comp_idx['compartment_id'] = range(1,len(B_comp_idx)+1)
B_comp_idx['compartment_id'] = -B_comp_idx['compartment_id']
A_comp['compartment_id'] = A_comp['compartment_key'].map(A_comp_idx.set_index('compartment_key')['compartment_id'])
B_comp['compartment_id'] = B_comp['compartment_key'].map(B_comp_idx.set_index('compartment_key')['compartment_id'])

A_id_bed = A_comp[['compartment_chrom','compartment_start','compartment_end','compartment_id']]
B_id_bed = B_comp[['compartment_chrom','compartment_start','compartment_end','compartment_id']]
A_id_bed.to_csv('A_compartments_with_ID.bed', sep="\t", index = False, header = False)
B_id_bed.to_csv('B_compartments_with_ID.bed', sep="\t", index = False, header = False)

A_comp['motif_ID'] = A_comp['peak_chrom'] + ':' + A_comp['peak_start'].astype(str) + '-' + A_comp['peak_end'].astype(str)
B_comp['motif_ID'] = B_comp['peak_chrom'] + ':' + B_comp['peak_start'].astype(str) + '-' + B_comp['peak_end'].astype(str)
A_comp.index = A_comp['motif_ID']
B_comp.index = B_comp['motif_ID']
A_dict = A_comp[['compartment_id']].to_dict('index')
B_dict = B_comp[['compartment_id']].to_dict('index')


# MSP peak codependency chromosome-wide ------------------------------
min_act = 0.2

# generate list of comparisons (limit by chromosome)
tot_peaks = 0
combos = []
for ch in act_df['chrom'].unique():
    ch_peaks = act_df[act_df['chrom'] == ch]
    ch_peaks = list(ch_peaks[ch_peaks['prop_acc'] >=  min_act]['ID']) # filter by percent actuation
    tot_peaks += len(ch_peaks)
    n_combo = 0
    for combo in combinations(ch_peaks, 2):
        if combo[0] != combo[1]:
            combos.append(combo)
            n_combo += 1

act_df.index = act_df['ID']
fire_dict = act_df[['chrom','start','end','prop_acc','bin','is_proximal','arm']].to_dict('index')

# initialize codep DF
codep_df_rows = []
for c in combos:
    m1,m2 = c
    p1 = fire_dict[m1]['prop_acc']
    p2 = fire_dict[m2]['prop_acc']
    counts = [0,0,0,0]
    for k,v in hs_calls.items():
        c1 = v[m1]
        c2 = v[m2]
        if c1 >= 0 and c2 >= 0:
            if c1 == 1 and c2 == 1: # both
                counts[3] += 1
            elif c1 == 1: # first only
                counts[1] += 1
            elif c2 == 1: # second only
                counts[2] += 1
            else:  # neither
                counts[0] += 1
    new_row = [m1, fire_dict[m1]['chrom'], fire_dict[m1]['start'], fire_dict[m1]['end'], fire_dict[m1]['is_proximal'], p1, m2, fire_dict[m2]['chrom'], fire_dict[m2]['start'], fire_dict[m2]['end'], fire_dict[m2]['is_proximal'], p2] + counts
    if sum(counts) >= 8:
        score = codep_score(counts)
        fs_score = codep_score_fsExp(counts, p1*p2)
        new_row.append(score)
        new_row.append(fs_score)
        arms = fire_dict[m1]['arm'] + fire_dict[m2]['arm']
        new_row.append(arms)
        # A/B compartments
        if m1 in A_dict.keys():
            comp1 = 'A'
            comp_id_1 = A_dict[m1]['compartment_id']
        elif m1 in B_dict.keys():
            comp1 = 'B'
            comp_id_1 = B_dict[m1]['compartment_id']
        else:
            comp1 = 'N'
            comp_id_1 = 0
        if m2 in A_dict.keys():
            comp2 = 'A'
            comp_id_2 = A_dict[m2]['compartment_id']
        elif m2 in B_dict.keys():
            comp2 = 'B'
            comp_id_2 = B_dict[m2]['compartment_id']
        else:
            comp2 = 'N'
            comp_id_2 = 0
        new_row.append(comp1)
        new_row.append(comp_id_1)
        new_row.append(comp2)
        new_row.append(comp_id_2)
        codep_df_rows.append(new_row)

codep_df = pd.DataFrame(codep_df_rows, columns=['motif_1','chrom_1','start_1','end_1','proximal_1','prop_1','motif_2','chrom_2','start_2','end_2','proximal_2','prop_2','neither','first','second','both','score','fs_score','arms','m1_compartment','m1_compartment_id','m2_compartment','m2_compartment_id'])


# bin footprint combos by distance
codep_df['dist'] = abs((codep_df['start_1'] + round((codep_df['end_1']-codep_df['start_1'])/2)) - (codep_df['start_2'] + round((codep_df['end_2']-codep_df['start_2'])/2)))
bin_len = 100
codep_df['dist_bin'] = np.floor(codep_df['dist'] / bin_len).astype(int)
# Merge < 200 bp distance into a single bin
codep_df.loc[codep_df['dist_bin'] < 1, 'dist_bin'] = 1

codep_df.to_csv('codep_DAF_exp_df.csv', index=False)


# DAF exp, mean scores
with_scores = codep_df[codep_df['score'].isna() == False]
bin_means = with_scores.groupby('dist_bin')['score'].mean()
bin_medians = with_scores.groupby('dist_bin')['score'].median()

binned_codep = pd.DataFrame(bin_means)
binned_codep['dist_bin'] = binned_codep.index
binned_codep['dist'] = binned_codep['dist_bin'] * bin_len
binned_codep['count'] = codep_df['dist_bin'].value_counts().sort_index()
binned_codep.to_csv('bin_codep_DAF_exp_avg_score.csv', index=False)

binned_codep = pd.DataFrame(bin_medians)
binned_codep['dist_bin'] = binned_codep.index
binned_codep['dist'] = binned_codep['dist_bin'] * bin_len
binned_codep['count'] = codep_df['dist_bin'].value_counts().sort_index()
binned_codep.to_csv('bin_codep_DAF_exp_Median_score.csv', index=False)


# codep score by distance for OPPOSITE HAPLOTYPES ------------------------------------------------------------------------

# initialize codep DF
opp_codep_df_rows = []
for c in combos:
    m1,m2 = c
    p1 = fire_dict[m1]['prop_acc']
    p2 = fire_dict[m2]['prop_acc']
    counts = [0,0,0,0]
    for samp in samples:
        c1 = hs_calls[f'{samp}_H1'][m1]
        c2 = hs_calls[f'{samp}_H2'][m2]
        if c1 >= 0 and c2 >= 0:
            if c1 == 1 and c2 == 1: # both
                counts[3] += 1
            elif c1 == 1: # first only
                counts[1] += 1
            elif c2 == 1: # second only
                counts[2] += 1
            else: # neither
                counts[0] += 1
        # Check both haplotype combinations
        c1 = hs_calls[f'{samp}_H2'][m1]
        c2 = hs_calls[f'{samp}_H1'][m2]
        if c1 >= 0 and c2 >= 0:
            if c1 == 1 and c2 == 1: # both
                counts[3] += 1
            elif c1 == 1: # first only
                counts[1] += 1
            elif c2 == 1: # second only
                counts[2] += 1
            else: # neither
                counts[0] += 1
    new_row = [m1, fire_dict[m1]['chrom'], fire_dict[m1]['start'], fire_dict[m1]['end'], fire_dict[m1]['is_proximal'], p1, m2, fire_dict[m2]['chrom'], fire_dict[m2]['start'], fire_dict[m2]['end'], fire_dict[m2]['is_proximal'], p2] + counts
    if sum(counts) >= 8:
        score = codep_score(counts)
        fs_score = codep_score_fsExp(counts, p1*p2)
        new_row.append(score)
        new_row.append(fs_score)
        arms = fire_dict[m1]['arm'] + fire_dict[m2]['arm']
        new_row.append(arms)
        # A/B compartments
        if m1 in A_dict.keys():
            comp1 = 'A'
            comp_id_1 = A_dict[m1]['compartment_id']
        elif m1 in B_dict.keys():
            comp1 = 'B'
            comp_id_1 = B_dict[m1]['compartment_id']
        else:
            comp1 = 'N'
            comp_id_1 = 0
        if m2 in A_dict.keys():
            comp2 = 'A'
            comp_id_2 = A_dict[m2]['compartment_id']
        elif m2 in B_dict.keys():
            comp2 = 'B'
            comp_id_2 = B_dict[m2]['compartment_id']
        else:
            comp2 = 'N'
            comp_id_2 = 0
        new_row.append(comp1)
        new_row.append(comp_id_1)
        new_row.append(comp2)
        new_row.append(comp_id_2)
        opp_codep_df_rows.append(new_row)

codep_df_oppHap = pd.DataFrame(opp_codep_df_rows, columns=['motif_1','chrom_1','start_1','end_1','proximal_1','prop_1','motif_2','chrom_2','start_2','end_2','proximal_2','prop_2','neither','first','second','both','score','fs_score','arms','m1_compartment','m1_compartment_id','m2_compartment','m2_compartment_id'])


# bin footprint combos by distance
codep_df_oppHap['dist'] = abs((codep_df_oppHap['start_1'] + round((codep_df_oppHap['end_1']-codep_df_oppHap['start_1'])/2)) - (codep_df_oppHap['start_2'] + round((codep_df_oppHap['end_2']-codep_df_oppHap['start_2'])/2)))
bin_len = 100
codep_df_oppHap['dist_bin'] = np.floor(codep_df_oppHap['dist'] / bin_len).astype(int)
# Merge < 200 bp distance into a single bin
codep_df_oppHap.loc[codep_df_oppHap['dist_bin'] < 1, 'dist_bin'] = 1

codep_df_oppHap.to_csv('codep_DAF_exp_OPPOSITE_HAP_df.csv', index=False)


# DAF exp, mean scores
with_scores = codep_df_oppHap[codep_df_oppHap['score'].isna() == False]
bin_means = with_scores.groupby('dist_bin')['score'].mean()
bin_medians = with_scores.groupby('dist_bin')['score'].median()

binned_codep_oppHap = pd.DataFrame(bin_means)
binned_codep_oppHap['dist_bin'] = binned_codep_oppHap.index
binned_codep_oppHap['dist'] = binned_codep_oppHap['dist_bin'] * bin_len
binned_codep_oppHap['count'] = codep_df_oppHap['dist_bin'].value_counts().sort_index()
binned_codep_oppHap.to_csv('bin_codep_DAF_exp_avg_score_OPPOSITE_HAP.csv', index=False)

binned_codep_oppHap = pd.DataFrame(bin_medians)
binned_codep_oppHap['dist_bin'] = binned_codep_oppHap.index
binned_codep_oppHap['dist'] = binned_codep_oppHap['dist_bin'] * bin_len
binned_codep_oppHap['count'] = codep_df_oppHap['dist_bin'].value_counts().sort_index()
binned_codep_oppHap.to_csv('bin_codep_DAF_exp_Median_score_OPPOSITE_HAP.csv', index=False)



# codep score by distance for DIFFERENT CELLS (Negative Control) ------------------------------------------------------------------------

# initialize codep DF
diff_codep_df_rows = []
for c in combos:
    m1,m2 = c
    p1 = fire_dict[m1]['prop_acc']
    p2 = fire_dict[m2]['prop_acc']
    counts = [0,0,0,0]
    for samp in samples:
        for samp2 in samples:
            if samp2 != samp:
                c1 = hs_calls[f'{samp}_H1'][m1]
                c2 = hs_calls[f'{samp2}_H2'][m2]
                if c1 >= 0 and c2 >= 0:
                    if c1 == 1 and c2 == 1: # both
                        counts[3] += 1
                    elif c1 == 1: # first only
                        counts[1] += 1
                    elif c2 == 1: # second only
                        counts[2] += 1
                    else: # neither
                        counts[0] += 1
                # Check both haplotype combinations
                c1 = hs_calls[f'{samp}_H2'][m1]
                c2 = hs_calls[f'{samp2}_H1'][m2]
                if c1 >= 0 and c2 >= 0:
                    if c1 == 1 and c2 == 1: # both
                        counts[3] += 1
                    elif c1 == 1: # first only
                        counts[1] += 1
                    elif c2 == 1: # second only
                        counts[2] += 1
                    else: # neither
                        counts[0] += 1
    new_row = [m1, fire_dict[m1]['chrom'], fire_dict[m1]['start'], fire_dict[m1]['end'], fire_dict[m1]['is_proximal'], p1, m2, fire_dict[m2]['chrom'], fire_dict[m2]['start'], fire_dict[m2]['end'], fire_dict[m2]['is_proximal'], p2] + counts
    if sum(counts) >= 8:
        score = codep_score(counts)
        fs_score = codep_score_fsExp(counts, p1*p2)
        new_row.append(score)
        new_row.append(fs_score)
        arms = fire_dict[m1]['arm'] + fire_dict[m2]['arm']
        new_row.append(arms)
        # A/B compartments
        if m1 in A_dict.keys():
            comp1 = 'A'
            comp_id_1 = A_dict[m1]['compartment_id']
        elif m1 in B_dict.keys():
            comp1 = 'B'
            comp_id_1 = B_dict[m1]['compartment_id']
        else:
            comp1 = 'N'
            comp_id_1 = 0
        if m2 in A_dict.keys():
            comp2 = 'A'
            comp_id_2 = A_dict[m2]['compartment_id']
        elif m2 in B_dict.keys():
            comp2 = 'B'
            comp_id_2 = B_dict[m2]['compartment_id']
        else:
            comp2 = 'N'
            comp_id_2 = 0
        new_row.append(comp1)
        new_row.append(comp_id_1)
        new_row.append(comp2)
        new_row.append(comp_id_2)
        diff_codep_df_rows.append(new_row)

codep_df_DiffCells = pd.DataFrame(diff_codep_df_rows, columns=['motif_1','chrom_1','start_1','end_1','proximal_1','prop_1','motif_2','chrom_2','start_2','end_2','proximal_2','prop_2','neither','first','second','both','score','fs_score','arms','m1_compartment','m1_compartment_id','m2_compartment','m2_compartment_id'])


# bin footprint combos by distance
codep_df_DiffCells['dist'] = abs((codep_df_DiffCells['start_1'] + round((codep_df_DiffCells['end_1']-codep_df_DiffCells['start_1'])/2)) - (codep_df_DiffCells['start_2'] + round((codep_df_DiffCells['end_2']-codep_df_DiffCells['start_2'])/2)))
bin_len = 100
codep_df_DiffCells['dist_bin'] = np.floor(codep_df_DiffCells['dist'] / bin_len).astype(int)
# Merge < 200 bp distance into a single bin
codep_df_DiffCells.loc[codep_df_DiffCells['dist_bin'] < 1, 'dist_bin'] = 1

codep_df_DiffCells.to_csv('codep_DAF_exp_DIFF_CELLS_df.csv', index=False)

