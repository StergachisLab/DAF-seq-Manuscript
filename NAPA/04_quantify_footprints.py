import csv
import os
from glob import glob
import pandas as pd
from collections import Counter

def rc(s):
    bases = {'A':'T','C':'G','G':'C','T':'A','Y':'R','R':'Y'}
    new = ''
    for i in range(len(s)):
        new = bases[s[i]] + new
    return(new)

ft_data_dir = 'ft_out'
motifs = {name.lstrip('ft_CT_').lstrip('ft_GA_').rstrip('.bed') for name in os.listdir(ft_data_dir)}

motif_zmws = dict() # CT & GT fibers [MSP, bound]
motif_dict = dict() # create dict of footprint stats and G/C content within motifs
for mo in motifs:
    ct_file = glob(f'{ft_data_dir}/*CT_{mo}*.bed')[0]
    ga_file = glob(f'{ft_data_dir}/*GA_{mo}*.bed')[0]
    ct_df = pd.read_csv(ct_file, delimiter="\t")
    ga_df = pd.read_csv(ga_file, delimiter="\t")
    for i in range(len(ct_df)):
        if ct_df.iloc[i,:]['n_spanning_fibers'] > 0 and ga_df.iloc[i,:]['n_spanning_fibers'] > 0:
            # CT reads
            row = ct_df.iloc[i,:]
            motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
            motif_zmws[motif_key] = {'CT':[set(),set()], 'GA':[set(),set()]}
            fp_codes = row['footprint_codes'].split(',')
            fibers = row['fiber_names'].split(',')
            for j in range(len(fp_codes)):
                if (int(fp_codes[j]) & 1) > 0:
                    motif_zmws[motif_key]['CT'][0].add(fibers[j])
                if (int(fp_codes[j]) & (1 << 1)):
                    motif_zmws[motif_key]['CT'][1].add(fibers[j])
            motif_dict[motif_key] = {'name':mo, 'chr':row['#chrom'], 'start':row['start'], 'end':row['end'], 'ct_nuc':row['n_overlapping_nucs'], 'ct_msp':row['n_spanning_msps'], 'ct_bound':ct_df.iloc[i,7]}
            tot_msp = motif_dict[motif_key]['ct_bound'] + motif_dict[motif_key]['ct_msp']
            if tot_msp > 0:
                motif_dict[motif_key]['ct_prop'] = (motif_dict[motif_key]['ct_bound']/tot_msp)
            else:
                motif_dict[motif_key]['ct_prop'] = 0
            # GA reads
            row = ga_df.iloc[i,:]
            motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
            fp_codes = row['footprint_codes'].split(',')
            fibers = row['fiber_names'].split(',')
            for j in range(len(fp_codes)):
                if (int(fp_codes[j]) & 1) > 0:
                    motif_zmws[motif_key]['GA'][0].add(fibers[j])
                if (int(fp_codes[j]) & (1 << 1)):
                    motif_zmws[motif_key]['GA'][1].add(fibers[j])
            motif_dict[motif_key]['ga_nuc'] = row['n_overlapping_nucs']
            motif_dict[motif_key]['ga_msp'] = row['n_spanning_msps']
            motif_dict[motif_key]['ga_bound'] = ga_df.iloc[i,7]
            tot_msp = motif_dict[motif_key]['ga_bound'] + motif_dict[motif_key]['ga_msp']
            if tot_msp > 0:
                motif_dict[motif_key]['ga_prop'] = (motif_dict[motif_key]['ga_bound']/tot_msp)
            else:
                motif_dict[motif_key]['ga_prop'] = 0


# Update CTCF with mod 2/3 values ------------------------
del motif_dict['chr19:47515578-47515612']
del motif_dict['chr19:47515578-47515611']
ct_file = glob(f'{ft_data_dir}/*CT_ModCTCF*.bed')[0]
ga_file = glob(f'{ft_data_dir}/*GA_ModCTCF*.bed')[0]
ct_df = pd.read_csv(ct_file, delimiter="\t")
ga_df = pd.read_csv(ga_file, delimiter="\t")
# use Module 2 only !!!!!
# CT reads
row = ct_df.iloc[0,:]
ct_codes = ct_df.iloc[0,:]['footprint_codes'].split(',')
ct_counts = Counter(ct_codes)
n_bound = 0
for k,v in ct_counts.items():
    fp_code = int(k)
    if (fp_code & (1 << 3)) > 0 and (fp_code & 1) > 0: # mod 2 bound and spanning MSP
        n_bound += v
motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
motif_zmws[motif_key] = {'CT':[set(),set()], 'GA':[set(),set()]}
fp_codes = row['footprint_codes'].split(',')
fibers = row['fiber_names'].split(',')
for j in range(len(fp_codes)):
    if (int(fp_codes[j]) & 1) > 0:
        motif_zmws[motif_key]['CT'][0].add(fibers[j])
    if (int(fp_codes[j]) & (1 << 3)):
        motif_zmws[motif_key]['CT'][1].add(fibers[j])
motif_dict[motif_key] = {'name':'CTCF_Mod2', 'chr':row['#chrom'], 'start':row['start'], 'end':row['end'], 'ct_nuc':row['n_overlapping_nucs'], 'ct_msp':row['n_spanning_msps'], 'ct_bound':n_bound, 'ct_prop':(n_bound/row['n_spanning_msps'])}
# GA reads
row = ga_df.iloc[0,:]
ga_codes = ga_df.iloc[0,:]['footprint_codes'].split(',')
ga_counts = Counter(ga_codes)
n_bound = 0
for k,v in ga_counts.items():
    fp_code = int(k)
    if (fp_code & (1 << 3)) > 0 and (fp_code & 1) > 0: # mod 2 bound and spanning MSP
        n_bound += v
motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
fp_codes = row['footprint_codes'].split(',')
fibers = row['fiber_names'].split(',')
for j in range(len(fp_codes)):
    if (int(fp_codes[j]) & 1) > 0:
        motif_zmws[motif_key]['GA'][0].add(fibers[j])
    if (int(fp_codes[j]) & (1 << 3)):
        motif_zmws[motif_key]['GA'][1].add(fibers[j])
motif_dict[motif_key]['ga_nuc'] = row['n_overlapping_nucs']
motif_dict[motif_key]['ga_msp'] = row['n_spanning_msps']
motif_dict[motif_key]['ga_bound'] = n_bound
motif_dict[motif_key]['ga_prop'] = (n_bound/row['n_spanning_msps'])
strand = '-'
ref_mod2 = 'GTCGGCC'
seq = rc(ref_mod2)
motif_dict[motif_key]['strand'] = strand
motif_dict[motif_key]['seq'] = seq
motif_dict[motif_key]['nC'] = seq.count('C')
motif_dict[motif_key]['nG'] = seq.count('G')

fimo_file = 'motifs/fimo_motifs_filtered.tsv'
with open(fimo_file) as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        key = f"{line[2]}:{line[3]}-{line[4]}"
        if key in motif_dict.keys():
            strand = line[5]
            seq = line[9].upper()
            motif_dict[key]['strand'] = strand
            if strand == '-': # fimo reports the sequence of the motif, convert to ref sequence as with DddA
                seq = rc(seq)
            motif_dict[key]['seq'] = seq
            motif_dict[key]['nC'] = seq.count('C')
            motif_dict[key]['nG'] = seq.count('G')
            
variants = pd.read_csv('../../haplotype_correct_NAPA/NAPA_PS00626_haplotype_corrected_hap_basecalls.tsv', delimiter="\t")


# filter to NAPA promoter region chr19:47,514,957-47,515,660 FIRE narrow peak #99860
# make start position the TSS 47,515,063

prom_start = 47515063
prom_end = 47515660

# filter motifs by number of Cs and Gs !!!!!!!!!!!!!!!!!!!!!!!!!!!
prop_cutoff = 0.05
n_base_cutoff = 2

ft_both = []
ft_ct = []
ft_ga = []
for k,v in motif_dict.items():
    if v['start'] > prom_start and v['end'] < prom_end:
        if v['ct_prop'] >= prop_cutoff and v['ga_prop'] >= prop_cutoff: # filter by % bound
            if v['nC'] >= n_base_cutoff and v['nG'] >= n_base_cutoff:
                ft_both.append(k)
        elif v['nC'] >= n_base_cutoff and v['ct_prop'] >= prop_cutoff:
            ft_ct.append(k)
        elif v['nG'] >= n_base_cutoff and v['ga_prop'] >= prop_cutoff:
            ft_ga.append(k)


# motifs in regulatory element
re2_both = []
start = 47518185
end = 47518479
for k,v in motif_dict.items():
    if v['start'] > start and v['end'] < end:
        # if v['ct_prop'] >= prop_cutoff and v['ga_prop'] >= prop_cutoff: # filter by % bound
        if v['nC'] >= n_base_cutoff and v['nG'] >= n_base_cutoff:
            re2_both.append(k)

# motif_dict['chr19:47518341-47518361'] CT only 12.2% MA1987.1.ZNF701


# write BED file of footprinted motifs

# MUST BE FOOTPRINTED ON BOTH STRANDS
out_bed_both = 'footprinted_motifs_both.bed'
with open(out_bed_both, 'w') as out_bed:
    writer = csv.writer(out_bed, delimiter="\t")
    for motif in ft_both:
        line = [motif_dict[motif]['chr'], motif_dict[motif]['start'], motif_dict[motif]['end'], motif_dict[motif]['name'], '.', motif_dict[motif]['strand']]
        writer.writerow(line)

out_bed_ct = 'footprinted_motifs_CT.bed'
with open(out_bed_ct, 'w') as out_bed:
    writer = csv.writer(out_bed, delimiter="\t")
    for motif in ft_both:
        line = [motif_dict[motif]['chr'], motif_dict[motif]['start'], motif_dict[motif]['end'], motif_dict[motif]['name'], motif_dict[motif]['ct_prop'], motif_dict[motif]['strand']]
        writer.writerow(line)
    for motif in ft_ct:
        line = [motif_dict[motif]['chr'], motif_dict[motif]['start'], motif_dict[motif]['end'], motif_dict[motif]['name'], motif_dict[motif]['ct_prop'], motif_dict[motif]['strand']]
        writer.writerow(line)

out_bed_ga = 'footprinted_motifs_GA.bed'
with open(out_bed_ga, 'w') as out_bed:
    writer = csv.writer(out_bed, delimiter="\t")
    for motif in ft_both:
        line = [motif_dict[motif]['chr'], motif_dict[motif]['start'], motif_dict[motif]['end'], motif_dict[motif]['name'], motif_dict[motif]['ga_prop'], motif_dict[motif]['strand']]
        writer.writerow(line)
    for motif in ft_ga:
        line = [motif_dict[motif]['chr'], motif_dict[motif]['start'], motif_dict[motif]['end'], motif_dict[motif]['name'], motif_dict[motif]['ga_prop'], motif_dict[motif]['strand']]
        writer.writerow(line)


# sort-bed footprinted_motifs_both.bed > footprinted_motifs_both_sorted.bed
# sort-bed footprinted_motifs_CT.bed > footprinted_motifs_CT_sorted.bed
# sort-bed footprinted_motifs_GA.bed > footprinted_motifs_GA_sorted.bed

# bedmap --max-element footprinted_motifs_CT_sorted.bed | sort | uniq | sort-bed - > max_motifs_CT.bed
# bedmap --max-element footprinted_motifs_GA_sorted.bed | sort | uniq | sort-bed - > max_motifs_GA.bed

# Merge motifs to find footprinted regions, not motifs ------------------------------------------------------------------------
# merge motifs which overlap either motif by 80%, then merge regions by overlap of either by 90% to merge regions contained within another region

# MUST BE FOOTPRINTED ON BOTH STRANDS
# bedops --everything footprinted_motifs_both_sorted.bed footprinted_motifs_both_sorted.bed | bedmap --ec --fraction-either 0.8 --echo-map-range - | sort-bed - | bedmap --ec --fraction-either 0.9 --echo-map-range - | sort-bed - | uniq > merged_ft_on_both_strands.bed


# bedops --everything footprinted_motifs_GA_sorted.bed footprinted_motifs_CT_sorted.bed | bedmap --ec --fraction-either 0.8 --echo-map-range - | sort-bed - | bedmap --ec --fraction-either 0.9 --echo-map-range - | sort-bed - | uniq > merged_both_strands.bed
# bedops --everything footprinted_motifs_CT_sorted.bed footprinted_motifs_CT_sorted.bed | bedmap --ec --fraction-either 0.8 --echo-map-range - | sort-bed - | bedmap --ec --fraction-either 0.9 --echo-map-range - | sort-bed - | uniq > merged_CT.bed
# bedops --everything footprinted_motifs_GA_sorted.bed footprinted_motifs_GA_sorted.bed | bedmap --ec --fraction-either 0.8 --echo-map-range - | sort-bed - | bedmap --ec --fraction-either 0.9 --echo-map-range - | sort-bed - | uniq > merged_GA.bed


# Export footprinting as DataFrame
c = 1
regions = dict()
with open('merged_ft_on_both_strands.bed') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        regions[c] = line
        c += 1
df = pd.DataFrame(motif_dict).T
# ONLY TAKE motifs with footprints on BOTH strands ----------------
df = df.filter(items=ft_both, axis=0)
# df = df.filter(items=ft_both+ft_ct+ft_ga, axis=0)
df['strand'] = None
df.loc[df.index.isin(ft_both), 'strand'] = 'both'
# df.loc[df.index.isin(ft_ct), 'strand'] = 'ct'
# df.loc[df.index.isin(ft_ga), 'strand'] = 'ga'

df_reg_numbers = []
df_regions = []
for i in range(len(df)):
    r = None
    st = df.iloc[i,:]['start']
    end = df.iloc[i,:]['end']
    for n,reg in regions.items():
        if st >= int(reg[1]) and end <= int(reg[2]):
            df_reg_numbers.append(n)
            df_regions.append(f'{reg[0]}:{reg[1]}-{reg[2]}')
df['reg_num'] = df_reg_numbers
df['region'] = df_regions

df.to_csv('footprint_stats.csv', index_label='motif')


# DF of zmws by region footprint status (not overlapped by MSP -> -1, MSP & unbound -> 0, MSP & bound -> 1)
# in cases where a motif can only be footprinted on CT or GA strand: only count from the relevant strand, code -1 otherwise (GA & CT, ZNF692 & ZNF263)
reg_zmw_dict = {r:{'msp':set(),'bound':set()} for r in set(df['reg_num'])}
for i in range(len(df)):
    key = df.index[i]
    r = df.iloc[i,:]['reg_num']
    if key in ft_both or key in ft_ct:
        reg_zmw_dict[r]['msp'] = reg_zmw_dict[r]['msp'].union(motif_zmws[key]['CT'][0])
        reg_zmw_dict[r]['bound'] = reg_zmw_dict[r]['bound'].union(motif_zmws[key]['CT'][1])
    if key in ft_both or key in ft_ga:
        reg_zmw_dict[r]['msp'] = reg_zmw_dict[r]['msp'].union(motif_zmws[key]['GA'][0])
        reg_zmw_dict[r]['bound'] = reg_zmw_dict[r]['bound'].union(motif_zmws[key]['GA'][1])

all_msp_zmws = set()
for k,v in reg_zmw_dict.items():
    all_msp_zmws = all_msp_zmws.union(reg_zmw_dict[k]['msp'])
out_rows = []
for z in all_msp_zmws:
    z_row = [z]
    for k in reg_zmw_dict.keys():
        if z in reg_zmw_dict[k]['msp']:
            if z in reg_zmw_dict[k]['bound']:
                z_row.append(1)
            else:
                z_row.append(0)
        else:
            z_row.append(-1)
    out_rows.append(z_row)
z_df = pd.DataFrame(out_rows)
z_df = z_df.set_index([0])

z_df.to_csv('zmw_footprint_regions.csv', index_label='zmw')


# Footprinting difference by strand (# of C or G in the motif region) -----------------------------
# for now, not correcting for variants identified previously

ft_diff_df = pd.DataFrame(motif_dict).T # 269 motifs
pass_list = []
for i in range(len(ft_diff_df)):
    if ft_diff_df.index[i] in df.index:
        pass_list.append(1)
    else:
        pass_list.append(0)
ft_diff_df['passing'] = pass_list


ft_diff_df.to_csv('diff_by_strand_df.csv', index_label='motif')

# ft_diff_df = ft_diff_df[(ft_diff_df['nC'] > 0) & (ft_diff_df['nG'] > 0)] # 224 motifs

# ft_diff_df = ft_diff_df[(ft_diff_df['nC'] >= 3) | (ft_diff_df['nG'] >= 3)] # 214 motifs

# ft_diff_df[['nC','nG']]



