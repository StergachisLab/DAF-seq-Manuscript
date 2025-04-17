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
files = os.listdir(ft_data_dir)
files.remove('ft_H1_ModCTCF.bed')
files.remove('ft_H2_ModCTCF.bed')
motifs = {name.split('_ft_')[1].rstrip('.bed') for name in files}

# Quantify occupancy separately for each haplotype ------------------------------------
motif_zmws = dict() # H1 & H2 fibers [MSP, bound]
motif_dict = dict() # create dict of footprint stats and G/C content within motifs
for mo in motifs:
    h1_file = glob(f'{ft_data_dir}/H1_GM12878_UBA1_PS00685_ft_{mo}*.bed')[0]
    h2_file = glob(f'{ft_data_dir}/H2_GM12878_UBA1_PS00685_ft_{mo}*.bed')[0]
    h1_df = pd.read_csv(h1_file, delimiter="\t")
    h2_df = pd.read_csv(h2_file, delimiter="\t")
    for i in range(len(h1_df)):
        if h1_df.iloc[i,:]['n_spanning_fibers'] > 0 and h2_df.iloc[i,:]['n_spanning_fibers'] > 0:
            # Hap 1
            row = h1_df.iloc[i,:]
            motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
            motif_zmws[motif_key] = {'H1':[set(),set()], 'H2':[set(),set()]}
            fp_codes = row['footprint_codes'].split(',')
            fibers = row['fiber_names'].split(',')
            for j in range(len(fp_codes)):
                if (int(fp_codes[j]) & 1) > 0:
                    motif_zmws[motif_key]['H1'][0].add(fibers[j])
                if (int(fp_codes[j]) & (1 << 1)):
                    motif_zmws[motif_key]['H1'][1].add(fibers[j])
            motif_dict[motif_key] = {'name':mo, 'chr':row['#chrom'], 'start':row['start'], 'end':row['end'], 'H1_nuc':row['n_overlapping_nucs'], 'H1_msp':row['n_spanning_msps'], 'H1_bound':h1_df.iloc[i,7]}
            tot_msp = motif_dict[motif_key]['H1_bound'] + motif_dict[motif_key]['H1_msp']
            if tot_msp > 0:
                motif_dict[motif_key]['H1_prop'] = (motif_dict[motif_key]['H1_bound']/tot_msp)
            else:
                motif_dict[motif_key]['H1_prop'] = 0
            # Hap 2
            row = h2_df.iloc[i,:]
            motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
            fp_codes = row['footprint_codes'].split(',')
            fibers = row['fiber_names'].split(',')
            for j in range(len(fp_codes)):
                if (int(fp_codes[j]) & 1) > 0:
                    motif_zmws[motif_key]['H2'][0].add(fibers[j])
                if (int(fp_codes[j]) & (1 << 1)):
                    motif_zmws[motif_key]['H2'][1].add(fibers[j])
            motif_dict[motif_key]['H2_nuc'] = row['n_overlapping_nucs']
            motif_dict[motif_key]['H2_msp'] = row['n_spanning_msps']
            motif_dict[motif_key]['H2_bound'] = h2_df.iloc[i,7]
            tot_msp = motif_dict[motif_key]['H2_bound'] + motif_dict[motif_key]['H2_msp']
            if tot_msp > 0:
                motif_dict[motif_key]['H2_prop'] = (motif_dict[motif_key]['H2_bound']/tot_msp)
            else:
                motif_dict[motif_key]['H2_prop'] = 0


# Update CTCF with mod 2/3 values ------------------------
del motif_dict['chrX:47193672-47193690']

h1_file = 'ft_out/ft_H1_ModCTCF.bed'
h1_df = pd.read_csv(h1_file, delimiter="\t")
# H1 reads
row = h1_df.iloc[0,:]
h1_codes = h1_df.iloc[0,:]['footprint_codes'].split(',')
h1_counts = Counter(h1_codes)
n_bound = 0
for k,v in h1_counts.items():
    fp_code = int(k)
    if (fp_code & (1 << 3)) > 0 and (fp_code & (1 << 4)) > 0 and (fp_code & 1) > 0: # mods 2&3 bound and spanning MSP
        n_bound += v
motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
motif_zmws[motif_key] = {'H1':[set(),set()], 'H2':[set(),set()]}
fp_codes = row['footprint_codes'].split(',')
fibers = row['fiber_names'].split(',')
for j in range(len(fp_codes)):
    if (int(fp_codes[j]) & 1) > 0:
        motif_zmws[motif_key]['H1'][0].add(fibers[j])
    if (int(fp_codes[j]) & (1 << 3)) and (int(fp_codes[j]) & (1 << 4)): # footprint CTCF Modules 2 AND 3
        motif_zmws[motif_key]['H1'][1].add(fibers[j])
motif_dict[motif_key] = {'name':'CTCF_Mod2-3', 'chr':row['#chrom'], 'start':row['start'], 'end':row['end'], 'H1_nuc':row['n_overlapping_nucs'], 'H1_msp':row['n_spanning_msps'], 'H1_bound':n_bound, 'H1_prop':(n_bound/row['n_spanning_msps'])}
# H2 reads
h2_file = 'ft_out/ft_H2_ModCTCF.bed'
h2_df = pd.read_csv(h2_file, delimiter="\t")
row = h2_df.iloc[0,:]
h2_codes = h2_df.iloc[0,:]['footprint_codes'].split(',')
h2_counts = Counter(h2_codes)
n_bound = 0
for k,v in h2_counts.items():
    fp_code = int(k)
    if (fp_code & (1 << 3)) > 0 and (fp_code & (1 << 4)) > 0 and (fp_code & 1) > 0: # mods 2&3 bound and spanning MSP
        n_bound += v
motif_key = f"{row['#chrom']}:{str(row['start'])}-{str(row['end'])}"
fp_codes = row['footprint_codes'].split(',')
fibers = row['fiber_names'].split(',')
for j in range(len(fp_codes)):
    if (int(fp_codes[j]) & 1) > 0:
        motif_zmws[motif_key]['H2'][0].add(fibers[j])
    if (int(fp_codes[j]) & (1 << 3)) and (int(fp_codes[j]) & (1 << 4)): # footprint CTCF Modules 2 AND 3
        motif_zmws[motif_key]['H2'][1].add(fibers[j])
motif_dict[motif_key]['H2_nuc'] = row['n_overlapping_nucs']
motif_dict[motif_key]['H2_msp'] = row['n_spanning_msps']
motif_dict[motif_key]['H2_bound'] = n_bound
motif_dict[motif_key]['H2_prop'] = (n_bound/row['n_spanning_msps'])
strand = '+'
ref_mod2 = 'GTCCACC'
ref_mod3 = 'AGGTGG' 
seq = ref_mod2+ref_mod3
motif_dict[motif_key]['strand'] = strand
motif_dict[motif_key]['seq'] = seq
motif_dict[motif_key]['nC'] = seq.count('C')
motif_dict[motif_key]['nG'] = seq.count('G')

# Add sequence content for motifs
fimo_file = 'motifs/filtered_fimo.tsv'
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
# filter motifs by number of Cs and Gs !!!!!!!!!!!!!!!!!!!!!!!!!!!
prop_cutoff = 0.05
n_base_cutoff = 2
ft_H1 = []
for k,v in motif_dict.items():
    if v['H1_prop'] >= prop_cutoff: # filter by % bound
        if v['nG'] >= n_base_cutoff:
            ft_H1.append(k)
# write BED file of footprinted motifs
out_bed_H1 = f'UBA1_GM12878_footprinted_motifs_H1.bed' # FOOTPRINTED ON EITHER HAPLOTYPE
with open(out_bed_H1, 'w') as out_bed:
    writer = csv.writer(out_bed, delimiter="\t")
    for motif in ft_H1:
        line = [motif_dict[motif]['chr'], motif_dict[motif]['start'], motif_dict[motif]['end'], motif_dict[motif]['name'], motif_dict[motif]['H1_prop'], motif_dict[motif]['strand']]
        writer.writerow(line)

