import pysam
import csv
import pandas as pd

# NAPA titration haplotype corrected BAMs
napa_dir = '../'
napa_bams = ['NAPA_PS00626_haplotype_corrected.bam',
             'NAPA_PS00627_haplotype_corrected.bam',
             'NAPA_PS00628_haplotype_corrected.bam',
             'NAPA_PS00629_haplotype_corrected.bam',
             'NAPA_PS00630_haplotype_corrected.bam',
             'NAPA_PS00631_haplotype_corrected.bam']

hap_basecalls = {47514700: ['T', 'C'],
 47515000: ['T', 'G'],
 47515603: ['G', 'G'],
 47516117: ['A', 'T'],
 47517160: ['C', 'T'],
 47517464: ['T', 'C'],
 47518032: ['T', 'T'],
 47518383: ['C', 'G'],
 47518977: ['C', 'C']}


da_dicts = []
for i in range(len(napa_bams)):
    da_dicts.append(dict())
    bam = pysam.AlignmentFile(f'{napa_dir}/{napa_bams[i]}', "rb")
    for read in bam.fetch():
        seq = read.query_sequence
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        strand = read.get_tag('ST')
        try:
            hap = read.get_tag('HP')
        except:
            hap = 'UNK'
        if hap in [1,2]:
            for pos in pair:
                if pos[0] == None or pos[1] == None: # indel, ignore
                    pass
                else:
                    qi = pos[0]
                    ref_coord = pos[1] + 1
                    base = seq[qi]
                    if strand == 'CT':
                        if ref_coord in hap_basecalls.keys():
                            if hap_basecalls[ref_coord][hap-1] in ['C']:
                                if ref_coord not in da_dicts[i].keys():
                                    da_dicts[i][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                                da_dicts[i][ref_coord][base] += 1
                        elif pos[2].upper() in ['C']:
                            if ref_coord not in da_dicts[i].keys():
                                da_dicts[i][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                            da_dicts[i][ref_coord][base] += 1
                    elif strand == 'GA':
                        if ref_coord in hap_basecalls.keys():
                            if hap_basecalls[ref_coord][hap-1] in ['G']:
                                if ref_coord not in da_dicts[i].keys():
                                    da_dicts[i][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                                da_dicts[i][ref_coord][base] += 1
                        elif pos[2].upper() in ['G']:
                            if ref_coord not in da_dicts[i].keys():
                                da_dicts[i][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                            da_dicts[i][ref_coord][base] += 1
    bam.close()


shared_pos = set.intersection(*[set(d.keys()) for d in da_dicts])

out_pos = []
for p in shared_pos:
    props = []
    for i in range(len(da_dicts)):
        da_dicts[i][p]
        if da_dicts[i][p]['C'] > da_dicts[i][p]['G']:
            prop = da_dicts[i][p]['Y']/(da_dicts[i][p]['Y']+da_dicts[i][p]['C'])
        elif da_dicts[i][p]['C'] < da_dicts[i][p]['G']:
            prop = da_dicts[i][p]['R']/(da_dicts[i][p]['R']+da_dicts[i][p]['G'])
        props.append(prop)
    out_pos.append([p] + props)

with open('prop_da_titration_NAPA.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(['position','PS00626','PS00627','PS00628','PS00629','PS00630','PS00631'])
    for row in out_pos:
        writer.writerow(row)

# average % deaminations by condition
tot_bases = [0 for i in range(len(da_dicts))]
deam_bases = [0 for i in range(len(da_dicts))]
for p in shared_pos:
    for i in range(len(da_dicts)):
        if da_dicts[i][p]['C'] > da_dicts[i][p]['G']:
            tot_bases[i] += da_dicts[i][p]['Y']+da_dicts[i][p]['C']
            deam_bases[i] += da_dicts[i][p]['Y']
        elif da_dicts[i][p]['C'] < da_dicts[i][p]['G']:
            tot_bases[i] += da_dicts[i][p]['R']+da_dicts[i][p]['G']
            deam_bases[i] += da_dicts[i][p]['R']
overall_props = []
for i in range(len(da_dicts)):
    overall_props.append(deam_bases[i] / tot_bases[i])

with open('overall_prop_da_titration_NAPA.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(['PS00626','PS00627','PS00628','PS00629','PS00630','PS00631'])
    writer.writerow(overall_props)


# NAPA Promoter % deaminations by condition

# filter to NAPA promoter region (FIRE narrow peak #99860 chr19:47,514,957-47,515,660)
# Made start position the NAPA TSS 47,515,063!
# chr19 47515063    47515660
# Consistent with NAPA footprinting script: /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/footprinting/NAPA/quantify_footprints.py

prom_start = 47515063
prom_end = 47515660

prom_pos = set()
for k,v in da_dicts[0].items():
    if k >= prom_start and k <= prom_end:
        prom_pos.add(k)

tot_bases = [0 for i in range(len(da_dicts))]
deam_bases = [0 for i in range(len(da_dicts))]
for p in prom_pos:
    for i in range(len(da_dicts)):
        if da_dicts[i][p]['C'] > da_dicts[i][p]['G']:
            tot_bases[i] += da_dicts[i][p]['Y']+da_dicts[i][p]['C']
            deam_bases[i] += da_dicts[i][p]['Y']
        elif da_dicts[i][p]['C'] < da_dicts[i][p]['G']:
            tot_bases[i] += da_dicts[i][p]['R']+da_dicts[i][p]['G']
            deam_bases[i] += da_dicts[i][p]['R']
overall_props = []
for i in range(len(da_dicts)):
    overall_props.append(deam_bases[i] / tot_bases[i])

with open('overall_prop_da_titration_NAPA_promoter.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(['PS00626','PS00627','PS00628','PS00629','PS00630','PS00631'])
    writer.writerow(overall_props)




# --------------------------------------------------------------------------------------------------------------------

# WASF1 titration BAMs
wasf1_dir = '../'
wasf1_bams = ['PS00626.m84046_240619_124816_s1.bc2072.ft.map-pb_corrected_realigned.bam',
'PS00627.m84046_240619_124816_s1.bc2073.ft.map-pb_corrected_realigned.bam',
'PS00628.m84046_240619_124816_s1.bc2074.ft.map-pb_corrected_realigned.bam',
'PS00629.m84046_240619_124816_s1.bc2075.ft.map-pb_corrected_realigned.bam',
'PS00630.m84046_240619_124816_s1.bc2076.ft.map-pb_corrected_realigned.bam',
'PS00631.m84046_240619_124816_s1.bc2077.ft.map-pb_corrected_realigned.bam']

# chr6:110,176,758-110,181,579 WASF1
reg_chrom = 'chr6'
reg_start = 110176758
reg_end = 110181579

da_dicts = []
for i in range(len(wasf1_bams)):
    da_dicts.append(dict())
    bam = pysam.AlignmentFile(f'{wasf1_dir}/{wasf1_bams[i]}', "rb")
    for read in bam.fetch(reg_chrom, reg_start, reg_end):
        seq = read.query_sequence
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        strand = read.get_tag('ST')
        for pos in pair:
            if pos[0] == None or pos[1] == None: # indel, ignore
                pass
            else:
                qi = pos[0]
                ref_coord = pos[1] + 1
                base = seq[qi]
                if strand == 'CT':
                    if pos[2].upper() in ['C']:
                        if ref_coord not in da_dicts[i].keys():
                            da_dicts[i][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                        da_dicts[i][ref_coord][base] += 1
                elif strand == 'GA':
                    if pos[2].upper() in ['G']:
                        if ref_coord not in da_dicts[i].keys():
                            da_dicts[i][ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0, 'Y':0}
                        da_dicts[i][ref_coord][base] += 1
    bam.close()

shared_pos = set.intersection(*[set(d.keys()) for d in da_dicts])

out_pos = []
for p in shared_pos:
    props = []
    for i in range(len(da_dicts)):
        da_dicts[i][p]
        if da_dicts[i][p]['C'] > da_dicts[i][p]['G']:
            prop = da_dicts[i][p]['Y']/(da_dicts[i][p]['Y']+da_dicts[i][p]['C'])
        elif da_dicts[i][p]['C'] < da_dicts[i][p]['G']:
            prop = da_dicts[i][p]['R']/(da_dicts[i][p]['R']+da_dicts[i][p]['G'])
        props.append(prop)
    out_pos.append([p] + props)

with open('prop_da_titration_WASF1.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(['position','PS00626','PS00627','PS00628','PS00629','PS00630','PS00631'])
    for row in out_pos:
        writer.writerow(row)

# average % deaminations by condition
tot_bases = [0 for i in range(len(da_dicts))]
deam_bases = [0 for i in range(len(da_dicts))]
for p in shared_pos:
    for i in range(len(da_dicts)):
        if da_dicts[i][p]['C'] > da_dicts[i][p]['G']:
            tot_bases[i] += da_dicts[i][p]['Y']+da_dicts[i][p]['C']
            deam_bases[i] += da_dicts[i][p]['Y']
        elif da_dicts[i][p]['C'] < da_dicts[i][p]['G']:
            tot_bases[i] += da_dicts[i][p]['R']+da_dicts[i][p]['G']
            deam_bases[i] += da_dicts[i][p]['R']
overall_props = []
for i in range(len(da_dicts)):
    overall_props.append(deam_bases[i] / tot_bases[i])

with open('overall_prop_da_titration_WASF1.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(['PS00626','PS00627','PS00628','PS00629','PS00630','PS00631'])
    writer.writerow(overall_props)

# --------------------------------------------------------------------------------------------------

# mean of the overall % deamination between NAPA & WASF1
df_n = pd.read_csv('overall_prop_da_titration_NAPA.tsv', sep="\t")
df_w = pd.read_csv('overall_prop_da_titration_WASF1.tsv', sep="\t")

means = []
for i in range(len(df_n.columns)):
    pn = df_n.iloc[0,i]
    pw = df_w.iloc[0,i]
    means.append((pn+pw)/2)

with open('average_overall_prop_da_titration_NAPA_and_WASF1.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    writer.writerow(['PS00626','PS00627','PS00628','PS00629','PS00630','PS00631'])
    writer.writerow(means)
