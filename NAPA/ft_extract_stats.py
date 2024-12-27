import pandas as pd
import numpy as np


# Will SOON no Longer necessary!!! Use ft pileup instead !!!!!!!!!!
# https://fiberseq.github.io/fibertools/help.html


# ft extract NAPA_nuc_CT.bam --all ft_extract_all_CT_NAPA.bed
# ft extract NAPA_nuc_GA.bam --all ft_extract_all_GA_NAPA.bed

df_ct = pd.read_csv('ft_extract_all_CT_NAPA.bed', sep = "\t")
df_ga = pd.read_csv('ft_extract_all_GA_NAPA.bed', sep = "\t")
merged = pd.concat([df_ct,df_ga])

nuc_reg = []
small_msps = [] # <75bp
large_msps = [] # >150 bp
for i in range(len(merged)):
    chrom = merged.iloc[i,:]['#ct']
    nucs = merged.iloc[i,:]['ref_nuc_starts'].rstrip(',').split(',')
    nucl = merged.iloc[i,:]['ref_nuc_lengths'].rstrip(',').split(',')
    if len(nucs) > 0 and '.' not in nucs:
        for j in range(len(nucs)):
            if int(nucl[j]) > 0:
                nuc_reg.append([chrom, int(nucs[j]), int(nucs[j])+int(nucl[j])+1, merged.iloc[i,:]['fiber']])
    msps = merged.iloc[i,:]['ref_msp_starts'].rstrip(',').split(',')
    mspl = merged.iloc[i,:]['ref_msp_lengths'].rstrip(',').split(',')
    if len(msps) > 0 and '.' not in msps:
        for j in range(len(msps)):
            if int(mspl[j]) < 75 and int(mspl[j]) > 0:
                small_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, merged.iloc[i,:]['fiber']])
            elif int(mspl[j]) > 150:
                large_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, merged.iloc[i,:]['fiber']])

nuc_reg = pd.DataFrame(nuc_reg, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
small_msps = pd.DataFrame(small_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
large_msps = pd.DataFrame(large_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
nuc_reg.to_csv('nuc_positions.bed', sep="\t", header=False, index=False)
small_msps.to_csv('small_msp_positions.bed', sep="\t", header=False, index=False)
large_msps.to_csv('large_msp_positions.bed', sep="\t", header=False, index=False)


# bedtools intersect -a nuc_positions.bed -b ../napa_region.bed | bedtools genomecov -bg -i stdin -g /gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > nuc_positions.bg
# bedtools intersect -a small_msp_positions.bed -b ../napa_region.bed | bedtools genomecov -bg -i stdin -g /gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > small_msp_positions.bg
# bedtools intersect -a large_msp_positions.bed -b ../napa_region.bed | bedtools genomecov -bg -i stdin -g /gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes > large_msp_positions.bg
# /gscratch/stergachislab/install_dir/bedGraphToBigWig nuc_positions.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes nuc_positions.bw
# /gscratch/stergachislab/install_dir/bedGraphToBigWig small_msp_positions.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes small_msp_positions.bw
# /gscratch/stergachislab/install_dir/bedGraphToBigWig large_msp_positions.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes large_msp_positions.bw




# losing HP tag after ft extract !!!!!!!!!!!!!!!!!!!!!!!!!!


# filter to phased and non-empty. Filter by length as well ?????

