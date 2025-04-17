import pysam
import pandas as pd
import csv
import subprocess
import sys


bam_dir = '../'
bams = ['UBA1_region_GA_phased.bam']

uba1_region = '../region_UBA1.bed'

for b in bams:
    bam_path = f"{bam_dir}/{b}"
    extract_H1 = 'GM12878_UBA1_PS00685_ft_extract_all_H1.bed'
    extract_H2 = 'GM12878_UBA1_PS00685_ft_extract_all_H2.bed'
    # Generate BigWigs ------------------------------------------------------------------
    name = 'GM12878_UBA1_PS00685_H1'
    df_H1 = pd.read_csv(extract_H1, sep = "\t")
    # nuc_reg = []
    small_msps = [] # <75bp
    large_msps = [] # >150 bp
    for i in range(len(df_H1)):
        chrom = df_H1.iloc[i,:]['#ct']
        nucs = df_H1.iloc[i,:]['ref_nuc_starts'].rstrip(',').split(',')
        nucl = df_H1.iloc[i,:]['ref_nuc_lengths'].rstrip(',').split(',')
        if len(nucs) > 0 and '.' not in nucs:
            for j in range(len(nucs)):
                if int(nucl[j]) > 0:
                    nuc_reg.append([chrom, int(nucs[j]), int(nucs[j])+int(nucl[j])+1, df_H1.iloc[i,:]['fiber']])
        msps = df_H1.iloc[i,:]['ref_msp_starts'].rstrip(',').split(',')
        mspl = df_H1.iloc[i,:]['ref_msp_lengths'].rstrip(',').split(',')
        if len(msps) > 0 and '.' not in msps:
            for j in range(len(msps)):
                if int(mspl[j]) < 75 and int(mspl[j]) > 0:
                    small_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H1.iloc[i,:]['fiber']])
                elif int(mspl[j]) > 150:
                    large_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H1.iloc[i,:]['fiber']])
    nuc_reg = pd.DataFrame(nuc_reg, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
    small_msps = pd.DataFrame(small_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
    large_msps = pd.DataFrame(large_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
    nuc_bed = f'{name}_nuc_positions.bed'
    small_bed = f'{name}_small_msp_positions.bed'
    large_bed = f'{name}_large_msp_positions.bed'
    nuc_reg.to_csv(nuc_bed, sep="\t", header=False, index=False)
    small_msps.to_csv(small_bed, sep="\t", header=False, index=False)
    large_msps.to_csv(large_bed, sep="\t", header=False, index=False)
    shell_command = f"bedtools intersect -a {nuc_bed} -b {uba1_region} | bedtools genomecov -bg -i stdin -g ../hg38.analysisSet.chrom.sizes > {nuc_bed.replace('.bed','.bg')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedtools intersect -a {small_bed} -b {uba1_region} | bedtools genomecov -bg -i stdin -g ../hg38.analysisSet.chrom.sizes > {small_bed.replace('.bed','.bg')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedtools intersect -a {large_bed} -b {uba1_region} | bedtools genomecov -bg -i stdin -g ../hg38.analysisSet.chrom.sizes > {large_bed.replace('.bed','.bg')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {nuc_bed.replace('.bed','.bg')} ../hg38.analysisSet.chrom.sizes {nuc_bed.replace('.bed','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {small_bed.replace('.bed','.bg')} ../hg38.analysisSet.chrom.sizes {small_bed.replace('.bed','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {large_bed.replace('.bed','.bg')} ../hg38.analysisSet.chrom.sizes {large_bed.replace('.bed','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    name = 'GM12878_UBA1_PS00685_H2'
    df_H2 = pd.read_csv(extract_H2, sep = "\t")
    nuc_reg = []
    small_msps = [] # <75bp
    large_msps = [] # >150 bp
    for i in range(len(df_H2)):
        chrom = df_H2.iloc[i,:]['#ct']
        nucs = df_H2.iloc[i,:]['ref_nuc_starts'].rstrip(',').split(',')
        nucl = df_H2.iloc[i,:]['ref_nuc_lengths'].rstrip(',').split(',')
        if len(nucs) > 0 and '.' not in nucs:
            for j in range(len(nucs)):
                if int(nucl[j]) > 0:
                    nuc_reg.append([chrom, int(nucs[j]), int(nucs[j])+int(nucl[j])+1, df_H2.iloc[i,:]['fiber']])
        msps = df_H2.iloc[i,:]['ref_msp_starts'].rstrip(',').split(',')
        mspl = df_H2.iloc[i,:]['ref_msp_lengths'].rstrip(',').split(',')
        if len(msps) > 0 and '.' not in msps:
            for j in range(len(msps)):
                if int(mspl[j]) < 75 and int(mspl[j]) > 0:
                    small_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H2.iloc[i,:]['fiber']])
                elif int(mspl[j]) > 150:
                    large_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H2.iloc[i,:]['fiber']])
    nuc_reg = pd.DataFrame(nuc_reg, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
    small_msps = pd.DataFrame(small_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
    large_msps = pd.DataFrame(large_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
    nuc_bed = f'{name}_nuc_positions.bed'
    small_bed = f'{name}_small_msp_positions.bed'
    large_bed = f'{name}_large_msp_positions.bed'
    nuc_reg.to_csv(nuc_bed, sep="\t", header=False, index=False)
    small_msps.to_csv(small_bed, sep="\t", header=False, index=False)
    large_msps.to_csv(large_bed, sep="\t", header=False, index=False)
    shell_command = f"bedtools intersect -a {nuc_bed} -b {uba1_region} | bedtools genomecov -bg -i stdin -g ../hg38.analysisSet.chrom.sizes > {nuc_bed.replace('.bed','.bg')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedtools intersect -a {small_bed} -b {uba1_region} | bedtools genomecov -bg -i stdin -g ../hg38.analysisSet.chrom.sizes > {small_bed.replace('.bed','.bg')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedtools intersect -a {large_bed} -b {uba1_region} | bedtools genomecov -bg -i stdin -g ../hg38.analysisSet.chrom.sizes > {large_bed.replace('.bed','.bg')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {nuc_bed.replace('.bed','.bg')} ../hg38.analysisSet.chrom.sizes {nuc_bed.replace('.bed','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {small_bed.replace('.bed','.bg')} ../hg38.analysisSet.chrom.sizes {small_bed.replace('.bed','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {large_bed.replace('.bed','.bg')} ../hg38.analysisSet.chrom.sizes {large_bed.replace('.bed','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    # Prop deaminated ------------------------------------------------------------------
    bam = pysam.AlignmentFile(bam_path, "rb")
    # chrX:47190561-47194939 UBA1
    reg_chrom = 'chrX'
    reg_start = 47190561
    reg_end = 47194939
    H1_pos_da = {}
    H2_pos_da = {}
    # Allow for non-ref C/G to be modified. % R/G|A & % Y/C|T
    for read in bam.fetch(reg_chrom, reg_start, reg_end):
        if read.is_secondary == False and read.is_supplementary == False:
            seq = read.seq
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            hap = read.get_tag('HP')
            for pos in pair:
                mcoord = pos[0]
                ref_coord = pos[1]
                if mcoord != None and pos[1] != None:
                    base = pos[2].upper()
                    if base == 'G':
                        if hap == 1:              
                            if ref_coord not in H1_pos_da.keys():
                                H1_pos_da[ref_coord] = {'total':0, 'da':0}
                            H1_pos_da[ref_coord]['total'] += 1
                            if seq[mcoord] == 'R':
                                H1_pos_da[ref_coord]['da'] += 1
                        elif hap == 2:              
                            if ref_coord not in H2_pos_da.keys():
                                H2_pos_da[ref_coord] = {'total':0, 'da':0}
                            H2_pos_da[ref_coord]['total'] += 1
                            if seq[mcoord] == 'R':
                                H2_pos_da[ref_coord]['da'] += 1
    # out bedgraph of DA proportions
    out_bg = f'GM12878_UBA1_PS00685_H1_DA_density.bg'
    with open(out_bg,'w') as fout:
        writer = csv.writer(fout, delimiter="\t")
        for k in H1_pos_da.keys():
            total = H1_pos_da[k]['total']
            da = H1_pos_da[k]['da']
            if da > 0: # filter 0 positions
                writer.writerow([reg_chrom, k, k+1, da/total])
    shell_command = f"sort -k1,1 -k2,2n {out_bg} > temp.bg"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"mv temp.bg {out_bg}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {out_bg} ../hg38.analysisSet.chrom.sizes {out_bg.replace('.bg','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    out_bg = f'GM12878_UBA1_PS00685_H2_DA_density.bg'
    with open(out_bg,'w') as fout:
        writer = csv.writer(fout, delimiter="\t")
        for k in H2_pos_da.keys():
            total = H2_pos_da[k]['total']
            da = H2_pos_da[k]['da']
            if da > 0: # filter 0 positions
                writer.writerow([reg_chrom, k, k+1, da/total])
    shell_command = f"sort -k1,1 -k2,2n {out_bg} > temp.bg"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"mv temp.bg {out_bg}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f"bedGraphToBigWig {out_bg} ../hg38.analysisSet.chrom.sizes {out_bg.replace('.bg','.bw')}"
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)

