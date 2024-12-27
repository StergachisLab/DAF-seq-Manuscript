import pysam
import csv


bam_dir='/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/consensus_bams'

samples = ['PS00718','PS00756','PS00757','PS00758']
CT_bams = [f'{bam_dir}/{s}_consensus_CT_HG38_corrected.bam' for s in samples]
GA_bams = [f'{bam_dir}/{s}_consensus_GA_HG38_corrected.bam' for s in samples]

out_rows = [['Cell','DA_bp','Tot_bp','prop_CT_da','prop_GA_da','prop_da_both_strands']]

def da_rate(bam_name):
    total_bp = 0
    da_bp = 0
    total_C = 0
    total_G = 0
    bam = pysam.AlignmentFile(bam_name, "rb")
    for read in bam.fetch():
        total_bp += len(read.seq)
        da_bp += read.seq.count('R')
        da_bp += read.seq.count('Y')
        total_G += read.seq.count('G')
        total_C += read.seq.count('C')
        total_G += read.seq.count('R')
        total_C += read.seq.count('Y')
    return(da_bp, total_bp, total_C, total_G)

for i in range(len(samples)):
    ct_da,ct_tot,ct_tot_c,ct_tot_g = da_rate(CT_bams[i])
    ga_da,ga_tot,ga_tot_c,ga_tot_g = da_rate(GA_bams[i])
    out_rows.append([samples[i], ct_da+ga_da, ct_tot+ga_tot, ct_da/ct_tot_c, ga_da/ga_tot_g, (ct_da+ga_da)/(ct_tot_c+ga_tot_g)])

with open('deamination_rate_by_cell.csv','w') as fw:
    writer = csv.writer(fw)
    for line in out_rows:
        writer.writerow(line)

