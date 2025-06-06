import csv
import subprocess
from glob import glob
import os


motif_bed = 'motifs/filtered_fimo.bed'

motifs = dict() # dict to group by motif type
with open(motif_bed) as fr:
    reader = csv.reader(fr, delimiter = "\t")
    for line in reader:
        if line[3] not in motifs.keys():
            motifs[line[3]] = [line]
        else:
            motifs[line[3]].append(line)
# write motif BEDs for ft footprint
for k,v in motifs.items():
    with open(f'motifs/ft_{k}_motifs.bed','w') as outBed:
        writer = csv.writer(outBed, delimiter = "\t")
        for i in range(len(v)):
            writer.writerow(v[i])
# write motif module YAMLs for ft footprint
for k,v in motifs.items():
    motif_length = int(v[0][2]) - int(v[0][1])
    with open(f'modules/ft_{k}_modules.yaml','w') as yaml:
        yaml.write('modules:\n')
        yaml.write(f'  - [0, {motif_length}]\n')


H1_BAM = '../pileups/GM12878_UBA1_PS00685_nuc_H1.bam'
H2_BAM = '../pileups/GM12878_UBA1_PS00685_nuc_H2.bam'

# run ft footprint
name = 'H1_GM12878_UBA1_PS00685'
for k in motifs.keys():
    bed = f'motifs/ft_{k}_motifs.bed'
    yaml = f'modules/ft_{k}_modules.yaml'
    shell_command = f'/mmfs1/gscratch/stergachislab/bin/ft footprint -x "len(msp)>150" -t 30 --bed {bed} --yaml {yaml} {H1_BAM} > ft_out/{name}_ft_{k}.bed'
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
name = 'H2_GM12878_UBA1_PS00685'
for k in motifs.keys():
    bed = f'motifs/ft_{k}_motifs.bed'
    yaml = f'modules/ft_{k}_modules.yaml'
    shell_command = f'/mmfs1/gscratch/stergachislab/bin/ft footprint -x "len(msp)>150" -t 30 --bed {bed} --yaml {yaml} {H2_BAM} > ft_out/{name}_ft_{k}.bed'
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)


# manually ran ft footprint for CTCF using individual modules
ft footprint -x "len(msp)>150" -t 30 --bed motifs/canonical_CTCF_motif.bed --yaml modules/ctcf_modules.yaml $H1_BAM > ft_out/ft_H1_ModCTCF.bed
ft footprint -x "len(msp)>150" -t 30 --bed motifs/canonical_CTCF_motif.bed --yaml modules/ctcf_modules.yaml $H2_BAM > ft_out/ft_H2_ModCTCF.bed
