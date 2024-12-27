import csv
import subprocess


motif_bed = 'motifs/fimo_motifs_filtered.bed'

CT_BAM = 'NAPA_nuc_CT.bam'
GA_BAM = 'NAPA_nuc_GA.bam'

# dict to group by motif type
motifs = dict()

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


# run ft footprint
for k in motifs.keys():
    bed = f'motifs/ft_{k}_motifs.bed'
    yaml = f'modules/ft_{k}_modules.yaml'
    shell_command = f'ft footprint {CT_BAM} {bed} {yaml} > ft_out/ft_CT_{k}.bed'
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)
    shell_command = f'ft footprint {GA_BAM} {bed} {yaml} > ft_out/ft_GA_{k}.bed'
    subprocess.run(shell_command, shell=True, check=True, capture_output=True, text=True)

# manually ran ft footprint for CTCF using individual modules
# ft footprint nuc_CT.bam motifs/ctcf_motif.bed modules/ctcf_modules.yaml > ft_out/ft_CT_ModCTCF.bed
# ft footprint nuc_GA.bam motifs/ctcf_motif.bed modules/ctcf_modules.yaml > ft_out/ft_GA_ModCTCF.bed

# ft footprint nuc_CT.bam motifs/ctcf_motif.bed modules/ctcf_modules.yaml > ft_out/ft_CT_Mod23CTCF.bed
# ft footprint nuc_GA.bam motifs/ctcf_motif.bed modules/ctcf_modules.yaml > ft_out/ft_GA_Mod23CTCF.bed
