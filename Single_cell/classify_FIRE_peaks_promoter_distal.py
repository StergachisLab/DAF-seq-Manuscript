""" Classify HG002 FIRE peaks as promoter-proximal or promoter-distal. """


# # filtering genocde to Ensembl Canonical TSSs
# cat gencode.v45.annotation_all_tss.bed | grep Ensembl_canonical | sort-bed - > gencodev45_Ensembl_canonical_TSS.bed

# # Within 1 Kb -- TSS set has been expanded from 1 bp to 20 bp for another project !!!!!
# zcat FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz | sort-bed - | bedmap --ec --echo-ref-name --echo-map --skip-unmapped --range 990 - gencodev45_Ensembl_canonical_TSS.bed > FIRE-peak_HG002_noUnrealiableCov_autosome_v45_TSS_bedmap.txt


tss_peaks = set()
with open('FIRE-peak_HG002_noUnrealiableCov_autosome_v45_TSS_bedmap.txt') as fr:
    for line in fr:
        peak = line.strip().split('|')[0]
        tss_peaks.add(peak)

with open('promoter-proximal_FIRE-peaks_HG002_noUnrealiableCov_autosome.txt','w') as fw:
    for p in tss_peaks:
        fw.write(f'{p}\n')

