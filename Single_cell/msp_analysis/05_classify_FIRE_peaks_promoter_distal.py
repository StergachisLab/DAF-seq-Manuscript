""" Classify HG002 FIRE peaks as promoter-proximal or promoter-distal. """

tss_peaks = set()
with open('FIRE-peak_HG002_noUnrealiableCov_autosome_v45_TSS_bedmap.txt') as fr:
    for line in fr:
        peak = line.strip().split('|')[0]
        tss_peaks.add(peak)

with open('promoter-proximal_FIRE-peaks_HG002_noUnrealiableCov_autosome.txt','w') as fw:
    for p in tss_peaks:
        fw.write(f'{p}\n')

