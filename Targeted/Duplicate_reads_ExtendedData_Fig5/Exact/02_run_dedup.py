
bams = ['../../data/PS00626_map-pb_corrected_realigned_NAPA_region.bam', # NAPA (4 uM 10 min, PS00626)
'../../data/PS00626_map-pb_corrected_realigned_WASF1_region.bam', # WASF1 (4 uM 10 min, PS00626),
'../../data/GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region.bam', # SLC39A4
'../../data/PS00685_GM12878_UBA1_map-pb_corrected_realigned_UBA1_region.bam', # UBA1
'../../data/PS00719_COLO_Region2_map-pb_corrected_realigned_COLO_region.bam', # COLO829 BL/T
'../../data/Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region.bam', # Liver
'../../data/Colon_SLC39A4_PS00682_map-pb_corrected_realigned_SLC39A4_region.bam', # Colon
'../../data/Heart_SLC39A4_PS00681_map-pb_corrected_realigned_SLC39A4_region.bam'] # Heart

for b in bams:
    print(f"python duplicate_ID.py --bam {b} &")
