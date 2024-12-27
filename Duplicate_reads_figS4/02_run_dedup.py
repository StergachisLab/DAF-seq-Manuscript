

# GM12878 --------------------------
bams = ['/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/PS00626_map-pb_corrected_realigned_NAPA_region.bam', # NAPA (4 uM 10 min, PS00626)
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/PS00626_map-pb_corrected_realigned_WASF1_region.bam', # WASF1 (4 uM 10 min, PS00626),
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region.bam', # SLC39A4
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/UBA1/PS00685_GM12878_UBA1_map-pb_corrected_realigned_UBA1_region.bam', # UBA1
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/COLO_Mix/PS00719_COLO_Region2_map-pb_corrected_realigned_COLO_region.bam', # COLO829 BL/T
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region.bam', # Liver
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/Colon_SLC39A4_PS00682_map-pb_corrected_realigned_SLC39A4_region.bam', # Colon
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/Heart_SLC39A4_PS00681_map-pb_corrected_realigned_SLC39A4_region.bam', # Heart
'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/PS00744_HepG2_SLC39A4_map-pb_corrected_realigned_SLC39A4_region.bam'] # HepG2

for b in bams:
    print(f"python duplicate_ID.py --bam {b} &")

# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/PS00626_map-pb_corrected_realigned_NAPA_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/PS00626_map-pb_corrected_realigned_WASF1_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/GM12878_SLC39A4_PS00686_map-pb_corrected_realigned_SLC39A4_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/UBA1/PS00685_GM12878_UBA1_map-pb_corrected_realigned_UBA1_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/COLO_Mix/PS00719_COLO_Region2_map-pb_corrected_realigned_COLO_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/Liver_SLC39A4_PS00680_map-pb_corrected_realigned_SLC39A4_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/Colon_SLC39A4_PS00682_map-pb_corrected_realigned_SLC39A4_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/Heart_SLC39A4_PS00681_map-pb_corrected_realigned_SLC39A4_region.bam &
# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/PS00744_HepG2_SLC39A4_map-pb_corrected_realigned_SLC39A4_region.bam &


# python duplicate_ID.py --bam /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/processed_bams/Colon_SLC39A4_PS00682_map-pb_corrected_realigned_SLC39A4_region.bam
