library(data.table)
library(readr)
library(dplyr)
library(stringr)



# get hg38 p-arm and q-arm coordiantes from UCSC file: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz

arms_df <- read_tsv("cytoBand.txt.gz", col_names = c("chrom","chromStart","chromEnd","name","gieStain")) |> 
                mutate(arm = substring(name, 1, 1)) |> 
                group_by(chrom, arm) |> 
                summarise(start_pos = min(chromStart),
                            end_pos = max(chromEnd),
                            length = end_pos - start_pos)

keep <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

main_chrom <- arms_df %>% filter(chrom %in% keep)
write_tsv(main_chrom, 'chrom_arm_coords_filtered.tsv')


bed_format <- main_chrom %>% select(chrom,start_pos,end_pos,arm)
write_tsv(bed_format, 'chrom_arm_coords_filtered.bed', col_names = FALSE)

