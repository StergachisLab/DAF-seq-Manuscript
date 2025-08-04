import csv

# Convert GFF to BED file of isoforms

out_lines = []
with open('GM12878_mas_collapsed.sorted.gff') as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        if len(line) > 1:
            if line[2] == "transcript":
                tsid = line[8].split(';')[-2].split(' ')[2].strip('"')
                if line [6] == "+":
                    start = int(line[3])
                    out_lines.append([line[0], start, start+1, tsid, '.', line[6]])
                elif line [6] == "-":
                    start = int(line[4])
                    out_lines.append([line[0], start, start+1, tsid, '.', line[6]])

with open('isoseq_counts_giab_hg002.bed','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    writer.writerows(out_lines)
