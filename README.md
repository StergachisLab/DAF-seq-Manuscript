# DAF-seq-Manuscript

Code used in data analysis as part of the DAF-seq manuscript: https://doi.org/10.1101/2024.11.06.622310.

Data processing and analyses are organized in three main folders: “General” for DAF-seq preprocessing such as strand assignment, “Targeted” and “Single_cell”. Environments used in these analyses are provided in “envs”. Scripts associated together for a particular analysis are numbered in the order they are to be executed (ex. 01, 02, etc.). Raw data for each analysis can be downloaded from SRA by the user into the applicable folder (“data” for targeted and single-cell analyses).

## General
The script “process_DddA_bam.py” takes an aligned BAM file as a command line argument and identifies putative deamination events. It then assigns read strands by adding a “ST” tag which designate top strand reads as “CT”, bottom strand reads as “GA”, and unassigned reads as “undetermined”. It then uses IUPAC degenerate base characters to mark deamination events by updating C-to-T bases to Y and G-to-A bases to R for top and bottom stranded reads, respectively. All primary alignments will be written to a new “corrected” BAM file.


## Targeted
Targeted DAF-seq analyses are divided my genomic region (ex. NAPA). Downloaded data is expected to be stored in the "data" folder.

The expected order of NAPA analyses is alignment and preprocessing, followed by scripts in the "footprinting", "qc_plots", "saturation" folders. Thermodynamic analysis of TF occupancy at the NAPA promoter is contained with the "Thermodynamic_Daf-seq_footprinting" folder within the "footprinting" section.

The expected order of UBA1 analyses is alignment and preprocessing, followed by cripts in the "pileups", "footprinting", "TSS_codependency" folders. The "FINAL_UBA1_plots.R" script creates the figure plots.


## Single Cell
Single cell analyses are organized by folder and are expected to be processed in the following order: data, collapse, phasing, msp_analysis, CTCF_occupancy, msp_codependency, scDAF_figures.

Fiber-seq data for HG002 and GM12878 can be downloaded from the links below. HG002FIRE peaks are required for scDAF-seq chromatin actuation analyses.

https://s3-us-west-1.amazonaws.com/stergachis-manuscript-data/index.html?prefix=2024/Vollger_et_al/FIRE/HG002/

https://s3-us-west-1.amazonaws.com/stergachis-manuscript-data/index.html?prefix=2024/Vollger_et_al/FIRE/GM12878/

https://s3-us-west-1.amazonaws.com/stergachis-manuscript-data/2024/Vollger_et_al/FIRE/HG002/FDR-peaks/FDR-FIRE-peaks.bed.gz


## Accessory Data
Data files used in analyses that were not included in the SRA upload are available at the following Zenodo repository: 10.5281/zenodo.15391983

