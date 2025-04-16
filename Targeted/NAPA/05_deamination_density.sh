#!/bin/bash
set -euo pipefail

sort -k1,1 -k2,2n NAPA_PS00626_DA_density.bg > temp.bg
mv temp.bg NAPA_PS00626_DA_density.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig NAPA_PS00626_DA_density.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes NAPA_PS00626_DA_density.bw

sort -k1,1 -k2,2n NAPA_PS00626_Hap1_DA_density.bg > temp.bg
mv temp.bg NAPA_PS00626_Hap1_DA_density.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig NAPA_PS00626_Hap1_DA_density.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes NAPA_PS00626_Hap1_DA_density.bw

sort -k1,1 -k2,2n NAPA_PS00626_Hap2_DA_density.bg > temp.bg
mv temp.bg NAPA_PS00626_Hap2_DA_density.bg
/gscratch/stergachislab/install_dir/bedGraphToBigWig NAPA_PS00626_Hap2_DA_density.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes NAPA_PS00626_Hap2_DA_density.bw
