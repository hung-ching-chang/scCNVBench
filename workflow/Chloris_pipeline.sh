#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#         Step 0: setting         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
sampleID=toy_data
rootDir=/your/local/path/to/repo/
refPath=$rootDir/data/reference_panel/He_et_al_2020/GSM4850590_Stomach_Counts.csv
codeDir=$rootDir/script/
outDir=$rootDir/Result/$sampleID/Chloris/
# create output folder
mkdir $outDir
cd $outDir

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#         Step 1: cellsnp         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
bam=$rootDir/data/$sampleID/possorted_genome_bam.bam
cb=$outDir/cell_barcodes.tsv
gunzip < $rootDir/data/$sampleID/filtered_feature_bc_matrix/barcodes.tsv.gz > $cb

# SNP vcf
# cd $rootDir/resources/download/
# wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
# gunzip genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
cellsnp-lite -s $rootDir/data/$sampleID/possorted_genome_bam.bam \
             -b $outDir/cell_barcodes.tsv \
             -O $outDir \
             -R $rootDir/resources/download/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
             -p 25 --minMAF 0.1 --minCOUNT 20 --UMItag Auto --cellTAG CB

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Step 2: run           #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
Rscript $codeDir/Chloris_run.R \
    --rootDir $rootDir \
    --sampleID $sampleID \
    --refPath $refPath

