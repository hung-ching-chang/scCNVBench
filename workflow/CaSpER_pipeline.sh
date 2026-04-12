#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#         Step 0: setting         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
sampleID=toy_data
rootDir=/your/local/path/to/repo/
refPath=$rootDir/data/reference_panel/He_et_al_2020/GSM4850590_Stomach_Counts.csv
codeDir=$rootDir/script/
bam_file=$rootDir/data/$sampleID/possorted_genome_bam.bam

outDir=$rootDir/Result/$sampleID/CaSpER/
mkdir $outDir
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#     Step 1: BAF generation      #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
module load samtools/1.22.1

# download "genome_list" from https://www.dropbox.com/s/rq7v67tiou1qwwg/hg38.list?dl=0
samtools view $bam_file | BAFExtract -generate_compressed_pileup_per_SAM stdin $genome_list $outDir 50 0
# download "genome_fa_dir" from https://www.dropbox.com/s/ysrcfcnk7z8gyit/hg38.zip?dl=0
BAFExtract -get_SNVs_per_pileup  $genome_list $outDir $genome_fa_dir 20 4 0.1 $outDir/snps.baf

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#      Step 2: Main function      #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
Rscript $codeDir/CaSpER_run.R \
    --rootDir $rootDir \
    --sampleID $sampleID \
    --refPath $refPath
