#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#         Step 0: setting         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
sampleID=toy_data
rootDir=/your/local/path/to/repo/
refPath=$rootDir/data/reference_panel/He_et_al_2020/GSM4850590_Stomach_Counts.csv
refannoPath=$rootDir/data/reference_panel/He_et_al_2020/Annotation_AHCA_alltissues_meta.data_84363_cell.txt
refOrgan=Stomach_cDNA
codeDir=$rootDir/script/

alleleDir=$rootDir/resources/download/
outDir=$rootDir/Result/$sampleID/Numbat/
# create output folder
mkdir $outDir
cd $outDir

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#          Step 1: pileup         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
bam=$rootDir/data/$sampleID/possorted_genome_bam.bam
cb=$outDir/cell_barcodes.tsv
gunzip < $rootDir/data/$sampleID/filtered_feature_bc_matrix/barcodes.tsv.gz > $cb

# download "1000G SNP VCF" and "1000G Reference Panel" to $alleleDir
# wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
# wget http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip
gmap=$alleleDir/genetic_map_hg38_withX.txt.gz
vcf=$alleleDir/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf

Rscript $codeDir/util/pileup_and_phase.R \
    --label $sampleID \
    --samples $sampleID \
    --bams $bam \
    --barcodes $cb \
    --snpvcf $vcf \
    --gmap $gmap \
    --paneldir $alleleDir/1000G_hg38 \
    --outdir $outDir \
    --ncores 25
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Step 2: run           #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
Rscript $codeDir/Numbat_run.R \
    --rootDir $rootDir \
    --sampleID $sampleID \
    --refPath $refPath \
    --refannoPath $refannoPath \
    --refOrgan $refOrgan

