#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#         Step 0: setting         #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
sampleID=toy_data
rootDir=/your/local/path/to/repo/
refPath=$rootDir/data/reference_panel/He_et_al_2020/GSM4850590_Stomach_Counts.csv
codeDir=$rootDir/script/

# %%%%%%%%%%%%%%%%%%%% #
#      Step 1: run     #
# %%%%%%%%%%%%%%%%%%%% #
Rscript $codeDir/CopyKAT_run.R \
    --rootDir $rootDir \
    --sampleID $sampleID \
    --refPath $refPath

