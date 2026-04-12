#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#        Step 1: setting        #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
sampleID=toy_data
rootDir=/your/local/path/to/repo/
CellCyclePath=$rootDir/resources/copyVAE_ref/Macosko_cell_cycle_genes.txt

Out_Dir=$rootDir/Result/$sampleID/copyVAE/
mkdir $Out_Dir
cd $Out_Dir
echo "Method: copyVAE ..."
echo "Process $sampleID ..."
echo "Output Directory: $Out_Dir"

# %%%%%%%%%%%%%%%%% #
#     Step 2: run   #
# %%%%%%%%%%%%%%%%% #
UMIPath=$rootDir/data/$sampleID/filtered_feature_bc_matrix.h5
copyvae $UMIPath $CellCyclePath
