# %%%%%%%%%%%%%%%%%%%%%%%% #
#        Clonalscope       # 
# %%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute clonalscope')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
parser = add_option(parser, '--normalCell', default = "", type = "character")
args <- parse_args(parser)

library(Clonalscope)
library(Seurat)
root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/Clonalscope/")
normal_Cell <- args$normalCell

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: clonalscope ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
input_expr <- as.matrix(Read10X(data.dir = input_Dir, strip.suffix = TRUE))
num.input.cell <- ncol(input_expr)
input.barcodes <- colnames(input_expr)
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
num.ref.cell <- ncol(ref_expr)
colnames(ref_expr) <- paste0("Ref_", colnames(ref_expr))
ref.barcodes <- colnames(ref_expr)

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
input_expr <- cbind(input_expr[overlapped.genes,], ref_expr[overlapped.genes,])
mtx <- as(as(as(input_expr, "dMatrix"), "generalMatrix"), "TsparseMatrix")
barcodes <- data.frame(V1 = c(input.barcodes, ref.barcodes))
features <- read.table(paste0(input_Dir, "/features.tsv.gz"), stringsAsFactors = F, sep='\t', header=F)
features <- features[!duplicated(features$V2),]
features <- features[features$V2 %in% rownames(mtx),]
mtx <- mtx[features$V2,]

celltype <- data.frame(barcode = barcodes[,1],
                       celltype = rep(c("unknown","normal"),c(num.input.cell, num.ref.cell)))

# Size of each chromosome
size <- read.table(paste0(root_Dir, "/resources/Clonalscope_ref/sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = F)

# List of cyclegenes retrieved from the "CopyKAT" package
cyclegenes <- readRDS(paste0(root_Dir, "/resources/Clonalscope_ref/cyclegenes.rds"))

# bed file indicating gene positions
bed <- read.table(paste0(root_Dir, "/resources/Clonalscope_ref/hg38.genes.bed"), sep='\t', header = T)

# Select barcodes for clustering only on epithelial/tumor cells
consensus.malignant <- read.csv(normal_Cell, row.names = 1)
malignant.bc <- rownames(consensus.malignant)[consensus.malignant$consensus.score >= 2]
malignant.bc <- gsub("-1", "", malignant.bc)
clustering_barcodes <- malignant.bc[malignant.bc %in% input.barcodes]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
## step 0: identify malignant cell (for B1)
chrarm <- read.table(paste0(root_Dir, "/resources/Clonalscope_ref/cytoarm_table_hg38.txt"), 
                     stringsAsFactors = F, sep='\t', header=T)
chrarm <- chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
bin_bed <- chrarm[,-2]
seg_table_filtered <- data.frame("chr" = bin_bed[,1], 'start' = as.numeric(bin_bed[,2]),
                                 'end' = as.numeric(bin_bed[,3]), 'states' = 1, 'length' = as.numeric(bin_bed[,3]) - as.numeric(bin_bed[,2]),
                                 'mean' = 0, 'var' = 0, 'Var1' = 1:nrow(bin_bed),'Freq' = 50000,
                                 'chrr' = paste0(bin_bed[,1],":", bin_bed[,2]), stringsAsFactors = F)

# Select barcodes for clustering only on epithelial/tumor cells
all_barcodes <- input.barcodes

# Filtering HLA and cell cycle genes
Input_filtered <- FilterFeatures(mtx = mtx, barcodes = barcodes, features = features, cyclegenes = cyclegenes)

# mincell parameter (default value set as 5% of all clustering cells)
mincell <- length(all_barcodes) * 0.05

Cov_obj <- RunCovCluster(mtx = Input_filtered$mtx, barcodes = Input_filtered$barcodes,
                         features = Input_filtered$features, bed = bed,
                         celltype0 = celltype, var_pt = 0.99, var_pt_ctrl = 0.99, include = 'all',
                         alpha_source = 'all', ctrl_region = NULL,
                         seg_table_filtered = seg_table_filtered, size = size,
                         dir_path = Out_Dir, breaks = 50, prep_mode = 'intersect', est_cap = 2,
                         clust_mode = 'all',clustering_barcodes = all_barcodes, mincell = mincell)
tumor_normal_df <- data.frame(cell.names = names(Cov_obj$result_final$result$Zest),
                              Clonalscope.pred = ifelse(Cov_obj$result_final$result$annot == "T", "tumor", "normal"))
write.table(tumor_normal_df, file = paste0(Out_Dir, "tumor_cell_pred.txt"),
            row.names = F, quote = F)


## step 1: First Round Clonalscope Estimation using chromosomal arms (on malignant cells)
chrarm <- read.table(paste0(root_Dir, "/resources/Clonalscope_ref/cytoarm_table_hg38.txt"), 
                     stringsAsFactors = F, sep='\t', header=T)
chrarm <- chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
bin_bed <- chrarm[,-2]
seg_table_filtered <- data.frame("chr" = bin_bed[,1], 'start' = as.numeric(bin_bed[,2]),
                                 'end' = as.numeric(bin_bed[,3]), 'states' = 1, 'length' = as.numeric(bin_bed[,3]) - as.numeric(bin_bed[,2]),
                                 'mean' = 0, 'var' = 0, 'Var1' = 1:nrow(bin_bed),'Freq' = 50000,
                                 'chrr' = paste0(bin_bed[,1],":", bin_bed[,2]), stringsAsFactors = F)

# Filtering HLA and cell cycle genes
Input_filtered <- FilterFeatures(mtx = mtx, barcodes = barcodes, features = features, cyclegenes = cyclegenes)

# mincell parameter (default value set as 5% of all clustering cells)
mincell <- length(clustering_barcodes) * 0.05

# Remove raw inputs
rm(mtx); rm(barcodes); rm(features)

Cov_obj <- RunCovCluster(mtx = Input_filtered$mtx, barcodes = Input_filtered$barcodes,
                         features = Input_filtered$features, bed = bed,
                         celltype0 = celltype, var_pt = 0.99, var_pt_ctrl = 0.99, include = 'all',
                         alpha_source = 'all', ctrl_region = NULL,
                         seg_table_filtered = seg_table_filtered, size = size,
                         dir_path = Out_Dir, breaks = 50, prep_mode = 'intersect', est_cap = 2,
                         clust_mode = 'all',clustering_barcodes = clustering_barcodes, mincell = mincell)
saveRDS(Cov_obj, file = paste0(Out_Dir,"Cov_obj_chrarm.rds"))


# 2. Identify segments with finer resolution from 1st round estimation of Clonalscope
clustering = Cov_obj$result_final$clustering
clustering2 = Cov_obj$result_final$clustering2
result = Cov_obj$result_final$result
Zest = Cov_obj$result_final$result$Zest
pdf(paste0(Out_Dir,"chrarm_cluster_heatmap.pdf"), width = 10, height = 7)
print(PlotClusters(df = Cov_obj$result_final$df_obj$df, celltype = celltype, 
                   Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode='annot'))
dev.off()

# load 200k bp bins for each cell
bin_bed = read.table(paste0(root_Dir, "/resources/Clonalscope_ref/hg38_200kb.windows.bed"))

# load standard gene annotation file
gtf <- readRDS((paste0(root_Dir, "/resources/Clonalscope_ref/gtf_filtered.rds")))
# filter out genes without chr:pos info (in annotation file) --- which may cause error in CreateSegtableNoWGS function
gtf_ensg <- sapply(strsplit(gtf$V9, ";"), "[", 1)
gtf_ensg <- sapply(strsplit(gtf_ensg, "gene_id "), "[", 2)
olap_ensg <- intersect(Input_filtered$features$V1, gtf_ensg)
gtf <- gtf[gtf_ensg %in% olap_ensg,]
Input_filtered$features <- Input_filtered$features[Input_filtered$features$V1 %in% olap_ensg,]
Input_filtered$mtx <- Input_filtered$mtx[rownames(Input_filtered$mtx) %in% Input_filtered$features$V2,]
# only keep autosomes (chr 1- 22)
bin_bed = bin_bed[bin_bed[,1] %in% paste0("chr",c(1:22)),]
seg_filtered <- CreateSegtableNoWGS(mtx = Input_filtered$mtx, barcodes = Input_filtered$barcodes, features=Input_filtered$features, 
                                    Cov_obj = Cov_obj, Zest = Zest, gtf = gtf, size = size,
                                    bin_bed = bin_bed, celltype0 = celltype, dir_path = Out_Dir,
                                    bin_mtx = NULL, # estimating bin x cell count matrix
                                    plot_seg = TRUE, hmm_states = c(0.5, 1.5, 2), max_qt = 0.95, 
                                    nmean = 500, rm_extreme = 2, adj = -0.5 # HMM parameters that can be adjusted by users
)
# 3. Running Clonalscope again with the new segments
# second round estimation with new segments
# clustering_barcodes <- celltype[which(grepl("unknown",celltype[,2])),1]  # only clustering on Tumor cells
second_seg_table <- data.frame("chr" = seg_filtered[,1], 'start' = as.numeric(seg_filtered[,2]),
                               'end' = as.numeric(seg_filtered[,3]), 'states' = 1, 'length' = as.numeric(seg_filtered[,3])-as.numeric(seg_filtered[,2]),
                               'mean' = 0, 'var' = 0, 'Var1' = 1:nrow(seg_filtered),'Freq' = 50000,
                               'chrr' = paste0(seg_filtered[,1],":", seg_filtered[,2]), stringsAsFactors = F)

# Wrapper function of Clonalscope
set.seed(2025)
Cov_obj2 <- RunCovCluster(mtx = Input_filtered$mtx, barcodes = Input_filtered$barcodes, 
                          features = Input_filtered$features, bed = bed, 
                          celltype0 = celltype, var_pt = 0.99, var_pt_ctrl = 0.99, include = 'all',
                          alpha_source = 'all', ctrl_region = NULL, 
                          seg_table_filtered = second_seg_table, size = size,
                          dir_path = Out_Dir, breaks = 50, prep_mode = 'intersect', est_cap = 2,
                          clust_mode = 'all',clustering_barcodes = clustering_barcodes, mincell = mincell) 

saveRDS(Cov_obj2, file = paste0(Out_Dir,"/Clonalscope_final.rds"))

result = Cov_obj2$result_final$result
pdf(paste0(Out_Dir,"seg_cluster_heatmap.pdf"), width = 10, height = 7)
print(PlotClusters(df = Cov_obj2$result_final$df_obj$df, celltype = celltype, 
                   Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode='annot'))
dev.off()