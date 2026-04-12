# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#   inferCNV with reference    # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute inferCNV')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

library(infercnv)
library(Seurat)
root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/inferCNV/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: inferCNV ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
input_expr <- as.matrix(Read10X(data.dir = input_Dir, strip.suffix = TRUE))
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
colnames(ref_expr) <- paste0("Ref_", colnames(ref_expr))

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
merge_expr <- cbind(input_expr[overlapped.genes,], ref_expr[overlapped.genes,])
countfile <- paste0(Out_Dir, "count.txt")
write.table(merge_expr, countfile, sep = "\t", 
            row.names = T, col.names = T, quote = F)

input_anno <- cbind(colnames(merge_expr), 
                    c(rep("unknown", ncol(input_expr)), 
                      rep("normal", ncol(ref_expr))))
annofile <- paste0(Out_Dir, "anno.txt")
write.table(input_anno, annofile, sep = "\t", 
            row.names = F, col.names = F, quote = F)
gene_order_path <-  paste0(root_Dir, "resources/inferCNV_ref/hg38_gencode_v27.txt") # download from https://data.broadinstitute.org/Trinity/CTAT/cnv/

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = countfile,
                                    annotations_file = annofile,
                                    delim="\t",
                                    gene_order_file = gene_order_path,
                                    ref_group_names = "normal")

# use hclust method (recommended for small sample size)
options(scipen = 100)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir = Out_Dir,  # dir is auto-created for storing outputs
                             cluster_by_groups = FALSE,   # T for multiple samples
                             HMM = TRUE,
                             HMM_type = "i6",
                             analysis_mode = "subclusters",
                             HMM_report_by = "subcluster",
                             tumor_subcluster_partition_method = "random_trees",
                             hclust_method = 'ward.D2',
                             tumor_subcluster_pval = 0.05, 
                             num_threads = 10, # don't set too many threads, which may require lots of memory
                             denoise = TRUE)

