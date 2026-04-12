# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#.     SCEVAN with reference      # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute SCEVAN')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

library(SCEVAN)
library(Seurat)
root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/SCEVAN/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: SCEVAN ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
input_expr <- as.matrix(Read10X(data.dir = input_Dir, strip.suffix = TRUE))
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
colnames(ref_expr) <- gsub("[.]", "-", colnames(ref_expr))
colnames(ref_expr) <- paste0("Ref_", colnames(ref_expr))

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
input_expr <- cbind(input_expr[overlapped.genes,], ref_expr[overlapped.genes,])
normal_cell_names <- colnames(ref_expr)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
results <- pipelineCNA(input_expr, 
                       plotTree = F,
                       sample = sample_id, 
                       par_cores = 1, # parallel computin may cause some issues (https://github.com/AntonioDeFalco/SCEVAN/issues/115)
                       SUBCLONES = TRUE,
                       norm_cell = normal_cell_names)
saveRDS(results, paste0(Out_Dir, "SCEVAN_res.RDS"))