# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#      CopyKAT with reference     # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(optparse)
parser = OptionParser(description = 'Execute CopyKAT')
parser = add_option(parser, '--rootDir', default = "/your/local/path/to/repo/", type = "character")
parser = add_option(parser, '--sampleID', default = "toy_data", type = "character")
parser = add_option(parser, '--refPath', default = "", type = "character")
args <- parse_args(parser)

library(copykat)
library(Seurat)
root_Dir <- args$rootDir
sample_id <-args$sampleID
ref_Path <- args$refPath
input_Dir <- paste0(root_Dir, "data/", sample_id, "/filtered_feature_bc_matrix")
Out_Dir <- paste0(root_Dir, "Result/", sample_id, "/CopyKAT_wRef/")

dir.create(Out_Dir, showWarnings = FALSE, recursive = TRUE)
setwd(Out_Dir)
message("Method: CopyKAT ...")
message("Process ", sample_id, " ...")
message("Output Directory: ", Out_Dir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#           Input data            # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
input_expr <- as.matrix(Read10X(data.dir = input_Dir, strip.suffix = TRUE))
colnames(input_expr) <- gsub("-1", "", colnames(input_expr))
ref_expr <- as.matrix(read.csv(ref_Path, row.names = 1))
colnames(ref_expr) <- paste0("Ref_", colnames(ref_expr))

overlapped.genes <- intersect(rownames(input_expr), rownames(ref_expr))
input_expr <- cbind(input_expr[overlapped.genes,], ref_expr[overlapped.genes,])
normal_cell_names <- colnames(ref_expr)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#            Execute              # 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
copykat_test <- copykat(rawmat = input_expr, id.type = "S", ngene.chr = 5, 
                        win.size = 25, KS.cut = 0.1, sam.name= sample_id, 
                        distance = "euclidean", norm.cell.names = normal_cell_names,
                        cell.line = "no", output.seg = "FLASE", 
                        plot.genes = "FLASE", genome = "hg20", n.cores = 20)
saveRDS(copykat_test, paste0(Out_Dir, "CopyKAT_res.RDS")) # copyKAT will automatically save the result, user can skip this step!
